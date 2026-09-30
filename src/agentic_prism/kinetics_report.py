"""Kinetic sensorgrams and standalone reports, consuming saved results only."""
from pathlib import Path
from datetime import datetime,timezone
import json
import html
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import plot_style as pstyle
from .report import THEMES,font_setup,format_num
from .report_shell import CSS,data_uri,style_select,switch_script,theme_figures
from .workflow import verify_run,dump,sha

LABELS={"estimated":"模型估计","estimated_with_diagnostics":"模型估计 · 需要复核","limited":"参数受限，无法可靠报告","failed":"拟合失败"}
DIAGNOSTICS={
    "fit_window_sensitive":"解离时间窗缩短后，至少一个速率或 KD 改变超过两倍，或参数触及边界；暂不报告精确值。",
    "fit_window_check_failed":"时间窗敏感性检查优化失败；暂不报告精确值。",
    "fit_window_check_insufficient_points":"时间窗敏感性检查点数不足；暂不报告精确值。",
    "serially_correlated_residuals_review_model_and_block_length":"残差存在较强时间相关性。分块重采样区间仍依赖块长和误差假设；请复核模型、漂移及时间窗，不能仅按区间宽度判断精度。",
    "little_dissociation_in_observation_window":"拟合窗口内模型预测的解离变化很小，koff 可能缺乏充分信息。",
    "numerical_boundary_hit":"参数命中数值边界；数值估计仅用于审计。",
    "rates_not_locally_identifiable":"速率参数局部不可辨识。",
    "intervals_conditional_on_subtracted_baseline":"区间以预先扣除的基线为条件，未传播其独立估计误差。",
    "reference_uncertainty_not_separately_propagated":"参考通道误差未单独建模；扣除后的残差不能代替对参考不确定性的完整传播。",
    "requires_two_distinct_positive_concentrations":"首版全局拟合要求至少两个不同的正浓度。",
    "insufficient_association_or_dissociation_points":"结合或解离阶段的有效点不足。",
    "flat_trace":"存在平坦曲线，本拟合组未报告参数。",
    "bootstrap_block_length_exceeds_half_phase_length":"分块长度超过某阶段有效点数的一半，未计算区间。",
    "irregular_sampling_block_duration_varies":"采样间隔不规则，同样的重采样块长对应不同时长，区间需要额外复核。"}
CI_LABELS={"withheld_reliability_check":"可靠性检查未通过，区间不作报告（仅在审计记录中保留）","conditional_block_bootstrap":"条件分块 bootstrap 百分位区间","not_requested":"未请求区间",
           "noise_scale_not_estimable":"噪声尺度无法估计，未生成零宽区间","insufficient_points_for_block_length":"有效点不足以支持指定块长",
           "insufficient_valid_bootstrap_refits":"有效 bootstrap 重拟合不足","not_computed":"未计算"}


def render_kinetics(run,style=None):
    run=Path(run).resolve();manifest=verify_run(run)
    cfg=json.loads((run/"config.resolved.json").read_text())
    results=json.loads((run/"results.json").read_text())["fits"]
    d=pd.read_csv(run/"normalized_data.csv",keep_default_na=False,dtype={k:str for k in ("fit_group_id","curve_id","observation_id")})
    pp=pd.read_csv(run/"predictions.csv",keep_default_na=False,dtype={k:str for k in ("fit_group_id","curve_id","observation_id")})
    for frame in (d,pp):
        frame["used_in_fit"]=frame.used_in_fit.astype(str).str.lower().eq("true")
    r=cfg["report"];style=style or r["plot_style"]
    if style not in THEMES: raise ValueError("Unsupported kinetic plot style")
    folder=run/"figures";folder.mkdir(exist_ok=True)
    latin,cjk=font_setup();meta=[];cards=[]
    for idx,fit in enumerate(results,1):
        gid=fit["fit_group_id"];g=d[d.fit_group_id==gid];pr=pp[pp.fit_group_id==gid]
        key=f"kinetics-{idx:03d}"
        for theme,token in THEMES.items():
            with plt.rc_context(pstyle.rc(token,latin,cjk)):
                width=r["figure_width_mm"]/25.4
                fig,axes=plt.subplots(2,1,figsize=(width,width*.78),gridspec_kw={"height_ratios":[2.2,1]},sharex=True,layout="constrained")
                # Prism-like: coloured raw traces ordered by concentration, black model lines on top.
                prism=theme=="prism_like"
                curves=list(g.groupby("curve_id",sort=False))
                rank={cid:i for i,cid in enumerate(sorted((c for c,_ in curves),key=lambda c:float(g[g.curve_id==c].concentration_M.max())))}
                for j,(cid,h) in enumerate(curves):
                    color=(plt.get_cmap("viridis")(.08+.78*rank[cid]/max(len(curves)-1,1)) if prism else plt.get_cmap("tab10")(j%10))
                    fit_color="black" if prism else color
                    t=h.time_s-h.association_start_s.min()
                    steps=h.concentration_M.nunique()
                    label=(f"{cid} (single-cycle, {steps} injections)" if steps>1 else f"{float(h.concentration_M.iloc[0])*1e9:g} nM")
                    axes[0].plot(t,h.response_processed,color=color,lw=2.2 if prism else .65,alpha=.5 if prism else .38,
                                 solid_capstyle="round",label=label if prism else None)
                    q=pr[pr.curve_id==cid]
                    if len(q):
                        axes[0].plot(q.time_since_association_s,q.predicted_response,color=fit_color,lw=.8 if prism else 1.,ls="--",alpha=.6)
                        q=q[q.used_in_fit]
                        axes[0].plot(q.time_since_association_s,q.predicted_response,color=fit_color,lw=1. if prism else 1.5,
                                     label=None if prism else label)
                        axes[1].plot(q.time_since_association_s,q.residual,color=color,lw=.65)
                    elif not prism:
                        axes[0].plot([],[],color=color,label=label)
                if prism:
                    axes[0].plot([],[],color="black",lw=1.,label="1:1 global fit")
                first=g.groupby("curve_id").association_start_s.transform("min")
                boundaries=sorted(set((g.dissociation_start_s-first).to_list())|(set((g.association_start_s-first).to_list())-{0.}))
                for ax in axes:
                    for boundary in boundaries: ax.axvline(boundary,color=pstyle.MUTED,ls=(0,(3,2)),lw=.8)
                axes[1].axhline(0,color="black" if theme=="prism_like" else pstyle.MUTED,lw=.8)
                axes[0].set_ylabel(f"Response ({g.response_unit.iloc[0]})")
                axes[1].set_ylabel("Residual")
                axes[1].set_xlabel("Time since association start (s)")
                axes[0].set_title(gid[:60],loc="left")
                axes[0].legend(ncol=min(3,g.curve_id.nunique()+(theme=="prism_like")),fontsize=7.5,loc="best")
                meta.append({"fit_group_id":gid,"theme":theme,"time_limits_s":list(axes[0].get_xlim()),"response_limits":list(axes[0].get_ylim()),
                             "residual_limits":list(axes[1].get_ylim()),"n_input":len(g),"n_fitted":int(g.used_in_fit.sum()),
                             "prediction_sha256":hashlib.sha256(pr.to_csv(index=False).encode()).hexdigest(),
                             "width_mm":r["figure_width_mm"],"height_mm":r["figure_width_mm"]*.78})
                pstyle.save(fig,folder/f"{key}__{theme}",token,r["png_dpi"])
        rows=[]
        for param,label,unit,factor in (("kon_M_inv_s_inv","kon","M⁻¹ s⁻¹",1),("koff_s_inv","koff","s⁻¹",1),("kd_M","KD","nM",1e9)):
            bounds=fit["intervals"].get(param)
            ci="—" if bounds is None or not fit["reportable"] else " – ".join(format_num(x,factor) for x in bounds)
            val=format_num(fit[param],factor) if fit["reportable"] else "不报告精确值"
            rows.append(f"<tr><td>{label}</td><td>{val}</td><td>{unit}</td><td>{ci}</td></tr>")
        figures=theme_figures(folder,key,style,f"{gid} 动力学曲线和残差",THEMES)
        notices="".join('<p class="notice">'+html.escape(DIAGNOSTICS.get(diag,diag))+'</p>' for diag in fit["diagnostics"])
        nuisance="".join(f'<tr><td>{html.escape(c["curve_id"])}</td><td>{format_num(c["concentration_M"],1e9) if c.get("concentration_M") is not None else " / ".join(format_num(x,1e9) for x in c["injection_concentrations_M"])}</td><td>{format_num(c["rmax"])}</td><td>{format_num(c["offset"])}</td><td>{c["association_residual_lag1"]:.3f} / {c["dissociation_residual_lag1"]:.3f}</td></tr>' for c in fit["curve_parameters"])
        cards.append(f'''<article class="kinetic-group"><div class="head"><h2>{html.escape(gid)}</h2><span class="badge {fit['status']}">{LABELS[fit['status']]}</span></div>
<p>{fit['n_curves']} 条曲线 · 纳入拟合 {fit['n_fit']} / {fit['n_input']} 个原始点 · {cfg['uncertainty']['level']:.0%} {CI_LABELS.get(fit['ci_status'],fit['ci_status'])}</p>
<div class="table-wrap"><table><tr><th>参数</th><th>模型估计</th><th>单位</th><th>条件区间</th></tr>{''.join(rows)}</table></div>
{notices}{figures}<p class="caption">彩色线：全部输入观测（Prism-like 风格中按浓度着色，模型线为黑色）；实线：拟合窗口内模型；虚线：窗口外模型外推；竖线：进样与解离起点（时间从第一次进样起算）。残差仅显示实际纳入拟合的点，定义为观测 − 预测。绘图风格不改变数据或轴范围。</p>
<details><summary>逐传感器参数和残差相关性</summary><div class="table-wrap"><table><tr><th>曲线</th><th>浓度 / nM</th><th>Rmax</th><th>offset</th><th>结合 / 解离 lag-1</th></tr>{nuisance}</table></div></details></article>''')
    if cfg.get("steady_state", {}).get("enabled"):
        saved=json.loads((run/"steady_state.json").read_text())
        cards.append('<section class="panel"><h2>稳态亲和力（可选）</h2><p>KD_ss / KD_kin 仅作诊断，不取平均。逐曲线窗口、排除理由和报告状态如下。</p><pre>'+html.escape(json.dumps(saved,ensure_ascii=False,indent=2))+'</pre></section>')
    downloads=" ".join(f'<a download="{p.name}" href="{data_uri(p)}">{p.name}</a>' for p in sorted(run.iterdir()) if p.is_file() and (p.name in manifest["scientific_artifacts_sha256"] or p.name=="manifest.json"))
    doc='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,"><title>AgenticPrism · BLI/SPR 动力学</title><style>@@CSS@@</style></head><body><header><div class="eyebrow">AGENTICPRISM / BINDING KINETICS</div><h1>BLI / SPR · 结合动力学</h1><p>多浓度 1:1 全局拟合（独立循环或单循环进样序列）。参数估计、时间窗和模型限制共同呈现；KD 来自同一次拟合中的 koff / kon。</p></header><main><div class="panel">@@CONTROLS@@<p class="caption">@@WINDOW@@</p></div>@@CARDS@@
<section class="panel"><h2>方法、来源与限制</h2><p>@@SOURCE@@</p><p>@@ASSAY@@</p><p>每个拟合组共享 kon/koff；@@NUISANCE@@。模型在结合与解离阶段连续，假设第一次进样起点的特异性结合量为零。@@DESIGN@@不支持多价/异质性模型或质量传输模型。</p><p>区间方法：@@UNCERTAINTY@@。每条曲线各阶段独立进行连续残差块重采样，并对每次联合重拟合计算 koff/kon；不把两个边际区间相除。这是依赖模型与块长的条件区间，不包含浓度误差或模型错设的不确定性。没有完成普适覆盖率或原厂软件等价验证。</p><p>已验证导入格式为公开 Octet RED384 数据中的 Results.txt（Time1/Data1）布局；不代表支持所有 Octet 版本、.frd 或 Biacore 原生项目。SPR 可在明确阶段、单位和预处理后使用标准 CSV；尚未完成 Biacore 专用导入验收。</p><details><summary>完整配置</summary><pre>@@CONFIG@@</pre></details></section>
<section class="panel"><h2>下载与复现</h2><p>图形和以下数据文件均已内嵌，报告可单文件离线使用。</p><div class="downloads">@@DOWNLOADS@@</div><p class="caption">输入 SHA-256：@@HASH@@</p></section><footer>AgenticPrism @@VERSION@@ · 数值分析与渲染分离 · 图形字体 @@FONT@@</footer></main>@@SCRIPT@@</body></html>'''
    design=("单循环设计：同一表面依次进样递增浓度、中间不再生，结合量在各步之间按解析 1:1 解连续传递；时间窗敏感性检查只截短最后一段解离。"
            if cfg["assay"].get("injection_design")=="single_cycle" else "独立循环设计：每条曲线一次进样，循环之间已再生。")
    reference={"none":"未做参比扣除。","reference_column":"已逐行扣除同循环参比通道（reference_response）。",
               "double_reference":"双参比：逐行计算 (样品 − 参比通道) − (缓冲液空白 − 空白参比通道)；各通道误差未单独传播。"}[cfg["preprocessing"]["reference_mode"]]
    values={"DESIGN":html.escape(design+reference),"WINDOW":html.escape(f"结合起点后跳过 {cfg['fit']['association_skip_s']} s；解离起点后跳过 {cfg['fit']['dissociation_skip_s']} s；解离拟合时长 {cfg['fit']['dissociation_duration_s'] if cfg['fit']['dissociation_duration_s'] is not None else '使用所给数据'}；显式抽样步长 {cfg['preprocessing']['stride']}。"),
            "CARDS":"".join(cards),"SOURCE":html.escape(cfg["source"]),"ASSAY":html.escape(cfg["assay"]["rationale"]),
            "NUISANCE":html.escape(f"Rmax：{cfg['fit']['rmax']}，offset：{cfg['fit']['offset']}"),
            "UNCERTAINTY":html.escape(f"{cfg['uncertainty']['method']}，{cfg['uncertainty']['replicates']} 次，块长 {cfg['uncertainty']['block_length']} 个纳入时间点，随机种子 {cfg['uncertainty']['seed']}"),
            "CONFIG":html.escape(json.dumps(cfg,ensure_ascii=False,indent=2)),"DOWNLOADS":downloads,
            "HASH":manifest["input_sha256"],"VERSION":manifest["package_version"],"FONT":html.escape(latin),"CONTROLS":style_select(THEMES,style),"SCRIPT":switch_script(style),"CSS":CSS}
    for k,v in values.items(): doc=doc.replace("@@"+k+"@@",v)
    (run/"report.html").write_text(doc)
    dump(run/"render_manifest.json",{"analysis_type":"binding_kinetics","created_utc":datetime.now(timezone.utc).isoformat(),"initial_style":style,
         "fonts":{"latin":latin,"cjk":cjk},"figures":meta,"scientific_manifest_sha256":sha(run/"manifest.json"),
         "renderer_sha256":sha(Path(__file__)),"report_sha256":sha(run/"report.html"),
         "figure_sha256":{p.name:sha(p) for p in sorted(folder.iterdir())}})
    verify_run(run)
    return run/"report.html"
