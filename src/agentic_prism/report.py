"""Self-contained HTML and physical-size figures from immutable saved results."""
from pathlib import Path
from datetime import datetime, timezone
import html
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from .workflow import verify_run, dump, sha
from . import plot_style as pstyle
from .plot_style import THEMES  # noqa: F401  (re-exported for the other renderers)
from .report_shell import CSS, data_uri, style_select, switch_script, theme_figures  # noqa: F401
STATUS = {"estimated": "可报告的模型估计", "limited": "无法可靠报告精确 KD", "failed": "拟合失败",
          "within_range": "量程内", "below_range": "低于量程", "above_range": "高于量程", "unknown": "未确定"}
CI_STATUS = {"not_computed": "未计算", "noise_scale_not_estimable": "无法估计噪声尺度",
             "lower_open": "下端未闭合", "upper_open": "上端未闭合", "both_open": "两端均未闭合",
             "profile_failed": "区间计算失败", "two_sided": "已计算；精确结果仅保留供审计"}
SUMMARY_STATUS = {"independent_units_unconfirmed": "独立实验关系未确认，未汇总",
                  "comparability_unconfirmed": "实验条件可比性未确认，未汇总",
                  "withheld_incomplete_estimates": "包含失败或受限曲线，未汇总",
                  "one_experiment_no_between_experiment_ci": "仅一次实验，无实验间区间",
                  "summarized_log_t": "已按独立实验汇总"}
DIAG = {"ci_conditional_on_estimated_zero_baseline": "区间以零浓度对照估计的基线为固定条件，未传播基线误差。",
        "out_of_range_no_automatic_statistical_bound": "最优解超出浓度量程；不自动将最低/最高浓度当作统计界限。",
        "numerical_boundary_hit": "优化结果命中数值边界。", "rank_deficient": "参数局部不可辨识。",
        "limited_upper_plateau_model_based": "按模型估计，最高浓度未充分覆盖上平台。",
        "flat_response": "响应平坦，不能估计 KD。", "nonpositive_response_for_log": "非正响应不能进行 log 响应拟合。",
        "insufficient_distinct_concentrations_or_df": "有效浓度或残差自由度不足。"}


def font_setup():
    names = {f.name for f in font_manager.fontManager.ttflist}
    latin = next((n for n in ("Arial", "DejaVu Sans") if n in names), "DejaVu Sans")
    cjk = next((n for n in ("PingFang SC", "Heiti TC", "Arial Unicode MS", "Noto Sans CJK SC", "WenQuanYi Zen Hei") if n in names), None)
    return latin, cjk


def format_num(value, factor=1):
    return "—" if value is None else f"{value * factor:.3g}"


def render_report(run, style=None):
    run = Path(run).resolve()
    manifest = verify_run(run)
    cfg = json.loads((run / "config.resolved.json").read_text())
    results = json.loads((run / "results.json").read_text())
    d = pd.read_csv(run / "normalized_data.csv", keep_default_na=False, dtype={"curve_id": str, "sample_id": str})
    pred = pd.read_csv(run / "predictions.csv", dtype={"curve_id": str}, keep_default_na=False)
    resid = pd.read_csv(run / "residuals.csv", dtype={"curve_id": str}, keep_default_na=False)
    for col in ("residual_linear", "residual_log10"):
        resid[col] = pd.to_numeric(resid[col], errors="coerce")
    d["exclude"] = d.exclude.astype(str).str.lower().eq("true")
    r = cfg["report"]
    style = style or r["plot_style"]
    if style not in THEMES:
        raise ValueError("Unsupported plot style")
    figures = run / "figures"
    figures.mkdir(exist_ok=True)
    latin, cjk = font_setup()
    width = r["figure_width_mm"] / 25.4
    height = width * 1.2
    records, figures_meta, cards = [], [], []
    for index, fit in enumerate(results["fits"]):
        cid = fit["curve_id"]
        g = d[d.curve_id == cid]
        pp, rr = pred[pred.curve_id == cid], resid[resid.curve_id == cid]
        key = f"curve-{index + 1:03d}"
        for theme, token in THEMES.items():
            with plt.rc_context(pstyle.rc(token, latin, cjk)):
                fig = plt.figure(figsize=(width, height), layout="constrained")
                gs = fig.add_gridspec(2, 2, width_ratios=[1, 5], height_ratios=[2, 1])
                ax0, ax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
                res0, res = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
                good = g[~g.exclude]
                yvals = list(g.response.astype(float)) + pp.predicted_response.to_list()
                low, high = min(yvals), max(yvals)
                pad = max((high - low) * .12, abs(high) * .05, 1e-12)
                for axis, iszero in ((ax0, True), (ax, False)):
                    subset = good[good.concentration_M.eq(0) if iszero else good.concentration_M.gt(0)]
                    x = np.zeros(len(subset)) if iszero else subset.concentration_M * 1e9
                    axis.plot(x, subset.response, alpha=.85, label="Observed", **pstyle.point_style(token, 0, 4.5))
                    if len(subset):
                        agg = subset.groupby("concentration_M").response.agg(["mean", "std", "count"])
                        agg = agg[agg["count"] > 1]
                        if len(agg):
                            xx = np.zeros(len(agg)) if iszero else agg.index.to_numpy() * 1e9
                            axis.errorbar(xx, agg["mean"], yerr=agg["std"], fmt="none", color="black", capsize=3, lw=1)
                    excluded = g[g.exclude & (g.concentration_M.eq(0) if iszero else g.concentration_M.gt(0))]
                    if len(excluded):
                        axis.scatter(np.zeros(len(excluded)) if iszero else excluded.concentration_M * 1e9, excluded.response, marker="x", color=pstyle.EXCLUDED, s=24)
                    q = pp[pp.concentration_M.eq(0) if iszero else pp.concentration_M.gt(0)]
                    if len(q):
                        axis.plot(np.zeros(len(q)) if iszero else q.concentration_M * 1e9, q.predicted_response,
                                  "_" if iszero else "-", color=token["color"], lw=token["line_width"], ms=10, mew=token["line_width"])
                    axis.set_ylim(low - pad, high + pad)
                ax0.set_ylabel(f"Response ({fit['response_unit']})")
                for axis in (ax0, res0):
                    axis.set_xlim(-.8, .8)
                    axis.set_xticks([0], ["0"])
                for axis in (ax, res):
                    axis.set_xscale("log")
                    pos = g.loc[g.concentration_M > 0, "concentration_M"] * 1e9
                    if len(pos):
                        axis.set_xlim(pos.min() / 1.3, pos.max() * 1.3)
                    # The log panel shares the zero panel's y scale; one visible y axis reads as a broken x axis.
                    axis.spines["left"].set_visible(False)
                    axis.tick_params(axis="y", left=False, labelleft=False)
                ax.tick_params(axis="x", labelbottom=False)
                ax0.tick_params(axis="x", labelbottom=False)
                rk = "residual_log10" if fit["residual_scale"] == "log" else "residual_linear"
                for axis, iszero in ((res0, True), (res, False)):
                    q = rr[rr.concentration_M.eq(0) if iszero else rr.concentration_M.gt(0)]
                    axis.axhline(0, color=pstyle.MUTED, lw=.8, ls=(0, (3, 2)))
                    if len(q):
                        used = q.used_in_objective.astype(str).str.lower().eq("true")
                        qq = q[used]
                        axis.plot(np.zeros(len(qq)) if iszero else qq.concentration_M * 1e9, qq[rk], **pstyle.residual_style(token))
                        qq = q[~used]
                        axis.plot(np.zeros(len(qq)) if iszero else qq.concentration_M * 1e9, qq[rk], "x", color="#999999", ms=4)
                if len(rr) and rr[rk].notna().any():
                    m = max(float(rr[rk].abs().max()) * 1.2, 1e-8)
                    res0.set_ylim(-m, m)
                    res.set_ylim(-m, m)
                res0.set_ylabel("Residual\n(log10)" if rk == "residual_log10" else "Residual")
                res.set_xlabel("Concentration (nM)")
                res0.set_xlabel("Control")
                title = cid if len(cid) < 30 else cid[:27] + "…"
                fig.suptitle(title, fontsize=10.5, fontweight="bold" if token["bold_labels"] else "normal", x=.02, ha="left")
                figures_meta.append({"curve_id": cid, "theme": theme, "width_mm": r["figure_width_mm"],
                                     "height_mm": r["figure_width_mm"] * 1.2,
                                     "response_ylim": list(ax.get_ylim()), "concentration_xlim_nM": list(ax.get_xlim()),
                                     "n_observations": len(g), "prediction_rows": len(pp), "prediction_sha256": hashlib.sha256(pp.to_csv(index=False).encode()).hexdigest()})
                pstyle.save(fig, figures / f"{key}__{theme}", token, r["png_dpi"], r["export_formats"])
        report_value = format_num(fit["kd_M"], 1e9) if fit["reportable"] else "不报告精确值"
        # Open or out-of-range numerical endpoints remain audit data, not assay limits.
        ci = f"{format_num(fit['ci_low_M'], 1e9)} – {format_num(fit['ci_high_M'], 1e9)}" if fit["reportable"] and fit["ci_status"] == "two_sided" else CI_STATUS.get(fit["ci_status"], fit["ci_status"])
        labels = [DIAG.get(x, "KD 区间：" + CI_STATUS.get(x.removeprefix("interval_"), x)) if x.startswith("interval_") else DIAG.get(x, x) for x in fit["diagnostics"]]
        if fit["status"] == "failed":
            labels.insert(0, "这条曲线没有生成拟合线。")
        imgs = theme_figures(figures, key, style, f"{cid} 结合曲线与残差", THEMES, r["export_formats"])
        card = f'''<article id="{key}" class="curve"><div class="curve-head"><h3>{html.escape(cid)}</h3><span class="badge {fit['status']}">{STATUS[fit['status']]}</span></div>
<p class="metric">{html.escape(fit['interpretation'])} <strong>{report_value}</strong> {'nM' if fit['reportable'] else ''}</p>
<p>{cfg['uncertainty']['level']:.0%} KD 区间：{html.escape(ci)} · {STATUS[fit['range_status']]}</p>
<p class="muted">基线 {format_num(fit['baseline'])} · 振幅 {format_num(fit['amplitude'])} {html.escape(fit['response_unit'])}；此版本不计算这两个参数的区间。纳入 {fit['n_included']}/{fit['n_total']} 点。</p>
{imgs}<p class="caption">观测点与模型线；同曲线同浓度若有多个纳入读数，显示其均值 ± SD，拟合仍使用各读数。零浓度独立展示。残差 = 观测 − 预测（log 模式为 log10 之差）；灰叉表示未参与目标函数的残差。红叉为明确排除的观测。不显示曲线置信带。</p>
{''.join('<p class="notice">'+html.escape(s)+'</p>' for s in labels)}</article>'''
        cards.append(card)
        records.append(f'<tr><td><a href="#{key}">{html.escape(cid)}</a></td><td>{html.escape(fit["sample_id"])}</td><td>{report_value}</td><td>{html.escape(ci)}</td><td>{STATUS[fit["range_status"]]}</td></tr>')
    # Snapshot all data downloads inside the HTML for portable single-file use.
    download_names = ["input.csv", "normalized_data.csv", "fit_results.csv", "sample_summary.csv", "predictions.csv", "residuals.csv", "sensitivity.csv", "results.json", "config.resolved.json", "manifest.json", "diagnostics.json", "preprocessing_log.json", "rerun.txt"]
    if (run / "interpretation_facts.json").exists():
        download_names.append("interpretation_facts.json")
    downloads = " ".join(f'<a download="{name}" href="{data_uri(run/name)}">{name}</a>' for name in download_names)
    summary_rows = "".join(f'<tr><td>{html.escape(s["sample_id"])}</td><td>{s["n_experiments"]}</td><td>{format_num(s["geometric_mean_kd_M"], 1e9)}</td><td>{format_num(s["ci_low_M"], 1e9)} – {format_num(s["ci_high_M"], 1e9)}</td><td>{SUMMARY_STATUS.get(s["status"], s["status"])}</td></tr>' for s in results["summaries"])
    controls = style_select(THEMES, style) if r["allow_style_switch"] else f'<div class="toolbar"><span>图形风格：{html.escape(style)}</span></div>'
    template = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><link rel="icon" href="data:,"><title>AgenticPrism · 平衡 KD 分析</title>
<style>@@CSS@@</style></head><body><header><div class="eyebrow">AGENTICPRISM / EQUILIBRIUM BINDING</div><h1>平衡结合 · KD 分析</h1><p>@@COUNT@@ 条曲线 · @@REPORTABLE@@ 条可报告模型估计。原始观测、模型结果与限制共同呈现。</p><nav><a href="#overview">结果概览</a><a href="#curves">曲线与诊断</a><a href="#methods">方法与来源</a><a href="#downloads">下载与复现</a></nav></header><main>
<div class="panel">@@CONTROLS@@</div>
<section id="overview" class="card"><h2>结果概览</h2><div class="table-wrap"><table><thead><tr><th>曲线</th><th>样本</th><th>KD / nM</th><th>区间 / nM</th><th>量程</th></tr></thead><tbody>@@ROWS@@</tbody></table></div>
<p class="caption">不可可靠报告的数值仅保留在下载的审计结果中。超出量程不自动构成统计上界或下界。</p>
<details><summary>独立实验汇总</summary><p>仅在独立实验关系、条件可比性得到确认且全部曲线可报告时汇总：技术重复曲线先在实验内平均 log10(KD)，随后等权汇总各实验。实验间区间为 log 尺度 Student t 区间，不是单曲线参数区间。</p><div class="table-wrap"><table><tr><th>样本</th><th>已标识实验数</th><th>几何均值 / nM</th><th>实验间区间 / nM</th><th>状态</th></tr>@@SUMMARY@@</table></div></details></section>
<section id="curves"><h2>曲线与诊断</h2><div class="grid">@@CARDS@@</div></section>
<section id="methods" class="panel"><h2>方法与证据边界</h2><p><b>模型：</b>response = baseline + amplitude × C / (KD + C)。正振幅的单点结合模型；浓度内部单位 M。@@METHOD@@</p><p>@@ASSAY@@</p><p><b>来源：</b>@@SOURCE@@</p><p>固定基线模式的区间以所给基线为条件。profile-F 是指定误差模型下的近似区间，未完成与 Prism 参考项目的数值基准。参数收敛或区间闭合均不证明实验机制正确。</p><p>本报告只分析平衡结合；其他实验类型应选用各自的模块。Prism-like 是本项目的视觉预设。</p><details><summary>完整配置</summary><pre>@@CONFIG@@</pre></details><details><summary>敏感性分析（不替代主分析）</summary><pre>@@SENS@@</pre></details></section>
<section id="downloads" class="panel"><h2>下载与复现</h2><p>所有下载内容已内嵌；单个 HTML 文件可离线查看和下载。图形下载位于各曲线下，并随当前风格切换。</p><div class="downloads">@@DOWNLOADS@@</div><p class="caption">输入 SHA-256：@@HASH@@</p></section><footer>AgenticPrism @@VERSION@@ · @@FONT@@ · 数值分析与图形渲染分离</footer></main>
@@SCRIPT@@</body></html>'''
    replacements = {"COUNT": str(len(results["fits"])), "REPORTABLE": str(sum(f["reportable"] for f in results["fits"])),
                    "CONTROLS": controls, "ROWS": "".join(records), "CARDS": "".join(cards), "SUMMARY": summary_rows,
                    "METHOD": html.escape(f"残差尺度 {cfg['fit']['residual_scale']}；权重 {cfg['fit']['weighting']}；基线 {cfg['fit']['baseline_mode']}；KD 区间方法 {cfg['uncertainty']['parameter_ci']}。"),
                    "ASSAY": html.escape(cfg["assay"]["rationale"]), "SOURCE": html.escape(cfg["source"]),
                    "CONFIG": html.escape(json.dumps(cfg, ensure_ascii=False, indent=2)), "SENS": html.escape((run / "sensitivity.csv").read_text()),
                    "DOWNLOADS": downloads, "HASH": manifest["input_sha256"], "VERSION": manifest["package_version"],
                    "FONT": html.escape(f"图形字体 {latin}; CJK {cjk or 'unavailable (plots use Latin labels)'}"), "SCRIPT": switch_script(style), "CSS": CSS}
    for key, value in replacements.items():
        template = template.replace("@@" + key + "@@", value)
    (run / "report.html").write_text(template)
    dump(run / "render_manifest.json", {"created_utc": datetime.now(timezone.utc).isoformat(), "initial_style": style,
         "actual_fonts": {"latin": latin, "cjk": cjk}, "themes": THEMES, "figures": figures_meta,
         "scientific_manifest_sha256": sha(run / "manifest.json"),
         "report_sha256": sha(run / "report.html"), "renderer_sha256": sha(Path(__file__)),
         "figure_sha256": {p.name: sha(p) for p in sorted(figures.iterdir()) if p.is_file()}})
    verify_run(run)
    return run / "report.html"
