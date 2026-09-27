"""Offline 4PL dose-response report from saved, hash-verified results."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import html
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .dose_schema import DOSE_UNITS
from .report import THEMES, data_uri, font_setup, format_num
from .workflow import dump, sha, verify_run

SECOND_COLOR = "#b5651d"  # test curve in comparison overlays; reference uses the theme color
DIAG = {
    "insufficient_distinct_doses_or_residual_df": "有效浓度或残差自由度不足；4PL 未拟合。",
    "flat_response": "响应基本平坦，无法估计半效浓度。",
    "optimizer_failed": "优化未收敛。",
    "numerical_boundary_hit": "参数触及数值边界，精确值仅供审计。",
    "midpoint_or_hill_not_locally_identifiable": "半效浓度或斜率局部不可辨识。",
    "half_response_outside_tested_range": "半效浓度位于实测正浓度范围外，不能报告精确外推值。",
    "low_response_range_relative_to_noise": "响应变化相对残差较小。",
    "plateau_not_well_covered": "实测浓度未充分覆盖两个平台；半效浓度依赖平台外推。",
    "few_transition_doses_for_variable_slope": "10%–90% 转换区的浓度点较少，变斜率参数需复核。",
    "zero_controls_dominate_unweighted_fit": "零剂量对照超过有效观测的四分之一；未加权拟合可能受其影响。",
    "profile_interval_open": "参数区间未在设置的搜索范围内闭合。",
    "profile_hill_nuisance_at_boundary": "区间端点处 Hill 斜率触及约束，需复核。",
    "observed_direction_opposite_to_declared": "实测响应随浓度的变化方向与配置声明的方向相反（Spearman 检验 p<0.05）；请核对 direction 设置或数据，不报告该曲线。",
    "nonpositive_response_for_relative_weighting": "相对加权（1/Y²）要求所有响应为正；该曲线未拟合。",
    "individual_fit_failed": "至少一条曲线单独拟合失败，无法比较。",
    "parallel_optimizer_failed": "平行模型优化未收敛。",
    "parallel_model_not_locally_identifiable": "平行模型参数局部不可辨识。",
    "nonparallel_by_f_test": "平行性 F 检验拒绝共享平台与斜率（p 低于设定 α）；不报告相对效价。",
    "parallelism_equivalence_not_demonstrated": "未能证明平行性等效：至少一个参数差异的区间超出预设等效界限；不报告相对效价。",
    "individual_fit_not_reportable": "至少一条曲线的单独拟合不可报告；相对效价仅供审计。",
    "curves_from_different_experiments": "两条曲线来自不同实验；相对效价混入了实验间差异。",
}
SUMMARY_STATUS = {"independent_units_unconfirmed": "独立实验关系未确认，未汇总",
                  "comparability_unconfirmed": "实验条件可比性未确认，未汇总",
                  "mixed_dose_units_not_summarized": "剂量单位族不一致，未汇总",
                  "withheld_incomplete_estimates": "包含失败或受限曲线，未汇总",
                  "cross_experiment_comparison_not_summarized": "含跨实验比较，未汇总",
                  "multiple_comparisons_per_experiment_not_summarized": "同一实验内多次比较，未汇总",
                  "one_experiment_no_between_experiment_ci": "仅一次实验，无实验间区间",
                  "summarized_log_t": "已按独立实验汇总（几何均值，log t 区间）"}


def _style(font, cjk, token):
    return {"font.family": [font] + ([cjk] if cjk else []), "font.size": 9, "axes.labelsize": 10,
            "axes.titlesize": 12, "axes.linewidth": token["axes_width"], "xtick.direction": "out",
            "ytick.direction": "out", "pdf.fonttype": 42, "svg.fonttype": "path"}


def _save(fig, folder, key, theme, dpi):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(folder / f"{key}__{theme}.{ext}", dpi=dpi, facecolor="white")
    plt.close(fig)


def _figure_block(folder, key, style, alt):
    blocks = []
    for theme in THEMES:
        exports = " ".join(f'<a download="{key}__{theme}.{ext}" href="{data_uri(folder / f"{key}__{theme}.{ext}")}">{ext.upper()}</a>'
                           for ext in ("svg", "pdf", "png"))
        blocks.append(f'<div class="theme-figure" data-theme="{theme}" {"hidden" if theme != style else ""}>'
                      f'<img alt="{html.escape(alt)}" src="{data_uri(folder / f"{key}__{theme}.svg")}">'
                      f'<div class="downloads">{exports}</div></div>')
    return "".join(blocks)


def _fmt_range(low, high, factor=1., unit=""):
    if low is None or high is None:
        return "—"
    return f"{format_num(low, factor)}–{format_num(high, factor)}{' ' + unit if unit else ''}"


def render_dose(run, style=None):
    run = Path(run).resolve()
    manifest = verify_run(run)
    cfg = json.loads((run / "config.resolved.json").read_text())
    results = json.loads((run / "results.json").read_text())
    fits = results["fits"]
    d = pd.read_csv(run / "normalized_data.csv", keep_default_na=False, dtype={"curve_id": str})
    pred = pd.read_csv(run / "predictions.csv", keep_default_na=False, dtype={"curve_id": str})
    grid = pd.read_csv(run / "curve_grid.csv", keep_default_na=False, dtype={"curve_id": str})
    comparison_grid = (pd.read_csv(run / "comparison_grid.csv", keep_default_na=False, dtype={"curve_id": str, "comparison_id": str})
                       if (run / "comparison_grid.csv").exists() else pd.DataFrame())
    for frame in (d, pred):
        if "exclude" in frame:
            frame["exclude"] = frame.exclude.astype(str).str.lower().eq("true")
        if "used_in_fit" in frame:
            frame["used_in_fit"] = frame.used_in_fit.astype(str).str.lower().eq("true")
    style = style or cfg["report"]["plot_style"]
    if style not in THEMES:
        raise ValueError("Unsupported 4PL plot style")
    font, cjk = font_setup()
    folder = run / "figures"
    folder.mkdir(exist_ok=True)
    width = cfg["report"]["figure_width_mm"] / 25.4
    dpi = cfg["report"]["png_dpi"]
    cards, figure_meta = [], []
    for number, fit in enumerate(fits, 1):
        curve_id = fit["curve_id"]
        obs = d[d.curve_id == curve_id]
        predicted = pred[pred.curve_id == curve_id]
        gg = grid[grid.curve_id == curve_id]
        unit = fit["input_unit"]
        factor = DOSE_UNITS[unit]  # canonical / factor = value in the input unit
        response_unit = obs.response_unit.iloc[0]
        dose = obs.concentration_canonical.to_numpy(float) / factor
        zero = dose == 0
        has_zero = bool(zero.any())
        positive = dose[~zero]
        x_min = float(positive.min()) if len(positive) else 1.
        x_max = float(positive.max()) if len(positive) else 10.
        x_limits = [x_min / 1.7, x_max * 1.7]
        y_values = np.r_[obs.response.to_numpy(float), gg.predicted_response.to_numpy(float)]
        y_min, y_max = float(np.min(y_values)), float(np.max(y_values))
        pad = max((y_max - y_min) * .12, 1e-6)
        y_limits = [y_min - pad, y_max + pad]
        residuals = predicted[predicted.used_in_fit]
        r_max = max(float(np.max(np.abs(residuals.residual.to_numpy(float)))) * 1.2, 1e-6) if len(residuals) else max(pad, 1.)
        key = f"dose-{number:03d}"
        included = ~obs.exclude.to_numpy(bool)
        y_obs = obs.response.to_numpy(float)
        for theme, token in THEMES.items():
            with plt.rc_context(_style(font, cjk, token)):
                fig = plt.figure(figsize=(width, width * .78), layout="constrained")
                if has_zero:
                    # Zero dose has no log coordinate: it gets its own linear panel.
                    gs = fig.add_gridspec(2, 2, width_ratios=[1, 7], height_ratios=[2.2, 1.])
                    z_ax, ax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
                    z_res, res = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1], sharex=ax)
                else:
                    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.])
                    ax = fig.add_subplot(gs[0, 0])
                    res = fig.add_subplot(gs[1, 0], sharex=ax)
                    z_ax = z_res = None
                panels = [(ax, res, ~zero)] + ([(z_ax, z_res, zero)] if has_zero else [])
                for main, resid_ax, mask in panels:
                    x_plot = np.zeros(mask.sum()) if main is z_ax else dose[mask]
                    keep = included[mask]
                    main.scatter(x_plot[keep], y_obs[mask][keep], s=20, color=token["color"], alpha=.75,
                                 label="Included observations")
                    if (~keep).any():
                        main.scatter(x_plot[~keep], y_obs[mask][~keep], s=25, facecolors="none",
                                     edgecolors="#999999", label="Excluded")
                    rows = residuals[(residuals.concentration_canonical == 0) == (main is z_ax)]
                    rx = np.zeros(len(rows)) if main is z_ax else rows.concentration_canonical.to_numpy(float) / factor
                    resid_ax.scatter(rx, rows.residual, s=13, color=token["color"], alpha=.7)
                    resid_ax.axhline(0, color="#777777", lw=.8)
                    main.set_ylim(*y_limits)
                    resid_ax.set_ylim(-r_max, r_max)
                    for a in (main, resid_ax):
                        a.spines[["top", "right"]].set_visible(False)
                        a.grid(True, alpha=.18) if token["grid"] else a.grid(False)
                if len(gg):
                    ax.plot(gg.concentration_canonical / factor, gg.predicted_response, color=token["color"],
                            lw=token["line_width"] + .2, label="4PL fit")
                for a in (ax, res):
                    a.set_xscale("log")
                    a.set_xlim(*x_limits)
                res.set_xlabel(f"Concentration ({unit}; log scale)")
                ax.set_title(curve_id[:55], loc="left")
                ax.legend(frameon=False, fontsize=7, loc="best")
                if has_zero:
                    at_zero = predicted[(predicted.concentration_canonical == 0)]
                    if len(at_zero):
                        z_ax.plot([0.], [at_zero.predicted_response.iloc[0]], "_", color=token["color"], ms=16, mew=1.7)
                    for a in (z_ax, z_res):
                        a.set_xlim(-.8, .8)
                        a.set_xticks([0], ["0"])
                    z_res.set_xlabel("Zero dose")
                    ax.tick_params(labelleft=False)
                    res.tick_params(labelleft=False)
                    z_ax.set_ylabel(f"Response ({response_unit})")
                    z_res.set_ylabel("Residual")
                else:
                    ax.set_ylabel(f"Response ({response_unit})")
                    res.set_ylabel("Residual")
                plt.setp(ax.get_xticklabels(), visible=False)
                _save(fig, folder, key, theme, dpi)
                figure_meta.append({"curve_id": curve_id, "theme": theme, "x_unit": unit, "zero_dose_panel": has_zero,
                                    "x_limits": x_limits, "y_limits": y_limits, "residual_limits": [-r_max, r_max],
                                    "prediction_sha256": hashlib.sha256(gg.to_csv(index=False).encode()).hexdigest()})
        midpoint = format_num(fit["half_response_input_unit"]) + " " + unit
        if not fit["reportable"]:
            midpoint = "不报告精确值"
        ci = fit["ci_canonical"]
        interval = _fmt_range(ci[0], ci[1], 1. / factor, unit) if fit["reportable"] else "—"
        notices = "".join(f'<p class="notice">{html.escape(DIAG.get(code, code))}</p>' for code in fit["diagnostics"])
        constraint = _constraint_text(fit.get("fixed_bottom"), fit.get("fixed_top"))
        weighting = "相对加权 1/Ŷ²" if fit.get("weighting") == "relative" else "未加权"
        cards.append(f'''<article class="dose-curve"><div class="head"><h2>{html.escape(curve_id)}</h2>
<span class="badge {fit['status']}">{'可报告的模型估计' if fit['reportable'] else '受限或失败'}</span></div>
<p>{html.escape(fit['sample_id'])} · {html.escape(fit['experiment_id'])} · 有效观测 {fit['n_fit']} / {fit['n_input']} · {fit['n_distinct_positive_doses']} 个正浓度 · {weighting} · {constraint} · {html.escape(fit['ci_status'])}</p>
<div class="table-wrap"><table><tr><th>参数</th><th>估计</th><th>区间 / 单位</th></tr>
<tr><td>{fit['endpoint']}（相对平台中点）</td><td>{midpoint}</td><td>{interval}</td></tr>
<tr><td>Bottom / Top</td><td>{format_num(fit['bottom'])} / {format_num(fit['top'])}</td><td>{html.escape(response_unit)}</td></tr>
<tr><td>Hill slope</td><td>{format_num(fit['hill_slope_signed'])}</td><td>无量纲</td></tr>
<tr><td>RMSE（响应尺度）</td><td>{format_num(fit['rmse'])}</td><td>{html.escape(response_unit)}</td></tr></table></div>
{notices}{_figure_block(folder, key, style, curve_id + " 4PL 曲线及残差")}<p class="caption">点为逐条原始观测，线为 4PL；下图残差为观测减预测。横轴使用输入单位 {html.escape(unit)}。零剂量对照如存在，放在左侧单独的线性小图中（短横线为模型在零剂量的取值），不画到虚构的对数浓度上。风格切换不改变数值或坐标范围。</p></article>''')

    summaries = results.get("summaries", [])
    summary_html = ""
    if summaries:
        rows = "".join(f'<tr><td>{html.escape(s["sample_id"])}</td><td>{s["n_curves"]}</td><td>{s["n_experiments"]}</td>'
                       f'<td>{format_num(s["geometric_mean_input_unit"])} {html.escape(s["input_unit"] or "") if s["geometric_mean_input_unit"] is not None else ""}</td>'
                       f'<td>{_fmt_range(s["ci_low_input_unit"], s["ci_high_input_unit"])}</td>'
                       f'<td>{html.escape(SUMMARY_STATUS.get(s["status"], s["status"]))}</td></tr>' for s in summaries)
        summary_html = (f'<section class="panel"><h2>样本汇总（{html.escape(cfg["assay"]["endpoint"])}）</h2><div class="table-wrap"><table>'
                        f'<tr><th>样本</th><th>曲线</th><th>独立实验</th><th>几何均值</th><th>{cfg["uncertainty"]["level"]:.0%} 区间</th><th>状态</th></tr>{rows}</table></div>'
                        '<p class="caption">同一实验内的技术重复曲线先在 log10 C50 上平均，再对独立实验等权；区间为实验间 log 尺度 Student-t 区间，不传播单条曲线的拟合区间。只要有一条曲线失败或受限，整个样本不汇总。独立实验和条件可比性须在配置中明确声明。</p></section>')

    comparison_html = ""
    comparisons = results.get("comparisons", [])
    for number, c in enumerate(comparisons, 1):
        key = f"compare-{number:03d}"
        ref_fit = next(f for f in fits if f["curve_id"] == c["reference_curve"])
        unit = ref_fit["input_unit"]
        factor = DOSE_UNITS[unit]
        cg = comparison_grid[comparison_grid.comparison_id == c["comparison_id"]] if len(comparison_grid) else comparison_grid
        for theme, token in THEMES.items():
            with plt.rc_context(_style(font, cjk, token)):
                fig, ax = plt.subplots(figsize=(width, width * .55), layout="constrained")
                for curve_id, color, label in ((c["reference_curve"], token["color"], "Reference"),
                                               (c["test_curve"], SECOND_COLOR, "Test")):
                    o = d[(d.curve_id == curve_id) & (d.concentration_canonical > 0) & ~d.exclude]
                    ax.scatter(o.concentration_canonical / factor, o.response, s=16, color=color, alpha=.7,
                               label=f"{label}: {curve_id[:30]}")
                    q = cg[cg.curve_id == curve_id] if len(cg) else cg
                    if len(q):
                        ax.plot(q.concentration_canonical / factor, q.predicted_response, color=color, lw=token["line_width"])
                ax.set_xscale("log")
                ax.set_xlabel(f"Concentration ({unit}; log scale)")
                ax.set_ylabel(f"Response ({d[d.curve_id == c['reference_curve']].response_unit.iloc[0]})")
                ax.set_title(c["comparison_id"][:55], loc="left")
                ax.spines[["top", "right"]].set_visible(False)
                ax.grid(True, alpha=.18) if token["grid"] else ax.grid(False)
                ax.legend(frameon=False, fontsize=7, loc="best")
                _save(fig, folder, key, theme, dpi)
                figure_meta.append({"comparison_id": c["comparison_id"], "theme": theme, "x_unit": unit})
        pm, par, shared = c.get("parallel_model") or {}, c.get("parallelism_f_test") or {}, c.get("shared_c50_f_test") or {}
        rp = format_num(c["relative_potency"]) if c["reportable"] else "不报告精确值"
        rp_ci = _fmt_range(*c["rp_ci"]) if c["reportable"] else "—"
        notices = "".join(f'<p class="notice">{html.escape(DIAG.get(code, code))}</p>' for code in c["diagnostics"])
        f_row = lambda t: (f'F({t["df_numerator"]}, {t["df_denominator"]}) = {t["F"]:.3g}, p = {t["p_value"]:.3g}' if t else "—")
        eq = c.get("parallelism_equivalence")
        equivalence_rows = ""
        if eq:
            label = {"hill_ratio_test_over_reference": "Hill 斜率比（待测/参比）", "bottom_difference_test_minus_reference": "Bottom 差（待测−参比）",
                     "top_difference_test_minus_reference": "Top 差（待测−参比）"}
            for t in eq["tests"]:
                verdict = "界限内" if t["within"] else "未检验" if t["within"] is None else "超出界限"
                equivalence_rows += (f'<tr><td>平行性等效：{label.get(t["parameter"], t["parameter"])}</td><td>{format_num(t["estimate"])}，'
                                     f'{eq["confidence_level"]:.0%} 区间 {_fmt_range(*t["ci"])}</td><td>预设界限 {_fmt_range(*t["limits"])}：{verdict}</td></tr>')
            equivalence_rows += (f'<tr><td>平行性等效结论</td><td>{"等效（全部区间在界限内）" if eq["equivalent"] else "未证明等效"}</td>'
                                 f'<td>各自独立拟合、合并残差方差的 Wald 区间；界限来源：{html.escape(eq["rationale"])}</td></tr>')
        acc = c.get("rp_acceptance")
        if acc:
            verdict = "区间整体在界限内" if acc["within"] else "无可报告区间，不判定为通过" if acc["within"] is None else "区间超出界限"
            equivalence_rows += (f'<tr><td>RP 接受界限</td><td>{_fmt_range(*acc["limits"])}</td>'
                                 f'<td>{verdict}；依据：{html.escape(acc["rationale"])}</td></tr>')
        comparison_html += f'''<article class="dose-curve"><div class="head"><h2>{html.escape(c["comparison_id"])}</h2>
<span class="badge {c['status']}">{'可报告的相对效价' if c['reportable'] else '受限或失败'}</span></div>
<p>参比 {html.escape(c["reference_curve"])}（{html.escape(c["reference_sample"])}）· 待测 {html.escape(c["test_curve"])}（{html.escape(c["test_sample"])}）</p>
<div class="table-wrap"><table><tr><th>量</th><th>估计</th><th>说明</th></tr>
<tr><td>相对效价 RP</td><td>{rp}</td><td>{cfg["uncertainty"]["level"]:.0%} profile-F 区间 {rp_ci}；RP = C50参比 / C50待测，&gt;1 表示待测更强</td></tr>
<tr><td>平行模型 C50（参比 / 待测）</td><td>{format_num(pm.get("reference_c50_canonical"), 1 / factor)} / {format_num(pm.get("test_c50_canonical"), 1 / factor)} {html.escape(unit)}</td><td>共享 Bottom {format_num(pm.get("bottom"))}、Top {format_num(pm.get("top"))}、Hill {format_num(pm.get("hill_slope_signed"))}</td></tr>
<tr><td>平行性 F 检验{"（仅供参考）" if eq else ""}</td><td>{f_row(par)}</td><td>平行模型 vs 各自独立拟合；{"判定依据为下方等效性界限" if eq else f'α = {cfg["comparison_settings"]["parallelism_alpha"]}'}</td></tr>
{equivalence_rows}<tr><td>C50 是否不同</td><td>{f_row(shared)}</td><td>共享 C50（平台、斜率各自）vs 各自独立拟合</td></tr>
<tr><td>单独拟合 C50 之比</td><td>{format_num(c["individual_c50_ratio"])}</td><td>描述性，不带区间</td></tr></table></div>
{notices}{_figure_block(folder, key, style, c["comparison_id"] + " 平行模型比较")}<p class="caption">点为两条曲线的正剂量观测，线为共享平台与斜率的平行 4PL。F 检验是统计显著性检验，不是等效性检验：精密的实验可能以很小的非平行拒绝平行性，粗糙的实验也可能无法发现真实的非平行。</p></article>'''
    potency = results.get("potency_summaries", [])
    if comparisons:
        rows = "".join(f'<tr><td>{html.escape(p["test_sample"])} vs {html.escape(p["reference_sample"])}</td><td>{p["n_experiments"]}</td>'
                       f'<td>{format_num(p["geometric_mean_rp"])}</td><td>{_fmt_range(p["ci_low"], p["ci_high"])}</td>'
                       f'<td>{html.escape(SUMMARY_STATUS.get(p["status"], p["status"]))}'
                       + ('' if "rp_acceptance_within" not in p else
                          f'；接受界限 {_fmt_range(*p["rp_acceptance_limits"])}：' + ("区间在界限内" if p["rp_acceptance_within"] else
                          "无区间，不判定" if p["rp_acceptance_within"] is None else "区间超出界限"))
                       + '</td></tr>' for p in potency)
        comparison_html = (f'<section class="panel"><h2>曲线比较与相对效价</h2><div class="table-wrap"><table><tr><th>比较</th><th>独立实验</th>'
                           f'<th>几何均值 RP</th><th>{cfg["uncertainty"]["level"]:.0%} 区间</th><th>状态</th></tr>{rows}</table></div>'
                           '<p class="caption">每次比较只用同一实验内的参比与待测曲线；跨实验汇总时对各实验的 log RP 等权取均值并给出 Student-t 区间。任一比较不可报告则整组不汇总。</p></section>'
                           + comparison_html)

    downloads = " ".join(f'<a download="{p.name}" href="{data_uri(p)}">{p.name}</a>' for p in sorted(run.iterdir())
                         if p.is_file() and (p.name in manifest["scientific_artifacts_sha256"] or p.name == "manifest.json"))
    fit_cfg = cfg["fit"]
    weighting = ("相对加权：最小化 Σ((Y−Ŷ)/Ŷ)²，Ŷ 为曲线值（Prism 的 1/Y² 相对加权），适用于噪声与信号成比例的读数"
                 if fit_cfg.get("weighting") == "relative" else "原始响应尺度的未加权最小二乘")
    doc = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,"><title>AgenticPrism · 4PL 剂量反应</title><style>
*{box-sizing:border-box}body{margin:0;background:#f3f6f6;color:#19383e;font:15px/1.7 Arial,"PingFang SC",sans-serif}header{background:#143e45;color:white;padding:40px max(5vw,20px)}header h1{font-size:32px;margin:10px 0}header p{max-width:950px;color:#d3e4e5}.eyebrow{letter-spacing:.15em;font-size:12px}main{max-width:1120px;margin:auto;padding:24px 18px}.panel,.dose-curve{background:white;border:1px solid #dbe5e4;border-radius:10px;padding:25px;margin:22px 0}.head{display:flex;gap:14px;justify-content:space-between;align-items:center;flex-wrap:wrap}h2{font-size:23px;margin:6px 0}p,h2,footer{overflow-wrap:anywhere}.badge{font-size:12px;background:#e5f3ee;padding:5px 9px;border-radius:4px}.limited,.failed,.estimated_with_diagnostics{background:#fff0dd;color:#805321}.notice{padding:10px 14px;border-left:3px solid #c59c58;background:#fff9f0;font-size:13px}.caption,footer{font-size:12px;color:#62767a}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:10px;border-bottom:1px solid #e3ebea;white-space:nowrap}th{background:#f4f8f7}img{width:100%;height:auto;display:block}a{color:#147b80;text-decoration:none}.downloads{display:flex;flex-wrap:wrap;gap:14px;font-size:13px}select{padding:8px;border:1px solid #bdcece;border-radius:4px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}details{margin:16px 0}[hidden]{display:none!important}@media(max-width:700px){main{padding:14px 10px}.panel,.dose-curve{padding:16px}header h1{font-size:26px}}@media print{select,.downloads{display:none}body{background:white}.dose-curve{break-inside:avoid}}
</style></head><body><header><div class="eyebrow">AGENTICPRISM / DOSE RESPONSE</div><h1>4PL · 剂量反应</h1><p>EC50 / IC50 是拟合 Bottom 与 Top 之间的相对半效浓度；标签取决于实验目的，而非仅由曲线上下方向决定。</p></header><main><div class="panel"><label>图形风格 <select id="style"><option value="prism_like">Prism-like</option><option value="standard">标准</option></select></label><p>@@ASSAY@@</p></div>@@SUMMARY@@@@COMPARISONS@@@@CARDS@@
<section class="panel"><h2>方法、来源与限制</h2><p>@@SOURCE@@</p><p>@@RATIONALE@@</p><p>模型：Y = Bottom + (Top−Bottom) / [1 + 10^(−signed Hill × (log10 C−log10 C50))]。Top 始终是较高的平台。各曲线独立拟合，保留逐孔观测；目标函数：@@WEIGHTING@@。平台约束：@@CONSTRAINT@@；固定的平台不计入参数个数。没有自动归一化、参考扣除或异常值删除。</p><p>区间：@@CI@@。若未覆盖平台、转换区点数不足、半效浓度超出实测浓度范围或参数区间未闭合，需结合诊断判断；曲线形状不能证明分子机制。相对 IC50 不等于固定响应值 50 的绝对 IC50（除非平台被固定为 0 与 100），也不等于 KD 或 Ki。</p><p>本模块未完成 Prism 数值等价、真实团队抗体实验或 5PL、双相模型验收。若单位为 source_unit，来源未提供物理浓度单位，结果仅作数值示例。</p><details><summary>完整配置</summary><pre>@@CONFIG@@</pre></details></section>
<section class="panel"><h2>下载与复现</h2><p>图形和数据已内嵌，可单文件离线查看。</p><div class="downloads">@@DOWNLOADS@@</div><p class="caption">输入 SHA-256：@@HASH@@</p></section><footer>AgenticPrism @@VERSION@@ · 数值分析与渲染分离</footer></main><script>const s=document.getElementById('style');function change(){document.querySelectorAll('[data-theme]').forEach(e=>e.hidden=e.dataset.theme!==s.value);}s.value=@@STYLE@@;s.addEventListener('change',change);change();</script></body></html>'''
    replacements = {"ASSAY": html.escape(cfg["assay"]["response_definition"]), "SUMMARY": summary_html,
                    "COMPARISONS": comparison_html, "CARDS": "".join(cards),
                    "SOURCE": html.escape(cfg["source"]), "RATIONALE": html.escape(cfg["assay"]["rationale"]),
                    "WEIGHTING": html.escape(weighting),
                    "CONSTRAINT": html.escape(_constraint_text(fit_cfg.get("fixed_bottom"), fit_cfg.get("fixed_top"))),
                    "CI": html.escape(f"{cfg['uncertainty']['method']}，水平 {cfg['uncertainty']['level']:.0%}；在每个候选 log10(C50) 上重新优化 Hill 斜率与平台，使用估计噪声尺度的 profile-F 阈值"),
                    "CONFIG": html.escape(json.dumps(cfg, ensure_ascii=False, indent=2)),
                    "DOWNLOADS": downloads, "HASH": manifest["input_sha256"],
                    "VERSION": manifest["package_version"], "STYLE": json.dumps(style)}
    for key, value in replacements.items():
        doc = doc.replace("@@" + key + "@@", value)
    (run / "report.html").write_text(doc)
    dump(run / "render_manifest.json", {"analysis_type": "dose_response_4pl",
         "created_utc": datetime.now(timezone.utc).isoformat(), "initial_style": style,
         "figures": figure_meta, "scientific_manifest_sha256": sha(run / "manifest.json"),
         "renderer_sha256": sha(Path(__file__)), "report_sha256": sha(run / "report.html"),
         "figure_sha256": {p.name: sha(p) for p in sorted(folder.iterdir())}})
    verify_run(run)
    return run / "report.html"


def _constraint_text(fixed_bottom, fixed_top):
    parts = [f"Bottom 固定为 {fixed_bottom:g}" if fixed_bottom is not None else "Bottom 拟合",
             f"Top 固定为 {fixed_top:g}" if fixed_top is not None else "Top 拟合"]
    return "，".join(parts)
