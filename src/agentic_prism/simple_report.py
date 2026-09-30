"""Self-contained reports for the scoped ELISA and group-comparison workflows."""
from datetime import datetime, timezone
from pathlib import Path
import html
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .calibration import predict
from .dose_schema import DOSE_UNITS
from . import plot_style as pstyle
from .report import THEMES, font_setup
from .report_shell import data_uri, fmt, json_block, page, theme_figures
from .workflow import dump, sha, verify_run

EYEBROW = {"elisa_quantification": "ELISA QUANTIFICATION", "group_comparison": "GROUP COMPARISON",
           "multi_group_comparison": "MULTI-GROUP COMPARISON", "nonparametric": "RANK-BASED COMPARISON",
           "mmrm": "MIXED MODEL FOR REPEATED MEASURES", "repeated_measures": "REPEATED MEASURES",
           "time_to_event": "TIME TO EVENT", "tumor_growth": "TUMOR GROWTH"}
SIGNIFICANCE_NOTE = ("图中括号与星号只转写已保存的校正后 p 值（* p<0.05，** p<0.01，*** p<0.001，**** p<0.0001，ns ≥0.05），"
                     "渲染时不做任何新检验；星号不表示效应大小或生物学重要性，请以差值及其区间为准。")


def _table(rows, columns):
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(fmt(row.get(key)))}</td>"
                                      for key, _ in columns) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _save(fig, figures, key, theme):
    pstyle.save(fig, figures / f"{key}__{theme}", THEMES[theme], 300)


def _annotate(ax, positions, pairs, values):
    """Brackets for saved adjusted p values (render option); widens only the upper y limit."""
    finite = np.asarray([x for v in values for x in v], float)
    if not len(pairs) or not len(finite):
        return
    top, low = float(np.nanmax(finite)), float(np.nanmin(finite))
    upper = pstyle.significance_brackets(ax, positions, pairs, top, max(top - low, abs(top) * .1, 1e-12))
    ax.set_ylim(ax.get_ylim()[0], max(ax.get_ylim()[1], upper))


def render_simple(run, style=None, annotate_significance=False):
    run = Path(run).resolve()
    manifest = verify_run(run)
    cfg = json.loads((run / "config.resolved.json").read_text())
    result = json.loads((run / "results.json").read_text())
    d = pd.read_csv(run / "normalized_data.csv", keep_default_na=False,
                    dtype={"independent_unit_id": str, "group": str, "observation_id": str})
    style = style or cfg["report"]["plot_style"]
    if style not in THEMES:
        raise ValueError("Unsupported plot style")
    figures = run / "figures"
    figures.mkdir(exist_ok=True)
    font, cjk = font_setup()
    kind = cfg["analysis_type"]
    figure_keys = []
    if kind == "elisa_quantification":
        for number, fit in enumerate(result["fits"], 1):
            plate = fit["plate_id"]
            g = d[d.plate_id == plate]
            standards = g[(g.role == "standard") & (g.exclude.astype(str).str.lower() != "true")]
            key = f"elisa-{number:03d}"
            figure_keys.append(key)
            unit = str(g.concentration_unit.iloc[0])
            factor = DOSE_UNITS[unit]
            for theme, token in THEMES.items():
                with plt.rc_context(pstyle.rc(token, font, cjk)):
                    fig = plt.figure(figsize=(6.2, 4.2), layout="constrained")
                    positive = standards[standards.concentration_canonical > 0]
                    zero = standards[standards.concentration_canonical == 0]
                    if len(zero):
                        grid = fig.add_gridspec(1, 2, width_ratios=[1, 5])
                        ax0 = fig.add_subplot(grid[0, 0])
                        ax = fig.add_subplot(grid[0, 1], sharey=ax0)
                        ax0.plot(np.zeros(len(zero)), zero.response, **pstyle.point_style(token, 0))
                        ax0.set_xlim(-.7, .7)
                        ax0.set_xticks([0], ["0"])
                        ax0.set_xlabel("Zero")
                        ax0.set_ylabel(f"Response ({g.response_unit.iloc[0]})")
                        ax.spines["left"].set_visible(False)
                        ax.tick_params(axis="y", left=False, labelleft=False)
                    else:
                        ax = fig.add_subplot(111)
                        ax.set_ylabel(f"Response ({g.response_unit.iloc[0]})")
                    ax.plot(positive.concentration, positive.response, label="Standards", **pstyle.point_style(token, 0))
                    qc = fit.get("qc")
                    if qc and qc["accepted"]:
                        ax.axvspan(qc["lloq_input_unit"], qc["uloq_input_unit"], color=pstyle.color(token, 1), alpha=.07, lw=0,
                                   label="Plate screening range")
                        for edge in (qc["lloq_input_unit"], qc["uloq_input_unit"]):
                            ax.axvline(edge, color=pstyle.MUTED, ls=(0, (3, 2)), lw=.8)
                    unknown = [w for w in result["wells"] if w["plate_id"] == plate and w["status"] == "quantified"
                               and w.get("measured_concentration_in_well") is not None]
                    if unknown:
                        ax.scatter([w["measured_concentration_in_well"] for w in unknown], [w["response"] for w in unknown],
                                   marker="D", s=26, facecolors="white", edgecolors=pstyle.color(token, 2), linewidths=1.2, zorder=4,
                                   label="Unknowns (in-well conc.)")
                    if fit["bottom"] is not None and len(positive):
                        x = np.geomspace(float(positive.concentration_canonical.min()), float(positive.concentration_canonical.max()), 200)
                        curve = {**fit, "log10_c_canonical": fit.get("log10_c_canonical", fit["log10_half_response_canonical"])}
                        ax.plot(x / factor, predict(x, curve, cfg["direction"]), color=token["color"],
                                label=fit.get("calibration_model", "4pl").upper())
                    ax.set_xscale("log")
                    ax.set_xlabel(f"Concentration in well ({unit})")
                    ax.set_title(f"Plate {plate}")
                    ax.legend(fontsize=8)
                    _save(fig, figures, key, theme)
        title = "ELISA 标准曲线与未知样本反算"
        details = (f"逐板原始响应 {cfg['fit'].get('model', '4pl').upper()}（配置中预先声明）。标准品回算回收率/CV及其配置阈值见下表。通过的最低、最高水平只界定本板标准品通过范围，"
                   "兼容字段 LLOQ/ULOQ 不代表经实验验证的检测方法定量限。独立 QC 不参与拟合；提供的 QC 未通过时未知孔不报告。"
                   "未提供独立 QC 的历史配置可给出探索性反算，但不能据此声称方法已验证。未知孔仅在本板允许范围内反算后乘稀释倍数。"
                   "不自动剔除标准品。图中菱形为孔内浓度。未知浓度区间为 delta 法近似，只反映本板标准曲线参数和孔间响应噪声（沿用标准品残差模型），"
                   "不包含板间、基质或稀释操作误差。稀释线性比较同一样本在本板各个可定量稀释度的校正浓度；未通过时可能提示 hook 效应或基质干扰，但不能据此判定原因。"
                   "跨板 QC 仅为描述性统计，不是完整的中间精密度验证。")
        rows = _table(result["wells"], [("plate_id", "板"), ("well_id", "孔"), ("sample_id", "样本"),
                                      ("response", "响应"), ("dilution_factor", "稀释倍数"),
                                      ("estimated_sample_concentration", "反算浓度"), ("ci_low", "区间下限"), ("ci_high", "区间上限"),
                                      ("concentration_unit", "单位"), ("status", "状态"), ("diagnostics", "诊断")])
        rows += _table(result["summaries"], [("plate_id", "板"), ("sample_id", "样本"),
                                           ("n_quantified", "可定量孔数"), ("n_wells_in_estimate", "计入孔数"), ("mean_concentration", "均值"),
                                           ("ci_low", "区间下限"), ("ci_high", "区间上限"),
                                           ("sd_concentration", "SD"), ("cv_percent", "CV%"), ("unit", "单位"), ("status", "状态"),
                                           ("linearity_status", "稀释线性"), ("linearity_diagnostics", "线性诊断")])
        linearity_rows = [{"plate_id": f["plate_id"], **row} for f in result["fits"] for row in f.get("dilution_linearity", [])]
        if any(s.get("n_dilutions", 1) > 1 for s in result["summaries"]):
            rows += "<h2>稀释线性（同板、同样本）</h2>" + _table(linearity_rows, [
                ("plate_id", "板"), ("sample_id", "样本"), ("dilution_factor", "稀释倍数"), ("n_wells", "孔数"),
                ("in_range", "全部孔可定量"), ("mean_corrected_concentration", "校正浓度均值"),
                ("recovery_vs_dilution_mean_percent", "相对各稀释度均值%")])
        comparisons = [{"plate_id": f["plate_id"], **f["model_comparison"]} for f in result["fits"] if f.get("model_comparison")]
        if comparisons:
            rows += "<h2>5PL 与 4PL 比较（仅供参考，不改变预设模型）</h2>" + _table(comparisons, [
                ("plate_id", "板"), ("f_statistic", "F"), ("p_value", "p"), ("aicc_4pl", "AICc 4PL"), ("aicc_5pl", "AICc 5PL")])
        qc_rows = [{"plate_id": f["plate_id"], **level, "limits": "–".join(f"{x:g}" for x in level["limits_percent"]) + "%",
                    "recovery": None if level["recovery_percent"] is None else f"{level['recovery_percent']:.1f}%",
                    "cv": None if level["cv_percent"] is None else f"{level['cv_percent']:.1f}%",
                    "verdict": "通过" if level["passed"] else "未通过"}
                   for f in result["fits"] if f.get("qc") for level in f["qc"]["levels"]]
        if qc_rows:
            rows += "<h2>标准品回算质控</h2>" + "".join(
                f"<p>{html.escape(f['plate_id'])}：{f['qc']['n_passing']}/{f['qc']['n_levels']} 个水平通过（需要 ≥{f['qc']['n_required']}）；"
                f"本板标准品通过范围 {f['qc']['lloq_input_unit'] if f['qc']['accepted'] else '—'} – {f['qc']['uloq_input_unit'] if f['qc']['accepted'] else '—'}；"
                f"{'标准曲线接受' if f['qc']['accepted'] else '标准曲线未通过质控，未知样本不报告'}</p>" for f in result["fits"] if f.get("qc"))
            rows += _table(qc_rows, [("plate_id", "板"), ("nominal_input_unit", "标称浓度"), ("n_wells", "孔数"),
                                     ("mean_back_calculated_input_unit", "回算均值"), ("recovery", "回收率"), ("cv", "CV"),
                                     ("limits", "接受范围"), ("role", "位置"), ("verdict", "判定")])
        rows += "<h2>独立质控（不参与拟合）</h2>" + _table(
            [{"plate_id":f["plate_id"],"status":f.get("independent_qc",{}).get("status","not_assessed"),**row}
             for f in result["fits"] for row in (f.get("independent_qc",{}).get("levels",[]) or [{}])],
            [("plate_id","板"),("status","QC状态"),("nominal","标称浓度"),("unit","单位"),
             ("n","复孔数"),("recovery_percent","回收率%"),("cv_percent","CV%"),("passed","通过")])
        if result.get("cross_plate_qc"):
            rows += "<h2>跨板质控描述（非完整方法验证）</h2>" + _table(result["cross_plate_qc"],
                [("nominal","标称浓度"),("unit","单位"),("n_plates","板数"),("status","证据状态"),
                 ("mean_recovery_percent","平均回收率%"),("between_plate_mean_cv_percent","板均值间CV%")])
        rows += _table(result["fits"], [("plate_id", "板"), ("reportable", "标准曲线可报告"),
                                      ("calibration_model", "模型"), ("bottom", "Bottom"), ("top", "Top"), ("half_response_input_unit", "半响应浓度"),
                                      ("asymmetry", "不对称参数 g"),
                                      ("rmse", "RMSE"), ("diagnostics", "诊断")])
    elif kind == "mmrm":
        title = "Marginal repeated-measures model"
        q = cfg["comparison"]
        details = (f"REML {q['covariance']}; {q['inference']}. Arm by categorical visit cell means. "
                   "Primary test: arm by visit interaction. Contrasts: arm minus control at each visit, Holm p and Bonferroni family intervals. "
                   "AR(1) refers to ordered visit lag. Missing observations require MAR conditional on the model. "
                   "See saved covariance, diagnostics and interpretation facts.")
        rows = _table(result["fits"], [(k,k) for k in ("status","n_units","n_total","f_statistic","df_numerator","df_denominator","p_value","diagnostics")])
        rows += _table(result["contrasts"], [(k,k) for k in ("arm","control_arm","condition","estimate","standard_error","df","ci_low","ci_high","p_adjusted")])
        key = "mmrm-observed-001"
        figure_keys.append(key)
        used = d[d.exclude.astype(str).str.lower() != "true"]
        for theme, token in THEMES.items():
            with plt.rc_context(pstyle.rc(token, font, cjk)):
                fig, ax = plt.subplots(figsize=(5.6, 3.9), layout="constrained")
                for i, arm in enumerate(q["arms"]):
                    summaries = used[used.arm == arm].groupby("group").value.agg(["mean", "std"]).reindex(q["groups"])
                    marks = pstyle.point_style(token, i)
                    marks.pop("linestyle")
                    ax.errorbar(np.arange(len(q["groups"])) + .05 * (i - (len(q["arms"]) - 1) / 2), summaries["mean"],
                                yerr=summaries["std"], color=pstyle.color(token, i), capsize=3, label=arm, **marks)
                ax.set_xticks(range(len(q["groups"])), q["groups"])
                ax.set_xlim(-.5, len(q["groups"]) - .5)
                ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
                ax.set_title("Observed means ± SD; available observations")
                ax.legend()
                _save(fig, figures, key, theme)
    elif kind == "nonparametric":
        title = "Nonparametric group comparison"
        q, fit = cfg["comparison"], result["fits"][0]
        details = (f"{q['design']}; unit-level data. Two-sided tests. Shift is B minus A. "
                   "HL concerns a common location shift or symmetric paired differences, not a general difference in medians. "
                   "Dunn tests compare pooled mean ranks with the declared multiplicity adjustment; no shift intervals are supplied for Dunn. "
                   "See interpretation facts for unavailable intervals, actual confidence level and approximation warnings.")
        rows = _table(result["fits"], [(k, k) for k in ("design", "n_units", "statistic", "df", "p_value", "hodges_lehmann", "estimate", "ci_low", "ci_high", "confidence_level", "achieved_confidence_level", "confidence_level_basis", "inference", "status", "diagnostics")])
        rows += _table(result["group_summaries"], [(k, k) for k in ("group", "n", "median")])
        rows += _table(result["contrasts"], [(k, k) for k in ("group_a", "group_b", "mean_rank_difference_b_minus_a", "statistic", "p_adjusted", "adjustment")])
        key = "ranks-001"
        figure_keys.append(key)
        used = d[d.exclude.astype(str).str.lower() != "true"]
        for theme, token in THEMES.items():
            with plt.rc_context(pstyle.rc(token, font, cjk)):
                fig, ax = plt.subplots(figsize=(max(3.2, .85 * len(q["groups"]) + 1.6), 3.9), layout="constrained")
                groups = [(g, used[used.group == g].value.astype(float).to_numpy()) for g in q["groups"]]
                positions = pstyle.dot_plot(ax, token, groups, center="median", spread=None)
                ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
                ax.set_title("Lines: medians")
                if annotate_significance:
                    index = {g: i for i, g in enumerate(q["groups"])}
                    pairs = [(index[c["group_a"]], index[c["group_b"]], c.get("p_adjusted")) for c in result["contrasts"]
                             if c.get("group_a") in index and c.get("group_b") in index]
                    if not pairs and len(q["groups"]) == 2:
                        pairs = [(0, 1, fit.get("p_value"))]
                    _annotate(ax, positions, pairs, [v for _, v in groups])
                _save(fig, figures, key, theme)
    elif kind == "repeated_measures":
        from .repeated_report import repeated_sections
        title, details, rows, figure_keys = repeated_sections(result, cfg, d, figures, font, cjk)
    elif kind == "time_to_event":
        from .survival_report import survival_sections
        title, details, rows, figure_keys = survival_sections(result, cfg, d, figures, font, cjk)
    elif kind == "tumor_growth":
        from .tumor_report import tumor_sections
        title, details, rows, figure_keys = tumor_sections(result, cfg, d, figures, font, cjk)
    elif kind == "multi_group_comparison":
        key = "multigroup-001"
        figure_keys.append(key)
        omnibus, q = result["fits"][0], cfg["comparison"]
        used = d[d.exclude.astype(str).str.lower() != "true"]
        for theme, token in THEMES.items():
            with plt.rc_context(pstyle.rc(token, font, cjk)):
                fig, ax = plt.subplots(figsize=(max(3.2, .85 * len(q["groups"]) + 1.6), 3.9), layout="constrained")
                groups = [(s["group"], used[used.group == s["group"]].value.astype(float).to_numpy())
                          for s in result["group_summaries"]]
                positions = pstyle.dot_plot(ax, token, groups, center="mean", spread="sd")
                ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
                if annotate_significance:
                    index = {g: i for i, (g, _) in enumerate(groups)}
                    pairs = [(index[c["group_a"]], index[c["group_b"]], c.get("p_adjusted")) for c in result["contrasts"]
                             if c.get("group_a") in index and c.get("group_b") in index]
                    _annotate(ax, positions, pairs, [v for _, v in groups])
                _save(fig, figures, key, theme)
        title = "多组统计比较"
        design = {"one_way_anova": "经典单因素 ANOVA（方差齐性）", "welch_anova": "Welch 单因素 ANOVA（不假定方差齐性）"}[omnibus["design"]]
        details = (f"{design}；多重比较族：{omnibus['contrast_method'] or '未预设（仅总体检验）'}。"
                   "差值为 B−A（对照在前或按声明顺序后组减前组）；区间为族内同时置信区间，p 值已按族校正。"
                   "两两比较按预设族直接报告，不以总体检验显著为前提。图中为每个独立单位的数值、均值及 ±SD。"
                   "输入必须已经是一单位一值；未检验正态性，也不代表生物学重要性或等效性。")
        if omnibus["diagnostics"]:
            details += " 诊断：" + "；".join(omnibus["diagnostics"])
        rows = _table([omnibus], [("design", "设计"), ("n_groups", "组数"), ("n_total", "总 n"),
                                  ("f_statistic", "F"), ("df_numerator", "df1"), ("df_denominator", "df2"),
                                  ("p_value", "总体 p"), ("max_min_sd_ratio", "SD 最大/最小")])
        rows += _table(result["group_summaries"], [("group", "组"), ("n", "n"), ("mean", "均值"), ("sd", "SD"), ("sem", "SEM")])
        if result["contrasts"]:
            rows += f"<h2>预设比较族（{html.escape(omnibus['contrast_method'])}）</h2>" + _table(
                result["contrasts"], [("group_a", "A"), ("group_b", "B"), ("difference_b_minus_a", "均值差 B−A"),
                                      ("ci_low", "同时区间下限"), ("ci_high", "同时区间上限"), ("p_adjusted", "校正 p")])
    else:
        key = "groups-001"
        figure_keys.append(key)
        comp = result["fits"][0]
        q = cfg["comparison"]
        used = d[d.exclude.astype(str).str.lower() != "true"]
        a = used[used.group == q["group_a"]].set_index("independent_unit_id").value.astype(float)
        b = used[used.group == q["group_b"]].set_index("independent_unit_id").value.astype(float)
        for theme, token in THEMES.items():
            with plt.rc_context(pstyle.rc(token, font, cjk)):
                fig, ax = plt.subplots(figsize=(3.0, 3.9), layout="constrained")
                if q["design"] == "paired_t":
                    # Before-after plot: each line is one independent unit.
                    for uid in a.index:
                        ax.plot([0, 1], [a[uid], b[uid]], color=pstyle.MUTED, lw=.8, alpha=.8, zorder=2)
                    for i, v in enumerate((a, b)):
                        ax.plot(np.full(len(v), float(i)), v, zorder=3, **pstyle.point_style(token, i, 5.5))
                    ax.set_xticks([0, 1], [q["group_a"], q["group_b"]])
                    ax.set_xlim(-.4, 1.4)
                    positions = np.array([0., 1.])
                else:
                    positions = pstyle.dot_plot(ax, token, [(q["group_a"], a.to_numpy()), (q["group_b"], b.to_numpy())])
                ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
                if annotate_significance:
                    _annotate(ax, positions, [(0, 1, comp.get("p_two_sided"))], [a.to_numpy(), b.to_numpy()])
                _save(fig, figures, key, theme)
        title = "两组统计比较"
        details = ("差值定义为 B−A；双侧 t 检验和区间。独立组使用 Welch 自由度，配对组按匹配的独立单位求差；输入必须已经是一单位一值。"
                   + ("图中连线连接同一独立单位的两次测量。" if q["design"] == "paired_t" else "图中每点为一个独立单位，横线与误差线为均值 ± SD（描述性）。"))
        rows = _table([comp], [("design", "设计"), ("group_a", "A"), ("group_b", "B"), ("n_a", "n A"),
                               ("n_b", "n B"), ("mean_a", "均值 A"), ("mean_b", "均值 B"),
                               ("mean_difference_b_minus_a", "均值差 B−A"), ("ci_difference", "差值区间"),
                               ("t_statistic", "t"), ("df", "自由度"), ("p_two_sided", "双侧 p")])
    annotated = annotate_significance and kind in ("group_comparison", "multi_group_comparison", "nonparametric")
    blocks = "".join(theme_figures(figures, key, style, f"{title} · {key}", THEMES) for key in figure_keys)
    if annotated:
        blocks += f'<p class="caption">{html.escape(SIGNIFICANCE_NOTE)}</p>'
    facts_section = ""
    facts_path = run / "interpretation_facts.json"
    if "interpretation_facts.json" in manifest["scientific_artifacts_sha256"]:
        facts_section = ('<section id="facts" class="card"><h2>解读依据</h2><p>可报告范围、区间方法、诊断、限制与结果来源已单独保存。'
                         '解读时须保留受限或不可报告的状态。</p><div class="downloads">'
                         f'<a download="interpretation_facts.json" href="{data_uri(facts_path)}">下载解读依据 JSON</a></div></section>')
    body = (f'<section id="figures" class="card"><h2>图形</h2>{blocks}</section>'
            f'<section id="methods" class="card"><h2>方法与说明</h2><p>{html.escape(details)}</p></section>'
            f'<section id="results" class="card results"><h2>结果</h2>{rows}</section>{facts_section}'
            '<section id="config" class="card"><h2>配置与验证</h2><p>结果哈希可用 agentic-prism verify 核验。图形重绘不改变计算结果。</p>'
            f'{json_block(cfg, "完整配置")}</section>')
    nav = [("figures", "图形"), ("methods", "方法与说明"), ("results", "结果")] + ([("facts", "解读依据")] if facts_section else []) + [("config", "配置与验证")]
    page_html = page(title, eyebrow="AgenticPrism / " + EYEBROW.get(kind, kind), heading=title,
                     lede=f"来源：{cfg['source']} · 版本 {manifest['package_version']} · 不宣称 Prism 数值等价",
                     nav=nav, body=body, style=style, themes=list(THEMES),
                     footer=f"AgenticPrism {html.escape(manifest['package_version'])} · 数值分析与图形渲染分离")
    (run / "report.html").write_text(page_html)
    dump(run / "render_manifest.json", {"rendered_utc": datetime.now(timezone.utc).isoformat(), "style": style,
                                         "significance_annotations": bool(annotated),
                                         "report_sha256": sha(run / "report.html"),
                                         "scientific_manifest_sha256": sha(run / "manifest.json")})
    return run / "report.html"
