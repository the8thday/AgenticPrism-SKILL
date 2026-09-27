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
from .report import THEMES, data_uri, font_setup
from .workflow import dump, sha, verify_run


def _table(rows, columns):
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(key, '—') if row.get(key) is not None else '—'))}</td>"
                                      for key, _ in columns) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def render_simple(run, style=None):
    run = Path(run).resolve()
    manifest = verify_run(run)
    cfg = json.loads((run / "config.resolved.json").read_text())
    result = json.loads((run / "results.json").read_text())
    d = pd.read_csv(run / "normalized_data.csv", keep_default_na=False)
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
                with plt.rc_context({"font.family": [font] + ([cjk] if cjk else []), "font.size": 9}):
                    fig = plt.figure(figsize=(7.1, 4.5), layout="constrained")
                    positive = standards[standards.concentration_canonical > 0]
                    zero = standards[standards.concentration_canonical == 0]
                    if len(zero):
                        grid = fig.add_gridspec(1, 2, width_ratios=[1, 5])
                        ax0 = fig.add_subplot(grid[0, 0])
                        ax = fig.add_subplot(grid[0, 1], sharey=ax0)
                        ax0.scatter(np.zeros(len(zero)), zero.response, color=token["color"])
                        ax0.set_xlim(-.7, .7)
                        ax0.set_xticks([0], ["0"])
                        ax0.set_xlabel("Zero")
                        ax0.set_ylabel(f"Response ({g.response_unit.iloc[0]})")
                        ax.tick_params(labelleft=False)
                    else:
                        ax = fig.add_subplot(111)
                        ax.set_ylabel(f"Response ({g.response_unit.iloc[0]})")
                    ax.scatter(positive.concentration, positive.response, color=token["color"], label="Standards")
                    qc = fit.get("qc")
                    if qc and qc["accepted"]:
                        ax.axvspan(qc["lloq_input_unit"], qc["uloq_input_unit"], color=token["color"], alpha=.07, lw=0,
                                   label="Plate screening range")
                        for edge in (qc["lloq_input_unit"], qc["uloq_input_unit"]):
                            ax.axvline(edge, color="#808080", ls="--", lw=.8)
                    unknown = [w for w in result["wells"] if w["plate_id"] == plate and w["status"] == "quantified"
                               and w.get("measured_concentration_in_well") is not None]
                    if unknown:
                        ax.scatter([w["measured_concentration_in_well"] for w in unknown], [w["response"] for w in unknown],
                                   marker="D", s=26, facecolors="white", edgecolors="#b5651d", linewidths=1.2, zorder=4,
                                   label="Unknowns (in-well conc.)")
                    if fit["bottom"] is not None and len(positive):
                        x = np.geomspace(float(positive.concentration_canonical.min()), float(positive.concentration_canonical.max()), 200)
                        curve = {**fit, "log10_c_canonical": fit.get("log10_c_canonical", fit["log10_half_response_canonical"])}
                        ax.plot(x / factor, predict(x, curve, cfg["direction"]), color=token["color"],
                                label=fit.get("calibration_model", "4pl").upper())
                    ax.set_xscale("log")
                    ax.set_xlabel(f"Concentration in well ({unit}; log scale)")
                    ax.set_title(f"Plate {plate}")
                    ax.legend(frameon=False, fontsize=8)
                    if token["grid"]:
                        ax.grid(True, alpha=.2)
                    for ext in ("svg", "pdf", "png"):
                        fig.savefig(figures / f"{key}__{theme}.{ext}", dpi=300)
                    plt.close(fig)
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
    elif kind == "multi_group_comparison":
        key = "multigroup-001"
        figure_keys.append(key)
        omnibus, q = result["fits"][0], cfg["comparison"]
        used = d[d.exclude.astype(str).str.lower() != "true"]
        for theme, token in THEMES.items():
            with plt.rc_context({"font.family": [font] + ([cjk] if cjk else []), "font.size": 9}):
                fig, ax = plt.subplots(figsize=(max(4.8, 1.3 * len(q["groups"]) + 2), 4.5), layout="constrained")
                for i, summary in enumerate(result["group_summaries"]):
                    v = used[used.group == summary["group"]].value.astype(float).to_numpy()
                    jitter = np.linspace(-.12, .12, len(v)) if len(v) > 1 else np.zeros(1)
                    ax.scatter(i + jitter, v, color=token["color"], s=22, zorder=3)
                    ax.hlines(summary["mean"], i - .25, i + .25, color="#202020", lw=1.4, zorder=4)
                    ax.errorbar(i, summary["mean"], yerr=summary["sd"], color="#202020", capsize=5, lw=1, zorder=4)
                ax.set_xticks(range(len(q["groups"])), q["groups"], rotation=20 if len(q["groups"]) > 4 else 0)
                ax.set_ylabel(f"{q['outcome']} ({q['unit']})")
                ax.set_xlim(-.6, len(q["groups"]) - .4)
                if token["grid"]:
                    ax.grid(True, axis="y", alpha=.2)
                for ext in ("svg", "pdf", "png"):
                    fig.savefig(figures / f"{key}__{theme}.{ext}", dpi=300)
                plt.close(fig)
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
            with plt.rc_context({"font.family": [font] + ([cjk] if cjk else []), "font.size": 9}):
                fig, ax = plt.subplots(figsize=(6.2, 4.5), layout="constrained")
                if q["design"] == "paired_t":
                    for uid in a.index:
                        ax.plot([0, 1], [a[uid], b[uid]], color="#999999", lw=.8, alpha=.7)
                ax.scatter(np.zeros(len(a)), a, color=token["color"], label=q["group_a"], zorder=3)
                ax.scatter(np.ones(len(b)), b, color="#b5651d", label=q["group_b"], zorder=3)
                ax.set_xticks([0, 1], [q["group_a"], q["group_b"]])
                ax.set_ylabel(f"{q['outcome']} ({q['unit']})")
                ax.set_xlim(-.4, 1.4)
                if token["grid"]:
                    ax.grid(True, axis="y", alpha=.2)
                for ext in ("svg", "pdf", "png"):
                    fig.savefig(figures / f"{key}__{theme}.{ext}", dpi=300)
                plt.close(fig)
        title = "两组统计比较"
        details = "差值定义为 B−A；双侧 t 检验和区间。独立组使用 Welch 自由度，配对组按匹配的独立单位求差；输入必须已经是一单位一值。"
        rows = _table([comp], [("design", "设计"), ("group_a", "A"), ("group_b", "B"), ("n_a", "n A"),
                               ("n_b", "n B"), ("mean_a", "均值 A"), ("mean_b", "均值 B"),
                               ("mean_difference_b_minus_a", "均值差 B−A"), ("ci_difference", "差值区间"),
                               ("t_statistic", "t"), ("df", "自由度"), ("p_two_sided", "双侧 p")])
    blocks = []
    for key in figure_keys:
        for theme in THEMES:
            links = " ".join(f'<a download="{key}__{theme}.{ext}" href="{data_uri(figures / f"{key}__{theme}.{ext}")}">{ext.upper()}</a>'
                             for ext in ("svg", "pdf", "png"))
            blocks.append(f'<div class="figure" data-theme="{theme}" {"" if theme == style else "hidden"}>'
                          f'<img alt="{html.escape(key)}" src="{data_uri(figures / f"{key}__{theme}.svg")}"><p>{links}</p></div>')
    page = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>body{{font:16px/1.6 system-ui,sans-serif;max-width:1050px;margin:auto;padding:22px;color:#203040;background:#f6f8fa}}
section{{background:white;padding:20px;margin:18px 0;border-radius:10px;overflow:auto}}table{{border-collapse:collapse;font-size:13px}}td,th{{padding:6px 9px;border-bottom:1px solid #ddd;white-space:nowrap;text-align:left}}
img{{max-width:100%}}button{{margin-right:8px;padding:6px 12px}}.figure[hidden]{{display:none}}</style>
<h1>{html.escape(title)}</h1><p>{html.escape(details)}</p><p>来源：{html.escape(cfg['source'])} · 版本 {html.escape(manifest['package_version'])} · 不宣称 Prism 数值等价</p>
<section><button onclick="selectTheme('prism_like')">Prism-like</button><button onclick="selectTheme('standard')">Standard</button>{''.join(blocks)}</section>
<section>{rows}</section><section><h2>配置与验证</h2><p>结果哈希可用 agentic-prism verify 核验。图形重绘不改变计算结果。</p>
<pre>{html.escape(json.dumps(cfg, ensure_ascii=False, indent=2))}</pre></section>
<script>function selectTheme(s){{document.querySelectorAll('.figure').forEach(e=>e.hidden=e.dataset.theme!==s)}}</script></html>'''
    (run / "report.html").write_text(page)
    dump(run / "render_manifest.json", {"rendered_utc": datetime.now(timezone.utc).isoformat(), "style": style,
                                         "report_sha256": sha(run / "report.html"),
                                         "scientific_manifest_sha256": sha(run / "manifest.json")})
    return run / "report.html"
