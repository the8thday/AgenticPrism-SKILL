"""Render saved time-to-event artifacts (Kaplan-Meier with number-at-risk table); never refit."""
import numpy as np
import matplotlib.pyplot as plt
from . import plot_style as pstyle
from .report import THEMES


def _step(curve, key):
    t = np.concatenate([[0.], [c["time"] for c in curve]])
    v = np.concatenate([[1.], [np.nan if c[key] is None else c[key] for c in curve]])
    return t, v


def survival_sections(result, cfg, d, figures, font, cjk):
    from .simple_report import _table
    q, study, fit = cfg["comparison"], cfg["study"], result["fits"][0]
    arms = q["arms"]
    horizon = max(c["time"] for c in result["km_curves"])
    ticks = np.linspace(0, horizon, 6)
    for theme, token in THEMES.items():
        with plt.rc_context(pstyle.rc(token, font, cjk)):
            fig = plt.figure(figsize=(6.2, 4.9), layout="constrained")
            grid = fig.add_gridspec(2, 1, height_ratios=[4, 1.1 + .18 * len(arms)])
            ax, table_ax = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
            for i, arm in enumerate(arms):
                curve = [c for c in result["km_curves"] if c["arm"] == arm]
                color = pstyle.color(token, i)
                t, s = _step(curve, "survival")
                ax.step(t, s, where="post", color=color, lw=token["line_width"], label=arm)
                if cfg["report"]["show_confidence_bands"]:
                    _, lo = _step(curve, "lower")
                    _, hi = _step(curve, "upper")
                    ax.fill_between(t, lo, hi, step="post", color=color, alpha=.12, lw=0)
                censored = [c for c in curve if c["n_censor"] > 0]
                ax.plot([c["time"] for c in censored], [c["survival"] for c in censored], "+", color=color, ms=8, mew=1.3)
            ax.set_ylim(-.02, 1.04)
            ax.set_xlim(0, horizon * 1.02)
            ax.set_xticks(ticks)
            ax.set_ylabel("Survival probability")
            ax.set_xlabel(f"Time since {study['time_origin']} ({study['time_unit']})")
            ax.legend(loc="lower left")
            table_ax.axis("off")
            table_ax.set_xlim(ax.get_xlim())
            table_ax.set_ylim(-.5, len(arms) + .3)
            table_ax.text(0, len(arms), "Number at risk", fontsize=8, va="center", fontweight="bold" if token["bold_labels"] else "normal")
            for i, arm in enumerate(arms):
                sub = d[(d.arm == arm) & (d.exclude.astype(str).str.lower() != "true")]
                times = sub.time.astype(float).to_numpy()
                table_ax.text(-.02, len(arms) - 1 - i, arm, transform=table_ax.get_yaxis_transform(), ha="right", va="center",
                              fontsize=8, color=pstyle.color(token, i))
                for tick in ticks:
                    table_ax.text(tick, len(arms) - 1 - i, str(int((times >= tick).sum())), ha="center", va="center",
                                  fontsize=8, color=pstyle.color(token, i))
            pstyle.save(fig, figures / f"survival-001__{theme}", token, 300)
    details = (f"终点：{study['endpoint']}；时间零点：{study['time_origin']}；删失：{study['censoring_rationale']}。"
               f"Kaplan–Meier 估计，{q['confidence_level']:.0%} 逐点置信限按 {q['conf_type']} 变换（Greenwood 方差）；"
               "中位生存时间及其区间取曲线（及其置信限）首次降至 0.5 的时间，未达到时标为 not reached。"
               f"组间比较用 log-rank 检验（{fit['logrank_inference']}）；预设两两比较为 log-rank，Holm 校正。"
               "Cox 比例风险模型（Efron 处理并列时间）给出风险比；比例风险假设用 cox.zph 同法的 score 检验评估，"
               "不成立时风险比只是随时间变化的效应的平均，需结合曲线解读。删失被假定为与预后无关（非信息性删失），"
               "这一点无法从数据中验证。与 R survival 的逐项对照见验证记录。")
    details += " 状态：" + fit["status"] + "；诊断：" + ("；".join(fit["diagnostics"]) or "无")
    landmark_cols = [(k, k.replace("survival_at_", "生存率@")) for k in result["arm_summaries"][0] if k.startswith("survival_at_")
                     and not k.endswith(("_lower", "_upper"))]
    rows = "<h2>各组摘要</h2>" + _table(result["arm_summaries"], [("arm", "组"), ("n", "n"), ("events", "事件"), ("censored", "删失"),
        ("median", "中位时间"), ("median_lower", "中位下限"), ("median_upper", "中位上限"), ("median_status", "中位状态")] + landmark_cols)
    lr = result["logrank"]
    rows += "<h2>Log-rank 检验（全部组）</h2>" + _table([lr], [("chisq", "χ²"), ("df", "df"), ("p_value", "p"), ("inference", "推断")])
    if result["contrasts"]:
        rows += "<h2>预设两两比较</h2>" + _table(result["contrasts"], [("arm_a", "A"), ("arm_b", "B"), ("logrank_chisq", "log-rank χ²"),
            ("p_adjusted", "Holm 校正 p"), ("hazard_ratio_b_vs_a", "HR (B/A)"), ("hr_lower", "同时区间下限"), ("hr_upper", "同时区间上限")])
    if result["cox_terms"]:
        rows += "<h2>Cox 模型</h2>" + _table(result["cox_terms"], [("term", "项"), ("hazard_ratio", "HR"), ("hr_lower", "下限"),
            ("hr_upper", "上限"), ("p_value", "Wald p"), ("adjusted_for", "调整协变量")])
        rows += "<h2>比例风险检验（cox.zph，KM 时间变换）</h2>" + _table(result["ph_test"], [("term", "项"), ("chisq", "χ²"), ("df", "df"), ("p_value", "p")])
    return "生存分析（Kaplan–Meier、log-rank 与 Cox）", details, rows, ["survival-001"]
