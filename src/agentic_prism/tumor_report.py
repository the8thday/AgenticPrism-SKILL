"""Render saved tumor-growth artifacts; never refit."""
import numpy as np
import matplotlib.pyplot as plt
from .report import THEMES
from .survival_report import COLORS


def tumor_sections(result, cfg, d, figures, font, cjk):
    from .simple_report import _table
    a, s, fit = cfg["analysis"], cfg["study"], result["fits"][0]
    arms = a["arms"]
    used = d[d.exclude.astype(str).str.lower() != "true"]
    beta = (result.get("fitted_coefficients") or {}).get("beta")
    for theme, token in THEMES.items():
        with plt.rc_context({"font.family": [font] + ([cjk] if cjk else []), "font.size": 9}):
            fig, ax = plt.subplots(figsize=(7.4, 4.9), layout="constrained")
            for i, arm in enumerate(arms):
                color = COLORS[i % len(COLORS)]
                sub = used[used.arm == arm]
                for _, animal in sub.groupby("animal_id"):
                    ax.plot(animal.day.astype(float), animal.volume.astype(float), color=color, alpha=.18, lw=.7)
                summ = [r for r in result["arm_day_summaries"] if r["arm"] == arm]
                ax.errorbar([r["day"] for r in summ], [r["mean"] for r in summ], yerr=[r["sem"] or 0 for r in summ],
                            fmt="o", color=color, ms=4, capsize=3, lw=1.2, label=f"{arm} (mean ± SEM)")
                if beta:
                    grid = np.linspace(a["baseline_day"], max(r["day"] for r in summ), 100)
                    k = len(arms)
                    ax.plot(grid, np.exp(beta[i] + beta[k + i] * grid) - a["log_offset"], color=color, lw=1.8, ls="--")
            ax.set_yscale("log")
            ax.axvline(a["analysis_day"], color="#888", ls=":", lw=.9)
            ax.set_xlabel(f"Time since {s['time_origin']} ({s['time_unit']})")
            ax.set_ylabel(f"Tumor volume ({s['volume_unit']}; log scale)")
            ax.set_title("Dashed: model-estimated geometric mean growth; dotted line: analysis day", fontsize=9)
            ax.legend(frameon=False, fontsize=8)
            if token["grid"]:
                ax.grid(True, which="both", alpha=.15)
            for ext in ("svg", "pdf", "png"):
                fig.savefig(figures / f"tumor-001__{theme}.{ext}", dpi=220)
            plt.close(fig)
    details = (f"移除规则：{s['removal_rule']}。生长模型：log(V + {a['log_offset']}) = 各组截距 + 各组斜率 × 时间 + 动物随机截距与随机斜率"
               "（非结构化协方差）+ 残差，REML 拟合；Satterthwaite 小样本 t/F（与 R lmerTest 同法）。"
               "生长速率为对数体积每单位时间的斜率，倍增时间 = ln2/斜率。模型 T/C 为分析日模型估计的几何均值之比。"
               f"动物因人道终点离开数据后，模型在随机缺失（MAR）假设下仍有效：{s['dropout_rationale']}。"
               "观测 TGI% = 100 ×（1 − 处理组平均体积增量 / 对照组平均体积增量），T/C% 为分析日平均体积之比，"
               "只用分析日仍被测量的动物，区间为 Fieller 法（按比较族做 Bonferroni 校正）；若有动物在分析日前被移除，这两个指标会有偏倚，"
               "应以模型结果和生存分析为主。")
    details += " 诊断：" + ("；".join(fit["diagnostics"]) or "无")
    rows = ""
    if result["growth_tests"]:
        rows += "<h2>各组生长速率相等的检验</h2>" + _table(result["growth_tests"], [("f_statistic", "F"), ("df_numerator", "df1"), ("df_denominator", "df2"), ("p_value", "p")])
    if result["growth_rates"]:
        rows += "<h2>各组生长速率</h2>" + _table(result["growth_rates"], [("arm", "组"), ("log_growth_rate_per_time", "对数生长速率"), ("rate_lower", "下限"),
            ("rate_upper", "上限"), ("doubling_time", "倍增时间"), ("doubling_time_lower", "倍增下限"), ("doubling_time_upper", "倍增上限")])
    if result["model_contrasts"]:
        rows += "<h2>与对照比较（模型）</h2>" + _table(result["model_contrasts"], [("arm", "组"), ("log_rate_difference", "速率差"),
            ("rate_difference_lower", "同时下限"), ("rate_difference_upper", "同时上限"), ("p_adjusted", "Holm p"),
            ("model_tc_geometric_ratio_at_analysis_day", "模型 T/C"), ("model_tc_lower", "T/C 下限"), ("model_tc_upper", "T/C 上限")])
    rows += f"<h2>分析日（{a['analysis_day']:g}）观测 TGI% 与 T/C%</h2>" + _table(result["observed_tgi"], [("arm", "组"), ("n_arm_on_day", "n"),
        ("n_control_on_day", "对照 n"), ("tgi_percent", "TGI%"), ("tgi_lower", "下限"), ("tgi_upper", "上限"), ("tc_percent", "T/C%"),
        ("tc_lower", "下限"), ("tc_upper", "上限")])
    rows += "<h2>各组各时间点观测摘要</h2>" + _table(result["arm_day_summaries"], [("arm", "组"), ("day", "时间"), ("n", "n"), ("mean", "均值"), ("sem", "SEM"), ("median", "中位数")])
    if result["dropout"]:
        rows += "<h2>分析日前离开的动物</h2>" + _table(result["dropout"], [("animal_id", "动物"), ("arm", "组"), ("last_day", "最后测量")])
    return "肿瘤生长曲线（随机斜率混合模型、TGI 与 T/C）", details, rows, ["tumor-001"]
