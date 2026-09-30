"""Render saved repeated-measures artifacts; never refit."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import probplot
from . import plot_style as pstyle
from .report import THEMES


def two_way_sections(result, cfg, d, figures, font, cjk):
    from .simple_report import _table
    q, fit = cfg["comparison"], result["fits"][0]
    used = d[d.exclude.astype(str).str.lower() != "true"].copy()
    used["group"], used["arm"] = used.group.astype(str), used.arm.astype(str)
    pos = {g: i for i, g in enumerate(q["groups"])}
    estimates = {(r["arm"], r["group"]): r for r in result["fixed_effects"]}
    for theme, token in THEMES.items():
        with plt.rc_context(pstyle.rc(token, font, cjk)):
            fig, ax = plt.subplots(figsize=(max(5.6, .8 * len(q["groups"]) + 2.6), 4.2), layout="constrained")
            for i, arm in enumerate(q["arms"]):
                color = pstyle.color(token, i)
                sub = used[used.arm == arm]
                for _, unit in sub.groupby("independent_unit_id"):
                    unit = unit.sort_values("group", key=lambda s: s.map(pos))
                    ax.plot(unit.group.map(pos), unit.value.astype(float), color=color, alpha=.18, lw=.7)
                means = [sub.value[sub.group == g].astype(float).mean() for g in q["groups"]]
                marks = pstyle.point_style(token, i, 5)
                marks.pop("linestyle")
                ax.plot(range(len(q["groups"])), means, "-", color=color, lw=token["line_width"], label=f"{arm} (observed mean)", **marks)
                if estimates:
                    est = np.array([estimates[(arm, g)]["estimated_mean"] for g in q["groups"]])
                    se = np.array([estimates[(arm, g)]["standard_error"] for g in q["groups"]])
                    ax.errorbar(np.arange(len(q["groups"])) + .06 * (i + 1), est, yerr=se, fmt="s", color=color, ms=3.5,
                                capsize=3, lw=1, alpha=.9)
            ax.set_xticks(range(len(q["groups"])), q["groups"], rotation=30 if len(q["groups"]) > 6 else 0)
            ax.set_xlabel("Within-unit condition (declared order)")
            ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
            ax.set_title(f"{fit['n_units']} units in {len(q['arms'])} arms; {fit['n_missing_cells']} missing cells"
                         + ("; squares: model-estimated mean ± SE" if estimates else ""))
            ax.legend(fontsize=8)
            pstyle.save(fig, figures / f"repeated-001__{theme}", token, 300)
    if q["design"] == "two_way_rm_anova":
        title = "两因素混合设计 ANOVA（组间 × 组内，完整数据）"
        details = ("组间因素（arm）用各单位跨条件均值的单因素 ANOVA 检验；条件及其交互作用在单位内正交对比上按 III 型（各组等权）"
                   "检验，自由度始终用合并组内协方差的 Greenhouse–Geisser ε 校正。预设比较：同一条件下组间用 Welch t，"
                   "同一组内条件间用配对 t；Holm 校正 p，Bonferroni 同时区间。")
    else:
        title = "两因素线性混合模型（组间 × 组内，单位随机截距）"
        details = ("REML 拟合：每个组 × 条件单元格一个固定均值，加单位随机截距与独立同方差残差（复合对称）。"
                   "三项检验为等权单元格均值的 III 型假设。"
                   + ("Satterthwaite 小样本 F/t（与 lmerTest 相同方法）。" if q["inference"] == "satterthwaite" else "Wald χ²/z 大样本近似，无小样本校正。")
                   + "缺失值不填补，依赖给定模型下的 MAR 假设。随机截距不能代替随机斜率或时间自相关。")
    details += " 交互作用是首要检验：交互显著时，组间差异随条件变化，主效应需结合各条件下的比较解读。"
    details += " 状态：" + fit["status"] + "；诊断：" + "；".join(fit["diagnostics"])
    effect_names = {"arm": "组间（arm）", "condition": "组内条件", "arm_x_condition": "交互作用"}
    tests = [{**t, "effect_label": effect_names[t["effect"]]} for t in result["tests"]]
    rows = "<h2>III 型检验</h2>" + _table(tests, [("effect_label", "效应"), ("f_statistic", "F"), ("wald_statistic", "Wald χ²"),
        ("df_numerator", "df1"), ("df_denominator", "df2"), ("epsilon_gg", "GG ε"), ("p_value", "p")])
    rows += "<h2>各单元格观测摘要</h2>" + _table(result["group_summaries"], [("arm", "组"), ("group", "条件"), ("n", "n"), ("mean", "均值"), ("sd", "SD")])
    if result["fixed_effects"]:
        rows += "<h2>模型估计单元格均值（SE 为 GLS 代入估计）</h2>" + _table(result["fixed_effects"], [("arm", "组"), ("group", "条件"), ("estimated_mean", "均值"), ("standard_error", "SE")])
    if result["contrasts"]:
        rows += "<h2>预设比较族：B − A</h2>" + _table(result["contrasts"], [("condition", "条件"), ("arm", "组"), ("group_a", "A"), ("group_b", "B"),
            ("difference_b_minus_a", "差值"), ("ci_low", "同时区间下限"), ("ci_high", "同时区间上限"), ("p_adjusted", "校正 p"), ("method", "方法")])
    rows += "<h2>缺失测量审计</h2>" + _table([r for r in result["missingness"] if not r["observed"]], [("independent_unit_id", "单位"), ("arm", "组"), ("group", "缺失条件")])
    return title, details, rows, ["repeated-001"]


def repeated_sections(result, cfg, d, figures, font, cjk):
    from .simple_report import _table
    q, fit = cfg["comparison"], result["fits"][0]
    if q["design"].startswith("two_way"):
        return two_way_sections(result, cfg, d, figures, font, cjk)
    used = d[d.exclude.astype(str).str.lower() != "true"].copy()
    # Reading normalized CSV must preserve literal subject/condition IDs like NA.
    used["group"] = used.group.astype(str)
    matrix = used.pivot(index="independent_unit_id", columns="group", values="value").reindex(columns=q["groups"])
    keys = ["repeated-001"]
    if result["residuals"]:
        keys.append("residuals-001")
    for theme, token in THEMES.items():
        with plt.rc_context(pstyle.rc(token,font,cjk)):
            fig, ax = plt.subplots(figsize=(max(5.2,.9*len(q["groups"])+1.8),4.1),layout="constrained")
            positions = np.arange(len(q["groups"]))
            for _, row in matrix.iterrows():
                ax.plot(positions,row.to_numpy(float),'-o',color=pstyle.MUTED if token["bold_labels"] else token["color"],alpha=.45,lw=.8,ms=3)
            ax.plot(positions,matrix.mean().to_numpy(),'-D',color="black",lw=1.8,ms=5,label="Observed mean")
            if result["fixed_effects"]:
                est=np.array([r["estimated_mean"] for r in result["fixed_effects"]])
                se=np.array([r["standard_error"] for r in result["fixed_effects"]])
                ax.errorbar(positions+.08,est,yerr=se,fmt='s',color=pstyle.color(token,2),ms=5,capsize=4,lw=1.2,
                            label="Model-estimated mean ± SE")
            ax.set_xticks(positions,q["groups"],rotation=20 if len(positions)>4 else 0)
            ax.set_ylabel(f"{pstyle.display_label(q['outcome'])} ({q['unit']})")
            ax.set_xlabel("Within-unit condition (declared order)")
            ax.set_title(f"{fit['n_units']} independent units; {fit['n_missing_cells']} missing cells")
            ax.legend()
            pstyle.save(fig,figures/f"repeated-001__{theme}",token,300)
            if result["residuals"]:
                r=result["residuals"]
                residual=np.array([x["conditional_residual"] for x in r])
                fitted=np.array([x["conditional_fitted"] for x in r])
                fig, axes=plt.subplots(1,2,figsize=(7.4,3.4),layout="constrained")
                axes[0].plot(fitted,residual,**pstyle.residual_style(token))
                axes[0].axhline(0,color=pstyle.MUTED,lw=.8,ls=(0,(3,2)))
                axes[0].set(xlabel="Conditional fitted value",ylabel="Conditional residual")
                (theoretical,ordered),_=probplot(residual,dist="norm")
                axes[1].plot(theoretical,ordered,**pstyle.residual_style(token))
                axes[1].set(xlabel="Normal theoretical quantile",ylabel="Ordered conditional residual")
                pstyle.save(fig,figures/f"residuals-001__{theme}",token,300)
    if q["design"]=="one_way_rm_anova":
        title="单因素重复测量 ANOVA"
        details=("同一独立单位在各条件各有一个观测。总体检验始终使用 Greenhouse–Geisser 校正自由度；"
                 "未校正 p 仅供审计。预设两两比较使用配对 t 检验、Holm 校正 p 与 Bonferroni 同时区间。"
                 "不根据总体检验结果筛选比较。图中连线为同一单位，黑线为观测均值。")
    else:
        title="线性混合效应模型：条件效应与单位随机截距"
        details=("REML 拟合：value = 条件均值 + 单位随机截距 + 残差。随机截距与残差独立、服从正态分布，"
                 "残差同方差且条件独立。固定效应协方差为代入估计方差分量的 GLS 协方差。"
                 "缺失值不填补；允许缺失时依赖给定模型下的 MAR 假设，不能由这些数据验证。"
                 "随机截距不能代替随机斜率或时间自相关。条件残差图仅供检查，不证明模型成立。")
        details+=(" 参数 bootstrap 重拟合整个单位的模拟观测，以中心化 Wald 统计量校准总体 p，"
                  "以 max-|t| 校准预设比较族的 p 和同时区间；属于依赖模型的近似，并非 KR/Satterthwaite。"
                  if q["inference"]=="parametric_bootstrap" else
                  " Satterthwaite 小样本推断（与 lmerTest 相同的方法）：单个比较用 Satterthwaite 自由度的 t 检验，"
                  "总体检验用 Fai–Cornelius 合并自由度的 F 检验；多重比较为 Holm p 与 Bonferroni 同时区间。"
                  "随机截距方差落在边界时不给出推断。"
                  if q["inference"]=="satterthwaite" else
                  " Wald z/χ² 为大样本近似，无小样本自由度校正；多重比较为 Holm p 与 Bonferroni 同时区间。")
        details+=" 图中棕色方块为模型估计均值 ± SE（描述性，不是置信区间）；数据不完整时它与观测均值不同。"
    details+=" 状态："+fit["status"]+"；诊断："+"；".join(fit["diagnostics"])
    rows="<h2>总体检验</h2>"+_table([fit],[("design","设计"),("inference","推断方法"),("reportable","推断可报告"),
        ("n_units","单位数"),("n_observations","观测数"),("f_statistic","F"),("wald_statistic","Wald χ²统计量"),
        ("epsilon_gg","GG ε"),("df_numerator","df1"),("df_denominator","df2"),("p_value","总体 p"),
        ("bootstrap_successful","成功重拟合"),("p_mc_standard_error","p 的 Monte Carlo SE"),
        ("random_intercept_variance","随机截距方差"),("residual_variance","残差方差")])
    rows+="<h2>观测数据摘要</h2>"+_table(result["group_summaries"],[("group","条件"),("n","观测单位数"),("mean","均值"),("sd","SD")])
    if result["fixed_effects"]:
        rows+="<h2>模型估计条件均值（SE 为 GLS 代入估计）</h2>"+_table(result["fixed_effects"],[("group","条件"),("estimated_mean","均值"),("standard_error","SE")])
    if result["contrasts"]:
        rows+="<h2>预设比较族：B − A</h2>"+_table(result["contrasts"],[("group_a","A"),("group_b","B"),
            ("difference_b_minus_a","差值"),("ci_low","同时区间下限"),("ci_high","同时区间上限"),
            ("p_adjusted","校正 p"),("reportable","推断可报告")])
    rows+="<h2>缺失测量审计</h2>"+_table([r for r in result["missingness"] if not r["observed"]],[("independent_unit_id","单位"),("group","缺失条件")])
    return title,details,rows,keys
