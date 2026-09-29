# 0.11.1–0.13.0 本地发布验证汇总

四个开发版本按顺序建立分支。工程检查通过不代表五个科学证据门槛全部通过；
新方法的未满足范围和实验性状态保留。没有 push、合并 main、发布或公开同步。
PK/PD 0.14.0 延后，未加入 NCA、房室 PK、TMDD、PD 或 population NLME。

## 每个版本的检查

| 版本 / 分支 | pytest | 数值比较通过 | 校准通过 / 失败判定 | legacy 配置 / 一致文件 | wheel smoke |
|---|---:|---:|---:|---:|---:|
| [0.11.1](RELEASE_0.11.1.md) / `dev/0.11.1` | 585 | 66/66 | 62 / 26 | 110 / 922 | 136 |
| [0.12.0](RELEASE_0.12.0.md) / `dev/0.12.0` | 602 | 16/19 | 4 / 44 | 136 / 1055 | 141 |
| [0.12.1](RELEASE_0.12.1.md) / `dev/0.12.1` | 612 | 3/3 | 10 / 0 | 141 / 1098 | 147 |
| [0.13.0](RELEASE_0.13.0.md) / `dev/0.13.0` | 652 | 16/20 | 13 / 6 | 147 / 1128 | 154 |

每个版本均完成 clean git archive 构建、独立 wheel-only 环境检查和 changed-Skill 校验。
每次全量 pytest 有8个既有 warnings。具体 archive commit、哈希、安装精确比较和
工作平台范围见各 `release_checks_<version>.json`。最终分支 head 在交付消息中列出；
archive commit 是构建输入，后续验证记录提交不改变该 wheel 的运行时代码。

## 实现范围、跳过范围与五个证据门槛

### 0.11.1

- 代码：独立两因素 ANOVA（SS II/III、交互及 EMM）、列联表/比例、相关与
  OLS/WLS、Deming、Passing–Bablok/CUSUM、正态 Bland–Altman LoA；补齐两组/多组 facts。
- 数值：66项通过；确定性比较最大相对差5.509263432e-5。
  RxC Monte Carlo p值相对差0.114267，在注册的4MCSE加离散分辨率容差内。
- 公开实例：factorial 与 mcr 肌酐实例复现；Agresti 列联表、匹配的 PEF/LoA、
  独立相关/回归 published worked example 为 UNMET。
- 校准：32,000条注册模拟记录，88个判定中62通过、26失败。
  相关/回归、RR、趋势及 McNemar 缺少独立校准行，此范围仍未满足。
- live-agent：6场景 PENDING。facts：新类型及两个旧模块提取/状态/字节检查通过。

### 0.12.0

- 代码：opt-in 漂移、异质配体、双价分析物、传质、off-rate screening；共享 solve_ivp
  原语；实测且有许可的 T200 XY 文本/Carterra XY 工作簿导入；动力学 facts。
  默认1:1保留。未猜测其他布局，native binary/.frd/未核实布局仍不支持。
- 数值：4项正向计算、4项 fixture refit、4项 public refit、7项 koff 比较中16/19通过。
  deSolve正向最大相对差5.604716907e-7；koff最大3.733090264e-7。
  失败：双价 fixture R未收敛（相对差2.997120263e-5）；公开异质/双价分别
  0.02327040303、0.02305965544。原始 refit 首次失败记录也保留。
- 公开实例：有真实导出文件，但匹配的复杂机制论文拟合复现仍 UNMET；共享 Rmax
  与作者逐曲线 Rmax 不同，不声称复制其生物学拟合。
- 校准：9,000原始加3,000补充记录；4个判定通过，44个覆盖率判定失败。
  真实1:1下误选：异质0/987（13拟合失败）、漂移1/1000、双价1/1000、传质0/1000。
  这些结果只针对预声明的 AICc+可辨识性规则及各备选模型，不支持任意择优选模。
- live-agent：5场景 PENDING。facts：默认/高级/部分失败/steady-state 状态保留，旧文件一致。

### 0.12.1

- 代码：tandem/premix/sandwich 有向阻断矩阵、声明对照与阈值、不对称诊断、
  平均连接聚类、条件 feature-bootstrap BP、双向阻断图分组。接受标准响应表；
  未猜测 vendor binning parser。没有 multiscale AU 或结构表位身份推断。
- 数值：3/3通过；1000个匹配R重采样下，4个 synthetic / 10个 public-lung 分支的
  BP差为0；最大相对 linkage 差3.095360458e-15，归一化最大相对差1.249926126e-14。
- 公开实例：Abdiche文章可读且CC BY，但缺少带匹配对照的数字响应矩阵，UNMET。
  肺数据只验证通用聚类算法，不代替表位实验。
- 校准：4,000条；10个可评估判定全通过，0 miss。两个人群的覆盖率分别
  n4=0.943、n8=0.951、n16=0.955911824（n=998）。图分组恢复率1、0.999、0.996。
  缺失对照1000/1000 withheld；对应两项coverage n=0，不算通过。
- live-agent：3场景 PENDING。facts：六种状态、禁用拟合器的重提取与30文件一致检查通过。

### 0.13.0

- 代码：Bliss/Loewe/HSA/ZIP、独立矩阵区间与模型分歧；HTS Z-prime、SSMD、
  median-polish B-score、板位置诊断与名义 predictive-t/BH 筛选；dose/ELISA facts，
  最早六个模块的 retrofit 完成。条件 ZIP 4PL参数与失败原因保存在结果中。
- 数值：16/20通过。两个公开 ONEIL 的 Bliss/HSA 接近机器精度；Loewe仅21/25、
  17/25单元可计算，gate失败。公开 ZIP最大绝对差0.299139007、76.329718890百分点，
  相对差0.023823519、0.782433046，超出0.05百分点容差；工作流withhold，不当作oracle通过。
  Synthetic Loewe绝对差0.466203（容差1.1，R使用100效应级网格），不是精确等价。
  公共 KcViab median-polish残差差0，B-score相对差1.37093e-14。
- 公开实例：ONEIL支持限于点计算，重复行无独立实验ID。KcViab缺少重复对照，
  完整HTS QC/命中 worked-example UNMET。cellHTS2当前环境无法安装，采用明确允许的
  base R medpolish oracle；归档包数据许可及哈希保留。
- 校准：4,000条，19个判定13通过、6失败，另1个条件withholding行n=0。
  HSA固定/自由平台覆盖率0.915/0.928；HTS全零假设假命中概率0.100，超过0.063784。
  因此要求的已验证FDR控制为 UNMET，输出仅为实验性的名义BH筛选。
  失败对照1000/1000 withheld。未调整阈值或替换原始失败行。
- live-agent：4场景 PENDING。facts：新类型全部fixture状态及 dose/ELISA布尔和status-only状态检查。

## Self-review：领域审阅重点

- **0.11.1**：重复单位不得充当独立样本；SS类型和比较家族必须声明；关联不等于一致。
  Deming误差方差比有来源；正态LoA及共同方差ANOVA的假设失败由保留的压力场景说明。
- **0.12.0**：双价信号为X1+X2而占位为X1+2X2，二者不混用；第一臂统计因子与第二步单位明确。
  共享响应尺度/Rmax、有效参比、valency和机制预声明；profiles、相关性、时间窗和残差可withhold。
  传感曲线/时间点不是生物学独立实验；没有单一双价KD。
- **0.12.1**：读出须对应实验格式；自阻断需实际抑制而非归一化恒等于1。
  对照必须同尺度；保留两个方向、缺失和不对称；连通分组不等于clique或同一结构表位。
  panel特征重采样不替代独立实验。
- **0.13.0**：记录检测终点说明或明确未知；同条件/同百分比尺度/同剂量网格；Loewe声明并检查效应范围。
  独立矩阵支持区间，技术孔不能制造独立性。HTS声明96/384孔尺寸并拒绝整行/整列缺失；
  位置校正的布局/多数无效应/正态假设必须成立，现有全零假设校准仍显示FDR控制不足。

## 偏差与限制的原因

- 0.11.1首次校准在JSON序列化处失败；仅修复类型序列化、保持种子/设计/边界并重跑，首轮stdout保留。
  Dunnett初始仅相对容差失败（绝对差2.7991e-8），最终记录加入1e-7绝对积分容差，首次失败不隐藏。
  Spearman仅n<=9无ties时枚举；其他用R exact=FALSE约定。缺少的独立校准范围如上。
- 0.11.1 wheel先于README Skill数量修正构建；后续clean-archive sdist修正文档，运行时字节相同。
  其validation标题遗漏的0.11.0版本字样在0.12.0修正，历史分支不重写。
  clean-environment先于扩展旧分组exact列表运行，额外独立wheel-only精确检查覆盖30配置。
- 0.12.0异质模型的原定adequate设计浓度范围足够，但快结合时间采样不足，1000/1000 withheld；
  这一设计要求未满足，原标签、结果和miss保留，不把它事后改名为stress或替换。
  注册使用50次bootstrap，不验证默认200次。原始每次模拟保留结果/状态，不序列化全部原始曲线或区间端点；固定seed重建输入。
- 0.12.0补齐其他模型的null-selection检查单独预注册于a3dc3f8；747c02f仅优化执行，
  对已确定拒绝的拟合跳过不能挽回该拒绝的bootstrap。等价性先测试，1104条原完整bootstrap记录保留复用。
  原9行覆盖率进程未改变，任何失败未被补充记录替代。
- 后续版本草稿曾在仓库外预检以利用计算等待时间；正式分支来自上一版本最终提交，所有校准均先提交后运行。
- 公开数据不可用、无法满足完整设计或oracle未一致的范围均标UNMET；未用synthetic冒充published。
  0.13.0的FDR控制要求未验证，原因是冻结设计中实测0.100假命中率；没有为过关调低阈值。
- 所有live-agent执行按用户要求延后；没有运行Claude CLI或用自身代理替代通过证据。

## 所有校准失败判定

以下逐项保留所有miss，包括同一场景在不同报告人群中的判定；rate/MCSE/bound均为比例。
原始每次模拟、种子及完整通过行见各release记录链接的JSON。

| 版本 | 场景 | 指标 | 人群 | n | rate | MCSE | bound |
|---|---|---|---|---:|---:|---:|---:|
| 0.11.1 | `anova_2_heteroscedastic` | `C(a, Sum)` | all_evaluable | 1000 | 0.216000000 | 0.013013224 | 0.063784049 |
| 0.11.1 | `anova_2_heteroscedastic` | `C(a, Sum)` | reportable_only | 1000 | 0.216000000 | 0.013013224 | 0.063784049 |
| 0.11.1 | `anova_2_heteroscedastic` | `C(a, Sum):C(b, Sum)` | all_evaluable | 1000 | 0.284000000 | 0.014259874 | 0.063784049 |
| 0.11.1 | `anova_2_heteroscedastic` | `C(a, Sum):C(b, Sum)` | reportable_only | 1000 | 0.284000000 | 0.014259874 | 0.063784049 |
| 0.11.1 | `anova_2_heteroscedastic` | `C(b, Sum)` | all_evaluable | 1000 | 0.188000000 | 0.012355404 | 0.063784049 |
| 0.11.1 | `anova_2_heteroscedastic` | `C(b, Sum)` | reportable_only | 1000 | 0.188000000 | 0.012355404 | 0.063784049 |
| 0.11.1 | `anova_3_equal` | `C(a, Sum):C(b, Sum)` | all_evaluable | 1000 | 0.073000000 | 0.008226239 | 0.063784049 |
| 0.11.1 | `anova_3_equal` | `C(a, Sum):C(b, Sum)` | reportable_only | 1000 | 0.073000000 | 0.008226239 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(a, Sum)` | all_evaluable | 1000 | 0.219000000 | 0.013078188 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(a, Sum)` | reportable_only | 1000 | 0.219000000 | 0.013078188 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(a, Sum):C(b, Sum)` | all_evaluable | 1000 | 0.265000000 | 0.013956181 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(a, Sum):C(b, Sum)` | reportable_only | 1000 | 0.265000000 | 0.013956181 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(b, Sum)` | all_evaluable | 1000 | 0.259000000 | 0.013853483 | 0.063784049 |
| 0.11.1 | `anova_3_heteroscedastic` | `C(b, Sum)` | reportable_only | 1000 | 0.259000000 | 0.013853483 | 0.063784049 |
| 0.11.1 | `wilson_n20_p0.05` | `coverage` | all_evaluable | 1000 | 0.902000000 | 0.009401915 | 0.936215951 |
| 0.11.1 | `wilson_n20_p0.05` | `coverage` | reportable_only | 1000 | 0.902000000 | 0.009401915 | 0.936215951 |
| 0.11.1 | `wilson_n100_p0.5` | `coverage` | all_evaluable | 1000 | 0.934000000 | 0.007851369 | 0.936215951 |
| 0.11.1 | `wilson_n100_p0.5` | `coverage` | reportable_only | 1000 | 0.934000000 | 0.007851369 | 0.936215951 |
| 0.11.1 | `loa_n30_skew` | `limit_0_coverage` | all_evaluable | 1000 | 0.015000000 | 0.003843826 | 0.936215951 |
| 0.11.1 | `loa_n30_skew` | `limit_0_coverage` | reportable_only | 1000 | 0.015000000 | 0.003843826 | 0.936215951 |
| 0.11.1 | `loa_n30_skew` | `limit_1_coverage` | all_evaluable | 1000 | 0.430000000 | 0.015655670 | 0.936215951 |
| 0.11.1 | `loa_n30_skew` | `limit_1_coverage` | reportable_only | 1000 | 0.430000000 | 0.015655670 | 0.936215951 |
| 0.11.1 | `loa_n100_skew` | `limit_0_coverage` | all_evaluable | 1000 | 0.000000000 | 0.000000000 | 0.936215951 |
| 0.11.1 | `loa_n100_skew` | `limit_0_coverage` | reportable_only | 1000 | 0.000000000 | 0.000000000 | 0.936215951 |
| 0.11.1 | `loa_n100_skew` | `limit_1_coverage` | all_evaluable | 1000 | 0.178000000 | 0.012096115 | 0.936215951 |
| 0.11.1 | `loa_n100_skew` | `limit_1_coverage` | reportable_only | 1000 | 0.178000000 | 0.012096115 | 0.936215951 |
| 0.12.0 | `one_to_one_drift` | `ka_coverage` | evaluable | 1000 | 0.899000000 | 0.009528851 | 0.936215951 |
| 0.12.0 | `one_to_one_drift` | `kd_coverage` | evaluable | 1000 | 0.877000000 | 0.010386096 | 0.936215951 |
| 0.12.0 | `one_to_one_drift` | `ka_coverage` | reportable | 879 | 0.902161547 | 0.010020803 | 0.935297799 |
| 0.12.0 | `one_to_one_drift` | `kd_coverage` | reportable | 879 | 0.882821388 | 0.010848412 | 0.935297799 |
| 0.12.0 | `heterogeneous_ligand` | `ka_coverage` | evaluable | 1000 | 0.902000000 | 0.009401915 | 0.936215951 |
| 0.12.0 | `heterogeneous_ligand` | `ka2_coverage` | evaluable | 1000 | 0.906000000 | 0.009228434 | 0.936215951 |
| 0.12.0 | `heterogeneous_ligand` | `kd_coverage` | evaluable | 1000 | 0.903000000 | 0.009359006 | 0.936215951 |
| 0.12.0 | `heterogeneous_ligand` | `kd2_coverage` | evaluable | 1000 | 0.920000000 | 0.008579044 | 0.936215951 |
| 0.12.0 | `bivalent_analyte` | `ka_coverage` | evaluable | 1000 | 0.889000000 | 0.009933730 | 0.936215951 |
| 0.12.0 | `bivalent_analyte` | `ka2_coverage` | evaluable | 1000 | 0.881000000 | 0.010239092 | 0.936215951 |
| 0.12.0 | `bivalent_analyte` | `kd_coverage` | evaluable | 1000 | 0.896000000 | 0.009653186 | 0.936215951 |
| 0.12.0 | `bivalent_analyte` | `kd2_coverage` | evaluable | 1000 | 0.885000000 | 0.010088360 | 0.936215951 |
| 0.12.0 | `bivalent_analyte` | `ka_coverage` | reportable | 883 | 0.886749717 | 0.010664483 | 0.935331137 |
| 0.12.0 | `bivalent_analyte` | `ka2_coverage` | reportable | 883 | 0.877689694 | 0.011026094 | 0.935331137 |
| 0.12.0 | `bivalent_analyte` | `kd_coverage` | reportable | 883 | 0.894677237 | 0.010330325 | 0.935331137 |
| 0.12.0 | `bivalent_analyte` | `kd2_coverage` | reportable | 883 | 0.886749717 | 0.010664483 | 0.935331137 |
| 0.12.0 | `mass_transport` | `ka_coverage` | evaluable | 1000 | 0.888000000 | 0.009972763 | 0.936215951 |
| 0.12.0 | `mass_transport` | `kd_coverage` | evaluable | 1000 | 0.882000000 | 0.010201765 | 0.936215951 |
| 0.12.0 | `mass_transport` | `km_coverage` | evaluable | 1000 | 0.873000000 | 0.010529530 | 0.936215951 |
| 0.12.0 | `mass_transport` | `ka_coverage` | reportable | 890 | 0.887640449 | 0.010585918 | 0.935388938 |
| 0.12.0 | `mass_transport` | `kd_coverage` | reportable | 890 | 0.882022472 | 0.010812957 | 0.935388938 |
| 0.12.0 | `mass_transport` | `km_coverage` | reportable | 890 | 0.874157303 | 0.011117671 | 0.935388938 |
| 0.12.0 | `off_rate_screening` | `kd_coverage` | evaluable | 1000 | 0.870000000 | 0.010634848 | 0.936215951 |
| 0.12.0 | `off_rate_screening` | `kd_coverage` | reportable | 979 | 0.872318693 | 0.010666204 | 0.936068899 |
| 0.12.0 | `transport_slow` | `ka_coverage` | evaluable | 1000 | 0.879000000 | 0.010313050 | 0.936215951 |
| 0.12.0 | `transport_slow` | `kd_coverage` | evaluable | 1000 | 0.881000000 | 0.010239092 | 0.936215951 |
| 0.12.0 | `transport_slow` | `km_coverage` | evaluable | 1000 | 0.882000000 | 0.010201765 | 0.936215951 |
| 0.12.0 | `transport_slow` | `ka_coverage` | reportable | 886 | 0.876975169 | 0.011035021 | 0.935355993 |
| 0.12.0 | `transport_slow` | `kd_coverage` | reportable | 886 | 0.886004515 | 0.010676894 | 0.935355993 |
| 0.12.0 | `transport_slow` | `km_coverage` | reportable | 886 | 0.874717833 | 0.011121459 | 0.935355993 |
| 0.12.0 | `transport_fast` | `ka_coverage` | evaluable | 705 | 0.887943262 | 0.011880016 | 0.933583437 |
| 0.12.0 | `transport_fast` | `kd_coverage` | evaluable | 705 | 0.880851064 | 0.012201178 | 0.933583437 |
| 0.12.0 | `transport_fast` | `km_coverage` | evaluable | 705 | 0.879432624 | 0.012263703 | 0.933583437 |
| 0.12.0 | `transport_fast` | `ka_coverage` | reportable | 222 | 0.765765766 | 0.028424751 | 0.920744986 |
| 0.12.0 | `transport_fast` | `kd_coverage` | reportable | 222 | 0.725225225 | 0.029960454 | 0.920744986 |
| 0.12.0 | `transport_fast` | `km_coverage` | reportable | 222 | 0.693693694 | 0.030937510 | 0.920744986 |
| 0.12.0 | `bivalent_short_association_stress` | `ka_coverage` | evaluable | 955 | 0.877486911 | 0.010609866 | 0.935894934 |
| 0.12.0 | `bivalent_short_association_stress` | `ka2_coverage` | evaluable | 955 | 0.910994764 | 0.009214336 | 0.935894934 |
| 0.12.0 | `bivalent_short_association_stress` | `kd_coverage` | evaluable | 955 | 0.878534031 | 0.010570728 | 0.935894934 |
| 0.12.0 | `bivalent_short_association_stress` | `kd2_coverage` | evaluable | 955 | 0.894240838 | 0.009951409 | 0.935894934 |
| 0.12.0 | `bivalent_short_association_stress` | `ka_coverage` | reportable | 554 | 0.848375451 | 0.015237858 | 0.931480813 |
| 0.12.0 | `bivalent_short_association_stress` | `ka2_coverage` | reportable | 554 | 0.880866426 | 0.013763131 | 0.931480813 |
| 0.12.0 | `bivalent_short_association_stress` | `kd_coverage` | reportable | 554 | 0.877256318 | 0.013941450 | 0.931480813 |
| 0.12.0 | `bivalent_short_association_stress` | `kd2_coverage` | reportable | 554 | 0.889891697 | 0.013299145 | 0.931480813 |
| 0.13.0 | `combination_fixed_plateaus` | `HSA_coverage` | evaluable | 1000 | 0.915000000 | 0.008819014 | 0.936215951 |
| 0.13.0 | `combination_fixed_plateaus` | `HSA_coverage` | reportable | 1000 | 0.915000000 | 0.008819014 | 0.936215951 |
| 0.13.0 | `combination_free_plateaus` | `HSA_coverage` | evaluable | 1000 | 0.928000000 | 0.008174105 | 0.936215951 |
| 0.13.0 | `combination_free_plateaus` | `HSA_coverage` | reportable | 1000 | 0.928000000 | 0.008174105 | 0.936215951 |
| 0.13.0 | `hts_null_hits` | `false_hit` | evaluable | 1000 | 0.100000000 | 0.009486833 | 0.063784049 |
| 0.13.0 | `hts_null_hits` | `false_hit` | reportable | 1000 | 0.100000000 | 0.009486833 | 0.063784049 |

合计 **76** 个失败判定。没有删除、事后重设门槛或以补充检查替换这些原始结果。

## PENDING live-agent reviewer commands

每条命令在其对应分支执行；`--release`选择场景，不会切换运行时代码。保留原始trace，由人审阅。

```sh
git switch dev/0.11.1
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0111 --release 0.11.1 --attempt reviewer-1
git switch dev/0.12.0
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0120 --release 0.12.0 --attempt reviewer-1
git switch dev/0.12.1
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0121 --release 0.12.1 --attempt reviewer-1
git switch dev/0.13.0
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0130 --release 0.13.0 --attempt reviewer-1
```

## 独立复核（2026-09-29，reviewer，不是实现者）

复核在各分支的独立 worktree 上进行，没有改动任何运行时代码、校准或数值记录。
`dev/0.13.0` 全量 pytest 652 通过。

**校准失败的原因核查**

- `anova_3_equal` 交互项 7.3% 属于蒙特卡洛偶然波动，不是实现错误。同一份数据上，
  statsmodels 的交互项 p 值与直接嵌套模型 F 检验相差不超过 1.4e-14。另用新种子重模拟
  10,000 次得到 5.59%，再用 numpy 直接做 F 检验模拟 40,000 次得到 4.82%，合计 51,000 次约 5.0%。
  注意：0.11.1 的 ANOVA、Newcombe 和 Fisher 校准直接调用第三方库，没有经过本包的封装层。
  ANOVA 的核心 p 值与本包用的是同一个计算，但数据处理和 EMM 不在校准范围内。
- Wilson 区间实现正确。精确覆盖率（对二项分布逐点求和，无需模拟）：n=20、p=0.05 为
  92.45%（Wilson 的已知性质）；n=100、p=0.5 为 94.31%，高于界限，模拟结果 93.4% 未达标
  只是蒙特卡洛误差。比例区间以后应直接报告精确覆盖率。
- 0.12.0 三项数值比较失败都不是实现错误。bivalent fixture 两侧 SSE 相对差约 5e-7，
  只是 R 报告未收敛。公开 T200 的异质/双价拟合两侧 SSE 一致到约 1e-8，参数差 2.3%，
  说明目标函数在该方向上很平（参数弱确定）。1% 的参数容差不适合这种情况，
  以后应以 SSE 一致、参数差相对区间宽度较小作为判据。
- 0.12.0 的 44 个覆盖率失败不能用来判断实际默认设置。注册校准用 50 次 bootstrap，
  默认是 200 次。0.6.0 的 1:1 记录在默认 200 次、独立噪声下覆盖率为 93–94.5%，
  因此块 bootstrap 本身略偏乐观，再加上 B=50，大致可以解释 87–90%。
  结论：**复杂动力学模型在默认设置下的区间没有经过校准**。另外，transport_fast 中
  通过全部可报告门槛的 222 个拟合覆盖率只有 69–77%，说明现有门槛拦不住这类不可靠区间。
- HTS 全零假设下假命中率偏高的根因已定位。同一设计下，不做 median polish 时
  P(任一假命中) = 4.3%，做 median polish 后为 9.3%（各 3,000 次）：polish 后的残差不再
  独立同方差（样本/阴性残差方差比中位数约 1.10），只用 16 个阴性孔估计尺度的
  predictive t 偏乐观。常见替代做法更差：全板 MAD 的 B-score 配正态参考为 26%，
  样本 MAD 为 53%。修复需要按实际板布局用模拟或置换构造参考分布；目前的实验性标注和披露是合适的。
- 预注册顺序全部正确：每个校准脚本都先于结果提交。0.12.0 补充检查在预注册后改过执行方式，
  但只跳过"已被拒绝"拟合的 bootstrap，不可能改变结论，部分原始记录也保留了。
  0.12.0、0.12.1、0.13.0 的校准 JSON 没有记录脚本哈希，可追溯性弱于 0.11.0；
  由 git 历史可确认这些脚本只提交过一次。

**live-agent（reviewer 用 Claude Code 运行）**：0.11.1 6/6、0.12.0 5/5、0.12.1 2/3、0.13.0 4/4。
0.12.1 `bins_by_eye` 失败且可复现（三次中两次 0 次工具调用，未读 Skill 就给出按肉眼合并的选项）。
详见各版本 `agent-scenarios/<ver>/RESULTS.md`。另有场景设计问题：多份输入与 fixture 逐字节相同，
0.13.0 部分代理因此沿用了 fixture 的声明。

**建议的后续工作（本次未实施）**

1. 为复杂动力学模型补做一次预注册的 B=200 校准，至少覆盖漂移、传质和双价模型；
   在此之前，复杂模型区间应在 facts 中明确标为"默认设置未校准"。
2. 传质模型在传质不限速（快速传质）时，即使通过门槛，区间也不可靠；应增加扣留条件，或把区间降为审计值。
3. HTS 命中：按板布局构造经过校准的参考分布后，才能去掉"实验性"标注。
4. 比例区间的校准改为精确覆盖率计算；校准脚本应调用本包的封装函数，并在 JSON 中记录脚本哈希。
5. 场景输入不要与 fixture 完全相同。
