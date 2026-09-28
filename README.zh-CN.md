# AgenticPrism — modular scientific analysis skills

[English](README.md) | **简体中文** | [版本历史（英文）](CHANGELOG.md)

面向抗体与生物药研发的统计分析 Agent Skill 集合，覆盖从结合表征、体外功能、
体内药效，到生物分析方法学验证、免疫原性和 CMC 质量（稳定性、效价、可比性、
质量标准）。AI Agent 读取 Skill，先确认分析是否适用于这个实验，再写出明确的
JSON 配置，调用同一个版本化的 Python 计算包。每次运行都保存输入、配置、哈希、
结果、诊断以及可离线打开的 HTML 报告。0.8 以后新增的模块还会生成
`interpretation_facts.json`，即 Agent 写结论必须依据的内容：主要结果、哪些可以
报告或被扣留及原因、必须说明的事项和限制。最早的六个模块（平衡结合、动力学、
剂量反应、ELISA、两组与多组比较）还没有这个文件，Agent 直接读取其结果文件。

本地开发版 0.10.1。

> **AI Agent：** 如需替用户安装本集合，请按 [INSTALL_FOR_AGENTS.md](INSTALL_FOR_AGENTS.md) 操作。

## 覆盖范围

| 研发阶段 | 要回答的问题 | Skill | 方法 |
|---|---|---|---|
| 结合表征 | 平衡亲和力（KD） | [equilibrium-binding](skills/equilibrium-binding/SKILL.md) | 单位点模型，profile-F 区间；跨独立实验汇总 |
| | 结合与解离速率 | [binding-kinetics](skills/binding-kinetics/SKILL.md) | 多循环或单循环 BLI/SPR 全局 1:1 拟合，参比或双参比扣除，带可靠性门控的块 bootstrap 区间，Octet `Results.txt` 导入 |
| 体外功能 | EC50/IC50 与相对效价 | [dose-response](skills/dose-response/SKILL.md) | 对称 4PL，相对中点及 profile 区间，平行线相对效价（F 检验或预设界限的等效性平行性） |
| | 由标准曲线求浓度 | [elisa-quantification](skills/elisa-quantification/SKILL.md) | 逐板 4PL/5PL，标准品回算质控，独立质控，delta 法未知浓度区间，稀释线性；读板仪网格导入 |
| 体内药效 | 肿瘤体积随时间变化 | [tumor-growth](skills/tumor-growth/SKILL.md) | 对数体积随机斜率模型：生长速率、倍增时间、速率差、模型 T/C；观测 TGI%、T/C% 的 Fieller 区间及脱落诊断 |
| | 生存、到达人道终点的时间 | [time-to-event](skills/time-to-event/SKILL.md) | Kaplan–Meier、log-rank（渐近或置换）、Cox 风险比、比例风险检验 |
| | 体重等按计划时间点的测量 | [repeated-measures](skills/repeated-measures/SKILL.md) | GG 校正的重复测量及裂区 ANOVA，Satterthwaite 随机截距模型，边际 US/AR(1) MMRM（可选 Kenward–Roger） |
| 生物分析与免疫原性 | 方法学验证（ICH M10 风格，配体结合法） | [method-validation](skills/method-validation/SKILL.md) | 准确度/精密度与总误差、稀释线性与钩状效应、带趋势检查的平行性、选择性、特异性、稳定性；accuracy profile |
| | ADA 切点、灵敏度、药物耐受 | [ada-cut-point](skills/ada-cut-point/SKILL.md) | 筛选/确证/滴度切点，固定或浮动，置信下限；阳性对照灵敏度及批间预测上限；药物耐受 |
| | 重复性与中间精密度 | [variance-components](skills/variance-components/SKILL.md) | 嵌套/交叉因素的 REML，不平衡数据，MLS/MOVER 区间 |
| CMC 与质量 | 长期稳定性与有效期 | [stability](skills/stability/SKILL.md) | Q1E 线性回归、先斜率后截距的可合并性检验、均值置信界限、声明支持条件后限制外推 |
| | 跨运行相对效价及验证 | [potency-assay](skills/potency-assay/SKILL.md) | 对数 RP 随机运行 REML、MLS/MOVER 中间精密度、偏倚、线性等效与实测范围；保留失败运行 |
| | 批次可比性与生物类似性（单个属性） | [comparability](skills/comparability/SKILL.md) | 批均值对声明界限的 TOST 等效性检验、声明 k 的质量范围或描述性比较；记录属性分级，不代为推断 |
| | 容忍区间与过程能力 | [specifications](skills/specifications/SKILL.md) | 精确正态及非参数容忍区间（批次不足时给出所需最少数量）；Pp/Ppk 及子组内 Cp/Cpk 与区间 |
| 各阶段通用 | 组间比较 | [group-comparison](skills/group-comparison/SKILL.md) | Welch 或配对 t 检验，单因素 ANOVA 加 Dunnett/Tukey/Games-Howell/Holm 比较族，Mann–Whitney 与 signed-rank（Hodges–Lehmann），Kruskal–Wallis + Dunn，Friedman |

[总入口 Skill](skills/agentic-prism/SKILL.md) 按实验目的选择专用 Skill；
[模块登记表](skills/agentic-prism/references/module-registry.md) 列出每个模块的确切能力和边界。

## 正确性如何核对

每种方法都在公开数据或 fixture 上与成熟实现逐值比对；区间和错误率用模拟检验，
通过/未通过的界限在运行前固定，未通过的结果如实记为未通过。详细范围和全部数字见
[validation/README.md](validation/README.md)。

| 模块 | 比对对象 | 最大差异 | 模拟校准 |
|---|---|---|---|
| 平衡结合 | BindCurve 0.2.0（25 条曲线） | KD 相对 1.3e-4 | 区间覆盖检查 |
| 动力学 | 已发表的 Octet RED384 数据 TitrationAnalysis 拟合 | 速率/KD 相对 2e-4 | bootstrap 覆盖；单循环 koff 91.3%，界限 91.4%（未通过） |
| 剂量反应 | 独立的 SciPy 四参数重拟合 | 相对 5.8e-7 | 4PL 与相对效价覆盖 |
| ELISA | 独立全参数拟合加求根 | 3.9e-10 | 未知浓度区间覆盖 |
| 组间比较 | SciPy；R `wilcox.test`、`kruskal.test`、`friedman.test`、`dunn.test` | F 相对 1e-10；秩统计量精确一致或 1e-14；HL 区间端点绝对 3e-4 以内（求根容差） | 族错误率；有 4 个单位的高方差组时 Welch 类方法未通过 |
| 重复测量、MMRM | R nlme、lmerTest、afex、emmeans、mmrm、pbkrtest | 相对 1.8e-6（自由度，MMRM） | Satterthwaite 通过；每组 12 个单位的非结构化 MMRM 未通过（7.0–7.1%） |
| 生存分析 | R survival（lung、veteran） | 相对 1.6e-14 | 置换 log-rank、KM log-log、Cox、比例风险检验均通过 |
| 肿瘤生长 | R lmerTest | 相对 1.9e-7 | 有脱落时速率差覆盖通过 |
| ADA 切点 | R base/car/lme4；rADA 示例 | 相对 2.1e-12 | 置信下限通过；点估计滴度切点未通过 |
| 精密度方差组件 | R VCA、lme4；CLSI EP05 示例 | 绝对 9.3e-7 | MLS 总方差通过；一个分量上界未通过（92.8%） |
| 方法学验证、ADA 灵敏度 | R VCA `anovaVCA`、`lm`、`binom.test`、`t.test`、`approx` | 相对 3.4e-12 | 18 行中 17 行通过；3 倍稀释间隔的灵敏度上限未通过 |
| CMC 稳定性 | R `lm`/`anova`/`emmeans`；Koleva 回归示例 | 相对 3.49e-13 | 10 行中 6 行通过；4 行覆盖不足均披露；已发表有效期数值未复现 |
| 跨运行效价 | R lme4；独立 MLS/MOVER 计算 | 相对 6.90e-7 | 12/12 行通过；未取得匹配的已发表完整算例 |
| 可比性、质量标准 | R `tolerance`（EXACT 因子）、`t.test`；NIST/SEMATECH 印出的容忍因子与过程能力数值 | 相对 3.3e-9 | 15 行中 13 行通过；n = 10 的两行过程能力因蒙特卡洛误差未通过（10 万次补充检查为 94.9%、95.4%） |

Agent 误用场景（例如没有接受标准就要求给出结论、要求省略钩状效应、把阳性对照
灵敏度说成患者的检测限）按版本记录在
[validation/agent-scenarios](validation/agent-scenarios/README.md)。
**不宣称与 GraphPad Prism 数值等价。**

## 安装

需要：git、首次安装时能联网，以及 [uv](https://docs.astral.sh/uv/)（推荐，会自动获取经过验证的 Python 3.13）或 Python 3.12 及以上。

```sh
git clone https://github.com/the8thday/AgenticPrism-SKILL.git AgenticPrism
cd AgenticPrism
python3 install.py          # Windows：py install.py
```

`install.py` 在克隆目录内建立 `.venv`，按 `requirements-lock.txt` 的版本安装依赖和本包，并做一次自检
（分析、出图和哈希校验各跑一遍）。`git pull` 更新后请重新运行；Skill 会检查运行程序与 Skill 的版本是否一致。
卸载时删除 `.venv` 即可。

**在 Agent 中使用。** 把 Agent 指向总入口 Skill：

> 使用 `<克隆目录>/skills/agentic-prism/SKILL.md` 分析我的实验文件。

如果 Agent 从某个目录自动发现 Skill（例如 Claude Code 的 `~/.claude/skills`），可以建立链接：

```sh
python3 install.py --link-skills ~/.claude/skills
```

它只创建符号链接，不会覆盖已有条目。不要复制 Skill 文件夹：复制出来的 Skill 找不到运行程序。
每个 Skill 会按自身的真实位置定位运行程序，并用 `agentic-prism doctor --collection <克隆目录>` 检查
（见[运行程序说明](skills/agentic-prism/references/runtime.md)）。十六个 Skill 与计算包属于同一版本，需一起使用。

**可选 R。** 只有 MMRM 的 Kenward–Roger 需要 R。先安装 R，再用 `python3 install.py --with-r`
安装固定版本的可选包到 `.r-lib/`；`doctor` 会显示它们。缺少 R 时只有这一方法拒绝运行，其余照常。

平台证据：安装与结果已在 macOS 上核对（uv + Python 3.13、pip + Python 3.14，结果完全一致）。
Linux 和 Windows 预期可用，但 CI 尚未实际运行。见 [0.7.1 验证记录](validation/RELEASE_0.7.1.md)。

## 使用

```sh
.venv/bin/agentic-prism analyze --config fixtures/public/config.json --output runs/my-public-run
.venv/bin/agentic-prism render --run runs/my-public-run --style standard
.venv/bin/agentic-prism verify --run runs/my-public-run
.venv/bin/agentic-prism analyze --config fixtures/kinetics_octet/config_300s.json --output runs/my-octet-300s
.venv/bin/agentic-prism import-octet --manifest fixtures/kinetics_octet/octet_manifest.json --output runs/my-octet-import
.venv/bin/agentic-prism analyze --config fixtures/dose_synthetic/config_ic50.json --output runs/my-ic50-run
.venv/bin/agentic-prism analyze --config fixtures/elisa_synthetic/config.json --output runs/my-elisa-run
.venv/bin/agentic-prism analyze --config fixtures/multigroup_synthetic/config_dunnett.json --output runs/my-multigroup-run
.venv/bin/agentic-prism analyze --config fixtures/method_validation/config_accuracy_precision.json --output runs/my-ap-run
.venv/bin/agentic-prism analyze --config fixtures/stability/config_pooled.json --output runs/my-shelf-life
.venv/bin/agentic-prism analyze --config fixtures/comparability/config_tost_absolute.json --output runs/my-comparability
.venv/bin/agentic-prism analyze --config fixtures/specification/config_normal_exact.json --output runs/my-tolerance-interval
.venv/bin/agentic-prism analyze --config fixtures/ada_performance/config_sensitivity.json --output runs/my-ada-sensitivity
.venv/bin/agentic-prism import-plate --manifest fixtures/plate_import_example/plate_manifest.json --output runs/my-plate-import
```

每次分析使用一个新输出目录，不覆盖旧结果。只改图用 `render`；它会验证结果哈希，
不重新拟合。曲线拟合报告内嵌图形与下载；ADA、方差组件和方法学验证报告可离线阅读结果，
下载链接指向同目录文件，请复制整个运行目录以保留下载。

新数据可以从对应的 fixture 配置开始，并按各 Skill 的输入契约填写适用性证据
（例如[平衡](skills/equilibrium-binding/references/input-schema.md)、
[动力学](skills/binding-kinetics/references/input-and-model.md)、
[剂量反应](skills/dose-response/references/input-and-model.md)、
[重复测量](skills/repeated-measures/references/input-and-model.md)、
[生存分析](skills/time-to-event/references/input-and-model.md)、
[肿瘤生长](skills/tumor-growth/references/input-and-model.md)、
[方法学验证](skills/method-validation/references/contract.md)、
[ADA](skills/ada-cut-point/references/contract.md)、
[方差组件](skills/variance-components/references/contract.md)、
[稳定性](skills/stability/references/input-contract.md)、
[跨运行效价](skills/potency-assay/references/input-contract.md)、
[可比性](skills/comparability/references/input-contract.md)、
[质量标准](skills/specifications/references/input-contract.md)）。
不能直接复制模拟数据中的 `true` 作为真实实验已验证的声明。

## 示例报告

上面每条命令都会生成一个新的运行目录，其中的 `report.html` 可离线阅读；请保留运行目录以保存完整产物。
仓库中不附带预先生成的示例报告。

## 输出和退出状态

运行包含输入快照、配置、输入/代码/环境哈希、规范化观测、结果与诊断、
两种主题的 SVG/PDF/PNG，以及 HTML。曲线拟合模块还保存预测、残差及适用的参数区间；
各模块另有自己的结果表（例如 ELISA 的标准品回算表、多组比较的比较族、方法学验证的逐水平结果）。
输出区分审计值与可报告估计。量程外值不会写成有统计保证的上/下界。

- 退出 0：工作流完成，没有硬失败曲线；仍可能有受限或量程外结果。
- 退出 3：含失败曲线，已保存可检查的结果与报告。模拟示例故意包含平坦曲线。
- 退出 2：配置、输入或执行错误；已开始的分析目录含 `failure.json`。

## 范围与限制

- **平衡结合：** 只计算 KD 的区间；不计算基线/振幅区间、曲线置信带或预测区间。由对照估计后固定的基线，其误差未传播，报告会注明。
- **动力学：** 公开示例是作者处理过的 Octet RED384 `Results.txt` 布局，参考结果来自 TitrationAnalysis，而非 Octet 原厂软件。不读取 `.frd`；SPR 可用标准 CSV 输入，但没有专用 Biacore 导入器。双参比只在配置声明时执行。
- **剂量反应：** 报告的是上下平台之间的**相对** EC50/IC50，平台固定为 0/100 时才等于响应值 50 处的浓度，也不等于 KD 或 Ki。默认的平行性 F 检验不是 USP <1032>/<1034> 推荐的等效性检验；等效性界限和 RP 接受界限必须来自实验室历史数据或已批准方案。
- **ELISA：** 标准品回算规则是软件筛选规则，不是完整的 ICH M10 验证；未知浓度区间只反映本板标准曲线和孔噪声，不含板间、基质或稀释误差。
- **组间比较：** 需要单位级数据和预先声明的比较族；Welch 类方法在组内少于 6 个单位时会提示可能偏宽松。
- **体内与纵向数据：** 重复测量用复合对称随机截距、裂区 ANOVA 或边际 US/AR(1) MMRM；肿瘤生长为对数尺度的指数生长，仅支持 MAR 脱落；生存分析不含竞争风险、脆弱性模型或区间删失。小样本未通过项见验证记录。
- **生物分析与免疫原性：** 方法学验证从回算浓度开始（不拟合标准曲线，不含残留和 ISR），接受标准必须由用户声明。ADA 切点需要单一试剂批次的完整平衡阴性样本组；动态切点未实现。阳性对照灵敏度不等于患者抗体的灵敏度。
- **CMC：** 稳定性采用 ICH Q1E 线性回归及 α = 0.25 的可合并性预检验，这种先检验再选模型的做法在部分设计下会降低覆盖率（已披露）；外推需要声明支持数据。效价合并要求每轮有重复测定。可比性一次只比较一个属性，不给出整体的生物类似性结论；容忍区间描述数据，本身不能作为质量标准。属性分级、界限、k 值和标准限都必须由用户声明。
- **数据来源：** 公开 fixture 保留作者数据与来源校验值，只复现原研究的特定设置，不是通用默认值。未确认第三方数据的单独授权时，不代为指定许可证。
