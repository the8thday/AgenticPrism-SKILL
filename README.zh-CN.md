# AgenticPrism — modular scientific analysis skills

[English](README.md) | **简体中文**

> **AI Agent：** 如需替用户安装本集合，请按 [INSTALL_FOR_AGENTS.md](INSTALL_FOR_AGENTS.md) 操作。

实际分析入口见下方各 Skill。

本地开发版 0.7.1。包含总入口、平衡结合、BLI/SPR 动力学、4PL 剂量反应、
ELISA 定量和组间统计比较六个 Skill，以及版本化的共享 Python 计算/绘图库。当前支持平衡 KD、独立
循环或单循环 1:1 动力学 kon/koff/KD（可显式参比或双参比扣除），以及相对 EC50/IC50（可固定平台、相对加权、
跨独立实验汇总，参比/待测曲线的相对效价，平行性 F 检验或预设界限的等效性检验，可选 RP 接受界限）、逐板 ELISA 标准曲线
（4PL 或预设 5PL，含标准品回算、独立质控、本板通过范围、稀释线性和以标准曲线为条件的未知浓度区间）反算未知浓度、
两组 Welch/配对 t 检验，以及三组及以上独立组的单因素 ANOVA（经典或 Welch）加预设多重比较族。
读板仪网格 + 板图可用 `import-plate` 转为长表。重复测量 ANOVA、混合效应模型和生存分析仍不在范围内。

## 从这里开始

- [总入口 Skill](skills/agentic-prism/SKILL.md)
- [平衡结合 Skill](skills/equilibrium-binding/SKILL.md)
- [BLI/SPR 动力学 Skill](skills/binding-kinetics/SKILL.md)
- [4PL 剂量反应 Skill](skills/dose-response/SKILL.md)
- [ELISA 定量 Skill](skills/elisa-quantification/SKILL.md)
- [组间比较 Skill（两组与多组）](skills/group-comparison/SKILL.md)
- [验证记录与边界](validation/README.md)

安装（见下文）后可直接执行：

```sh
.venv/bin/agentic-prism analyze --config fixtures/public/config.json --output runs/my-public-run
.venv/bin/agentic-prism render --run runs/my-public-run --style standard
.venv/bin/agentic-prism verify --run runs/my-public-run
.venv/bin/agentic-prism analyze --config fixtures/kinetics_octet/config_300s.json --output runs/my-octet-300s
.venv/bin/agentic-prism import-octet --manifest fixtures/kinetics_octet/octet_manifest.json --output runs/my-octet-import
.venv/bin/agentic-prism analyze --config fixtures/dose_synthetic/config_ic50.json --output runs/my-ic50-run
.venv/bin/agentic-prism analyze --config fixtures/dose_synthetic/config_potency.json --output runs/my-potency-run
.venv/bin/agentic-prism analyze --config fixtures/elisa_synthetic/config.json --output runs/my-elisa-run
.venv/bin/agentic-prism analyze --config fixtures/groups_synthetic/paired_config.json --output runs/my-paired-run
.venv/bin/agentic-prism analyze --config fixtures/multigroup_synthetic/config_dunnett.json --output runs/my-multigroup-run
.venv/bin/agentic-prism analyze --config fixtures/elisa_dilution_5pl/config.json --output runs/my-5pl-run
.venv/bin/agentic-prism analyze --config fixtures/kinetics_single_cycle/config.json --output runs/my-sck-run
.venv/bin/agentic-prism import-plate --manifest fixtures/plate_import_example/plate_manifest.json --output runs/my-plate-import
```

每次分析使用一个新输出目录，不覆盖旧结果。只改图用 `render`；它会验证结果哈希，
不重新拟合。报告是自包含 HTML，图形和文件下载均内嵌；可以单独拷贝。

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
（见[运行程序说明](skills/agentic-prism/references/runtime.md)）。六个 Skill 与计算包属于同一版本，需一起使用。

平台证据：安装与结果已在 macOS 上核对（uv + Python 3.13、pip + Python 3.14，结果完全一致）。
Linux 和 Windows 预期可用，但 CI 尚未实际运行。见 [0.7.1 验证记录](validation/RELEASE_0.7.1.md)。

## 示例报告

上面每条命令都会生成一个新的运行目录，其中的 `report.html` 是自包含的离线报告（图形与下载内嵌）。
仓库中不附带预先生成的示例报告。

平衡数据从 `fixtures/synthetic/config.json` 开始；动力学数据从
`fixtures/kinetics_synthetic/config.json` 开始，按照对应的
[平衡输入契约](skills/equilibrium-binding/references/input-schema.md)或
[动力学输入契约](skills/binding-kinetics/references/input-and-model.md)填写适用性证据。
4PL 使用[剂量反应输入契约](skills/dose-response/references/input-and-model.md)。
不能直接复制模拟数据中的 `true` 作为生物实验已验证的声明。

## 输出和退出状态

运行包含输入快照、配置、输入/代码/环境哈希、规范化观测、结果与诊断、
两种主题的 SVG/PDF/PNG，以及 HTML。曲线拟合模块还保存预测、残差及适用的参数区间。平衡模块另有实验汇总；
动力学模块保存逐传感器参数与联合 bootstrap 样本。
ELISA 模块保存标准品回算质控表、未知孔反算值（含区间）、样本汇总与稀释线性表；两组统计保存均值差、区间和检验结果；
多组统计保存总体检验、分组描述和预设比较族；相对效价等效性检验另存 `parallelism_equivalence.csv`。
输出区分数值审计值与可报告估计。量程外值不会自动写成有统计保证的上/下界。

- 退出 0：工作流完成，没有硬失败曲线；仍可能有受限/量程外结果。
- 退出 3：含失败曲线，已保存可检查的结果与报告。模拟示例故意包含平坦曲线。
- 退出 2：配置、输入或执行错误；已开始的分析目录含 `failure.json`。

公开 fixture 保留作者数据与来源校验值，复现原研究的特定设置。它不是所有抗体
实验的通用默认值，也不是独立实验关系的证明。仓库数据的单独授权未从保留的
README 确认，因此没有擅自给第三方数据重新指定许可证。

**不宣称 Prism 数值等价。** 平衡模块只计算 KD 的区间；不计算基线/振幅区间、曲线
置信带或预测区间。对照估计后固定的基线，其误差未传播，报告明确标注这一条件。

动力学公开示例是经过作者处理的 Octet RED384 `Results.txt` 布局，参考结果来自
TitrationAnalysis，而非 Octet 原厂软件。300 s 与 600 s 是同一次实验的两种拟合窗；
现有适配器不读取 `.frd`，SPR 可通过标准 CSV 输入，但暂无专用 Biacore 导入器。
单循环（连续进样、不再生）数据按每一行所属进样步骤声明起止时间和浓度；双参比按
(样品 − 参比通道) − (缓冲液空白 − 空白参比通道) 逐行计算，只在配置声明时执行。

4PL 拟合原始逐孔响应值，输出拟合（或固定）上下平台之间的**相对** EC50/IC50 与
profile-F 区间。IC50 不自动等于固定响应值 50 的浓度（平台固定为 0/100 时才等于），
也不等于 KD 或 Ki。相对效价 RP = C50参比 / C50待测，来自共享平台与斜率的平行模型；
默认的平行性 F 检验不是 USP <1032>/<1034> 推荐的等效性检验；0.7.0 可改用预设界限的等效性检验
（Hill 斜率比及可选的平台差 90% 区间须全部落在界限内），并可设定 RP 接受界限。界限必须来自实验室历史数据或已批准方案，软件不会也不应代为选择。
公开示例没有物理浓度单位，仅用于数值和流程演示；团队真实实验仍需明确剂量
单位、读数含义、重复结构及归一化历史。

ELISA 定量按板拟合已知浓度标准品，并对每个标准品水平回算：默认中间水平回收率 80–120%、
LLOQ/ULOQ 处 75–125%、复孔 CV ≤ 20%，至少 75% 且不少于 6 个水平通过才接受该板
（软件筛选规则，均可配置，并非完整 ICH M10 验证）。旧字段 LLOQ/ULOQ 表示本板通过范围。独立质控还可进一步限制未知孔是否可报告，
并乘以输入的稀释倍数；未通过的标准品只报告、不自动剔除。0.7.0 起可预设 5PL，未知孔与样本给出 delta 法区间：
它只反映本板标准曲线参数和孔响应噪声，不包含板间、基质或稀释操作误差。同一样本多个稀释度时检查稀释线性，
并提示可能的 hook 效应或基质干扰模式（仅为模式提示，不是诊断）。
两组比较支持预先声明的 Welch 独立组或完整配对 t 检验；三组及以上独立组支持经典或 Welch 单因素 ANOVA，
并预设 Dunnett、Tukey-Kramer、Games-Howell 或 Holm-Welch 比较族。均要求每个独立单位每条件一个值；
不支持重复测量 ANOVA、混合效应、协变量或非参数检验。

## 0.6.0 可靠性增强

- 动力学增加解离时间窗敏感性检查；强残差相关、解离信息不足等问题触发受限状态。审计区间保留，但不作为可报告置信区间。
- ELISA 增加不参与拟合的独立 QC、回收率/CV 检查和跨板描述统计。旧字段 LLOQ/ULOQ 只表示本板通过范围，不能直接声称方法定量限已验证。
- 公开数据参考库重新执行当前版本，并同时检查数值和应触发的受限状态。
- 增加 Agent 误用场景、后端拒绝测试、独立环境 wheel 安装验证和三平台 CI 配置。

详细已验证范围和尚缺证据见 [0.6.0 验证记录](validation/RELIABILITY_0.6.0.md)。真实 SPR/原厂软件对照、完整 ELISA 方法验证，以及 CI 未实际运行的平台均不标为已验证。

## 0.7.0 新增

- 多组比较：经典/Welch 单因素 ANOVA 与预设比较族，已对照独立公式和 SciPy 核对，并做族错误率与同时区间覆盖模拟。
- 读板仪网格 + 板图导入（`import-plate`），逐格原样复制并记录来源哈希。
- ELISA：可预设 5PL、未知浓度 delta 法区间、同板稀释线性。
- 相对效价：等效性平行性检验与 RP 接受界限。
- 动力学：单循环（连续进样）模型与显式双参比。
- 真实 Agent 会话场景记录（见 [agent 场景](validation/agent-scenarios/README.md)）。

实际观察到的数值、模拟结果（包括未达到预设界限的情形）和证据边界见 [0.7.0 验证记录](validation/RELEASE_0.7.0.md)。
仓库未配置 git remote，三平台 CI 尚未实际运行。
