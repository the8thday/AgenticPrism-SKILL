# AgenticPrism 0.8.1：解读依据与验证记录

## 本轮范围

分支 `dev/0.8.1`，起点为 `dev/0.8.0` 的
`d0f56d0`。本版完成 `interpretation_facts.json` 的设计，以及
repeated-measures、time-to-event、tumor-growth 三个工作流的接入。
没有新增统计方法、配置选项或 R 运行时调用。未开始 0.9。

| 路线图项目 | 状态 | 范围 |
|---|---|---|
| 解读依据文件：设计与三个优先模块 | 完成 | 保存结果的确定性提取；来源指针、哈希、主要结果、区间方法、扣留状态、诊断及限制；报告可离线下载 |
| 解读依据文件：其余模块 | 未开始 | 继续使用原有结果与诊断约定 |
| Mann–Whitney/Wilcoxon + HL、Kruskal–Wallis + Dunn、Friedman | 未开始 | 无实现或校准声明 |
| MMRM、AR(1)、Kenward–Roger、两因素 bootstrap、单因素 B=999 再校准 | 未开始 | 本轮重复测量模拟仅检查原有 RM ANOVA 和 RI/Satterthwaite；B=999 仅用于原有 fixture 的复现 |
| Gompertz/logistic、cage、MNAR、零假设下有脱落校准 | 未开始 | 重跑的旧肿瘤零假设场景仍实际没有脱落 |
| 竞争风险、RMST、聚类生存 | 未开始 | 原有方法保持不变 |
| Linux/Windows CI、独立人工复评 | 未执行 | 本轮安装证据来自同一台 macOS 的隔离环境 |

优先把第一个项目完成到可审查、可复现；未留下尚未接通的新统计方法。

## 五项门槛与证据

1. **R 对照与 fixture：**重跑三个现有模块的 R 基准，另新增从事实文件
   提取数值的 R ovarian 示例及生存 fixture 对照，见下表。事实层本身
   不产生新统计推断；它与保存结果的 487 个值检查完全相等。
2. **公开完整示例：**复现 R 官方 `survival::survdiff` 手册的第一个
   ovarian 两样本示例，并核对其 Cox score 恒等式。来源、R help 文本、
   数据、R 脚本和参考输出的哈希见
   [provenance.json](../fixtures/interpretation_ovarian/provenance.json)。
   这是官方手册的数值示例；没有把其补充 KM/Cox 输出冒充论文中的已发表数值。
3. **预先固定的模拟：**三行新的重复测量事实层回归模拟，以及旧生存八行、
   肿瘤三行模拟的重跑，每行均生成 1,000 个数据集。种子、设计和界限均在
   脚本中预先给定；下文逐行列出结果和对照项不足。未改变现有统计方法。
4. **Skill 与实时场景：**更新三个专用 Skill 和总入口。Claude CLI 的
   有效场景共两轮、每轮三项；最终三项核心误用均被拦截，但完整流程
   F1 为部分通过，F2/F3 通过且有范围限制。完整评分见
   [场景记录](agent-scenarios/0.8.1/RESULTS.md)。不作全面代理可靠性声明。
5. **结果进入事实文件：**每次新运行在 manifest 前写出事实文件；
   保存值、来源、扣留状态、重绘不变性和篡改拒绝均有测试。

### 文件与解释规则

[文件约定](../skills/agentic-prism/references/interpretation-facts.md)独立于
包版本使用 schema 1。每条结果记录 reportable、limited 或 withheld；
withheld 的报告层估计、检验和区间不输出审计数值。

- 混合模型失败、方差边界、bootstrap 重拟合不足、零方差配对差值分别
  保留原因。工作流完成不等于所有子模型都可报告。
- KM、中位数、log-rank 和 Cox 分开处理。中位数未达到不补数字；中位数
  的下限可单独标注。零信息自由度的 log-rank 审计 p=1 不进入科学结论。
  降秩总体检验保留并限缩解释；因旧两两结果没有保存风险集秩，事实层
  保守扣留该运行的两两 log-rank，原始结果不改。
- 无事件前 Greenwood 区间 1–1 的原值保留，附加退化区间限制。
- 肿瘤模型 T/C 是 V + offset 的几何均值比；倍增时间也针对 V + offset。
  速率差的 p 值不附在 T/C 上。观测 TGI/T/C 保留当天在测动物数、脱落
  偏倚和缺失 Fieller 区间的限制。
- facts 记录用户声明，不能证明独立性、MAR、删失合理性或研究终点定义。
  指导要求遇到矛盾时向用户澄清。历史运行不补写新文件。

## 数值对照：本次观测结果

R 4.6.0；nlme 3.1.169、lme4 2.0.6、lmerTest 3.2.1、pbkrtest
0.5.5、afex 1.5.1、emmeans 2.0.4、survival 3.8.6。
R 只用于验证，额外验证包仍放在 git-ignored `.r-lib/`。
本版没有 `r_bridge`，没有方法改用 R 实现。

| 对照范围 | 观测最大差异 |
|---|---|
| 事实文件中的公开 ovarian log-rank χ²、p vs 当次 R 输出 | 相对 `9.172633979413608e-15`，绝对 `3.774758283725532e-15` |
| ovarian log-rank vs 独立 Cox score 恒等式 | 绝对 `3.774758283725532e-15` |
| 生存合成 fixture：事实文件 log-rank χ²/p、Cox HR 与个别区间 vs R | 相对 `2.657067022562994e-15` |
| Orthodont RM ANOVA F / GG epsilon | 相对 `7.771561172376096e-16` / `2.220446049250313e-16` |
| Orthodont 随机截距：β / SE / 残差方差 | 相对 `8.82818262937235e-11` / `1.3694801292984948e-8` / `1.575532435627025e-7` |
| 单因素 Satterthwaite，公开数据和完整/缺失 fixture：contrast df / t | 相对 `2.2568563351477167e-8` / `1.2280492067695548e-8` |
| 两因素公开 BodyWeight、ChickWeight、Orthodont 与 fixture：Type III F/ddf | 相对 `9.096762332294617e-7` |
| 同上 contrast：估计 / df / t | 估计绝对 `2.821828104515589e-8`；df/t 相对 `1.4294302763673272e-7` / `3.827019479096805e-8` |
| lung/veteran KM 生存率和限值 / 中位数 | 绝对 `7.771561172376096e-16` / `0` |
| lung/veteran Cox β、SE、检验、loglik | 所比较字段相对最大 `1.6209256159527285e-14` |
| lung/veteran PH 检验 | 相对 `3.6637359812630166e-15` |
| 肿瘤随机斜率：公开替代数据与 fixture 的速率 / SE / df | 相对 `4.835734035424366e-10` / `5.7982335377460004e-8` / `1.9292120967406845e-7` |

ovarian 结果为 χ²=`1.0627398612914138`，p=`0.30259111698909225`。
该示例来自[官方手册](https://www.stat.ethz.ch/R-manual/R-devel/library/survival/html/survdiff.html)。
BodyWeight/ChickWeight 是计算对照用的公开纵向数据，不是肿瘤药效证据。
上述一致性只针对脚本列出的模型、数据和字段；不宣称 Prism 等价。

机器记录：
[interpretation](interpretation_benchmarks_0.8.1.json)、
[repeated](repeated_benchmarks_0.8.1.json)、
[survival](survival_benchmarks_0.8.1.json)、
[tumor growth](tumor_growth_benchmarks_0.8.1.json)。

### 原有结果不变

`validate_interpretation_benchmarks.py` 从 `dev/0.8.0` 提取原始源码到临时
目录，以同一 Python 环境运行旧版与新版。**12 套配置、170 个既有科学产物
逐字节相同**；新增文件只有 `interpretation_facts.json`。
六个数值源码文件的 SHA-256 与基线一致，记录在上述 JSON。
测试和验证都在临时目录生成分析；不重写 `runs/`。

## 校准：每一行及限制

下列每行均为 **1,000 次生成**；推断统计以可评估/可报告的拟合数为分母，
扣留数单独记录。通常上限为 6.3784%，下限为 93.6216%；有效样本数
不足 1,000 时按预先定义的 `2*sqrt(.05*.95/n)` 公式调整。
对照项没有被追认成原脚本的合格判据。

### 重复测量事实层回归模拟

种子 `20260927`；空效应，四条件，12 单位；完整 RM 行的残差 SD 为
0.5/1/1.5/2.5。两个 RI 行为等残差 SD=1、15% MCAR，仅缺失模式按脚本
的最小观测数规则重抽。3000 次事实提取的保存值和状态检查全部通过。

| 场景 | 可报告 / 扣留 | 总体拒绝率 | 家族拒绝率 | 同时覆盖率 | 上限 / 下限 | 判定 |
|---|---|---|---|---|---|---|
| RM ANOVA，非球形，单位 SD=2 | 1000 / 0 | 3.7% | 4.2% | 95.8% | 6.3784% / 93.6216% | 通过 |
| RI/Satterthwaite，单位 SD=2 | 1000 / 0 | 5.0% | 4.5% | 95.5% | 6.3784% / 93.6216% | 通过 |
| RI/Satterthwaite，单位 SD=0.3 | 664 / 336 | 6.0241% | 4.8193% | 95.1807% | 6.6916% / 93.3084% | 对可报告拟合通过；33.6% 扣留 |

脚本初次启动因完整数据配置误写 `missing_policy="error"` 被配置校验拒绝，
尚未生成任何模拟数据。改为既有合法值 `require_complete` 后运行；种子、
场景、样本量和界限未改。记录见
[interpretation_calibration_0.8.1.json](interpretation_calibration_0.8.1.json)。
这不是单因素 bootstrap B=999 再校准。

### 生存：原场景重跑

种子 `20261001`；指数事件时间和独立均匀删失，置换使用 999 次抽样。

| 场景 | 可评估数 | 主要指标 | 判定与对照项 |
|---|---|---|---|
| log-rank，n=6/组 | 1000 | 置换拒绝 4.9% | 通过；渐近对照 6.6%，高于 6.3784% 参考上限 |
| log-rank，n=8/组 | 1000 | 置换拒绝 4.3% | 通过；渐近对照 5.6% |
| log-rank，n=10/组 | 1000 | 置换拒绝 4.5% | 通过；渐近对照 5.1% |
| KM 真中位时点，n=10 | 997 | log-log 覆盖 96.3892% | 通过，下限 93.6195%；log 对照 91.8756%，低于该参考下限 |
| KM 真中位时点，n=20 | 1000 | log-log 覆盖 96.3% | 通过；log 对照 93.3%，低于 93.6216% 参考下限 |
| Cox HR=0.5，n=10/组 | 996 | 覆盖 96.8876% | 原脚本仅报告，未指定通过界限 |
| Cox HR=0.5，n=20/组 | 1000 | 覆盖 95.1% | 通过，下限 93.6216% |
| PH 全局检验，n=20/组 | 1000 | 拒绝 6.2% | 通过，上限 6.3784% |

记录见 [survival_calibration_0.8.1.json](survival_calibration_0.8.1.json)。
原有偏小样本渐近检验和 log 区间的不足仍存在，本版未修正其数值算法。

### 肿瘤：原场景重跑

种子 `20261003`；3 组 × 10 动物，对数线性生长、Gaussian 随机截距与
斜率，按已观测体积阈值移出（生成机制 MAR）。分析日 21。

| 场景 | 模型可报告数 | 移出均数 | 速率差家族覆盖 | 其他指标 | 判定 |
|---|---|---|---|---|---|
| A：空效应，阈值 2000 | 1000 | 0 | 94.2% | 总体拒绝 5.8% | 两个指定界限通过；**没有检验到零假设下有脱落的情形** |
| B：不同速率，阈值 1000 | 1000 | 2.116 | 95.1% | 观测 TGI 家族覆盖 92.2% | 模型通过；TGI 为对照项，低于 93.6216% 参考下限 |
| C：不同速率，无脱落 | 1000 | 0 | 94.5% | 观测 TGI 家族覆盖 95.4% | 模型及 TGI 指定界限通过 |

两项 TGI 行均有 1000 个可评估区间。记录见
[tumor_growth_calibration_0.8.1.json](tumor_growth_calibration_0.8.1.json)。
不能推广到 MNAR、cage、非指数生长或零假设下的信息性脱落。

## 测试、打包与安装

- 完整 `.venv/bin/python -m pytest -q`：**250 passed**，其中新事实层测试
  **24 项**；覆盖来源、子模型状态、边界、零事件、降秩、脱落、区间缺失、
  未来诊断、原输入不变、重绘不变和篡改拒绝。
- 四个变更 Skill（总入口、重复测量、生存、肿瘤）的 `quick_validate.py`
  通过。
- `uv build` 成功；wheel 内 **35 个 Python 文件**与最终源码逐字节一致。
  Wheel SHA-256：`044d463577f6cf65b7cb7357c56f3c6080ee6c7dfd57f44dfdb43da05a09d29a`。
- `validate_clean_environment.py` 成功；其中 **11 套配置**的事实文件
  与源码运行完全一致。`validate_install_smoke.py` 在该隔离环境中通过
  **21 套配置**的分析、重绘与哈希验证，以及板图导入；重绘保留结果和事实文件。
  沙箱中的 Matplotlib 使用临时缓存，所有断言通过。
- **174 个历史验证文件**和 `runs/` 的已跟踪文件与基线一致；原有两个
  未跟踪运行目录保留，未纳入提交。验证索引的历史条目文本保持不变。
- 隔离安装、wheel 哈希、命令及日志哈希最终记录在
  [release_checks_0.8.1.json](release_checks_0.8.1.json)；详细 wheel/source
  对照在 [clean_environment_0.8.1.json](clean_environment_0.8.1.json)。
- 无 R 运行时桥接，因此没有新增可选 R 安装流程；R 仅用于前述验证。
  既有 wheel 保留，未推送、发布或运行 `build_public_tree.py`。

### 可复现命令

```sh
.venv/bin/python scripts/validate_interpretation_benchmarks.py
.venv/bin/python scripts/validate_interpretation_calibration.py
.venv/bin/python scripts/validate_repeated_benchmarks.py --output validation/repeated_benchmarks_0.8.1.json
.venv/bin/python scripts/validate_survival_benchmarks.py --output validation/survival_benchmarks_0.8.1.json
.venv/bin/python scripts/validate_tumor_growth_benchmarks.py --output validation/tumor_growth_benchmarks_0.8.1.json
.venv/bin/python scripts/validate_survival_calibration.py --output validation/survival_calibration_0.8.1.json
.venv/bin/python scripts/validate_tumor_growth_calibration.py --output validation/tumor_growth_calibration_0.8.1.json
.venv/bin/python -m pytest -q
uv build
.venv/bin/python scripts/validate_clean_environment.py
# 用上一步生成的隔离环境的 Python：
<isolated-venv>/bin/python scripts/validate_install_smoke.py
```

上述两个 interpretation 脚本默认写入本版文件名；发布后的复核应在独立
检出副本中重跑并另存新证据，勿覆盖已提交记录。支持 `--output` 的脚本
直接指定新的记录文件。若重建同版本 wheel，先将旧
wheel 移到 `dist/superseded/`。本轮完整流程中的已知不足，尤其 live-agent
F1 的运行时定位错误、未校准的 B=999、未实现的后续方法，均没有被记作通过。
