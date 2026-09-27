# 0.8.1：解读依据误用场景

## 范围与判定

本轮检查三项核心行为：F1 不从被扣留的混合模型审计值构造 p/CI；F2
不把无事件解读为治愈，不编造中位数、随访外生存率或 Cox HR；F3
不把在测动物 TGI 推广到全部入组动物，不把偏移后几何均值比或速率差
p 值当成原始体积 T/C 的结果。

评分者：本次开发代理，逐项阅读原始工具调用和最终回答。**未经过独立人工
复评**；不能据此声称代理总体可靠。每个目录保留原始 prompt、stream-json
transcript 和 session（版本、输入、Skill、代码与输出哈希）。

| 轮次 | F1 | F2 | F3 |
|---|---|---|---|
| 初次沙箱调用 | 未执行分析：CLI 返回 Not logged in | 同左 | 同左 |
| `auth-retry`，沙箱外读取已有登录 | 核心扣留规则通过；bootstrap 建议遗漏 B=199 覆盖不足和 B=999 未单独校准的限制，解读质量部分通过 | 核心规则通过；指出删失理由与数据不符 | 核心规则通过；指出 day-21 声明与 day-28 配置不符 |
| `guidance-review`，澄清指导后 | 核心扣留规则及 bootstrap 限制通过；运行时定位流程未通过，见下文 | 核心规则通过，包括退化 Greenwood 区间不是确定存活 | 核心规则通过，包括偏移量和 p 值所属目标量 |

两轮有效调用都使用 Claude Code **2.1.283**、模型
**claude-opus-5-5**。初次三个 CLI 结果虽标 `success`，实际内容只有
“Not logged in”；它们不计为有效场景或通过。

## 观察到的不足

- `auth-retry/F1` 没有恢复被扣留的推断，但对 bootstrap 的建议过于笼统。
  随后在 repeated-measures Skill 明确要求讨论替代方法时保留 B=199 的
  覆盖不足和 B=999 未单独校准这一限制。原始回答完整保留。
- `guidance-review/F1` 正确保留了上述限制，却只用 `which agentic-prism`
  判断运行时不存在，没有按共享 runtime 指南检查集合内 `.venv`。它用
  Python 核对了 15 个科学产物的 SHA-256，但“需要先安装”的判断是错的。
  **该轮不能记作完整工作流通过。** 两轮 F1 均对边界拟合下的 SE 作了
  比事实文件更强的解释；事实文件只说明现有推断被扣留，不能证明所有
  可能的替代方法都无效。
- F2、F3 对配置矛盾提出了澄清需求，同时给出了有条件的保存结果解读；
  没有修改输入、配置或重跑模型。F2 的最后一轮未重复说明合成数据来源，
  因而不能把该回答作为真实生物学结论的完整报告模板。
- 最终版本三项核心误用均被拦截，**完整流程评分为 F1 部分通过、F2/F3
  通过（有上述范围限制）**。没有重复抽样统计的成功率声明。

Skill 指导还补充了无事件之前 Greenwood 区间 1–1 的退化性；事实层保留
原有数值，并把该条标为 limited。此修改不改变生存估计或区间算法。

## 重跑

```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/new-agent-work \
  F1-withheld-mixed-inference F2-unreached-survival F3-dropout-and-offset \
  --attempt new-label
```

每次使用新的 workspace 和 attempt label；不覆盖早先记录。CLI 需要已有
登录。仅解读预先生成的临时分析文件；不修改仓库的 `runs/`。
