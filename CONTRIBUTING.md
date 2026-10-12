# 参与贡献

欢迎人类和 AI 代理贡献代码。所有贡献都遵守 [`AGENTS.md`](AGENTS.md) 中的规则。

- **来源声明**：每个 PR 都要填写模板中的代码来源声明，并注明生成方式（人工，或模型与代理名称）。
- **签署**：提交需带 `Signed-off-by`（`git commit -s`），表示你有权以 Apache-2.0 提交这些代码（[DCO](https://developercertificate.org/)）。
- **验证器变更**：`verifier/`、`gensolid-check` 和容差模块属于验证层，需要 CODEOWNERS 审核。新增测试用例，尤其是让现有内核失败的用例，是最有价值的贡献之一。
- **设计变更**：影响数据模型、API 或流程的决定，先在 `docs/adr/` 新增一篇 ADR 再动手。
