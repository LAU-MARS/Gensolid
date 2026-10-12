## 变更内容

<!-- 做了什么，为什么。 -->

## 代码来源声明（必填）

- [ ] 本 PR 不包含从 OCCT 源码翻译或改写的代码（只参考了文档、论文或黑盒行为）
- [ ] 若移植了 truck 等 Apache-2.0 项目的代码，已在文件头注明来源并更新 `NOTICE`
- [ ] 未粘贴 Parasolid 文档原文
- 生成方式：<!-- 人工 / 某模型 + 代理名称；RSI 运行编号（如有） -->

## 验证

- [ ] `cargo fmt --all --check`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo test` 通过
- [ ] 新的建模操作输出都通过 `gensolid-check`，并有属性测试
- [ ] 未修改 `verifier/` 或测试断言来让测试通过
