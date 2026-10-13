# 验证器（L0）

验证器是人类掌控的“尺子”：它决定一段生成代码是否被接受。代理可以读，不能改（见 `.github/CODEOWNERS`）。

| 目录 | 状态 | 内容 |
|---|---|---|
| `policy/` | 已有 | 仓库规则测试，随 `cargo test` 运行 |
| `corpus/` | 已有 `blocks.json`，`abc/` 可放入 STEP | 参数化长方体总会跑。把 ABC 的 `.step`/`.stp` 放进 `corpus/abc/` 后，对拍会把它们算进一致率 |
| `oracle/` | 已有 | 长方体几何对拍；STEP 有效性与 OCCT `BRepCheck` 的一致率；解析求交完整率。未安装 pythonocc 时 OCCT 部分跳过 |
| `private/` | 不入库 | 若存在 `private/blocks.json`，只回传通过率 |
| `bench/` | 计划（M0） | 阶段验收基准：布尔成功率、圆角成功率、XT 往返率、性能 |
| `private/` | 不入库 | 私有保留测试集；只把分数回传给 RSI 循环 |

`private/` 已列入 `.gitignore`。私有集应存放在独立的、代理无权访问的存储中，由评分服务读取。
