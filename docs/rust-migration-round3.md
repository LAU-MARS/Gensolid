# Rust 迁移第三轮：有理斜切与多段圆柱面偏置

本文是第三轮历史记录。PR #5 的七项 CI 全部通过后已合入 `main`，合并提交为 `618efeb`；最新能力见[第四轮：多段有理斜切](rust-migration-round4.md)。

PR #4 的 Linux、Windows、macOS、C 接口、Wasm、Python 和依赖许可七项 CI 全部成功后，通过正常合并流程合入 `main`，合并提交为 `2d7e856`。本轮三个 agent 均从该主线提交开始，分别实现两项能力和独立反例测试，Hilora 负责集成、联合示例与验收。

| Agent | 分支 | 职责 |
|---|---|---|
| Eulon / 欧纶 | `eulon/rust-rational-oblique-cuts-20261023` | 有理斜切的完整域证明、公开选项和真实共享拓扑 |
| Gaurion / 高里昂 | `gaurion/rust-multispan-cylinder-offset-20261023` | 多段圆柱偏置、原始 NURBS 区间抽取与全局 graph 证明 |
| Noetra / 诺特拉 | `noetra/rust-rational-offset-validation-20261023` | 独立几何方程、恶意 provider、缝线、全局折返和事务回滚测试 |
| Hilora / 希洛拉 | `hilora/rust-rational-cuts-multispan-offset-20261023` | 集成审查、联合示例、文档和整轮回归 |

## 有理斜切：从候选到共享拓扑

新增 `SheetCutOptions { allow_rational_oblique: true }`，传给 `Session::body_split_face_by_sheet_with_options(body, face, tool, backend, options)`。本轮支持正权重、次数不超过 16 的单段 clamped Bézier 截面的直线挤出面，被有界凸平面完整横截；需要严格内部交支，源曲面须通过原有 regular graph 验证。结果是保留原 carrier 的两张裁剪面，共享真实三维交线及各自的 UV 曲线；操作失败时不发布部分模型或占用标签。

旧 `body_split_face_by_sheet` 使用默认选项 `allow_rational_oblique: false`，保持其已公开能力契约；原“未实现”测试断言保留，新能力由显式入口和新的正向、反向测试验证。新选项不把有理多段斜切、局部交支、相切或任意曲面裁剪宣称为已实现。

验收从原始控制网、权重和参数域出发，证明完整横截及交支的空间误差。候选中的有理参数恒等式通过齐次乘积构造，使用存储后 UV 的 `U′W − UW′` 的整个 Bernstein 控制包证明单调性，而不是要求 UV 控制点坐标逐个递增。再通过 `(U 升阶) − tW` 约束参数误差，结合曲面导数、平面残差和横截斜率下界得到整支位置误差。UV 域量化及坐标舍入同样参与判定。

外部 backend 的成功、XYZ、两侧 UV 或 `reported_deviation=0` 都不构成证明。新增路径重建候选后独立验证原始曲面关系；错误候选不能把错误拓扑写进模型，空候选仍按原契约明确拒绝。

## 多段圆柱偏置：逐段证明与全局一致性

保留 `body_offset_sheet_backend` 和已有偏置接口。扩展的对象是同轴、同半径、正权重的圆弧直线挤出面：次数 `(2,1)`，V 方向为两列，U 方向 clamped 且连续。内部结点重数 1 或 2 可表示连续细分，重数 3 的不连续输入明确拒绝。偏置距离仍按 `du × dv` 定向，半径塌缩或越轴反转拒绝。

候选对原控制网应用统一的径向仿射变换，保留结点和权重。接受证明对原始 NURBS 做向外舍入的齐次区间抽取，不能把浮点 `trim` 的结果当作无误差的源曲面。每个跨度验证实际 provider 输出的法向距离、切向误差与取向，还必须属于同一圆柱，接缝保持连续。

所有跨度使用整面共同的投影坐标框架，以及统一归一化参数下的 Jacobian。局部跨度各自可投影，不足以证明整面可建模；折返、完整闭合圆周或不能通过该全局证明的输入仍拒绝。原尺寸、连续性、容差、4096 次通用区间细分预算和最终 `check_body` 不放宽。

## 复验与联合示例

本轮在 Windows 使用以下低内存配置；设置后再执行下方命令。其他平台可设置同名环境变量。

```powershell
$env:CARGO_BUILD_JOBS = '1'
$env:CARGO_INCREMENTAL = '0'
$env:CARGO_PROFILE_DEV_DEBUG = '0'
$env:CARGO_PROFILE_TEST_DEBUG = '0'
$env:CARGO_PROFILE_DEV_CODEGEN_UNITS = '1'
$env:CARGO_PROFILE_TEST_CODEGEN_UNITS = '1'
$env:RUST_TEST_THREADS = '1'
```

```sh
cargo run -p gensolid-api --example rust_rational_modeling
cargo test -p gensolid-api --test rational_oblique_cut
cargo test -p gensolid-api --test rust_multispan_cylindrical_offset
cargo test -p gensolid-api --test rust_offset_cut_independent
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
cargo check -p gensolid-api --target wasm32-unknown-unknown
```

联合示例连续执行两条实际建模路径：单段圆柱 → 偏置 → 有理斜切，以及多段圆柱 → 偏置 → 等参切分，并重新检查源实体。多段有理斜切仍在后续范围，不能将第二条路径的等参切分写成多段斜切成功。

| 验收项 | 最终集成结果 |
|---|---|
| 完整 `cargo test` | **312 passed，0 failed，0 ignored**；48 个测试结果组，本机完整运行约 396 秒 |
| 本轮新增测试 | 25 项：有理斜切 6、多段偏置 10、独立组合与反例 9；其中包含每轮 96 个配置的随机案例 |
| 格式与静态检查 | `cargo fmt --all --check`、整个工作区及所有 target 的 Clippy `-D warnings` 均通过 |
| 联合示例 | 两条切分路径均输出 2 面、7 边、6 顶点；结果及保留的源实体全部 0 fault |
| Wasm | `cargo check -p gensolid-api --target wasm32-unknown-unknown` 通过 |
| 独立复验 | Noetra 在正式实现的同一内容上另行完整运行 312 项测试，全部通过 |

本轮使用三个 agent 的独立工作树，每个编译任务数与测试线程数为 1，关闭增量和调试信息控制内存。实现与独立测试按 agent 分别提交，集成审查通过后推送各分支及 Hilora 集成分支。上述结果为本地验收，GitHub 跨平台 CI 状态以本轮 PR 的检查页为准。此过程由当前会话直接调度 agent 完成，不声称无人值守 RSI 控制器已完成自动闭环。

冻结的检查器、容差、verifier、CI 和既有测试断言均保留。原生 SMLib 差分验证不包含在本机纯 Rust 回归中；C/Python/Wasm 二进制发布也不由 Rust API 验收自动完成。后续优先扩展多段有理斜切、一般扫掠框架和自由曲面偏置，每一项继续以可证明的实际建模能力计入迁移进度。
