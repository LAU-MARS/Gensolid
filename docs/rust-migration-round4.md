# Rust 迁移第四轮：多段有理斜切

PR #5 的 Linux、Windows、macOS、C 接口、Wasm、Python 和依赖许可七项 CI 全部成功后，已正常合入 `main`，合并提交为 `618efeb`。本轮从该提交继续，将有理斜切从单段扩展到连续多段截面的直线挤出面。

| Agent | 分支 | 职责 |
|---|---|---|
| Eulon / 欧纶 | `eulon/rust-multispan-oblique-cuts-20261023` | 原始 NURBS 的区间抽取、多段交线构造与全域证明 |
| Gaurion / 高里昂 | `gaurion/rust-multispan-cut-topology-20261023` | 公开 API、真实共享拓扑、正反向实体面的回归 |
| Noetra / 诺特拉 | `noetra/rust-multispan-cut-validation-20261023` | 独立几何、跨缝、错误候选与事务回滚测试草稿 |
| Hilora / 希洛拉 | `hilora/rust-multispan-oblique-cuts-20261023` | 集成审查、接手独立验证、联合示例和完整验收 |

Noetra 在保存独立测试后遇到账户额度限制，由 Hilora 接手审查、执行和收尾。该过程由当前会话调度，不代表无人值守 RSI 控制器已经完成自动闭环。

## 使用入口与范围

新增 `Session::body_split_face_by_sheet_multispan(body, face, tool, backend)`，返回 `FaceSplit`。目标仍须是经过检查的完整 NURBS 矩形面，工具须是有界凸仿射平面，交支必须完整横截、严格处于内部并可建立空间误差界。新增有理斜切路径支持次数不超过 16、2–256 个非空跨度的连续 clamped 正权重截面；内部结点重数不能超过次数。满足这些表示条件仍须通过完整证明，不能据此保证任意病态输入成功。

原 `body_split_face_by_sheet`、`body_split_face_by_sheet_with_options` 和 `SheetCutOptions` 保留既有能力契约；旧测试中的多段拒绝仍成立。新入口通过新增正向、反向和组合测试验证，不修改既有断言。

源面创建补充了有理直纹面的 Bernstein 区间证明路径，解决部分正常多段曲面在通用区间细分中耗尽预算的问题。它使用原来的整面投影 Frame、全域归一化参数和对称 Jacobian 正定门槛；尺寸、连续性检查及 4096 次通用回退预算保持原样。

实体切分要求源实体已经通过构造和检查。原有实体缝合仍可能拒绝部分带非二进制精确结点的细分截面，错误为 `Unimplemented`，发生在切分之前；本轮保留该限制，并测试构造失败的完整回滚。独立验证另用可创建的真实多段实体检查负向面的切分与邻面边传播。

## 跨段交线与参数

两张结果面保留原曲面 carrier，共享一条跨段三维交线及对应的 fin 参数曲线。内部 NURBS 接缝由曲线内部结点表达，不因增加跨度就增加拓扑边或顶点；实体面的边界端点切分仍传播到邻面。

新的交线采用共同的整数分段参数：每个非空原始跨度依次映射到 `[0,1]`、`[1,2]` 等区间，由 UV 曲线映射回原始非均匀参数域。在支持的段数范围内，整数结点反向及双反向可以精确表示，避免浮点反射导致结点漂移。**交线参数不保证等于原曲面的 U/V 参数**；调用方应通过 fin 的参数曲线获取曲面坐标。原曲面的参数、控制网和权重保留。

证明从原始存储的控制点、权重和结点做向外舍入的齐次区间抽取，检查最终存储后的 XYZ/UV；浮点拆段结果、后端候选及其报告的误差不能代替证明。对每段核对横截条件、参数覆盖、单调性、曲面位置和空间误差，并验证跨段及反向使用时的一致性。

## 复验

Windows 低内存配置与[第三轮](rust-migration-round3.md)相同：一个编译任务、一个测试线程、关闭增量与调试信息、一个 codegen unit。

```sh
cargo run -p gensolid-api --example rust_multispan_cut
cargo test -p gensolid-api --test multispan_cut_topology
cargo test -p gensolid-api --test rust_multispan_cut_independent
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
cargo check -p gensolid-api --target wasm32-unknown-unknown
```

最终集成快照在 Windows 上串行验收通过，测试耗时约 404 秒（本次缓存构建环境）。正式提交逐文件核对 Git blob，与被验证的生产代码、测试和示例一致。

| 检查 | 结果 |
|---|---|
| `cargo fmt --all --check` | 通过 |
| `cargo clippy --workspace --all-targets -- -D warnings` | 通过 |
| `cargo test`，含文档测试 | **334 passed，0 failed，0 ignored** |
| 五段圆柱 → 偏置 → 斜切示例 | **2 个面、7 条边、6 个顶点，0 故障**；保留的源圆柱和偏置体也为 0 故障 |
| `cargo check -p gensolid-api --target wasm32-unknown-unknown` | 通过；这是目标编译检查，不是浏览器运行验证 |
| `git diff --check` | 通过 |

本轮新增 22 项测试，其中两组 property test 各配置 16 个随机样例，另保留触发过图证明预算问题的回归 seed。覆盖不同直纹方向、正反向实体面、非均匀接缝、偏置与切分组合、错误或空后端候选、事务及编号回滚、极端公共权重缩放、近相切与过粗参数域等情况。证明模块还直接验证最终 XYZ/UV 被破坏以及一个 ULP 的接缝漂移均不能混过证书。

联合示例使用具有原双重结点和新增单重结点的五段圆柱，连续执行法向偏置与多段斜切，检查结果及保留的源实体。以上为本地验收；推送后的多平台 CI 以 GitHub PR 的实际检查状态为准。

冻结的检查器、容差、verifier、CI 和既有测试断言均不修改。相切、局部裁剪、多交支和任意自由曲面求交仍属于后续能力；本轮本地验收不包含原生 SMLib 差分或二进制发行。
