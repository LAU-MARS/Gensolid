# SMLib 核心调用路径的 Rust 迁移

本文记录已通过 PR #3 合入 `main` 的首轮实现，合并提交为 `a90e4a2`。下文“本轮”均指首轮；最新能力及第二轮验证见[曲面求交、曲线路径扫掠与非平面偏置](rust-migration-round2.md)。

2026-10-23 开始分阶段用 Gensolid 自有 Rust 算法替代目前借助 SMLib 的核心几何服务。本轮提供无需 C++ 的 `RustGeometryBackend`，通过原有 `GeometryBackend` → `Session` → 建模检查器进入真实 B-rep 操作。默认构建不需要 `SMLIB_SOURCE_DIR`、原生 SMLib 静态库或浏览器 companion Wasm。

这是一轮核心路径迁移，不是把约 53.4 万行上游 C++ 全部翻译完成，也不是宣称与整个 SMLib facade 功能等价。原来的可选 `smlib` feature、原生后端及其测试完整保留，供尚未迁移的能力使用和后续差分验证；本轮没有把未实现返回值计为成功。

## 本轮分工与能力边界

| Agent | 本轮工作 | 可用范围 |
|---|---|---|
| Bezira / 贝齐拉 | 有理曲线核心 | 线段、XY 圆/椭圆/圆弧；解析二阶导数；clamped NURBS 插结、裁剪、升阶 |
| Gaurion / 高里昂 | 自有 Rust 曲面后端 | 有理挤出、兼容截面放样/直纹面、直线路径扫掠、任意轴旋转面、受限双线性平面偏置 |
| Noetra / 诺特拉 | 独立审核与内核验证 | 原 SMLib 三个建模场景的等价测试、随机尺寸/坐标框架、退化输入拒绝、失败回滚、源实体保持不变 |
| Hilora / 希洛拉 | 求交与集成 | 多项式单段 Bézier 挤出面 × 仿射平面完整交支；XYZ 与两侧 UV 成对返回；原有拓扑证明保持不变 |

名字表示职责；本轮由当前会话直接调度三个子 agent 并由主 agent 集成，不将其冒称为无人值守 RSI 控制器的一次自动闭环。

### 首轮 Rust 实现明确拒绝的范围

- 任意曲线路径扫掠、非平面法向偏置、通用曲面/曲面求交，以及需要裁剪的局部交支。
- 不兼容的放样截面需要先完成几何兼容化，不能用重采样隐藏误差。
- 求交候选只覆盖受限多项式挤出；有理挤出构造可用，不代表有理挤出的所有相交情况已实现。
- 库中的曲面构造能力不等于建模层能为所有这些曲面建立实体。现有 `Session` 的 regular graph、闭合性、容差与 `check_body` 限制仍然生效；例如完整旋转面可能超出当前 sheet/prism 建模范围。
- 新曲线编辑会保留必要或额外的内部结点，不把结点删除/最简表示当作已完成。

## 使用与复验

Rust 调用方可从 `gensolid_api` 导入 `RustGeometryBackend` 和 `SurfaceConstruction`，传给已有的 `body_create_sheet_backend`、`body_create_solid_backend_prism`、`body_split_face_by_sheet`、`body_offset_sheet_backend`。所有结果仍由原检查器验证，并使用 Session 的事务边界。

```sh
cargo run -p gensolid-api --example rust_geometry
cargo test -p gensolid-api --test rust_geometry_backend
cargo test -p gensolid-geom
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
cargo check -p gensolid-api --target wasm32-unknown-unknown
```

这些命令没有打开 `smlib` feature。原 SMLib 原生测试须在其支持的 C++ 构建环境单独执行；本轮 Windows 的纯 Rust 测试通过不能代替那套原生差分测试。Wasm 的 Rust 编译检查也不等于新能力已经发布到 npm 或完成浏览器端验收。

### 本轮已取得的验证结果

| 验证 | 结果 |
|---|---|
| 集成生产提交 `baa32a5` 的完整 `cargo test` | 255 passed，0 failed，0 ignored |
| 仅新增测试的 `c5dc5d3` 独立 API 验证 | 7 passed，0 failed，0 ignored；含 72 个随机案例 |
| 新曲线/曲面性质测试 | 另外覆盖 128 / 192 个随机输入 |
| 格式与 Clippy | 全工作区生产代码通过；最终新增 API 测试全目标复验通过 |
| 真实纯 Rust 示例 | sheet / solid 都是 0 faults，solid 为 6 faces |
| Wasm 目标 | `gensolid-api` 在 `wasm32-unknown-unknown` 编译检查通过 |
| 冻结范围审计 | 检查器、容差、verifier、原 SMLib 代码与测试均未修改 |

以上是 255 项完整回归与 7 项新增独立测试的组合验收，共 262 项不同测试，不把未开启的原 SMLib 原生用例计入。新增性质测试共 392 个随机案例，不把它们重复计算为 392 个测试函数。

独立审查实际修复了大起始角的圆弧相位问题、区间运算遗漏的舍入误差，以及大偏移 UV 域的控制点/内部求值量化问题。后者曾被外部客户端复现为报告偏差约 `5.45e-15`、实际配对残差达 `0.125`；最终同一输入明确返回 `Unimplemented`。误差检查覆盖两侧 UV，原有建模层的完整交支证明仍然保留。

新代码基于 NURBS 的齐次坐标、张量积、仿射投影等数学公式和项目已有 Rust 实现原创编写，没有搬运 OCCT 源代码，也没有声称逐行翻译 SMLib。原 SMLib Apache-2.0 来源、固定 revision 和 NOTICE 继续保留在 [SMLib 接入文档](smlib.md) 中。冻结的检查器、容差政策和既有测试断言没有修改。

## 后续全 Rust 路线

1. 继续迁移实际调用的几何服务：曲面解析二阶偏导、行列插结/升阶、曲线兼容化、投影和受限精确求交；每项都需原有 fixture 和独立几何契约。
2. 补齐通用算法及明确的完整性证明：通用曲面相交、扫掠框架、非平面偏置；数值候选不能直接当作拓扑结论。
3. 为公开 C/Python/Wasm 接口接入 Rust 能力，并在可运行 SMLib 的环境完成逐能力差分验证。旧后端仍可显式选择，已迁移路径无需依赖它。
4. 只有所有生产调用及验收都由 Rust 覆盖后，才删除 C++ adapter、外部构建脚本和 companion Wasm。移除依赖由证据决定，不以改名或禁用测试代替迁移。
