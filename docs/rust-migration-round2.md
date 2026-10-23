# Rust 迁移第二轮：求交、曲线路径扫掠和非平面偏置

本轮从 PR #3 的主线合并提交 `a90e4a2` 开始，按实际 `GeometryBackend → Session → check_body` 调用路径扩展 Rust 能力。三个 agent 在独立工作树并行实现，互相交叉审查，由 Hilora 集成。原检查器、容差政策、verifier 和既有测试断言保持不变。

## 分工与实际新增能力

| Agent / 分支 | 实现 | 可进入的真实建模流程 |
|---|---|---|
| Eulon / `eulon/rust-spline-intersections-20261023` | 有理及连续多段 NURBS 挤出面 × 有界仿射平面；每段返回 XYZ 和两侧 UV | 有理单段／多段挤出面的等参切分，原独立完整交支证明继续生效 |
| Bezira / `bezira/rust-curved-path-sweep-20261023` | `SurfaceConstruction::TranslationSweep { profile, path }`，固定截面方向沿空间 NURBS 路径平移 | 满足 regular graph 条件时创建 checked sheet 和 prism solid |
| Gaurion / `gaurion/rust-cylindrical-offset-20261023` | 单段有理二次圆弧沿直线挤出形成的圆柱面法向偏置 | 独立认证实际 provider 输出后创建 checked sheet，支持正、负及零距离 |
| Hilora / `hilora/rust-advanced-surfaces-20261023` | 默认后端接线、联合示例、整轮验收与文档 | 圆柱面 → 偏置 → 有理求交切分的连续操作，以及空间扫掠曲面／实体 |

这是当前会话直接调度的 agent 实现与提交过程，不将其描述为无人值守 RSI 控制器已经完成自动闭环。

## 数学契约和拒绝边界

### 曲线路径扫掠

新增显式平移扫掠的定义是 `S(u,v) = C(u) + P(v) − P(v_min)`。张量积控制点为 `C_i + P_j − P(v_min)`，权重为分别归一化后的 `w_i q_j`。非 clamped 路径使用自然参数域起点的实际求值，不以首个控制点代替。相对坐标减小路径公共大平移造成的消减误差；不可表示或次正规数的权重乘积、控制点组合舍入过大及整个控制包坐标间距无法满足精度的输入明确拒绝。

截面方向保持固定；这不涵盖 Frenet／平行移动框架、拐角连接、封闭路径接缝或 SMLib 任意扫掠语义。旧 `SurfaceConstruction::Sweep` 的直线路径契约保持原样。几何公式能构造曲面，并不意味着所有曲面都能成为有效 B-rep；奇异、折叠、不连续及低于会话精度的建模请求仍须拒绝并回滚。

### 有理和多段求交

支持正权重、clamped、连续的挤出截面，挤出方向的成对权重必须相等。按非空结点区间提取 Bézier 段，保持原参数域和每段独立交支。两种输入顺序及 U/V 挤出方向均可用。只接受完整且严格位于源曲面及有界工具内部的交支；部分裁剪、相切、重叠和无法确定的条件返回明确错误。

有理曲线的参数恒等式 `u(t)=t` 不能靠把 Greville 坐标填进有理控制点实现；这里构造齐次乘积 `t W(t)`，把对应分子及分母升到 `p+1` 次。权重条件数、UV 存储及参数量化参与拒绝判断。`reported_deviation` 只是候选诊断值，拓扑接受仍由建模层独立证明。

旧的公开 `intersect_extrusion_plane` 函数保留其多项式单段边界；默认 Rust 后端改用新增的 `intersect_spline_extrusion_plane`。有理斜切的几何候选可计算，但现有建模层尚不能证明通用有理斜切拓扑，因此 `Session` 明确拒绝该类操作；不会把“能算候选”写成“能完成切体”。通用 NURBS 曲面／曲面求交仍未实现。

### 非平面偏置

首个支持族是次数 `(2,1)`、控制网 `(3,2)`、单个 clamped 跨度的圆弧直线挤出面，圆弧不超过四分之一圈，沿轴权重重复。距离符号跟随 `du × dv`。半径塌陷或越轴翻转必须拒绝；球面、圆锥面、一般自由曲面及多段圆柱面偏置仍未覆盖。

接受条件独立检查 provider 真正返回的曲面。由原始坐标开始进行向外舍入的区间运算，使用整个 Bernstein 控制包约束位移长度、位移与两条切向的正交误差，以及源／输出的法向取向。正则性的严格下界把切向残差转为位置误差，误差须低于 `default_local − solver`。不是用角点或中心点采样证明整面正确，也不是仅将输出与同一个生成函数再比较一遍。错距离、内部扭曲、修改有理基底和越轴输出均有拒绝测试。

对支持的圆柱面还新增了更紧的 regular graph 证明：在原中点坐标框架和归一化参数下，用 Bernstein 齐次分子直接证明原有 `j00 > linear`、`j11 > linear`、`j00·j11 − off² > linear²`。原尺寸限制、连续性检查、通用区间细分路径及 4096 次预算保留，检查器和精度阈值没有放宽。

### 交叉审查实际找到的问题

Eulon 用独立客户端发现扫掠在 `2^50` 量级的截面坐标下会返回 `0.25` 的位置偏差；仅检查控制点舍入不足以保护内部求值。新增整个控制包的坐标间距保护后，同一反例明确返回错误。Gaurion 的随机测试还发现半径 `0.2`、圆弧角 `0.3`、高度 `0.3` 的旋转圆柱面在偏置前就耗尽原 graph 证明预算；上述等价 Bernstein 证明解决了该输入，失败种子与直接回归保留。Bezira 独立检查了偏置及新增 graph 证明，Gaurion 独立检查了有理求交的齐次 UV 和空交支判据。

## 复验

以下命令均可在不启用 `smlib` 的环境执行：

```sh
cargo run -p gensolid-api --example rust_advanced_geometry
cargo test -p gensolid-geom --test translation_sweep
cargo test -p gensolid-geom --test rust_spline_surface_intersection
cargo test -p gensolid-api --test translation_sweep
cargo test -p gensolid-api --test rust_rational_intersection
cargo test -p gensolid-api --test rust_cylindrical_offset
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
cargo check -p gensolid-api --target wasm32-unknown-unknown
```

可选 facade 的条件编译另用 `cargo check -p gensolid-ffi --target wasm32-unknown-unknown --features smlib` 检查；这只检查 Rust 接口，不链接或运行原生 SMLib。

本轮在 Windows 使用单编译任务、单测试线程、关闭增量和调试信息控制内存峰值。并行 agent 数量为 3，编译任务数每个工作树为 1；agent 并行不等于无限增加同机编译进程。

### 集成验收结果（2026-10-23）

| 检查 | 实际结果 |
|---|---|
| 完整 `cargo test` | **287 passed，0 failed，0 ignored** |
| 本轮新增测试 | 25 项：求交 10、扫掠 9、圆柱偏置 6；含 224 个配置的随机案例及保存的失败回归 |
| `cargo fmt --all --check` | 通过 |
| `cargo clippy --workspace --all-targets -- -D warnings` | 通过，包含 bindings 和新示例 |
| `gensolid-api` 的 `wasm32-unknown-unknown` 编译检查 | 通过 |
| `gensolid-ffi` 的 Wasm + `smlib` feature 编译检查 | 通过，新增枚举分支在可选 facade 中可编译 |
| 冻结范围审计 | `AGENTS.md`、`verifier/`、检查器、容差文件、`.github/` 均未修改；既有测试断言未修改 |

联合示例 `cargo run -p gensolid-api --example rust_advanced_geometry` 实际输出：

```text
Fixed-orientation spatial sweep sheet: 1 faces, 4 edges, 4 vertices, 0 faults
Spatial sweep prism: 6 faces, 12 edges, 8 vertices, 0 faults
Normal offset of cylindrical sheet: 1 faces, 4 edges, 4 vertices, 0 faults
Rational cylinder / plane cut: 2 faces, 7 edges, 6 vertices, 0 faults
Original cylindrical sheet retained: 1 faces, 4 edges, 4 vertices, 0 faults
Original offset sheet retained: 1 faces, 4 edges, 4 vertices, 0 faults
```

切分直接作用于新偏置结果，并在完成后重新检查两个源曲面，验证功能组合和源实体保持有效。

本地纯 Rust 验证不代替原生 SMLib 差分测试，Wasm 编译检查也不代表新能力已发布为浏览器／Python／C 的完整公开接口。新增 `TranslationSweep` 是 Rust 枚举变体，下游穷尽匹配该枚举时需要增加分支；C ABI 本轮未改。

## 下一步

后续优先把已具备几何候选的有理斜切接入独立拓扑证明，再扩展连续多段圆柱偏置、曲面导数及兼容化。移动框架扫掠和一般自由曲面偏置需单独定义误差和奇异性契约。每一族均先覆盖实际调用和反例、通过整轮验证，再提交推送；不以删除 C++、放宽检查器或把 `Unimplemented` 当成功来计算迁移进度。
