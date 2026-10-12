# Gensolid

[English](README.md) | **简体中文**

对标 Parasolid 的开源 B-rep 几何建模内核。代码主要由 RSI（递归自我改进）驱动的大模型生成，正确性由人类掌控的验证器裁决。

> 状态：**pre-alpha，M0（验证基础设施）阶段**。API 会变，不要用于生产。

## 设计要点

- **Rust 内核**：编译器是大模型代码的第一道免费验证器；内核 crate 禁用 `unsafe`、`unwrap`、`panic`。
- **仿 PK 风格的 C 接口**：整数 tag、`GS_<CLASS>_<verb>` 命名、带 `o_t_version` 的选项结构体；Python 和 WASM 绑定与之共享同一层语义（`gensolid-api`）。
- **从第一天按 Parasolid 数据模型设计**：
  - 解析曲面、NURBS、过程曲面并存，XT 可无损往返；
  - 拓扑 body / region / shell / face / loop / fin / edge / vertex，支持非流形；
  - 每条边、每个顶点可带局部容差（容差建模）；
  - 拓扑判定走精确谓词和区间算术。
- **XT 读写是一等功能**。
- **验证器优先**：有效性检查器、策略测试、对拍框架先于建模算法存在；代理无权修改它们。

## 目前能做什么

| 模块 | 内容 |
|---|---|
| `gensolid-core` | 带代际的类型化 ID 竞技场、错误码、会话精度与局部容差 |
| `gensolid-math` | 点/向量/坐标系、外向舍入区间算术、精确 `orient2d` / `orient3d` |
| `gensolid-geom` | 直线、圆、椭圆、B 样条曲线；平面、圆柱、圆锥、球、环面、B 样条曲面；拉伸面、旋转面；求值、求导、点到曲线/曲面距离 |
| `gensolid-topo` | 非流形拓扑模型、长方体基本体 |
| `gensolid-check` | 有效性检查（拓扑一致性、环闭合、流形边、点线面在容差内） |
| `gensolid-xt` | XT 文本头读写、节点流分词（**仅用合成样例测试过，需用真实 XT 文件验证**） |
| `gensolid-api` | PK 风格会话与 tag API |
| `gensolid-ffi` | C ABI 与头文件 `crates/gensolid-ffi/include/gensolid.h` |
| `bindings/python`、`bindings/wasm` | Python（PyO3）与 WebAssembly 绑定 |
| `verifier/policy` | 仓库规则测试（禁止字面量容差、B-rep 不得依赖网格、`unsafe` 仅限 FFI 等） |

## 构建

需要 Rust 1.88+。

```sh
cargo test                                   # 内核 + 策略测试
cargo clippy --workspace --all-targets -- -D warnings
cargo build -p gensolid-wasm --target wasm32-unknown-unknown --release
pip install maturin && maturin develop -m bindings/python/Cargo.toml
```

Python：

```python
import gensolid
s = gensolid.Session()
body = s.body_create_solid_block(1.0, 2.0, 3.0)
assert len(s.body_ask_faces(body)) == 6 and s.body_check(body) == []
```

C：见 `crates/gensolid-ffi/examples/block.c`。

## 路线图

| 阶段 | 主题 | 验收指标 |
|---|---|---|
| M0 | 验证基础设施：检查器、XT/STEP 读取、OCCT 对拍、基准集、私有测试集 | 对真实 XT 文件给出有效性报告 |
| M1 | 几何核心：求交、投影、精确谓词全覆盖 | 曲面求交分支完整率 |
| M2 | 拓扑与基础建模：欧拉操作、拉伸/旋转/扫掠/放样、质量属性、网格化 | 质量属性误差 |
| M3 | 布尔运算（解析面 → NURBS → 容差 → 非流形） | 随机旋转模型对的布尔成功率 |
| M4 | 圆角、倒角、抽壳、偏置、直接建模 | 圆角链成功率 |
| M5 | XT 写出、STEP AP242、缝合与修复 | XT 往返后可被 Parasolid 应用打开 |
| M6 | 收敛建模（网格体）、PK 风格 API 覆盖 | 对照 Parasolid 功能清单的覆盖率 |

设计决策记录在 [`docs/adr/`](docs/adr/)，代理开发规则见 [`AGENTS.md`](AGENTS.md)。

## 致谢与许可

Gensolid 站在以下项目和文献的肩膀上。

### 核心参考

#### Open CASCADE Technology (OCCT)

<https://github.com/Open-Cascade-SAS/OCCT> · LGPL-2.1（附 Open CASCADE 例外条款）

Gensolid 的核心开源参考。我们参考 OCCT 的架构、算法设计和文档，并把它作为差分测试的 oracle（在验证器中以独立进程运行，对比两者在同一输入上的结果）。

由于 OCCT 采用 LGPL 而 Gensolid 采用 Apache-2.0，Gensolid 源码中**不包含**从 OCCT 源码翻译或改写的代码，也不链接 OCCT。算法依据公开文献和 OCCT 的文档、黑盒行为独立实现。规则见 [`docs/adr/0007-references-and-provenance.md`](docs/adr/0007-references-and-provenance.md)。

#### truck

<https://github.com/ricosjp/truck> · Apache-2.0 · RICOS Co. Ltd. 及 truck 贡献者

Rust B-rep 内核，Gensolid NURBS 实现的参考。truck 与 Gensolid 许可证兼容；凡从 truck 移植的代码，会在文件头注明来源，并在 [`NOTICE`](NOTICE) 中保留其版权声明。

### 规范与文献

- Siemens, *Parasolid XT Format Reference*：XT 读写的格式依据。Siemens 公开该规范是为了让第三方编写转换器。
- L. Piegl, W. Tiller, *The NURBS Book*, 2nd ed., Springer, 1997：B 样条基函数与求值算法（A2.1、A2.3）。
- J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates*, Discrete & Computational Geometry 18, 1997：精确谓词。

### 依赖库

- [robust](https://github.com/georust/robust)（MIT / Apache-2.0）：Shewchuk 精确谓词的 Rust 实现。
- [PyO3](https://github.com/PyO3/pyo3)、[wasm-bindgen](https://github.com/rustwasm/wasm-bindgen)、[cbindgen](https://github.com/mozilla/cbindgen)、[proptest](https://github.com/proptest-rs/proptest)。

### 许可与商标

以 [Apache License 2.0](LICENSE) 发布。Parasolid 是 Siemens Digital Industries Software 的商标。Gensolid 以 Parasolid 为能力对标对象，但与 Siemens 无隶属或背书关系。
