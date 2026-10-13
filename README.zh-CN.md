# Gensolid

[English](README.md) | **简体中文**

对标 Parasolid 的开源 B-rep 几何建模内核。代码主要由 RSI（递归自我改进）驱动的大模型生成，正确性由人类掌控的验证器裁决。

> 状态：**pre-alpha，已建立 M0 公开开发基准，M1 开发目录 1088/1088，M2 基础建模已启动**。一般曲面求交尚未完整实现；具体证据与范围见 [M0/M1 与 M2 起步状态](docs/m0-m1-status.md)。API 会变，不要用于生产。

## 交互式几何 Viewer

运行 `python tools/viewer.py`，打开 `http://127.0.0.1:8765`。可旋转/缩放、修改参数、查看检查结果并导出 JSON。五个场景直接调用 Rust 内核，无需 npm 或外部服务。

新增凸轮廓离轴旋转、多截面凸多边形实体放样、兼容有理曲线的光滑曲面放样，以及任意次数 NURBS 曲线—曲面的区间根认证。一般曲面—曲面求交提供数值交线和全局可能区域，**分支连接完整性尚未认证**。预览网格没有弦高保证。具体范围见 [新增实现与公开证据](docs/advanced-modeling.md)；`python tools/benchmark-advanced.py` 可复现新增用例。

## 设计要点

- **Rust 内核**：编译器是大模型代码的第一道免费验证器；内核 crate 禁用 `unsafe`、`unwrap`、`panic`。
- **仿 PK 风格的 C 接口**：整数 tag、`GS_<CLASS>_<verb>` 命名、带 `o_t_version` 的选项结构体；Python 和 WASM 绑定与之共享同一层语义（`gensolid-api`）。
- **从第一天按 Parasolid 数据模型设计**：
  - 解析曲面、NURBS、过程曲面并存，已支持的 XT 导入保留几何；
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
| `gensolid-geom` | 解析与 NURBS 求值求导、投影、支持配对的解析求交、一/二/三次有理曲线与平面求交、双线性曲面与平面求交 |
| `gensolid-topo` | 非流形拓扑、凸多面体、多边形拉伸、直线扫掠、相似截面放样、直边分裂 |
| `gensolid-modeling` | 经检查器验证的建模入口；平面 B-rep 的体积、面积、质心、惯量，包含孔和空腔 |
| `gensolid-mesh` | 输出网格：平面面片和具有弦高上界的解析曲面参数区域 |
| `gensolid-check` | 有效性检查（拓扑一致性、环闭合、流形边、点线面在容差内） |
| `gensolid-xt` | 按 schema 读取文本/中性二进制并给出报告；已支持三个真实 ABC XT 件的精确嵌入 schema；明确报告不支持的 schema 与几何 |
| `gensolid-step` | Part21 流形实体导入，解析/NURBS 几何、经认证的凸平面空腔与逐体检查；平面 AP203 写出 |
| `gensolid-cli` | JSON 文件诊断、公开开发基准、OCCT 一致率与覆盖率、ABC 配套格式读取核对、仅接收汇总分数的私有接口 |
| `gensolid-api` | PK 风格会话与 tag API |
| `gensolid-ffi` | C ABI 与头文件 `crates/gensolid-ffi/include/gensolid.h` |
| `bindings/python`、`bindings/wasm` | Python（PyO3）与 WebAssembly 绑定 |
| `verifier/policy` | 仓库规则测试（禁止字面量容差、B-rep 不得依赖网格、`unsafe` 仅限 FFI 等） |
| `verifier/oracle` | 公开长方体语料与 OCCT 进程外对拍；安装 pythonocc 后才会真正比较 |

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
assert s.body_ask_volume(body) == 6.0
```

C：见 `crates/gensolid-ffi/examples/block.c`。

公开开发基准：

```sh
python tools/fetch-abc-smoke.py
cargo run -p gensolid-cli -- check corpus/abc-smoke
cargo run -p gensolid-cli -- benchmark --output target/m0-m1-m2-report.json
```

ABC 示例来自一个官方公开模型。默认基准还包含独立 OCCT 实体与负例；冻结 manifest 保留缺失或被改动输入的诊断。报告将未实现配对计入求交总数，同时输出 OCCT 对比覆盖率与一致率；`benchmark` 成功退出表示报告已生成且输入完整性检查通过，覆盖情况以报告字段为准。预留的 `--private-score <receipt.json>` 只接收 `suite_id`、`total`、`passed`，私有案例留在独立评分服务。

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

- [NVIDIA USD BRep / SMLib](https://github.com/NVIDIA-Omniverse/usd-brep)（Apache-2.0），锁定 `c9c5979f04cea8c38213004d4f04df0dcfb40d05`：感谢 NVIDIA 开放此内核；本分支通过 C ABI 和安全 Rust 接口封装其曲线/曲面能力，作为可选 C++ 依赖（`gensolid-ffi` 的 `smlib` feature）。构建、范围及限制见 [SMLib 集成说明](docs/smlib.md)。
- [robust](https://github.com/georust/robust)（MIT / Apache-2.0）：Shewchuk 精确谓词的 Rust 实现。
- [parasolid-core](https://github.com/monozukuri-ai/parasolid-kit)（MIT 和 Apache-2.0），锁定 0.3.7：XT 节点/schema 解码；拓扑和几何映射由 Gensolid 独立实现。
- [sha2](https://github.com/RustCrypto/hashes)（MIT / Apache-2.0）：CLI 冻结开发语料的完整性检查。
- [PyO3](https://github.com/PyO3/pyo3)、[wasm-bindgen](https://github.com/rustwasm/wasm-bindgen)、[cbindgen](https://github.com/mozilla/cbindgen)、[proptest](https://github.com/proptest-rs/proptest)。

### 许可与商标

以 [Apache License 2.0](LICENSE) 发布。Parasolid 是 Siemens Digital Industries Software 的商标。Gensolid 以 Parasolid 为能力对标对象，但与 Siemens 无隶属或背书关系。

包含 SMLib 的浏览器发布包：配置锁定的源码和 Emscripten 后，运行 `python3 tools/package-smlib-wasm.py`。脚本一起打包 Rust、SMLib 两个 Wasm 模块及 JS 加载器；见 [Wasm 构建与验证](docs/smlib.md#browser-wasm-release)。
