# ADR 0001：内核用 Rust 编写

状态：已采纳（2026-10-12）

## 背景

Gensolid 的代码主要由大模型生成。商业内核（Parasolid、ACIS）和 OCCT 都用 C/C++，但大模型写出的 C++ 中，未定义行为、悬垂指针和数据竞争很难被自动发现，而自动发现错误正是 RSI 循环的前提。

## 决定

- 内核全部用 Rust（edition 2024，最低 1.88）。
- 工作区 lint：`unsafe_code`、`clippy::unwrap_used`、`expect_used`、`panic`、`todo`、`unimplemented` 都是 deny；只有 `gensolid-ffi` 允许 `unsafe`。
- 实体用带代际的类型化 ID（`Id<Face>`、`Id<Edge>`），不用指针，类型系统阻止把面 ID 当边 ID 传。

## 后果

- 编译器和 clippy 成为验证器的第一层，零成本拦截一大类错误。
- 原生支持 WASM，与 truck、vcad、brepkit 等 Rust 生态可以相互借鉴。
- 与 C/C++ 世界的集成通过 C ABI 完成（见 ADR 0002）。
