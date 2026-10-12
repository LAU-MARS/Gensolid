# ADR 0002：仿 PK 风格的 C 接口，绑定共享一层语义

状态：已采纳（2026-10-12）

## 背景

Parasolid 的 PK 接口是业内事实标准，350 多款应用基于它开发。降低这些开发者的迁移成本，是对标 Parasolid 的关键。

## 决定

- `gensolid-api` 是唯一的语义层：一个会话、整数 tag（0 为空 tag，同一会话内不复用）、`<class>_<verb>_<object>` 命名。
- `gensolid-ffi` 把它导出为 C ABI，命名风格对齐 PK：
  - 函数名为 `GS_BODY_create_solid_block` 这种形式；
  - 每个函数返回 `GS_ERROR_code_t`；
  - 选项结构体带 `o_t_version`，用 `*_o_init` 初始化；
  - 内核分配的数组用 `GS_MEMORY_free` 释放；
  - 头文件由 cbindgen 生成，CI 校验它与源码一致。
- Python（PyO3）和 WASM（wasm-bindgen）直接包装 `gensolid-api`，不经过 C ABI，避免两次 FFI。
- 不追求与 PK 二进制兼容，也不复制 PK 的函数清单；只对齐风格和概念。日后若要做 PK 兼容层，可在 C ABI 之上单独实现。

## 后果

- 三种语言绑定行为一致，测试只需在 `gensolid-api` 层写一次。
- C ABI 有全局会话，与 PK 一致；多会话并发留待以后单独设计。
