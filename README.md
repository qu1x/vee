# vee

Vector Expression Emitter (VEE): Geometric Algebra Code Generator

[![Build][]](https://github.com/qu1x/vee/actions/workflows/build.yml)
[![Documentation][]](https://docs.rs/vee)
[![Downloads][]](https://crates.io/crates/vee)
[![Version][]](https://crates.io/crates/vee)
[![Rust][]](https://www.rust-lang.org)
[![License][]](https://mozilla.org/MPL)

[Build]: https://github.com/qu1x/vee/actions/workflows/build.yml/badge.svg
[Documentation]: https://docs.rs/vee/badge.svg
[Downloads]: https://img.shields.io/crates/d/vee.svg
[Version]: https://img.shields.io/crates/v/vee.svg
[Rust]: https://img.shields.io/badge/rust-v1.85.0-brightgreen.svg
[License]: https://img.shields.io/crates/l/vee

The goal of this crate is to generate optimized code for geometric algebra flavors. Currently, this
crate implements the symbolic reduction of multivector expressions up to polynomials with rational
coefficients.

## Features

  * Uniquely reduce symbolic multivector expressions.
  * Generate expressions in text form.
  * Define the metric-agnostic basis (i.e., elliptic, parabolic, hyperbolic) and the multivector
    entities of the projective geometric algebra (PGA) for dimensions D = N + 1 < 8.

## Roadmap

  * Eliminate orthonormalization conditions form expressions using reflection/projection operator by
    factoring pinned symbols and GCD coefficients in symbolic storage rather than during formatting.
  * Generate expressions in code form.
  * Define other geometric algebra flavors.

See the [release history](RELEASES.md) to keep track of the development.

## Pseudo-local Documentation Builds

Build the documentation with [cargo-tex](cargo-tex). Note that navigating the documentation requires
web access as KaTeX is embedded via remote CDN.

```sh
cargo tex --open
```

## License

Copyright © 2025 Rouven Spreckels <rs@qu1x.dev>

Licensed under the terms of the [`MPL-2.0`](LICENSES-MPL).

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted for inclusion in the
Covered Software by You shall be licensed as above, without any additional terms or conditions.
