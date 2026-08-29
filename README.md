# vee

Vector Expression Emitter (VEE): Geometric Algebra Code Generator

[![Build][]](https://github.com/qu1x/vee/actions/workflows/build.yml)
[![Documentation][]](https://docs.rs/vee)
[![Downloads][]](https://crates.io/crates/vee)
[![Version][]](https://crates.io/crates/vee)
[![Rust][]](https://www.rust-lang.org)
[![License][]](https://mozilla.org/MPL)
[![DOI][]](https://doi.org/10.5281/zenodo.21433093)

[Build]: https://github.com/qu1x/vee/actions/workflows/build.yml/badge.svg
[Documentation]: https://docs.rs/vee/badge.svg
[Downloads]: https://img.shields.io/crates/d/vee.svg
[Version]: https://img.shields.io/crates/v/vee.svg
[Rust]: https://img.shields.io/badge/rust-v1.94.1-brightgreen.svg
[License]: https://img.shields.io/crates/l/vee
[DOI]: https://img.shields.io/badge/DOI-5281/zenodo.21433093-blue.svg

The goal of this crate is to generate optimized code for geometric algebra flavors.

## Features

  * Zero non-optional dependencies.
  * Uniquely reduce symbolic multivector expressions for algebraic and structural equivalence to
    coincide.
  * Generate text form in Unicode/ASCII/LaTeX mode.
  * Generate code form in generic/Rust mode.
  * Generate tree form (i.e., DOT as in [`text/vnd.graphviz`]) in Unicode/ASCII mode.
  * Eliminate orthonormalization conditions from expressions using reflection/projection
    operator by factoring pinned symbols, GCD coefficients, and predominant signs.
  * Evaluate symbols as rationals.
  * Count operations (i.e., multiplications and additions).
  * Define the metric-agnostic basis, i.e., elliptic, hyperbolic, and parabolic (Euclidean)
    along with the multivector entities for dimensions D = N + 1 < 8 of the
    plane-based pistachio flavor, i.e., projective geometric algebra (PGA).

[`text/vnd.graphviz`]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)

## Roadmap

  * Simplify emitter by flattening expression tree into token stream and perform defer logic
    during stream iteration rather than tree traversal.
  * Further optimize expressions to reduce operation count by domain-specific common
    subexpression elimination (CSE) targeting exterior products. Reduce search space by
    leveraging applicable geometric decomposition, e.g., Euclidean decomposition for parabolic
    reflection operator.
  * Generate expressions in SIMD mode.
  * Define other geometric algebra flavors.

See the [release history](RELEASES.md) to keep track of the development.

## Pseudo-Local Documentation Builds

Build the documentation with [cargo-tex](cargo-tex). Note that navigating the documentation requires
web access as KaTeX is embedded via remote CDN.

```sh
cargo tex --open
```

## License

Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>

Licensed under the terms of the [`MPL-2.0`](LICENSES-MPL).

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted for inclusion in the
Covered Software by You shall be licensed as above, without any additional terms or conditions.
