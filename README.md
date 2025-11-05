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

The goal of this crate is to generate optimized code for geometric algebra flavors.

## Features

  * Uniquely reduce symbolic multivector expressions.
  * Generate expressions in text form.
      * Use symbols with Unicode *combining diacritical marks* (`"{}"`).
      * Use Symbols labelled after basis blades (`"{:#}"`).
      * Expand (`"{:+}"`) or reduce (i.e., factorize) expressions (`"{}"`).
      * Factor predominant sign (`"{:-}"`).
      * Omit newlines (`"{:0}"`).
      * Omit plus signs (`"{:<}"`).
  * Generate expressions in DOT form (i.e., [`text/vnd.graphviz`]).
      * Use symbols with Unicode diacritical marks (`"{:o}"`).
      * Use Symbols labelled after basis blades (`"{:#o}"`).
      * Expand (`"{:+o}"`) or reduce (i.e., factorize) expressions (`"{:o}"`).
      * Factor predominant sign (`"{:-o}"`).
      * Change rank direction from top-to-bottom (TB) to left-to-right (LR) (`"{:0o}"`).
  * Eliminate orthonormalization conditions from expressions using reflection/projection operator by
    factoring pinned symbols, GCD coefficients, and predominant signs.
  * Evaluate symbols as rationals.
  * Count operations (i.e., multiplications and additions).
  * Define the metric-agnostic basis (i.e., elliptic, parabolic, hyperbolic) and the multivector
    entities of the projective geometric algebra (PGA) for dimensions D = N + 1 < 8.

[`text/vnd.graphviz`]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)

## Roadmap

  * Further optimize expressions to reduce operation count by domain-specific common subexpression
    elimination (CSE) targeting exterior products.
  * Generate expressions in LaTeX and code form.
  * Define other geometric algebra flavors.

See the [release history](RELEASES.md) to keep track of the development.

## Pseudo-Local Documentation Builds

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
