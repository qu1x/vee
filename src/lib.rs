// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! $`V_{EE}`$ -- Vector Expression Emitter: Geometric Algebra Code Generator
//!
//! The goal of this crate is to generate optimized code for geometric algebra flavors. Currently,
//! this crate implements the symbolic reduction of multivector expressions up to polynomials with
//! rational coefficients. In contrast, rational polynomials and hence polynomial division is not
//! required for lower dimensional geometric algebra flavors as the inverse of a multivector is
//! given by multiplying it with the inverse of its mixed-grade norm realizing a Study number for
//! dimensions $`D < 6`$, i.e., a generalized complex number.[^1] See the [contents](#contents)
//! below where the symbolic expressions are generated in text, code, and tree form. For evaluating
//! symbols as rationals, see the [code form in Rust mode](#code-form-in-rust-mode) example.
//!
//! Currently, the plane-based pistachio flavor -- Projective Geometric Algebra (PGA) -- is
//! implemented for $`D \equiv N + 1 \le 8`$ in all three metrics, i.e., elliptic, hyperbolic, and
//! parabolic (Euclidean).[^2]
//!
//!   * The PGAs with $`N < 2`$ are gated behind the `rudimentary` feature as both lack rotations.
//!   * The PGAs with $`N > 4`$ are gated behind the `exploratory` feature as there are no inverses
//!     based on Study numbers.
//!
//! Both features provide dimension-agnostic insights regarding duality and the choice of basis
//! blades. Additionally, the latter feature explores grade-preserving conditions among
//! orthonormalization conditions. The PGA is especially of interest for computer graphics (e.g.,
//! game and physics engines) as it is the most compact flavor (i.e., a one-up flavor) unifying the
//! established but scattered frameworks, e.g., homogeneous coordinates, Plücker coordinates, (dual)
//! quaternions, and screw theory. Even without any knowledge of geometric algebra, an API can be
//! more intuitive as it unifies the positional and directional aspects of geometric entities (e.g.,
//! planes, lines, points) and the linear and angular aspects of rigid-body dynamics in a
//! dimension-agnostic as well as metric-agnostic way with closed-form (i.e., non-iterative)
//! solutions up to 4D (e.g., [`PgaP2`], [`PgaP3`], [`PgaP4`]).[^3]
//!
//! [^1]: S. De Keninck and M. Roelfs, “Normalization, square roots, and the exponential and
//! logarithmic maps in geometric algebras of less than 6D”, [Mathematical Methods in the Applied
//! Sciences 47, 1425–1441](https://doi.org/10.1002/mma.8639).
//! [^2]: M. Roelfs and S. De Keninck, “Graded Symmetry Groups: Plane and Simple”, [Advances in
//! Applied Clifford Algebras 33](https://doi.org/10.1007/s00006-023-01269-9).
//! [^3]: L. Dorst and S. De Keninck, “Physical Geometry by Plane-Based Geometric Algebra”,
//! [Advanced Computational Applications of Geometric Algebra,
//! 43–76](https://doi.org/10.1007/978-3-031-55985-3_2).
//!
//! # Contents
//!
//!   * [Symbolic Engine](#symbolic-engine)
//!   * [Multivector Operators](#multivector-operators)
//!   * [Text Form Arguments](struct.Multivector.html#text-form-in-unicodeasciilatex-mode)
//!       * [Unicode Mode Examples](#text-form-in-unicode-mode)
//!       * [ASCII Mode Examples](#text-form-in-ascii-mode)
//!       * $`\LaTeX`$ [Mode Examples](#text-form-in-latex-mode)
//!   * [Code Form Arguments](struct.Multivector.html#code-form-in-genericrust-mode)
//!       * [Generic Mode Examples](#code-form-in-generic-mode)
//!       * [Rust Mode Examples](#code-form-in-rust-mode)
//!   * [Tree Form Arguments](struct.Multivector.html#tree-form-in-unicodeascii-mode)
//!       * [Unicode Mode Examples](#tree-form-in-unicode-mode)
//!       * [ASCII Mode Examples](#tree-form-in-ascii-mode)
//!
//! # Symbolic Engine
//!
//! The symbolic engine uses [`BTreeMap`] as storage container such that the canonical form of the
//! [`Multivector`]/[`Polynomial`]/[`Monomial`] hierarchy over [`Symbol`] is structurally enforced
//! by the choice of keys and values. Therefore, the non-volatile [`BTreeMap`] fields are public.
//! The canonical forms of non-zero [`Rational`] coefficients and non-zero [`Integer`] exponents are
//! behaviourally enforced by internally finalizing operators with [`Rational::reduce()`] and
//! excluding zero with the [`None`] variant of <code>[Option]<[Rational]></code> and
//! <code>[Option]<[Integer]></code>. Therefore, the volatile numerator/denominator field of
//! [`Rational`] and the volatile [`i64`] field of [`Integer`] are private.
//!
//! [`BTreeMap`]: `std::collections::BTreeMap`
//!
//! ## Features
//!
//!   * Zero non-optional dependencies.
//!   * Uniquely reduce symbolic multivector expressions for algebraic and structural equivalence to
//!     coincide.
//!   * Generate text form in Unicode/ASCII/$`\LaTeX`$ mode.
//!   * Generate code form in generic/Rust mode.
//!   * Generate tree form (i.e., DOT graphs as in [`text/vnd.graphviz`]) in Unicode/ASCII mode.
//!   * Eliminate orthonormalization conditions from expressions using reflection/projection
//!     operator by factoring pinned symbols, GCD coefficients, and predominant signs.
//!   * Evaluate symbols as rationals.
//!   * Count operations (i.e., multiplications and additions).
//!   * Define the metric-agnostic basis, i.e., elliptic, hyperbolic, and parabolic (Euclidean)
//!     along with the multivector entities for dimensions $`D \equiv N + 1 \le 8`$ of the
//!     plane-based pistachio flavor, i.e., projective geometric algebra (PGA).
//!
//! [`text/vnd.graphviz`]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
//!
//! ## Roadmap
//!
//!   * Simplify emitter by flattening expression tree into token stream and perform defer logic
//!     during stream iteration rather than tree traversal.
//!   * Further optimize expressions to reduce operation count by domain-specific common
//!     subexpression elimination (CSE) targeting exterior products. Reduce search space by
//!     leveraging applicable geometric decomposition, e.g., Euclidean decomposition for parabolic
//!     reflection operator.
//!   * Generate expressions in SIMD mode.
//!   * Define other geometric algebra flavors.
//!
//! # Multivector Operators
//!
//! Following table lists some common operators shared between flavors. The first three and the
//! omitted even-grade operators (e.g., square root, exponential, logarithm, power) must be manually
//! implemented based on Study numbers (i.e., generalized complex numbers) whereas the expressions
//! for the remaining ones are generated based on [`Multivector`] leveraging its respective symbolic
//! operator implementations. Otherwise, [`Multivector::norm_squared()`] is the closest offered and
//! [`Multivector::unit()`] is a marker, affecting only the multi-plane reflection operator (i.e,
//! the group conjugation or the sandwich product) and the projection operator inclusive rejection.
//! The table suggests the mixed-grade selection `B::from(a)` which is applicable whenever the
//! target implementation supports distinct types `B` for blades and versors rather than one
//! multivector type for all. Symbolic grade selection is supported with [`Multivector::grade()`]
//! and [`Multivector::vector()`] or with their mixed-grade complements [`Multivector::grades()`]
//! and [`Multivector::vectors()`].
//!
//! ```gdef
//! \gdef\norm{
//!   \| a \| \equiv \sqrt{a \tilde a}
//! }
//! \gdef\unit{
//!   \hat a \equiv \dfrac{a}{\| a \|}
//! }
//! \gdef\inv{
//!   a^{-1} \equiv \dfrac{\tilde a}{\| a \|^2}
//! }
//! \gdef\rev{
//!   \tilde a \equiv \sum_s (-1)^{s \choose 2} \lang a \rang_s
//! }
//! \gdef\pol{
//!   a^{\perp} \equiv a\I
//! }
//! \gdef\not{
//!   a^* \equiv \sum_s \lang a \rang_s^*
//!     : \lang a \rang_s^* = \sum_i \alpha_i a_i^*
//!       : a_i a_i^* = \I
//! }
//! \gdef\unnot{
//!   a_* \equiv a^{***} \therefore (a^*)_* = a^{****} = a
//! }
//! \gdef\neg{
//!   -a \equiv (-1)a
//! }
//! \gdef\add{
//!   a + b \equiv \sum_s \lang a \rang_s + \sum_t \lang b \rang_t
//! }
//! \gdef\sub{
//!   a - b \equiv \sum_s \lang a \rang_s - \sum_t \lang b \rang_t
//! }
//! \gdef\mul{
//!   ab \equiv \sum_{s,t} \lang a \rang_s \lang b \rang_t
//! }
//! \gdef\div{
//!   \dfrac{a}{b} \equiv ab^{-1}
//! }
//! \gdef\rem{
//!   a \times b \equiv \frac{1}{2}(ab - ba)
//! }
//! \gdef\bitor{
//!   a \mid b \equiv \sum_{s,t}
//!     \lang
//!       \lang a \rang_s
//!       \lang b \rang_t
//!     \rang_{|s - t|}
//! }
//! \gdef\bitxor{
//!   a \wedge b \equiv \sum_{s,t}
//!     \lang
//!       \lang a \rang_s
//!       \lang b \rang_t
//!     \rang_{s + t}
//! }
//! \gdef\bitand{
//!   a \vee b \equiv {(a^* \wedge b^*)}_*
//! }
//! \gdef\shl{
//!   a \looparrowleft b
//!     \equiv (-1)^{\prod_{s,t} st} ba \tilde b
//! }
//! \gdef\shr{
//!   a \curvearrowright b \equiv (a \mid b) \tilde b
//! }
//! \gdef\from{
//!   \lang a \rang_b \equiv \sum_{s \in \{t|b = \sum_t \lang b \rang_t\}} \lang a \rang_s
//! }
//! ```
//!
//! Operator            | Name                          | Formula
//! ------------------- | ----------------------------- | --------
//! `a.norm()`          | Norm (mixed grade)            | $`\norm`$
//! `a.unit()`          | Unit (orthonormal)            | $`\unit`$
//! `a.inv()`           | Inverse                       | $`\inv`$
//! `a.rev()`           | Reverse                       | $`\rev`$
//! `a.pol()`           | Polarity                      | $`\pol`$
//! `!a`                | Dual (right complement)       | $`\not`$
//! `!!!a`              | Undual (left complement)      | $`\unnot`$
//! `-a`                | Negation (orientation)        | $`\neg`$
//! `B::from(a)`        | Selection (mixed grade)       | $`\from`$
//! `a + b`, `a += b`   | Sum                           | $`\add`$
//! `a - b`, `a -= b`   | Difference                    | $`\sub`$
//! `a * b`, `a *= b`   | Product (geometric)           | $`\mul`$
//! `a / b`, `a /= b`   | Quotient (geometric)          | $`\div`$
//! `a << b`, `a <<= b` | Reflection ($`a`$ by $`b`$)   | $`\shl`$
//! `a >> b`, `a >>= b` | Projection ($`a`$ onto $`b`$) | $`\shr`$
//! `a % b`             | Commutator                    | $`\rem`$
//! `a \| b`            | Contraction (symmetric)       | $`\bitor`$
//! `a ^ b`             | Meet (progressive)            | $`\bitxor`$
//! `a & b`             | Join (regressive)             | $`\bitand`$
//!
//! # Text Form in Unicode Mode
//!
//! Generates the expression for rotating and/or translating a point in [`PgaP3`], i.e., the type
//! alias of [`Multivector`] parameterized for the Parabolic (Euclidean) 3D PGA. The
//! [`PgaP3::pin()`] method pins symbols of [`PgaP3::point()`] with the *combining x below* (i.e.,
//! the Unicode *combining diacritical mark* `"◌͓"`) to distinguish them from the symbols of
//! [`PgaP3::motor()`]. This isometry (i.e., up to a screw motion) is isomorphic to the
//! transformation of a homogeneous point by a dual quaternion. When importing type aliases, the
//! examples rename the multivector type as `Vee` and its basis blade type as `Bee`, with the former
//! being parameterized by the latter.
//!
//! [`PgaP3::motor()`]: struct.Multivector.html#method.motor-1
#![cfg_attr(
    not(feature = "rudimentary"),
    doc = "[`PgaP3::point()`]: struct.Multivector.html#method.point-1"
)]
#![cfg_attr(
    feature = "rudimentary",
    doc = "[`PgaP3::point()`]: struct.Multivector.html#method.point-2"
)]
//!
//! ```
//! use vee::{PgaP3 as Vee, Symbol, format_eq, pga::PgaP3 as Bee};
//!
//! // Assumes motor is not orthonormalized.
//! format_eq!(Vee::point().pin() << Vee::motor(), [
//!     "+(+vv+xx+yy+zz)w͓e123",
//!     "+(+(+vv+xx-yy-zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓+2(-Vx-Xv-Yz+Zy)w͓)e032",
//!     "+(+2(-vz+xy)X͓+(+vv-xx+yy-zz)Y͓+2(+vx+yz)Z͓+2(-Vy+Xz-Yv-Zx)w͓)e013",
//!     "+(+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+vv-xx-yy+zz)Z͓+2(-Vz-Xy+Yx-Zv)w͓)e021",
//! ]);
//!
//! // Assumes motor is orthonormalized.
//! format_eq!(Vee::point().pin() << Vee::motor().unit(), [
//!     "+w͓e123",
//!     "+(+(+1-2yy-2zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓+2(-Vx-Xv-Yz+Zy)w͓)e032",
//!     "+(+2(-vz+xy)X͓+(+1-2xx-2zz)Y͓+2(+vx+yz)Z͓+2(-Vy+Xz-Yv-Zx)w͓)e013",
//!     "+(+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+1-2xx-2yy)Z͓+2(-Vz-Xy+Yx-Zv)w͓)e021",
//! ]);
//!
//! // Assumes motor and point are (ortho)normalized where point has positive orientation.
//! format_eq!(Vee::point().eval([(Bee::e123(), 1)]).pin() << Vee::motor().unit(), [
//!     "+e123",
//!     "+(+2(-Vx-Xv-Yz+Zy)+(+1-2yy-2zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓)e032",
//!     "+(+2(-Vy+Xz-Yv-Zx)+2(-vz+xy)X͓+(+1-2xx-2zz)Y͓+2(+vx+yz)Z͓)e013",
//!     "+(+2(-Vz-Xy+Yx-Zv)+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+1-2xx-2yy)Z͓)e021",
//! ]);
//!
//! // The `eval()` method accepts `Into<Symbol>` which is implemented for `(Symbol, Bee)`.
//! assert_eq!(
//!     const { Bee::e123() },
//!     const { (Symbol::new('w', "e123"), Bee::new("e123").unwrap()) },
//! );
//! ```
//!
//! The symbols are assigned to basis blades such that lowercase symbols are dual to their
//! corresponding uppercase symbols. For blades containing $`\e_0`$, uppercase symbols are used. The
//! [`PgaP3::swp()`] method swaps lowercase and uppercase symbols. This is useful for testing
//! duality equivalences.
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!(Vee::plane(), ["+We0", "+xe1", "+ye2", "+ze3"]);
//! format_eq!(Vee::point(), ["+we123", "+Xe032", "+Ye013", "+Ze021"]);
//!
//! assert_ne!(!Vee::plane(), Vee::point());
//! assert_eq!(!Vee::plane(), Vee::point().swp());
//! ```
//!
//! Optional plus signs are skipped with `"{:<}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:<}", Vee::point().pin() << Vee::motor().unit(), [
//!     "w͓e123",
//!     "+((1-2yy-2zz)X͓+2(vz+xy)Y͓+2(-vy+xz)Z͓+2(-Vx-Xv-Yz+Zy)w͓)e032",
//!     "+(2(-vz+xy)X͓+(1-2xx-2zz)Y͓+2(vx+yz)Z͓+2(-Vy+Xz-Yv-Zx)w͓)e013",
//!     "+(2(vy+xz)X͓+2(-vx+yz)Y͓+(1-2xx-2yy)Z͓+2(-Vz-Xy+Yx-Zv)w͓)e021",
//! ]);
//! ```
//!
//! The predominant sign is factored as well with `"{:-}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! // Unfactored predominant sign.
//! format_eq!(Vee::point().pin() << Vee::motor(), [
//!    "+(+vv+xx+yy+zz)w͓e123",
//!    "+(+(+vv+xx-yy-zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓+2(-Vx-Xv-Yz+Zy)w͓)e032",
//!    "+(+2(-vz+xy)X͓+(+vv-xx+yy-zz)Y͓+2(+vx+yz)Z͓+2(-Vy+Xz-Yv-Zx)w͓)e013",
//!    "+(+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+vv-xx-yy+zz)Z͓+2(-Vz-Xy+Yx-Zv)w͓)e021",
//! //                                          ^^^^^^^^^^^^^^^^^^
//! ]);
//!
//! // Factored predominant sign.
//! format_eq!("{:-}", Vee::point().pin() << Vee::motor(), [
//!    "+(+vv+xx+yy+zz)w͓e123",
//!    "+(+(+vv+xx-yy-zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓-2(+Vx+Xv+Yz-Zy)w͓)e032",
//!    "+(+2(-vz+xy)X͓+(+vv-xx+yy-zz)Y͓+2(+vx+yz)Z͓-2(+Vy-Xz+Yv+Zx)w͓)e013",
//!    "+(+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+vv-xx-yy+zz)Z͓-2(+Vz+Xy-Yx+Zv)w͓)e021",
//! //                                          ^^^^^^^^^^^^^^^^^^
//! ]);
//! ```
//!
//! The factorization is skipped with `"{:+}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:+}", Vee::point().pin() << Vee::motor(), [
//!     "+(+vvw͓+w͓xx+w͓yy+w͓zz)e123",
//!     "+(-2Vw͓x-2Xvw͓+X͓vv+X͓xx-X͓yy-X͓zz-2Yw͓z+2Y͓vz+2Y͓xy+2Zw͓y-2Z͓vy+2Z͓xz)e032",
//!     "+(-2Vw͓y+2Xw͓z-2X͓vz+2X͓xy-2Yvw͓+Y͓vv-Y͓xx+Y͓yy-Y͓zz-2Zw͓x+2Z͓vx+2Z͓yz)e013",
//!     "+(-2Vw͓z-2Xw͓y+2X͓vy+2X͓xz+2Yw͓x-2Y͓vx+2Y͓yz-2Zvw͓+Z͓vv-Z͓xx-Z͓yy+Z͓zz)e021",
//! ]);
//! ```
//!
//! # Text Form in ASCII Mode
//!
//! Alternatively, symbols are labelled after their initially assigned basis blades starting with:
//!
//!   * `'p'` if pinned with [`PgaP3::pin()`],
//!   * `'l'` if left-hand side as in [`PgaP3::lhs()`],
//!   * `'r'` if right-hand side as in [`PgaP3::rhs()`],
//!   * `'v'` otherwise.
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:#}", Vee::point().pin() << Vee::motor().unit(), [
//!     "+p123*e123",
//!     "+(+(+1-2*v31*v31-2*v12*v12)*p032+2*(+v*v12+v23*v31)*p013+2*(-v*v31+v23*v12)*p021\
//!        +2*(-v0123*v23-v01*v-v02*v12+v03*v31)*p123)*e032",
//!     "+(+2*(-v*v12+v23*v31)*p032+(+1-2*v23*v23-2*v12*v12)*p013+2*(+v*v23+v31*v12)*p021\
//!        +2*(-v0123*v31+v01*v12-v02*v-v03*v23)*p123)*e013",
//!     "+(+2*(+v*v31+v23*v12)*p032+2*(-v*v23+v31*v12)*p013+(+1-2*v23*v23-2*v31*v31)*p021\
//!        +2*(-v0123*v12-v01*v31+v02*v23-v03*v)*p123)*e021",
//! ]);
//!
//! format_eq!("{:#}", Vee::line().lhs() * Vee::line().rhs(), [
//!     "-l23*r23-l31*r31-l12*r12",
//!     "+(-l02*r12+r02*l12+l03*r31-r03*l31)*e01",
//!     "+(+l01*r12-r01*l12-l03*r23+r03*l23)*e02",
//!     "+(-l01*r31+r01*l31+l02*r23-r02*l23)*e03",
//!     "+(-l31*r12+r31*l12)*e23",
//!     "+(+l23*r12-r23*l12)*e31",
//!     "+(-l23*r31+r23*l31)*e12",
//!     "+(+l01*r23+r01*l23+l02*r31+r02*l31+l03*r12+r03*l12)*I",
//! ]);
//! ```
//!
//! # Text Form in $`\LaTeX`$ Mode
//!
//! Generate $`\LaTeX`$ documentation as in
//!
//! ```math
//! \begin{aligned}
//!   (-r_{01} l_1 - r_{02} l_2 - r_{03} l_3) & \boldsymbol{e}_0 \\
//!   + (-l_2 r_{12} + r_{31} l_3) & \boldsymbol{e}_1 \\
//!   + (l_1 r_{12} - r_{23} l_3) & \boldsymbol{e}_2 \\
//!   + (-l_1 r_{31} + r_{23} l_2) & \boldsymbol{e}_3 \\
//!   + (l_1 r_{23} + l_2 r_{31} + l_3 r_{12}) & \boldsymbol{e}_{123} \\
//!   + (-l_0 r_{23} - r_{02} l_3 + r_{03} l_2) & \boldsymbol{e}_{032} \\
//!   + (-l_0 r_{31} + r_{01} l_3 - r_{03} l_1) & \boldsymbol{e}_{013} \\
//!   + (-l_0 r_{12} - r_{01} l_2 + r_{02} l_1) & \boldsymbol{e}_{021}
//! \end{aligned}
//! ```
//!
//! with `"{:$>}"`, using only the standard `amsmath` package:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:$>}", Vee::plane().lhs() * Vee::line().rhs(), [
//!     r"\begin{aligned}",
//!     r"  (-r_{01} l_1 - r_{02} l_2 - r_{03} l_3) & \boldsymbol{e}_0 \\",
//!     r"  + (-l_2 r_{12} + r_{31} l_3) & \boldsymbol{e}_1 \\",
//!     r"  + (l_1 r_{12} - r_{23} l_3) & \boldsymbol{e}_2 \\",
//!     r"  + (-l_1 r_{31} + r_{23} l_2) & \boldsymbol{e}_3 \\",
//!     r"  + (l_1 r_{23} + l_2 r_{31} + l_3 r_{12}) & \boldsymbol{e}_{123} \\",
//!     r"  + (-l_0 r_{23} - r_{02} l_3 + r_{03} l_2) & \boldsymbol{e}_{032} \\",
//!     r"  + (-l_0 r_{31} + r_{01} l_3 - r_{03} l_1) & \boldsymbol{e}_{013} \\",
//!     r"  + (-l_0 r_{12} - r_{01} l_2 + r_{02} l_1) & \boldsymbol{e}_{021}",
//!     r"\end{aligned}",
//! ]);
//! ```
//!
//! # Code Form in Generic Mode
//!
//! Generate generic statements with `"{:x}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:x}", Vee::rotator().lhs() * Vee::rotator().rhs(), [
//!     "e=l*r-l23*r23-l31*r31-l12*r12",
//!     "e23=l*r23+r*l23-l31*r12+r31*l12",
//!     "e31=l*r31+r*l31+l23*r12-r23*l12",
//!     "e12=l*r12+r*l12-l23*r31+r23*l31",
//! ]);
//! ```
//!
//! # Code Form in Rust Mode
//!
//! Generate Rust code with `"{:#x}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:#x}", Vee::rotator().lhs() * Vee::rotator().rhs(), [
//!     "let e = l * r - l23 * r23 - l31 * r31 - l12 * r12;",
//!     "let e23 = l * r23 + r * l23 - l31 * r12 + r31 * l12;",
//!     "let e31 = l * r31 + r * l31 + l23 * r12 - r23 * l12;",
//!     "let e12 = l * r12 + r * l12 - l23 * r31 + r23 * l31;",
//! ]);
//! ```
//!
//! Either emit or omit zero vectors:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq, pga::PgaP3 as Bee};
//!
//! // Emit zero vector as numeric zero.
//! let plane = Vee::plane()
//!     .eval([(Bee::e0(), 0), (Bee::e1(), 1)])
//!     .eval([(Bee::e2(), (1, 2))]);
//!
//! format_eq!("{:#x}", plane, [
//!     "let e0 = 0.0;",
//!     "let e1 = 1.0;",
//!     "let e2 = 1.0 / 2.0;",
//!     "let e3 = v3;",
//! ]);
//!
//! // Omit zero vector as symbolic zero.
//! let normal = plane.omit();
//!
//! format_eq!("{:#x}", normal, [
//!     "let e1 = 1.0;",
//!     "let e2 = 1.0 / 2.0;",
//!     "let e3 = v3;",
//! ]);
//! ```
//!
//! Generate Rust code dereferencing fields with `"{:^x}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:^x}", Vee::rotator().lhs() * Vee::rotator().rhs(), [
//!     "o.e = l.e * r.e - l.e23 * r.e23 - l.e31 * r.e31 - l.e12 * r.e12;",
//!     "o.e23 = l.e * r.e23 + r.e * l.e23 - l.e31 * r.e12 + r.e31 * l.e12;",
//!     "o.e31 = l.e * r.e31 + r.e * l.e31 + l.e23 * r.e12 - r.e23 * l.e12;",
//!     "o.e12 = l.e * r.e12 + r.e * l.e12 - l.e23 * r.e31 + r.e23 * l.e31;",
//! ]);
//! ```
//!
//! # Tree Form in Unicode Mode
//!
//! Generate DOT graphs (i.e., [`text/vnd.graphviz`]) in Unicode mode with `"{:o}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:o}", Vee::plane(), [
//!     r#"digraph vee {"#,
//!     r#"  n0 [label="∑" shape=box]"#,
//!     r#"  n1 [label="∏" shape=box]"#,
//!     r#"  n0 -> n1"#,
//!     r#"  n2 [label="W" shape=ellipse]"#,
//!     r#"  n1 -> n2"#,
//!     r#"  n3 [label="e0" shape=diamond]"#,
//!     r#"  n1 -> n3"#,
//!     r#"  n4 [label="∏" shape=box]"#,
//!     r#"  n0 -> n4"#,
//!     r#"  n5 [label="x" shape=ellipse]"#,
//!     r#"  n4 -> n5"#,
//!     r#"  n6 [label="e1" shape=diamond]"#,
//!     r#"  n4 -> n6"#,
//!     r#"  n7 [label="∏" shape=box]"#,
//!     r#"  n0 -> n7"#,
//!     r#"  n8 [label="y" shape=ellipse]"#,
//!     r#"  n7 -> n8"#,
//!     r#"  n9 [label="e2" shape=diamond]"#,
//!     r#"  n7 -> n9"#,
//!     r#"  n10 [label="∏" shape=box]"#,
//!     r#"  n0 -> n10"#,
//!     r#"  n11 [label="z" shape=ellipse]"#,
//!     r#"  n10 -> n11"#,
//!     r#"  n12 [label="e3" shape=diamond]"#,
//!     r#"  n10 -> n12"#,
//!     r#"}"#,
//! ]);
//! ```
//!
//! # Tree Form in ASCII Mode
//!
//! Generate DOT graphs (i.e., [`text/vnd.graphviz`]) in ASCII mode with `"{:#o}"`:
//!
//! ```
//! use vee::{PgaP3 as Vee, format_eq};
//!
//! format_eq!("{:#o}", Vee::plane(), [
//!     r#"digraph vee {"#,
//!     r#"  n0 [label="+" shape=box]"#,
//!     r#"  n1 [label="*" shape=box]"#,
//!     r#"  n0 -> n1"#,
//!     r#"  n2 [label="v0" shape=ellipse]"#,
//!     r#"  n1 -> n2"#,
//!     r#"  n3 [label="e0" shape=diamond]"#,
//!     r#"  n1 -> n3"#,
//!     r#"  n4 [label="*" shape=box]"#,
//!     r#"  n0 -> n4"#,
//!     r#"  n5 [label="v1" shape=ellipse]"#,
//!     r#"  n4 -> n5"#,
//!     r#"  n6 [label="e1" shape=diamond]"#,
//!     r#"  n4 -> n6"#,
//!     r#"  n7 [label="*" shape=box]"#,
//!     r#"  n0 -> n7"#,
//!     r#"  n8 [label="v2" shape=ellipse]"#,
//!     r#"  n7 -> n8"#,
//!     r#"  n9 [label="e2" shape=diamond]"#,
//!     r#"  n7 -> n9"#,
//!     r#"  n10 [label="*" shape=box]"#,
//!     r#"  n0 -> n10"#,
//!     r#"  n11 [label="v3" shape=ellipse]"#,
//!     r#"  n10 -> n11"#,
//!     r#"  n12 [label="e3" shape=diamond]"#,
//!     r#"  n10 -> n12"#,
//!     r#"}"#,
//! ]);
//! ```

#![cfg_attr(docsrs, feature(doc_cfg))]

mod factor;
use factor::Factor;

mod algebra;
mod choose;
mod factorization;
mod integer;
mod monomial;
mod multivector;
mod polynomial;
mod rational;
mod symbol;
mod tree;

pub use algebra::Algebra;
pub use choose::choose;
pub use factorization::Factorization;
pub use integer::Integer;
pub use monomial::Monomial;
pub use multivector::Multivector;
pub use polynomial::Polynomial;
pub use rational::Rational;
pub use symbol::Symbol;
pub use tree::Tree;

#[cfg(not(feature = "pretty_assertions"))]
#[doc(hidden)]
pub use core::assert_eq as __assert_eq;
#[cfg(feature = "pretty_assertions")]
#[doc(hidden)]
pub use pretty_assertions::assert_eq as __assert_eq;

/// Formats the `$lhs` expression using [`Display`] and asserts the `$rhs` string literals.
///
/// Passes `$fmt` to [`Display`] with `{}` as default if omitted. Appends `"\n"` to each `$rhs`
/// literal and asserts the concatenation thereof.
///
/// With the `pretty_assertions` feature, the respective [`assert_eq!`] macro is used. In this way,
/// the Unicode *combining diacritical marks* are rendered as in the examples using [`format_eq!`].
///
/// [`Display`]: `core::fmt::Display`
#[macro_export]
macro_rules! format_eq {
    ($lhs:expr, [$($rhs:literal),* $(,)?]) => {{
        $crate::format_eq!("{}", $lhs, [$($rhs),*]);
    }};
    ($fmt:literal, $lhs:expr, [$($rhs:literal),* $(,)?]) => {{
        let lhs = format!($fmt, $lhs);
        let mut rhs = String::with_capacity(lhs.len());
        $(
            rhs.push_str($rhs);
            rhs.push_str("\n");
        )*
        $crate::__assert_eq!(lhs, rhs);
    }};
}

/// The unary logical negation assignment operator.
pub trait NotAssign {
    /// Performs the unary logical negation assignment operation.
    fn not_assign(&mut self);
}

/// The unary negation assignment operator.
pub trait NegAssign {
    /// Performs the unary negation assignment operation.
    fn neg_assign(&mut self);
}

/// The unary reversion operator.
pub trait Rev {
    /// The resulting type after applying the unary reversion operator.
    type Output;

    /// Performs the unary reversion operation.
    #[must_use]
    fn rev(self) -> Self::Output;
}

/// The unary reversion assignment operator.
pub trait RevAssign {
    /// Performs the unary reversion assignment operation.
    fn rev_assign(&mut self);
}

/// The unary inversion operator.
pub trait Inv {
    /// The resulting type after applying the unary inversion operator.
    type Output;

    /// Performs the unary inversion operation.
    #[must_use]
    fn inv(self) -> Self::Output;
}

/// The unary inversion assignment operator.
pub trait InvAssign {
    /// Performs the unary inversion assignment operation.
    fn inv_assign(&mut self);
}

pub mod pga;

/// Multivector for Elliptic 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaE0 = Multivector<pga::PgaE0>;
/// Multivector for Elliptic 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaE1 = Multivector<pga::PgaE1>;
/// Multivector for Elliptic 2D PGA.
pub type PgaE2 = Multivector<pga::PgaE2>;
/// Multivector for Elliptic 3D PGA.
pub type PgaE3 = Multivector<pga::PgaE3>;
/// Multivector for Elliptic 4D PGA.
pub type PgaE4 = Multivector<pga::PgaE4>;
/// Multivector for Elliptic 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE5 = Multivector<pga::PgaE5>;
/// Multivector for Elliptic 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE6 = Multivector<pga::PgaE6>;
/// Multivector for Elliptic 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE7 = Multivector<pga::PgaE7>;

/// Multivector for Hyperbolic 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaH0 = Multivector<pga::PgaH0>;
/// Multivector for Hyperbolic 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaH1 = Multivector<pga::PgaH1>;
/// Multivector for Hyperbolic 2D PGA.
pub type PgaH2 = Multivector<pga::PgaH2>;
/// Multivector for Hyperbolic 3D PGA.
pub type PgaH3 = Multivector<pga::PgaH3>;
/// Multivector for Hyperbolic 4D PGA.
pub type PgaH4 = Multivector<pga::PgaH4>;
/// Multivector for Hyperbolic 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH5 = Multivector<pga::PgaH5>;
/// Multivector for Hyperbolic 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH6 = Multivector<pga::PgaH6>;
/// Multivector for Hyperbolic 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH7 = Multivector<pga::PgaH7>;

/// Multivector for Parabolic (Euclidean) 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaP0 = Multivector<pga::PgaP0>;
/// Multivector for Parabolic (Euclidean) 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaP1 = Multivector<pga::PgaP1>;
/// Multivector for Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Multivector<pga::PgaP2>;
/// Multivector for Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Multivector<pga::PgaP3>;
/// Multivector for Parabolic (Euclidean) 4D PGA.
pub type PgaP4 = Multivector<pga::PgaP4>;
/// Multivector for Parabolic (Euclidean) 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP5 = Multivector<pga::PgaP5>;
/// Multivector for Parabolic (Euclidean) 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP6 = Multivector<pga::PgaP6>;
/// Multivector for Parabolic (Euclidean) 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP7 = Multivector<pga::PgaP7>;
