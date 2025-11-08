// Copyright © 2025 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! $`V_{EE}`$ -- Vector Expression Emitter: Geometric Algebra Code Generator
//!
//! The goal of this crate is to generate optimized code for geometric algebra flavors. Currently,
//! this crate implements the symbolic reduction of multivector expressions up to polynomials with
//! rational coefficients. In contrast, rational polynomials and hence polynomial division is not
//! required for lower dimensional geometric algebra flavors as the inverse of a multivector is
//! given by multiplying it with the inverse of its mixed-grade norm, i.e., a Study number for
//! dimensions $`D < 6`$.[^1] See the [examples](#examples) below where the symbolic expressions are
//! generated in text form. The next releases will implement code forms (e.g., Rust code in various
//! profiles based on SIMD using [`lav`] with and without generics or arbitrary precision types
//! using [`rug`]). The pre-generated code forms will be provided along with the code generator
//! behind respective feature gates. When [`packages_as_namespaces`] is stable, each code form will
//! become a crate. Currently, the plane-based pistachio flavor -- Projective Geometric Algebra
//! (PGA) -- is implemented for $`D \equiv N + 1 \le 8`$ in all three metrics, i.e., elliptic,
//! hyperbolic, and parabolic (Euclidean).[^2] The 5D, 6D, and 7D PGAs (i.e., $`N = 5`$, $`N = 6`$,
//! and $`N = 7`$) are exploratory as there are no inverses based on Study numbers. They provide
//! dimension-agnostic insights regarding duality, the choice of basis blades, and grade-preserving
//! conditions among orthonormalization conditions. The PGA is especially of interest for computer
//! graphics (e.g., game and physics engines) as it is the most compact flavor (i.e., a one-up
//! flavor) unifying the established but scattered frameworks, e.g., homogeneous coordinates,
//! Plücker coordinates, (dual) quaternions, and screw theory. Even without any knowledge of
//! geometric algebra, an API can be more intuitive as it unifies the positional and directional
//! aspects of geometric entities (e.g., planes, lines, points) and the linear and angular aspects
//! of rigid-body dynamics in a dimension-agnostic way with closed-form (i.e., non-iterative)
//! solutions up to 4D (e.g., [`PgaP2`], [`PgaP3`], [`PgaP4`]).[^3]
//!
//! [`packages_as_namespaces`]:
//! https://rust-lang.github.io/rfcs/3243-packages-as-optional-namespaces.html
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
//! [`lav`]: https://docs.rs/lav
//! [`rug`]: https://docs.rs/rug
//!
//! # Operators
//!
//! Following table lists the common operators shared between flavors. The code for the first three
//! will be manually written based on Study numbers whereas the code for the remaining ones will be
//! automatically generated based on [`Multivector`].
//!
//! ```gdef
//! \gdef\e{
//!   \boldsymbol e
//! }
//! \gdef\I{
//!   \boldsymbol I
//! }
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
//!     \rang_{\|s - t\|}
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
//!     \equiv \sum_{s,t} (-1)^{st} \lang b \rang_t \lang a \rang_s \lang \tilde b \rang_t
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
//! # Examples
//!
//! Generates the expression for rotating and/or translating a point in [`PgaP3`], i.e., the type
//! alias of [`Multivector`] parameterized for the Parabolic (Euclidean) 3D PGA. The
//! [`PgaP3::pin()`] method pins symbols of [`PgaP3::point()`] with the *combining x below* (i.e.,
//! the Unicode *combining diacritical mark* `"◌͓"`) to distinguish them from the symbols of
//! [`PgaP3::motor()`](struct.Multivector.html#method.rotator-1). This isometry (i.e., up to a screw
//! motion) is isomorphic to the transformation of a homogeneoous point by a dual quaternion.
//!
//! ```
//! use vee::{format_eq, PgaP3 as Vee};
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
//! format_eq!(Vee::point().eval([(('w', "e123"), 1)]).pin() << Vee::motor().unit(), [
//!     "+e123",
//!     "+(+2(-Vx-Xv-Yz+Zy)+(+1-2yy-2zz)X͓+2(+vz+xy)Y͓+2(-vy+xz)Z͓)e032",
//!     "+(+2(-Vy+Xz-Yv-Zx)+2(-vz+xy)X͓+(+1-2xx-2zz)Y͓+2(+vx+yz)Z͓)e013",
//!     "+(+2(-Vz-Xy+Yx-Zv)+2(+vy+xz)X͓+2(-vx+yz)Y͓+(+1-2xx-2yy)Z͓)e021",
//! ]);
//! ```
//!
//! The symbols are assigned to basis blades such that lowercase symbols are dual to their
//! corresponding uppercase symbols. For blades containing $`\e_0`$, uppercase symbols are used. The
//! [`PgaP3::swp()`] method swaps lowercase and uppercase symbols. This is useful for testing
//! duality equivalences.
//!
//! ```
//! use vee::{format_eq, PgaP3 as Vee};
//!
//! format_eq!(Vee::plane(), ["+We0", "+xe1", "+ye2", "+ze3"]);
//! format_eq!(Vee::point(), ["+we123", "+Xe032", "+Ye013", "+Ze021"]);
//!
//! assert_ne!(!Vee::plane(), Vee::point());
//! assert_eq!(!Vee::plane(), Vee::point().swp());
//! ```
//!
//! Alternatively, symbols are labelled after their initially assigned basis blades starting with:
//!
//!   * `'p'` if pinned with [`PgaP3::pin()`],
//!   * `'l'` if left-hand side as in [`PgaP3::lhs()`],
//!   * `'r'` if right-hand side as in [`PgaP3::rhs()`],
//!   * `'v'` otherwise.
//!
//! ```
//! use vee::{format_eq, PgaP3 as Vee};
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
//! Optional plus signs are skipped with `"{:<}"`:
//!
//! ```
//! use vee::{format_eq, PgaP3 as Vee};
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
//! use vee::{format_eq, PgaP3 as Vee};
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
//! use vee::{format_eq, PgaP3 as Vee};
//!
//! format_eq!("{:+}", Vee::point().pin() << Vee::motor(), [
//!     "+(+vvw͓+w͓xx+w͓yy+w͓zz)e123",
//!     "+(-2Vw͓x-2Xvw͓+X͓vv+X͓xx-X͓yy-X͓zz-2Yw͓z+2Y͓vz+2Y͓xy+2Zw͓y-2Z͓vy+2Z͓xz)e032",
//!     "+(-2Vw͓y+2Xw͓z-2X͓vz+2X͓xy-2Yvw͓+Y͓vv-Y͓xx+Y͓yy-Y͓zz-2Zw͓x+2Z͓vx+2Z͓yz)e013",
//!     "+(-2Vw͓z-2Xw͓y+2X͓vy+2X͓xz+2Yw͓x-2Y͓vx+2Y͓yz-2Zvw͓+Z͓vv-Z͓xx-Z͓yy+Z͓zz)e021",
//! ]);
//! ```

/// Formats the `$lhs` expression using [`Display`] and asserts the `$rhs` string literals.
///
/// Passes `$fmt` to [`Display`] with `{}` as default if omitted. Appends `"\n"` to each `$rhs`
/// literal and asserts the concatenation thereof.
///
/// With the `pretty_assertions` feature, the respective [`assert_eq!`] macro is used. In this way,
/// the Unicode *combining diacritical marks* are rendered as in the examples using [`format_eq!`].
#[macro_export]
macro_rules! format_eq {
    ($lhs:expr, [$($rhs:literal),* $(,)?]) => {{
        format_eq!("{}", $lhs, [$($rhs),*]);
    }};
    ($fmt:literal, $lhs:expr, [$($rhs:literal),* $(,)?]) => {{
        #[cfg(feature = "pretty_assertions")]
        use pretty_assertions::assert_eq;
        let lhs = format!($fmt, $lhs);
        let mut rhs = String::with_capacity(lhs.len());
        $(
            rhs.push_str($rhs);
            rhs.push_str("\n");
        )*
        assert_eq!(lhs, rhs);
    }};
}

use core::{
    cmp::{Ordering, min},
    fmt::{self, Debug, Display, Octal},
    iter::{FromIterator, repeat},
    mem::{swap, take},
    num::{NonZero, NonZeroI32},
    ops::{
        Add, AddAssign, BitAnd, BitOr, BitXor, Div, DivAssign, Mul, MulAssign, Neg, Not, Rem, Shl,
        Shr, Sub, SubAssign,
    },
};
use std::collections::{BTreeMap, BTreeSet};

/// Finds the binomial coefficient.
pub trait Choose {
    /// The output type.
    type Output;

    /// Finds the binomial coefficient as in `self` over `other`.
    #[must_use]
    fn choose(self, other: Self) -> Self::Output;
}

impl Choose for u32 {
    type Output = Self;

    fn choose(self, other: Self) -> Self::Output {
        let (n, k) = (self, other);
        if n < k {
            0
        } else {
            (0..k).fold(1, |r, i| r * (n - i) / (i + 1))
        }
    }
}

/// Finds the greatest common divisor (GCD) and the least common multiple (LCM).
pub trait Factor
where
    Self: Copy + Mul<Output = Self> + Div<Output = Self> + PartialEq + Eq + Default,
{
    /// The zero constant.
    const ZERO: Self;

    /// Finds the greatest common divisor (GCD) of `self` and `other`.
    #[must_use]
    fn gcd(self, other: Self) -> Self;
    /// Finds the least common multiple (LCM) of `self` and `other`.
    ///
    /// Calls [`Self::gcd_lcm()`] and discards GCD.
    #[inline]
    #[must_use]
    fn lcm(self, other: Self) -> Self {
        self.gcd_lcm(other).1
    }
    /// Returns <code>([Self::gcd()], [Self::lcm()])</code> of `self` and `other`.
    ///
    /// Calls [`Self::gcd()`] to find LCM.
    #[inline]
    #[must_use]
    fn gcd_lcm(self, other: Self) -> (Self, Self) {
        if self == Self::ZERO && other == Self::ZERO {
            (Self::ZERO, Self::ZERO)
        } else {
            let gcd = self.gcd(other);
            let lcm = self * (other / gcd);
            (gcd, lcm)
        }
    }
    /// Finds the GCD of iterator over `Self`.
    #[inline]
    #[must_use]
    fn gcd_bulk(r: impl IntoIterator<Item = Self>) -> Self {
        r.into_iter().reduce(Self::gcd).unwrap_or_default()
    }
    /// Finds the LCM of iterator over `Self`.
    #[inline]
    #[must_use]
    fn lcm_bulk(r: impl IntoIterator<Item = Self>) -> Self {
        r.into_iter().reduce(Self::lcm).unwrap_or_default()
    }
    /// Finds the [`Self::gcd`] or [`Self::lcm`] and the predominant sign of iterator over `Self`.
    #[must_use]
    fn signed(f: impl Fn(Self, Self) -> Self, r: impl IntoIterator<Item = Self>) -> Self
    where
        Self: Neg<Output = Self> + PartialOrd + Ord,
    {
        let (acc, neg, len) = r
            .into_iter()
            .fold((Self::ZERO, 0, 0), |(acc, neg, len), r| {
                (
                    f(acc, r),
                    if r < Self::ZERO { neg + 1 } else { neg },
                    len + 1,
                )
            });
        if neg > len / 2 { -acc } else { acc }
    }
}

impl Factor for Rational {
    const ZERO: Self = Self::ZERO;

    #[inline]
    fn gcd(self, other: Self) -> Self {
        Self {
            p: self.p.gcd(other.p),
            q: self.q.lcm(other.q),
        }
    }
    #[inline]
    fn lcm(self, other: Self) -> Self {
        Self {
            p: self.p.lcm(other.p),
            q: self.q.gcd(other.q),
        }
    }
    #[inline]
    fn gcd_lcm(self, other: Self) -> (Self, Self) {
        let (g_p, l_p) = self.p.gcd_lcm(other.p);
        let (g_q, l_q) = self.q.gcd_lcm(other.q);
        (Self { p: g_p, q: l_q }, Self { p: l_p, q: g_q })
    }
}

impl Factor for i32 {
    const ZERO: Self = 0;

    #[inline]
    fn gcd(self, other: Self) -> Self {
        self.unsigned_abs()
            .gcd(other.unsigned_abs())
            .try_into()
            .expect("GCD is `i32::MIN.unsigned_abs() == i32::MAX.unsigned_abs() + 1`.")
    }
}

impl Factor for u32 {
    const ZERO: Self = 0;

    #[allow(clippy::many_single_char_names, clippy::debug_assert_with_mut_call)]
    fn gcd(self, other: Self) -> Self {
        let mut a = self;
        let mut b = other;
        let g = if a == 0 || b == 0 {
            a | b
        } else {
            let mut u = a.trailing_zeros();
            let v = b.trailing_zeros();
            let x = min(u, v);
            b >>= v;
            while a != 0 {
                a >>= u;
                let d = b.abs_diff(a);
                u = d.trailing_zeros();
                b = min(a, b);
                a = d;
            }
            b << x
        };
        debug_assert_eq!(g, {
            let mut a = self;
            let mut b = other;
            while b != 0 {
                a %= b;
                swap(&mut a, &mut b);
            }
            a
        });
        g
    }
}

/// A geometric algebra defined by a flavor's basis (i.e., all its basis blades).
///
/// Implementations for this trait (e.g., [`pga::Pga`]) define the basis of a particular flavor from
/// which the Cayley table can be constructed. The generators of a basis (i.e., the basis blades
/// of <code>[Self::grade()] == 1</code>) define the whole algebra whereas the required [`Mul`]
/// operator implements the signature-aware anti-commutative multiplication for the chosen storage.
pub trait Algebra
where
    Self: Copy
        + Clone
        + Eq
        + PartialEq
        + Ord
        + PartialOrd
        + Default
        + Into<Symbol>
        + TryFrom<Symbol, Error = Symbol>
        + Debug
        + Display
        + Mul<Output = (i8, Self)>
        + Not<Output = (i8, Self)>,
{
    /// The embedded dimension.
    ///
    /// Not to be confused with the embedding dimension, e.g., `N + 1` is the embedding dimension
    /// for one-up flavors of embedded dimension `N`.
    const N: u32;

    /// The ordered basis (i.e., all basis blades).
    #[must_use]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self>;
    /// The scalar.
    #[must_use]
    fn scalar() -> Self;
    /// The pseudoscalar.
    #[must_use]
    fn pseudoscalar() -> Self;
    /// The grade.
    #[must_use]
    fn grade(&self) -> u32;
    /// The number $`n`$ of basis blades with the same grade $`g`$ of this basis blade.
    ///
    /// ```math
    /// n = { N + 1  \choose g }
    /// ```
    #[must_use]
    fn blade_len(&self) -> usize;
    /// The reverse.
    #[must_use]
    #[inline]
    fn rev(self) -> (i8, Self) {
        (1 - (self.grade().choose(2) & 1) as i8 * 2, self)
    }
}

/// Uniquely reduced form of a symbolic multivector expression.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
///
/// A multivector $`V`$ is realized as the sum of products of a Laurent [`Polynomial`] $`P_b`$ and a
/// basis blade $`\e_b`$ of an ordered basis `B`.
///
/// ```math
/// V \equiv \sum_b P_b \e_b
/// ```
///
/// All operators (e.g., [`Add`], [`Mul`]) implemented for [`Multivector`] reduce an arbitrary
/// expression into this unique form.
///
/// Generate text form with:
///
///   * `"{}"` for factorization of pinned symbols and GCDs,
///   * `"{:-}"` for factorization of pinned symbols and GCDs inclusive the predominant sign,
///   * `"{:+}"` for expanded form (i.e., no factorization),
///   * `"{:#}"` for alternative symbols labelled after basis blades,
///   * `"{:0}"` for zero newlines,
///   * `"{:<}"` for omitting plus signs.
///
/// Generate DOT form (i.e., [`text/vnd.graphviz`]) with:
///
///   * `"{:o}"` for factorization of pinned symbols and GCDs,
///   * `"{:-o}"` for factorization of pinned symbols and GCDs inclusive the predominant sign,
///   * `"{:+o}"` for expanded form (i.e., no factorization),
///   * `"{:#o}"` for alternative symbols labelled after basis blades,
///   * `"{:0o}"` for left-to-right rank direction.
///
/// [`text/vnd.graphviz`]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Multivector<B: Algebra> {
    /// Symbolic storage.
    pub map: BTreeMap<B, Polynomial>,
    /// Whether to leverage orthonormalization conditions (ONC).
    pub onc: bool,
}

impl<B: Algebra> TryFrom<Tree> for Multivector<B> {
    type Error = Symbol;

    fn try_from(tree: Tree) -> Result<Self, Self::Error> {
        match tree {
            Tree::Add(sib) => sib
                .into_iter()
                .try_fold(Self::zero(), |add, sib| Ok(add + Self::try_from(sib)?)),
            Tree::Mul(sib) => sib
                .into_iter()
                .try_fold(Self::one(), |mul, sib| Ok(mul * Self::try_from(sib)?)),
            Tree::Num(num) => Ok(num.into()),
            Tree::Sym(sym) => sym.try_into(),
        }
    }
}

impl<B: Algebra> From<Rational> for Multivector<B> {
    fn from(r: Rational) -> Self {
        if r.is_zero() {
            Self::zero()
        } else {
            let mut p = Polynomial::default();
            p.map.insert(Monomial::one(), r);
            let mut v = Self::default();
            v.map.insert(B::scalar(), p);
            v
        }
    }
}

impl<B: Algebra> TryFrom<Symbol> for Multivector<B> {
    type Error = Symbol;

    fn try_from(s: Symbol) -> Result<Self, Self::Error> {
        if s.is_vec() {
            Ok(Self::new([(Symbol::one(), s.try_into()?)]))
        } else {
            Ok(Self::new([(s, B::scalar())]))
        }
    }
}

impl<B: Algebra> Multivector<B> {
    /// Creates a new multivector from an iterator over tuples of symbols and basis blades.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee, pga::Pga};
    ///
    /// let plane = Vee::new([
    ///     (("W", "e0"), Pga::new("e0")),
    ///     (("x", "e1"), Pga::new("e1")),
    ///     (("y", "e2"), Pga::new("e2")),
    ///     (("z", "e3"), Pga::new("e3")),
    /// ]);
    ///
    /// assert_eq!(plane, Vee::plane());
    /// format_eq!(plane, ["+We0", "+xe1", "+ye2", "+ze3"]);
    /// ```
    #[must_use]
    #[inline]
    pub fn new<E, S>(iter: E) -> Self
    where
        E: IntoIterator<Item = (S, B)>,
        S: Into<Symbol>,
    {
        iter.into_iter().map(|(s, b)| ([[s]], b)).collect()
    }
    /// Appends Unicode *combining dot above* (i.e., `"◌̇"`) to all symbols.
    ///
    /// This is orthogonal to [`Self::cdm()`] extending the symbol space.
    #[must_use]
    pub fn alt(mut self) -> Self {
        self.map.values_mut().for_each(|p| *p = take(p).alt());
        self
    }
    /// Appends Unicode *combining x below* (i.e., `"◌͓"`) to all symbols.
    ///
    /// Pins this multivector as being sandwiched by the reflection or projection operator.
    ///
    /// Calls <code>[Self::cdm]\([Symbol::PIN]\)</code>.
    #[must_use]
    #[inline]
    pub fn pin(self) -> Self {
        self.cdm(Symbol::PIN)
    }
    /// Appends Unicode *combining left arrowhead below* (i.e., `"◌͔"`) to all symbols.
    ///
    /// Pins this multivector as left-hand side.
    ///
    /// Calls <code>[Self::cdm]\([Symbol::LHS]\)</code>.
    #[must_use]
    #[inline]
    pub fn lhs(self) -> Self {
        self.cdm(Symbol::LHS)
    }
    /// Appends Unicode *combining right arrowhead below* (i.e., `"◌͕"`) to all symbols.
    ///
    /// Pins this multivector as right-hand side.
    ///
    /// Calls <code>[Self::cdm]\([Symbol::RHS]\)</code>.
    #[must_use]
    #[inline]
    pub fn rhs(self) -> Self {
        self.cdm(Symbol::RHS)
    }
    /// Appends Unicode *combining diacritical mark* to all symbols.
    ///
    /// This example appends *combining double breve below* (i.e., `"◌͜◌"`) to plane $`p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// format_eq!(Vee::plane(), ["+We0", "+xe1", "+ye2", "+ze3"]);
    /// format_eq!(Vee::plane().cdm('\u{035c}'), ["+W͜e0", "+x͜e1", "+y͜e2", "+z͜e3"]);
    /// ```
    #[must_use]
    pub fn cdm(mut self, mark: char) -> Self {
        self.map.values_mut().for_each(|p| *p = take(p).cdm(mark));
        self
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(mut self) -> Self {
        self.map.values_mut().for_each(|p| *p = take(p).swp());
        self
    }
    /// Collects all grades.
    #[must_use]
    pub fn grades(&self) -> BTreeSet<u32> {
        self.map.keys().map(Algebra::grade).collect()
    }
    /// Returns the grade or `None` if it is of mixed-grade.
    #[must_use]
    pub fn grade(&self) -> Option<u32> {
        let grades = self.grades();
        (grades.len() == 1)
            .then(|| grades.first().copied())
            .flatten()
    }
    /// Whether being an entity (i.e., having unique symbols and exactly one per basis blade).
    #[must_use]
    pub fn is_entity(&self) -> bool {
        let mut set = BTreeSet::new();
        for p in self.map.values() {
            if let Some(m) = p.map.keys().next().filter(|_| p.map.len() == 1) {
                if let Some(s) = m.map.keys().next().filter(|_| m.map.len() == 1) {
                    if set.insert(s) {
                        continue;
                    }
                }
            }
            return false;
        }
        true
    }
    /// Collects the vectors per grade.
    #[must_use]
    pub fn vectors(/*mut*/ self) -> BTreeMap<u32, Self> {
        let mut vectors = BTreeMap::new();
        for grade in self.grades() {
            // let map = self
            // 	.map
            // 	.extract_if(.., |b, _p| vec.grade() == grade)
            // 	.collect();
            let mut map = self.map.clone();
            map.retain(|b, _p| b.grade() == grade);
            vectors.insert(grade, Self { map, onc: false });
        }
        vectors
    }
    /// Returns the vector of `grade`. The vector is empty if there is no `grade`.
    #[must_use]
    pub fn vector(mut self, grade: u32) -> Self {
        self.map.retain(|b, _p| grade == b.grade());
        self
    }
    /// Collects the basis blades.
    #[must_use]
    pub fn basis_blades(&self) -> BTreeSet<B> {
        self.map.keys().copied().collect()
    }
    /// The reverse.
    #[must_use]
    pub fn rev(mut self) -> Self {
        self.map.iter_mut().for_each(|(b, p)| {
            let (s, _b) = b.rev();
            if s < 0 {
                *p = -take(p);
            }
        });
        self
    }
    /// The zero.
    #[must_use]
    #[inline]
    pub fn zero() -> Self {
        Self::default()
    }
    /// The one.
    #[must_use]
    pub fn one() -> Self {
        Self::new([(Symbol::one(), B::scalar())])
    }
    /// Evaluates each symbol `S` of map `M` as respective rational `R`.
    #[must_use]
    #[allow(clippy::missing_panics_doc)]
    pub fn eval<M, S, R>(mut self, map: M) -> Self
    where
        M: IntoIterator<Item = (S, R)>,
        S: Into<Symbol>,
        R: Into<Rational>,
    {
        // Using inner non-generic function as monomorphization barrier.
        fn eval(vec: Vec<&mut Polynomial>, map: &BTreeMap<Symbol, Rational>) {
            for old_p in vec {
                let mut new_p = Polynomial::default();
                for (old_m, old_r) in &old_p.map {
                    let mut new_m = old_m.clone();
                    let mut new_r = *old_r;
                    for (old_s, old_e) in &old_m.map {
                        if let Some(map_r) = map.get(old_s) {
                            new_m.map.remove(old_s);
                            new_r *= map_r.pow(old_e.get());
                        }
                    }
                    if new_m.map.is_empty() {
                        new_m = Monomial::one();
                    }
                    new_p.map.insert(new_m, new_r);
                }
                *old_p = new_p;
            }
        }
        eval(
            self.map.values_mut().collect(),
            &map.into_iter().map(|(s, r)| (s.into(), r.into())).collect(),
        );
        self.map.retain(|_b, p| !p.map.is_empty());
        self
    }
    /// The polarity.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert_eq!(Vee::plane().pol(), Vee::direction().swp());
    /// ```
    #[must_use]
    #[inline]
    pub fn pol(self) -> Self {
        self * !Self::one()
    }
    /// The mixed-grade squared norm (i.e., a Study number).
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// format_eq!(Vee::plane().squared_norm(), ["+xx+yy+zz"]);
    /// format_eq!(Vee::point().squared_norm(), ["+ww"]);
    /// format_eq!(Vee::line().squared_norm(), ["+xx+yy+zz", "+2(-Xx-Yy-Zz)I"]);
    /// format_eq!(Vee::displacement().squared_norm(), ["+xx+yy+zz"]);
    /// format_eq!(Vee::moment().squared_norm(), []);
    /// ```
    #[must_use]
    pub fn squared_norm(self) -> Self {
        self.clone() * self.rev()
    }
    /// Leverages orthonormalization conditions.
    ///
    /// Assumes <code>[Self::squared_norm]\(self\) == [Self::one()]</code>.
    #[must_use]
    pub const fn unit(mut self) -> Self {
        self.onc = true;
        self
    }
    /// Applies `lhs == rhs` condition to `self`.
    ///
    /// Factors the GCD and the predominant sign of `lhs`. The remaining polynomial is matched with
    /// each remaining polynomial of <code>self.map.into_values().map([Factorization::from])</code>.
    /// The matched polynomials are replaced with the `rhs` vector of the respective `lhs` vector.
    #[must_use]
    pub fn cond(mut self, lhs: &Self, rhs: &Self) -> Self {
        for (lhs_b, lhs_p) in lhs.map.clone() {
            let rhs_p = rhs.map.get(&lhs_b).cloned().unwrap_or_default();
            let lhs_g = lhs_p.signed_gcd();
            let lhs_p = lhs_p / lhs_g;
            for p in self.map.values_mut() {
                let mut f = Factorization::from(p.clone());
                for (p, _r) in f.map.values_mut() {
                    if p == &lhs_p {
                        *p = rhs_p.clone();
                    }
                }
                *p = f.into();
            }
        }
        self.map.retain(|_b, p| !p.map.is_empty());
        self
    }
    /// Returns the number of `(multiplications, additions)`.
    #[must_use]
    pub fn ops(&self) -> (usize, usize) {
        self.map.values().fold((0, 0), |(v_muls, v_adds), p| {
            let f = Factorization::from(p.clone());
            let (f_muls, f_adds) =
                f.map
                    .into_iter()
                    .fold((0, 0), |(mut f_muls, f_adds), (m, (p, r))| {
                        if !m.is_one() {
                            f_muls += 1;
                        }
                        if !r.abs().is_one() {
                            f_muls += 1;
                        }
                        let (p_muls, p_adds) = p.ops();
                        (f_muls + p_muls, f_adds + p_adds)
                    });
            (v_muls + f_muls, v_adds + f_adds)
        })
    }
}

impl<B: Algebra, P, M, S> FromIterator<(P, B)> for Multivector<B>
where
    P: IntoIterator<Item = M>,
    M: IntoIterator<Item = S>,
    S: Into<Symbol>,
{
    fn from_iter<V: IntoIterator<Item = (P, B)>>(iter: V) -> Self {
        Self {
            map: iter
                .into_iter()
                .map(|(p, b)| (b, Polynomial::from_iter(p)))
                .collect(),
            onc: false,
        }
    }
}

impl<B: Algebra> Add for Multivector<B> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: Self) -> Self::Output {
        self += other;
        self
    }
}

impl<B: Algebra> AddAssign for Multivector<B> {
    fn add_assign(&mut self, other: Self) {
        for (b, p) in other.map {
            *self.map.entry(b).or_default() += p;
        }
        self.map.retain(|_b, p| !p.map.is_empty());
    }
}

impl<B: Algebra> Sub for Multivector<B> {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: Self) -> Self::Output {
        self -= other;
        self
    }
}

impl<B: Algebra> SubAssign for Multivector<B> {
    fn sub_assign(&mut self, other: Self) {
        for (b, p) in other.map {
            *self.map.entry(b).or_default() -= p;
        }
        self.map.retain(|_b, p| !p.map.is_empty());
    }
}

impl<B: Algebra> Neg for Multivector<B> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.map.values_mut().for_each(|p| *p = -take(p));
        self
    }
}

impl<B: Algebra> Mul for Multivector<B> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        let mut map = BTreeMap::<B, Polynomial>::new();
        for (&lhs_b, lhs_p) in &self.map {
            for (&rhs_b, rhs_p) in &other.map {
                let (s, b) = lhs_b * rhs_b;
                let p = lhs_p.clone() * rhs_p.clone();
                *map.entry(b).or_default() += p * i32::from(s);
            }
        }
        map.retain(|_b, p| !p.map.is_empty());
        Self { map, onc: false }
    }
}

impl<B: Algebra> MulAssign for Multivector<B> {
    fn mul_assign(&mut self, other: Self) {
        *self = take(self) * other;
    }
}

impl<B: Algebra> Mul<i32> for Multivector<B> {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: i32) -> Self::Output {
        self *= other;
        self
    }
}

impl<B: Algebra> MulAssign<i32> for Multivector<B> {
    fn mul_assign(&mut self, other: i32) {
        if other == 0 {
            self.map = BTreeMap::default();
        } else {
            self.map.values_mut().for_each(|p| *p *= other);
        }
    }
}

impl<B: Algebra> Div<i32> for Multivector<B> {
    type Output = Self;

    #[inline]
    fn div(mut self, other: i32) -> Self::Output {
        self /= other;
        self
    }
}

impl<B: Algebra> DivAssign<i32> for Multivector<B> {
    fn div_assign(&mut self, other: i32) {
        assert_ne!(other, 0, "division by zero");
        self.map.values_mut().for_each(|p| *p /= other);
    }
}

impl<B: Algebra> BitOr for Multivector<B> {
    type Output = Self;

    fn bitor(self, other: Self) -> Self::Output {
        let mut mv = Self::default();
        for (lhs_grade, lhs_vector) in self.vectors() {
            for (rhs_grade, rhs_vector) in other.clone().vectors() {
                mv += (lhs_vector.clone() * rhs_vector).vector(rhs_grade.abs_diff(lhs_grade));
            }
        }
        mv
    }
}

impl<B: Algebra> BitXor for Multivector<B> {
    type Output = Self;

    fn bitxor(self, other: Self) -> Self::Output {
        let mut mv = Self::default();
        for (lhs_grade, lhs_vector) in self.vectors() {
            for (rhs_grade, rhs_vector) in other.clone().vectors() {
                mv += (lhs_vector.clone() * rhs_vector).vector(lhs_grade + rhs_grade);
            }
        }
        mv
    }
}

impl<B: Algebra> Not for Multivector<B> {
    type Output = Self;

    fn not(self) -> Self::Output {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (b, p)| {
            let (s, b) = !b;
            map.insert(b, p * i32::from(s));
            map
        });
        Self { map, onc: false }
    }
}

impl<B: Algebra> BitAnd for Multivector<B> {
    type Output = Self;

    #[inline]
    fn bitand(self, other: Self) -> Self::Output {
        !!!(!self ^ !other)
    }
}

impl<B: Algebra> Rem for Multivector<B> {
    type Output = Self;

    fn rem(self, other: Self) -> Self::Output {
        (self.clone() * other.clone() - other * self) / 2
    }
}

impl<B: Algebra> Shl for Multivector<B> {
    type Output = Self;

    fn shl(mut self, other: Self) -> Self::Output {
        let onc = other.onc.then(|| {
            let lhs = other.clone().squared_norm();
            (self.clone() | (Self::one() - lhs.clone()), lhs, Self::one())
        });
        if self
            .grade()
            .zip(other.grade())
            .map(|(lhs_grade, rhs_grade)| (lhs_grade * rhs_grade) % 2)
            == Some(1)
        {
            self = -self;
        }
        let shl = other.clone() * self * other.rev();
        if let Some((onc, lhs, rhs)) = onc {
            (shl + onc).cond(&lhs, &rhs)
        } else {
            shl
        }
    }
}

impl<B: Algebra> Shr for Multivector<B> {
    type Output = Self;

    fn shr(self, other: Self) -> Self::Output {
        let shr = (self | other.clone()) * other.clone().rev();
        if other.onc {
            shr.cond(&other.squared_norm(), &Self::one())
        } else {
            shr
        }
    }
}

impl<B: Algebra> Display for Multivector<B> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fn traverse<'a>(
            fmt: &mut fmt::Formatter,
            tree: &Tree,
            depth: usize,
            grasp: bool,
            mut defer: &'a str,
        ) -> Result<&'a str, fmt::Error> {
            match tree {
                Tree::Add(siblings) => {
                    let grasp = grasp || !(defer.is_empty() || defer == "+");
                    if grasp {
                        write!(fmt, "{defer}")?;
                        defer = "";
                        write!(fmt, "(")?;
                    }
                    for (index, sibling) in siblings.iter().enumerate() {
                        if fmt.align().is_none() || index != 0 {
                            defer = "+";
                        }
                        defer = traverse(fmt, sibling, depth + 1, grasp, defer)?;
                        write!(fmt, "{defer}")?;
                    }
                    if grasp {
                        write!(fmt, ")")?;
                    }
                    defer = "";
                }
                Tree::Mul(siblings) => {
                    let mut is_one = false;
                    if let Some(Tree::Sym(sym)) = siblings.last() {
                        if (sym.is_one() || sym.is_scalar()) && depth <= 1 && siblings.len() == 2 {
                            is_one = true;
                        }
                    }
                    for (index, sibling) in siblings.iter().enumerate() {
                        defer = if index == 0 {
                            defer
                        } else if fmt.alternate() && defer.is_empty() && !is_one {
                            "*"
                        } else {
                            ""
                        };
                        let grasp = !(index == 0 && is_one);
                        defer = traverse(fmt, sibling, depth + 1, grasp, defer)?;
                    }
                    defer = "";
                }
                Tree::Num(num) => {
                    if num.abs().is_one() {
                        if num.is_negative() {
                            write!(fmt, "-")?;
                        } else if fmt.align().is_none() {
                            write!(fmt, "+")?;
                        };
                        defer = "1";
                    } else if num.is_negative() {
                        write!(fmt, "{num}")?;
                        defer = "";
                    } else if !num.is_zero() {
                        write!(fmt, "{defer}{num}")?;
                        defer = "";
                    }
                }
                Tree::Sym(sym) => {
                    write!(fmt, "{defer}")?;
                    if !(sym.is_one() || sym.is_scalar()) || depth == 0 {
                        Display::fmt(sym, fmt)?;
                        defer = "";
                    }
                    if !fmt.sign_aware_zero_pad() && sym.is_vec() {
                        writeln!(fmt)?;
                    }
                }
            }
            if depth == 0 {
                if defer == "1" {
                    write!(fmt, "{defer}")?;
                }
                defer = "";
            }
            Ok(defer)
        }
        let tree = if fmt.sign_plus() {
            Tree::from(self.clone())
        } else {
            Tree::with_factorization(self.clone(), fmt.sign_minus())
        };
        let defer = if fmt.align().is_none() { "+" } else { "" };
        traverse(fmt, &tree, 0, false, defer)?;
        Ok(())
    }
}

impl<B: Algebra> Octal for Multivector<B> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fn traverse(
            fmt: &mut fmt::Formatter,
            tree: &Tree,
            inode: usize,
            mut index: usize,
        ) -> Result<usize, fmt::Error> {
            let [add, mul] = if fmt.alternate() {
                ["+", "*"]
            } else {
                ["∑", "∏"]
            };
            match tree {
                Tree::Add(siblings) => {
                    writeln!(fmt, "  n{index} [label=\"{add}\" shape=box]")?;
                    if inode != index {
                        writeln!(fmt, "  n{inode} -> n{index}")?;
                    }
                    let inode = index;
                    for sibling in siblings {
                        index = traverse(fmt, sibling, inode, index + 1)?;
                    }
                    Ok(index)
                }
                Tree::Mul(siblings) => {
                    writeln!(fmt, "  n{index} [label=\"{mul}\" shape=box]")?;
                    if inode != index {
                        writeln!(fmt, "  n{inode} -> n{index}")?;
                    }
                    let inode = index;
                    for sibling in siblings {
                        index = traverse(fmt, sibling, inode, index + 1)?;
                    }
                    Ok(index)
                }
                Tree::Num(num) => {
                    writeln!(fmt, "  n{index} [label=\"{num}\" shape=circle]")?;
                    if inode != index {
                        writeln!(fmt, "  n{inode} -> n{index}")?;
                    }
                    Ok(index)
                }
                Tree::Sym(sym) => {
                    let shape = if sym.is_vec() {
                        "diamond"
                    } else {
                        match sym.cdm {
                            Symbol::NIL => "ellipse",
                            Symbol::PIN => "hexagon",
                            Symbol::LHS => "larrow",
                            Symbol::RHS => "rarrow",
                            label => panic!("unknown symbol label `{label}`"),
                        }
                    };
                    if fmt.alternate() {
                        writeln!(fmt, "  n{index} [label=\"{sym:#}\" shape={shape}]")?;
                    } else {
                        writeln!(fmt, "  n{index} [label=\"{sym}\" shape={shape}]")?;
                    }
                    if inode != index {
                        writeln!(fmt, "  n{inode} -> n{index}")?;
                    }
                    Ok(index)
                }
            }
        }
        writeln!(fmt, "digraph vee {{")?;
        if fmt.sign_aware_zero_pad() {
            writeln!(fmt, "  rankdir=LR")?;
        }
        let t = if fmt.sign_plus() {
            Tree::from(self.clone())
        } else {
            Tree::with_factorization(self.clone(), fmt.sign_minus())
        };
        traverse(fmt, &t, 0, 0)?;
        writeln!(fmt, "}}")?;
        Ok(())
    }
}

/// Uniquely reduced but **volatile** form of symbolic polynomial factorization.
///
/// Factors pinned symbols as monomials and factors the greatest common divisors (GCDs) of the
/// remaining polynomials and the GCD among them. Optionally, the GCDs are
/// [`signed`](`Factor::signed()`) comprising the factored predominant sign.
///
/// Initially, a factorization is uniquely reduced but in contrast to [`Polynomial`], the invariants
/// are no longer enforced by storage, making this form volatile. As the members are public, the
/// invariants are unguarded. For instance, manipulating a polynomial in [`Self::map`], potentially
/// invalidates the GCD.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Factorization {
    /// Symbolic storage mapping factored monomials to remaining polynomials and their GCDs.
    pub map: BTreeMap<Monomial, (Polynomial, Rational)>,
    /// GCD among remaining polynomials' GCDs.
    pub gcd: Rational,
}

impl Factorization {
    /// Performs the factorization.
    ///
    /// Optionally, the GCDs are `signed` comprising the factored predominant sign.
    #[must_use]
    pub fn new(p: Polynomial, signed: bool) -> Self {
        let mut f = p
            .map
            .into_iter()
            .fold(Self::default(), |mut f, (mut m, c)| {
                let mut g = m.clone();
                g.map.retain(|s, _e| s.is_pin());
                m.map.retain(|s, _e| !s.is_pin());
                // TODO
                // let g = Monomial {
                //     map: m.map.extract_if(|s, _e| s.is_pin()).collect(),
                // };
                if g.map.is_empty() {
                    g = Monomial::one();
                }
                if m.map.is_empty() {
                    m = Monomial::one();
                }
                let (p, r) = f.map.entry(g).or_default();
                *p.map.entry(m).or_default() += c;
                *r = Rational::ONE;
                f
            });
        let gcd = if signed {
            Polynomial::signed_gcd
        } else {
            Polynomial::gcd
        };
        f.map.values_mut().for_each(|(p, r)| {
            *r = gcd(p);
            *p /= *r;
        });
        f.gcd = if signed {
            Rational::signed(Rational::gcd, f.map.values().map(|(_p, r)| *r))
        } else {
            Rational::gcd_bulk(f.map.values().map(|(_p, r)| *r))
        };
        f.map.values_mut().for_each(|(_p, r)| *r /= f.gcd);
        f
    }
    /// Whether this factorization is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty()
            || self.gcd.is_zero()
            || self
                .map
                .iter()
                .all(|(m, (p, c))| m.is_zero() || p.is_zero() || c.is_zero())
    }
    /// Whether this factorization is one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        let mut sum = 0;
        self.gcd.is_one()
            && self
                .map
                .iter()
                .filter(|(m, (p, c))| !(m.is_zero() || p.is_zero() || c.is_zero()))
                .all(|(m, (p, c))| {
                    sum += 1;
                    m.is_one() && p.is_one() && c.is_one()
                })
            && sum == 1
    }
}

impl Default for Factorization {
    fn default() -> Self {
        Self {
            map: BTreeMap::new(),
            gcd: Rational::ONE,
        }
    }
}

impl From<Polynomial> for Factorization {
    #[inline]
    fn from(p: Polynomial) -> Self {
        Self::new(p, true)
    }
}

impl From<Factorization> for Polynomial {
    fn from(f: Factorization) -> Self {
        let mut p = Self::default();
        for (f_m, (f_p, f_r)) in f.map {
            for (m, c) in f_p.map {
                *p.map.entry(f_m.clone() * m).or_default() += c * f_r * f.gcd;
            }
        }
        p
    }
}

/// Uniquely reduced form of a symbolic polynomial expression.
///
/// A Laurent polynomial $`P_b`$ is realized as the sum of products of a [`Rational`] coefficient
/// $`C_m`$ and a primitive Laurent [`Monomial`] $`M_m`$ (i.e., an element of an ordered polynomial
/// basis).
///
/// ```math
/// P_b \equiv \sum_m C_m M_m
/// ```
///
/// All operators (e.g., [`Add`], [`Mul`]) implemented for [`Polynomial`] reduce an arbitrary
/// expression into this unique form.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Polynomial {
    /// Symbolic storage.
    pub map: BTreeMap<Monomial, Rational>,
}

impl Polynomial {
    /// The one.
    #[must_use]
    pub fn one() -> Self {
        let mut map = BTreeMap::new();
        map.insert(Monomial::one(), Rational::ONE);
        Self { map }
    }
    /// Whether this polynomial is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty() || self.map.iter().all(|(m, c)| m.is_zero() || c.is_zero())
    }
    /// Whether this polynomial is one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.map
            .first_key_value()
            .is_some_and(|(m, c)| self.map.len() == 1 && m.is_one() && c.is_one())
    }
    /// Extends the symbol space.
    #[must_use]
    pub fn alt(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, c)| {
            map.insert(s.alt(), c);
            map
        });
        Self { map }
    }
    /// Appends combining diacritical `mark` to all symbols.
    #[must_use]
    pub fn cdm(self, mark: char) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, c)| {
            map.insert(s.cdm(mark), c);
            map
        });
        Self { map }
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, c)| {
            map.insert(s.swp(), c);
            map
        });
        Self { map }
    }
    /// Returns the number of `(multiplications, additions)`.
    #[must_use]
    pub fn ops(&self) -> (usize, usize) {
        (
            self.map
                .iter()
                .map(|(m, r)| {
                    (m.map
                        .values()
                        .filter_map(|e| usize::try_from(e.get()).ok())
                        .sum::<usize>()
                        + usize::from(!r.abs().is_one()))
                    .saturating_sub(1)
                })
                .sum::<usize>(),
            self.map.len().saturating_sub(1),
        )
    }
    /// GCD of coefficients.
    #[must_use]
    pub fn gcd(&self) -> Rational {
        Rational::gcd_bulk(self.map.values().copied())
    }
    /// LCM of coefficients.
    #[must_use]
    pub fn lcm(&self) -> Rational {
        Rational::lcm_bulk(self.map.values().copied())
    }
    /// GCD and predominant sign of coefficients.
    pub fn signed_gcd(&self) -> Rational {
        Rational::signed(Rational::gcd, self.map.values().copied())
    }
    /// LCD and predominant sign of coefficients.
    pub fn signed_lcd(&self) -> Rational {
        Rational::signed(Rational::lcm, self.map.values().copied())
    }
}

impl<M, S> FromIterator<M> for Polynomial
where
    M: IntoIterator<Item = S>,
    S: Into<Symbol>,
{
    fn from_iter<P: IntoIterator<Item = M>>(iter: P) -> Self {
        Self {
            map: iter
                .into_iter()
                .map(|sym| (Monomial::from_iter(sym), Rational::ONE))
                .collect(),
        }
    }
}

impl Add for Polynomial {
    type Output = Self;

    #[inline]
    fn add(mut self, other: Self) -> Self::Output {
        self += other;
        self
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, other: Self) {
        for (m, c) in other.map {
            *self.map.entry(m).or_default() += c;
        }
        self.map.retain(|_m, c| !c.is_zero());
    }
}

impl Sub for Polynomial {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: Self) -> Self::Output {
        self -= other;
        self
    }
}

impl SubAssign for Polynomial {
    fn sub_assign(&mut self, other: Self) {
        for (m, c) in other.map {
            *self.map.entry(m).or_default() -= c;
        }
        self.map.retain(|_m, c| !c.is_zero());
    }
}

impl Neg for Polynomial {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.map.values_mut().for_each(|c| *c = -*c);
        self
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        let mut map = BTreeMap::new();
        for (lhs_m, lhs_c) in &self.map {
            for (rhs_s, rhs_c) in &other.map {
                let m = lhs_m.clone() * rhs_s.clone();
                let c = *lhs_c * *rhs_c;
                *map.entry(m).or_default() += c;
            }
        }
        map.retain(|_m, c: &mut Rational| !c.is_zero());
        Self { map }
    }
}

impl MulAssign for Polynomial {
    fn mul_assign(&mut self, other: Self) {
        *self = take(self) * other;
    }
}

impl Mul<Rational> for Polynomial {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: Rational) -> Self::Output {
        self *= other;
        self
    }
}

impl MulAssign<Rational> for Polynomial {
    fn mul_assign(&mut self, other: Rational) {
        if other.is_zero() {
            self.map = BTreeMap::default();
        } else {
            self.map.values_mut().for_each(|c| *c *= other);
        }
    }
}

impl Mul<i32> for Polynomial {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: i32) -> Self::Output {
        self *= other;
        self
    }
}

impl MulAssign<i32> for Polynomial {
    #[inline]
    fn mul_assign(&mut self, other: i32) {
        *self *= Rational::from(other);
    }
}

impl Div<Rational> for Polynomial {
    type Output = Self;

    #[inline]
    fn div(mut self, other: Rational) -> Self::Output {
        self /= other;
        self
    }
}

impl DivAssign<Rational> for Polynomial {
    fn div_assign(&mut self, other: Rational) {
        assert!(!other.is_zero(), "division by zero");
        self.map.values_mut().for_each(|c| *c /= other);
    }
}

impl Div<i32> for Polynomial {
    type Output = Self;

    #[inline]
    fn div(mut self, other: i32) -> Self::Output {
        self /= other;
        self
    }
}

impl DivAssign<i32> for Polynomial {
    #[inline]
    fn div_assign(&mut self, other: i32) {
        *self /= Rational::from(other);
    }
}

/// Rational number in canonical form.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rational {
    p: i32,
    q: i32,
}

impl Rational {
    /// The zero constant.
    pub const ZERO: Self = Self { p: 0, q: 1 };
    /// The one constant.
    pub const ONE: Self = Self { p: 1, q: 1 };

    /// Finds the irreducible fraction of numerator `p` and denominator `q`.
    ///
    /// # Panics
    ///
    /// Panics if denominator is zero.
    #[must_use]
    pub fn new(mut p: i32, mut q: i32) -> Self {
        if p == 0 {
            Self::ZERO
        } else if p == q {
            Self::ONE
        } else {
            let mut g = p.gcd(q);
            assert_ne!(g, 0, "division by zero");
            if q < 0 {
                g = -g;
            }
            p /= g;
            q /= g;
            Self { p, q }
        }
    }
    /// The numerator.
    #[must_use]
    #[inline]
    pub const fn p(&self) -> i32 {
        self.p
    }
    /// The denominator.
    #[must_use]
    #[inline]
    pub const fn q(&self) -> i32 {
        self.q
    }
    /// Inverts the fraction.
    ///
    /// # Panics
    ///
    /// Panics if numerator is zero.
    #[must_use]
    pub fn inv(&self) -> Self {
        assert_ne!(self.p, 0, "division by zero");
        if self.p < 0 {
            Self {
                p: -self.p,
                q: -self.q,
            }
        } else {
            Self {
                p: self.q,
                q: self.p,
            }
        }
    }
    /// The absolute.
    #[must_use]
    #[inline]
    pub const fn abs(&self) -> Self {
        Self {
            p: self.p.abs(),
            q: self.q,
        }
    }
    /// Whether this rational number is negative.
    #[must_use]
    #[inline]
    pub const fn is_negative(&self) -> bool {
        self.p.is_negative()
    }
    /// Whether this rational number is positive.
    #[must_use]
    #[inline]
    pub const fn is_positive(&self) -> bool {
        self.p.is_positive()
    }
    /// Whether this rational number is one.
    #[must_use]
    #[inline]
    pub const fn is_one(&self) -> bool {
        self.p == 1 && self.q == 1
    }
    /// Whether this rational number is zero.
    #[must_use]
    #[inline]
    pub const fn is_zero(&self) -> bool {
        self.p == 0
    }
    /// Power.
    #[must_use]
    pub fn pow(&self, exp: i32) -> Self {
        let abs = exp.unsigned_abs();
        let pow = Self::new(self.p().pow(abs), self.q().pow(abs));
        if exp.is_negative() { pow.inv() } else { pow }
    }
}

impl Default for Rational {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl From<i32> for Rational {
    #[inline]
    fn from(p: i32) -> Self {
        Self { p, q: 1 }
    }
}

impl From<(i32, i32)> for Rational {
    #[inline]
    fn from((p, q): (i32, i32)) -> Self {
        Self::new(p, q)
    }
}

impl Add for Rational {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(self.p * other.q + self.q * other.p, self.q * other.q)
    }
}

impl AddAssign for Rational {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Add<i32> for Rational {
    type Output = Self;

    #[inline]
    fn add(self, other: i32) -> Self {
        Self::new(self.p + self.q * other, self.q)
    }
}

impl AddAssign<i32> for Rational {
    #[inline]
    fn add_assign(&mut self, other: i32) {
        *self = *self + other;
    }
}

impl Add<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn add(self, other: Rational) -> Rational {
        Rational::new(self * other.q + other.p, other.q)
    }
}

impl Sub for Rational {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(self.p * other.q - self.q * other.p, self.q * other.q)
    }
}

impl SubAssign for Rational {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Rational {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            p: -self.p,
            q: self.q,
        }
    }
}

impl Sub<i32> for Rational {
    type Output = Self;

    #[inline]
    fn sub(self, other: i32) -> Self {
        Self::new(self.p - self.q * other, self.q)
    }
}

impl SubAssign<i32> for Rational {
    #[inline]
    fn sub_assign(&mut self, other: i32) {
        *self = *self - other;
    }
}

impl Sub<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn sub(self, other: Rational) -> Rational {
        Rational::new(self * other.q - other.p, other.q)
    }
}

impl Mul for Rational {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self {
        Self::new(self.p * other.p, self.q * other.q)
    }
}

impl MulAssign for Rational {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Mul<i32> for Rational {
    type Output = Self;

    #[inline]
    fn mul(self, other: i32) -> Self {
        Self::new(self.p * other, self.q)
    }
}

impl MulAssign<i32> for Rational {
    #[inline]
    fn mul_assign(&mut self, other: i32) {
        *self = *self * other;
    }
}

impl Mul<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn mul(self, other: Rational) -> Rational {
        Rational::new(self * other.p, other.q)
    }
}

impl Div for Rational {
    type Output = Self;

    #[inline]
    fn div(self, other: Self) -> Self {
        Self::new(self.p * other.q, self.q * other.p)
    }
}

impl DivAssign for Rational {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}

impl Div<i32> for Rational {
    type Output = Self;

    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: i32) -> Self {
        Self::new(self.p, self.q * other)
    }
}

impl DivAssign<i32> for Rational {
    #[inline]
    fn div_assign(&mut self, other: i32) {
        *self = *self / other;
    }
}

impl Div<Rational> for i32 {
    type Output = Rational;

    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: Rational) -> Rational {
        Rational::new(self * other.q, other.p)
    }
}

impl PartialOrd for Rational {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        (self.p * other.q).cmp(&(other.p * self.q))
    }
}

impl Display for Rational {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.p, fmt)?;
        if self.q != 1 {
            write!(fmt, "/{}", self.q)?;
        }
        Ok(())
    }
}

/// Uniquely reduced form of a symbolic monomial expression.
///
/// A primitive Laurent monomial $`M_m`$ is realized as the product of <code>[Symbol]s</code>
/// $`S_s`$ with individual non-zero exponents $`E_s`$.
///
/// ```math
/// M_m \equiv \prod_s S_s^{E_s}
/// ```
///
/// All operators (e.g., [`Add`], [`Mul`]) implemented for [`Monomial`] reduce an arbitrary
/// expression into this unique form.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Monomial {
    /// Symbolic storage.
    pub map: BTreeMap<Symbol, NonZeroI32>,
}

impl Monomial {
    /// The one.
    #[must_use]
    #[allow(clippy::missing_panics_doc)]
    pub fn one() -> Self {
        let mut map = BTreeMap::new();
        map.insert(Symbol::one(), NonZero::new(1).unwrap());
        Self { map }
    }
    /// Whether this monomial is zero.
    #[must_use]
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty()
    }
    /// Whether this monomial is one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.map
            .iter()
            .all(|(s, e)| s.is_one() && e.get().abs() == 1)
    }
    /// Extends the symbol space.
    #[must_use]
    pub fn alt(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, e)| {
            map.insert(s.alt(), e);
            map
        });
        Self { map }
    }
    /// Appends combining diacritical `mark` to all symbols.
    #[must_use]
    pub fn cdm(self, mark: char) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, e)| {
            map.insert(s.cdm(mark), e);
            map
        });
        Self { map }
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, e)| {
            map.insert(!s, e);
            map
        });
        Self { map }
    }
}

impl<S> FromIterator<S> for Monomial
where
    S: Into<Symbol>,
{
    fn from_iter<M: IntoIterator<Item = S>>(iter: M) -> Self {
        Self {
            map: iter
                .into_iter()
                .map(|s| (s.into(), NonZeroI32::new(1).unwrap()))
                .collect(),
        }
    }
}

impl Mul for Monomial {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: Self) -> Self::Output {
        self *= other;
        self
    }
}

impl MulAssign for Monomial {
    fn mul_assign(&mut self, other: Self) {
        for (s, rhs_e) in other.map {
            if let Some(lhs_e) = self.map.get(&s) {
                if let Some(e) = NonZeroI32::new(lhs_e.get() + rhs_e.get()) {
                    assert!(self.map.insert(s, e).is_some());
                } else {
                    assert!(self.map.remove(&s).is_some());
                }
            } else {
                assert!(self.map.insert(s, rhs_e).is_none());
            }
        }
        self.map.retain(|s, _e| !s.is_one());
    }
}

impl Div for Monomial {
    type Output = Self;

    #[inline]
    fn div(mut self, other: Self) -> Self::Output {
        self /= other;
        self
    }
}

impl DivAssign for Monomial {
    fn div_assign(&mut self, other: Self) {
        for (s, rhs_e) in other.map {
            if let Some(lhs_e) = self.map.get(&s) {
                if let Some(e) = NonZeroI32::new(lhs_e.get() - rhs_e.get()) {
                    assert!(self.map.insert(s, e).is_some());
                } else {
                    assert!(self.map.remove(&s).is_some());
                }
            } else {
                assert!(self.map.insert(s, rhs_e).is_none());
            }
        }
        self.map.retain(|s, _e| !s.is_one());
    }
}

/// Symbol as Unicode character with optional *combining diacritical mark*.
#[derive(Debug, Clone, Copy, Default)]
pub struct Symbol {
    var: char,
    alt: char,
    cdm: char,
    lab: &'static str,
}

impl PartialEq for Symbol {
    fn eq(&self, other: &Self) -> bool {
        self.var == other.var && self.alt == other.alt && self.cdm == other.cdm
    }
}

impl Eq for Symbol {}

impl Ord for Symbol {
    fn cmp(&self, other: &Self) -> Ordering {
        self.var
            .cmp(&other.var)
            .then(self.alt.cmp(&other.alt))
            .then(self.cdm.cmp(&other.cdm))
    }
}

impl PartialOrd for Symbol {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Symbol {
    /// Unicode null (i.e., `'\0'`).
    pub const NIL: char = '\0';
    /// Unicode zero-width space.
    pub const VEC: char = '\u{200b}';

    /// Unicode *combining dot above* (i.e., `"◌̇"`).
    pub const ALT: char = '\u{0307}';
    /// Unicode *combining x below* (i.e., `"◌͓"`).
    pub const PIN: char = '\u{0353}';
    /// Unicode *combining left arrowhead below* (i.e., `"◌͔"`).
    pub const LHS: char = '\u{0354}';
    /// Unicode *combining right arrowhead below* (i.e., `"◌͕"`).
    pub const RHS: char = '\u{0355}';

    /// Creates empty [`Default`] symbol representing scalar one.
    #[must_use]
    #[inline]
    pub const fn one() -> Self {
        Self {
            var: Self::NIL,
            alt: Self::NIL,
            cdm: Self::NIL,
            lab: "",
        }
    }
    /// Whether this symbol is a basis blade.
    #[must_use]
    #[inline]
    pub const fn is_vec(&self) -> bool {
        self.var == Self::VEC
    }
    /// Whether this symbol is the scalar.
    #[must_use]
    #[inline]
    pub const fn is_scalar(&self) -> bool {
        self.is_vec() && self.lab.len() == 1
    }
    /// Whether this symbol is the pseudoscalar.
    #[must_use]
    #[inline]
    pub const fn is_pseudoscalar(&self) -> bool {
        self.is_vec() && self.alt == Self::ALT
    }
    /// Whether this symbol is [`Self::one()`].
    #[must_use]
    #[inline]
    pub const fn is_one(&self) -> bool {
        self.var == Self::NIL
    }
    /// Whether this symbol is pinned.
    #[must_use]
    #[inline]
    pub const fn is_pin(&self) -> bool {
        self.cdm == Self::PIN
    }
    /// Whether this symbol is alternative.
    #[must_use]
    #[inline]
    pub const fn is_alt(&self) -> bool {
        self.alt == Self::ALT
    }
    /// Creates symbol for variable `var` with symbol name `sym`.
    #[must_use]
    #[inline]
    pub const fn new(var: char, sym: &'static str) -> Self {
        Self {
            var,
            alt: Self::NIL,
            cdm: Self::NIL,
            lab: sym,
        }
    }
    /// Marks this symbol with [`Self::ALT`].
    #[must_use]
    #[inline]
    pub const fn alt(mut self) -> Self {
        self.alt = Self::ALT;
        self
    }
    /// Marks this symbol with Unicode *combining diacritical mark*.
    #[must_use]
    #[inline]
    pub(crate) const fn cdm(mut self, cdm: char) -> Self {
        self.cdm = cdm;
        self
    }
    /// Marks this symbol with [`Self::PIN`].
    #[must_use]
    #[inline]
    pub const fn pin(self) -> Self {
        self.cdm(Self::PIN)
    }
    /// Marks this symbol with [`Self::LHS`] as left-hand side.
    #[must_use]
    #[inline]
    pub const fn lhs(self) -> Self {
        self.cdm(Self::LHS)
    }
    /// Marks this symbol with [`Self::RHS`] as right-hand side.
    #[must_use]
    #[inline]
    pub const fn rhs(self) -> Self {
        self.cdm(Self::RHS)
    }
}

impl From<(&str, &'static str)> for Symbol {
    #[inline]
    fn from((vars, sym): (&str, &'static str)) -> Self {
        let mut vars = vars.chars();
        let var = vars.next().unwrap_or_default();
        assert_eq!(vars.next(), None, "multi-character symbol");
        Self::new(var, sym)
    }
}

impl From<(char, &'static str)> for Symbol {
    #[inline]
    fn from((var, sym): (char, &'static str)) -> Self {
        Self::new(var, sym)
    }
}

impl Not for Symbol {
    type Output = Self;

    /// Swaps lowercase and uppercase character.
    fn not(self) -> Self {
        let var = if self.var.is_lowercase() {
            let mut iter = self.var.to_uppercase();
            assert_eq!(iter.len(), 1, "no uppercase for {}", self.var);
            iter.next().unwrap()
        } else {
            let mut iter = self.var.to_lowercase();
            assert_eq!(iter.len(), 1, "no lowercase for {}", self.var);
            iter.next().unwrap()
        };
        Self {
            var,
            alt: self.alt,
            cdm: self.cdm,
            lab: self.lab,
        }
    }
}

impl Display for Symbol {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        if self.is_vec() {
            if self.is_scalar() {
                write!(fmt, "1")?;
            } else if self.is_pseudoscalar() {
                write!(fmt, "I")?;
            } else {
                write!(fmt, "{}", &self.lab)?;
            }
        } else if fmt.alternate() {
            if !self.is_one() {
                write!(
                    fmt,
                    "{}{}",
                    match self.cdm {
                        Self::PIN => 'p',
                        Self::LHS => 'l',
                        Self::RHS => 'r',
                        _ => 'v',
                    },
                    &self.lab[1..]
                )?;
            }
        } else if !self.is_one() {
            write!(fmt, "{}", self.var)?;
            if self.alt != Self::NIL {
                write!(fmt, "{}", self.alt)?;
            }
            if self.cdm != Self::NIL {
                write!(fmt, "{}", self.cdm)?;
            }
        }

        Ok(())
    }
}

/// Non-binary algebraic expression tree up to symbolic [`Multivector`] expressions.
///
/// Unifies and simplifies following implementations due to being an intermediate and recursive data
/// structure:
///
///   * <code>impl [Display] for [Multivector]</code>
///   * <code>impl [Octal] for [Multivector]</code>
///
/// Convert from:
///
///   * <code>Tree::from([Multivector])</code> or <code>[Tree::with_factorization](v, true)</code>
///   * <code>Tree::from([Factorization])</code>
///   * <code>Tree::from([Polynomial])</code>
///   * <code>Tree::from([Monomial])</code>
///   * <code>Tree::from([Rational])</code>
///   * <code>Tree::from([Symbol])</code>
///
/// Convert to:
///
///   * <code>[Multivector]::try_from(tree)</code> where <code>Error = [Symbol]</code> indicates
///     that [`Symbol`] is not part of the [`Algebra`]. This is due [`Tree`] and [`Symbol`] being
///     non-generic and hence agnostic of the [`Algebra`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tree {
    /// Sum of subtrees.
    Add(Vec<Self>),
    /// Product of subtrees.
    Mul(Vec<Self>),
    /// Rational leaf node.
    Num(Rational),
    /// Symbolic leaf node.
    Sym(Symbol),
}

impl Tree {
    /// The zero constant.
    pub const ZERO: Self = Self::Num(Rational::ZERO);
    /// The one constant.
    pub const ONE: Self = Self::Num(Rational::ONE);

    /// Performs factorization on `v` and creates tree.
    ///
    /// Optionally, the GCDs are `signed` comprising the factored predominant sign.
    #[must_use]
    #[allow(clippy::missing_panics_doc)]
    pub fn with_factorization<B: Algebra>(v: Multivector<B>, signed: bool) -> Self {
        let mut add = Vec::with_capacity(v.map.len());
        for (b, p) in v.map {
            let f = Factorization::new(p, signed);
            if f.is_one() {
                add.push(b.into().into());
            } else {
                add.push(Self::Mul(vec![f.into(), b.into().into()]));
            }
        }
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.pop().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
}

impl Default for Tree {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl<B: Algebra> From<Multivector<B>> for Tree {
    fn from(v: Multivector<B>) -> Self {
        let add = v
            .map
            .into_iter()
            .filter(|(_b, p)| !p.is_zero())
            .map(|(b, p)| {
                if p.is_one() {
                    b.into().into()
                } else {
                    Self::Mul(vec![p.into(), b.into().into()])
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.into_iter().next().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
}

impl From<Factorization> for Tree {
    fn from(f: Factorization) -> Self {
        let add = f
            .map
            .into_iter()
            .filter(|(m, (p, c))| !(m.is_zero() || p.is_zero() || c.is_zero()))
            .map(|(m, (p, c))| {
                let is_mul = [c.is_one(), p.is_one(), m.is_one()].map(bool::not);
                let len = is_mul.into_iter().map(usize::from).sum::<usize>();
                let mut mul = Vec::with_capacity(len);
                if is_mul[0] {
                    mul.push(c.into());
                }
                if is_mul[1] {
                    mul.push(p.into());
                }
                if is_mul[2] {
                    mul.push(m.into());
                }
                if len == 1 {
                    mul.pop().expect("unreachable")
                } else {
                    Self::Mul(mul)
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            let any = add.into_iter().next().expect("unreachable");
            if f.gcd.is_one() {
                any
            } else if let Self::Num(num) = any {
                Self::Num(f.gcd * num)
            } else {
                Self::Mul(vec![f.gcd.into(), any])
            }
        } else if f.gcd.is_one() {
            Self::Add(add)
        } else {
            Self::Mul(vec![f.gcd.into(), Self::Add(add)])
        }
    }
}

impl From<Polynomial> for Tree {
    fn from(p: Polynomial) -> Self {
        let add = p
            .map
            .into_iter()
            .filter(|(p, c)| !(p.is_zero() || c.is_zero()))
            .map(|(m, c)| {
                if c.is_one() {
                    m.into()
                } else {
                    let m = Self::from(m);
                    if let Self::Mul(mul) = m {
                        let mut vec = vec![Self::from(c)];
                        vec.extend(mul);
                        Self::Mul(vec)
                    } else if let Self::Num(num) = m {
                        Self::Num(c * num)
                    } else {
                        Self::Mul(vec![c.into(), m])
                    }
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.into_iter().next().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
}

impl From<Monomial> for Tree {
    fn from(m: Monomial) -> Self {
        let mul = m
            .map
            .into_iter()
            .filter(|(s, _e)| !s.is_one())
            .flat_map(|(s, e)| {
                repeat(s)
                    .take(e.get().try_into().expect("negative exponent"))
                    .map(From::from)
            })
            .collect::<Vec<Self>>();
        if mul.is_empty() {
            Self::ONE
        } else if mul.len() == 1 {
            mul.into_iter().next().expect("unreachable")
        } else {
            Self::Mul(mul)
        }
    }
}

impl From<Rational> for Tree {
    #[inline]
    fn from(r: Rational) -> Self {
        Self::Num(r)
    }
}

impl From<Symbol> for Tree {
    #[inline]
    fn from(s: Symbol) -> Self {
        if s.is_one() { Self::ONE } else { Self::Sym(s) }
    }
}

pub mod pga;

/// Multivector for Elliptic 0D PGA.
pub type PgaE0 = Multivector<pga::PgaE0>;
/// Multivector for Elliptic 1D PGA.
pub type PgaE1 = Multivector<pga::PgaE1>;
/// Multivector for Elliptic 2D PGA.
pub type PgaE2 = Multivector<pga::PgaE2>;
/// Multivector for Elliptic 3D PGA.
pub type PgaE3 = Multivector<pga::PgaE3>;
/// Multivector for Elliptic 4D PGA.
pub type PgaE4 = Multivector<pga::PgaE4>;
/// Multivector for Elliptic 5D PGA (exploratory).
pub type PgaE5 = Multivector<pga::PgaE5>;
/// Multivector for Elliptic 6D PGA (exploratory).
pub type PgaE6 = Multivector<pga::PgaE6>;
/// Multivector for Elliptic 7D PGA (exploratory).
pub type PgaE7 = Multivector<pga::PgaE7>;

/// Multivector for Hyperbolic 0D PGA.
pub type PgaH0 = Multivector<pga::PgaH0>;
/// Multivector for Hyperbolic 1D PGA.
pub type PgaH1 = Multivector<pga::PgaH1>;
/// Multivector for Hyperbolic 2D PGA.
pub type PgaH2 = Multivector<pga::PgaH2>;
/// Multivector for Hyperbolic 3D PGA.
pub type PgaH3 = Multivector<pga::PgaH3>;
/// Multivector for Hyperbolic 4D PGA.
pub type PgaH4 = Multivector<pga::PgaH4>;
/// Multivector for Hyperbolic 5D PGA (exploratory).
pub type PgaH5 = Multivector<pga::PgaH5>;
/// Multivector for Hyperbolic 6D PGA (exploratory).
pub type PgaH6 = Multivector<pga::PgaH6>;
/// Multivector for Hyperbolic 7D PGA (exploratory).
pub type PgaH7 = Multivector<pga::PgaH7>;

/// Multivector for Parabolic (Euclidean) 0D PGA.
pub type PgaP0 = Multivector<pga::PgaP0>;
/// Multivector for Parabolic (Euclidean) 1D PGA.
pub type PgaP1 = Multivector<pga::PgaP1>;
/// Multivector for Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Multivector<pga::PgaP2>;
/// Multivector for Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Multivector<pga::PgaP3>;
/// Multivector for Parabolic (Euclidean) 4D PGA.
pub type PgaP4 = Multivector<pga::PgaP4>;
/// Multivector for Parabolic (Euclidean) 5D PGA (exploratory).
pub type PgaP5 = Multivector<pga::PgaP5>;
/// Multivector for Parabolic (Euclidean) 6D PGA (exploratory).
pub type PgaP6 = Multivector<pga::PgaP6>;
/// Multivector for Parabolic (Euclidean) 7D PGA (exploratory).
pub type PgaP7 = Multivector<pga::PgaP7>;
