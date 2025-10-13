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
//! become a crate. Currently, the planed-based pistachio flavor -- Projective Geometric Algebra
//! (PGA) -- is implemented for $`D \equiv N + 1 \le 8`$ in all three metrics, i.e., elliptic,
//! hyperbolic, and parabolic (Euclidean).[^2] The 5D, 6D, and 7D PGAs (i.e., $`N = 5`$, $`N = 6`$,
//! and $`N = 7`$) are incomplete as there are no inverses based on Study numbers but they provide
//! dimension-agnostic insights regarding duality and the choice of basis blades. The PGA is
//! especially of interest for computer graphics, game engines, and physics simulations as it is the
//! most compact flavor (i.e., a one-up flavor) unifying the established but scattered frameworks,
//! e.g., homogeneous coordinates, Plücker coordinates, (dual) quaternions, and screw theory. Even
//! without any knowledge of geometric algebra, an API can be more intuitive as it unifies the
//! positional and directional aspects of geometric entities (e.g., planes, lines, points) and the
//! linear and angular aspects of rigid-body dynamics in a dimension-agnostic way with closed-form
//! (i.e., non-iterative) solutions up to 4D (e.g., [`PgaP2`], [`PgaP3`], [`PgaP4`]).[^3]
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
//! Generates the expression for rotating a plane in [`PgaP3`], i.e., the type alias of
//! [`Multivector`] parameterized for the Parabolic (Euclidean) 3D PGA. The [`PgaP3::pin()`] method
//! pins symbols of [`PgaP3::plane()`] with the *combining x below* (i.e., the Unicode *combining
//! diacritical mark* `"◌͓"`) to distinguish them from the symbols of
//! [`PgaP3::rotator()`](struct.Multivector.html#method.rotator-1).
//!
//! ```
//! use vee::{format_eq, PgaP3 as Vee};
//!
//! format_eq!(Vee::plane().pin() << Vee::rotator(), [
//!     "+(+[+vv+xx+yy+zz]W͓)e0",
//!     "+(+[+2vz+2xy]y͓+[-2vy+2xz]z͓+[+vv+xx-yy-zz]x͓)e1",
//!     "+(+[+2vx+2yz]z͓+[-2vz+2xy]x͓+[+vv-xx+yy-zz]y͓)e2",
//!     "+(+[+2vy+2xz]x͓+[-2vx+2yz]y͓+[+vv-xx-yy+zz]z͓)e3",
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
//! format_eq!(Vee::plane(), "We0+xe1+ye2+ze3");
//! format_eq!(Vee::point(), "we123+Xe032+Ye013+Ze021");
//!
//! assert_ne!(!Vee::plane(), Vee::point());
//! assert_eq!(!Vee::plane(), Vee::point().swp());
//! ```

/// Formats the `$left` expression using [`Display`] and asserts the `$right` string literal(s).
///
/// If `$right` is an array expression of string literals, formats `$left` in alternate form (i.e.,
/// `"{:#}"`) instead of default form (i.e., `"{}"`), appends `"\n"` to each `$right` literal, and
/// asserts the concatenation thereof.
///
/// With the `pretty_assertions` feature, the respective [`assert_eq!`] macro is used. In this way,
/// the Unicode *combining diacritical marks* are rendered as in the examples using [`format_eq!`].
#[macro_export]
macro_rules! format_eq {
    ($emitted:expr, $literal:literal) => {{
        #[cfg(feature = "pretty_assertions")]
        use pretty_assertions::assert_eq;
        let emitted = format!("{}", $emitted);
        assert_eq!(emitted, $literal);
    }};
    ($emitted:expr, [$($literal:literal),* $(,)?]) => {{
        #[cfg(feature = "pretty_assertions")]
        use pretty_assertions::assert_eq;
        let emitted = format!("{:#}", $emitted);
        let mut literal = String::with_capacity(emitted.len());
        $(
            literal.push_str($literal);
            literal.push_str("\n");
        )*
        assert_eq!(emitted, literal);
    }};
}

use core::{
    cmp::{Ordering, min},
    fmt::{self, Debug, Display},
    iter::FromIterator,
    mem::{swap, take},
    num::NonZeroI32,
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
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Multivector<B: Algebra> {
    /// Symbolic storage.
    pub map: BTreeMap<B, Polynomial>,
}

impl<B: Algebra> Multivector<B> {
    /// Creates a new multivector from an iterator over tuples of symbols and basis blades.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee, pga::Pga};
    ///
    /// let plane = Vee::new([
    ///     ("W", Pga::new("e0")),
    ///     ("x", Pga::new("e1")),
    ///     ("y", Pga::new("e2")),
    ///     ("z", Pga::new("e3")),
    /// ]);
    ///
    /// assert_eq!(plane, Vee::plane());
    /// format_eq!(plane, "We0+xe1+ye2+ze3");
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
    /// format_eq!(Vee::plane(), "We0+xe1+ye2+ze3");
    /// format_eq!(Vee::plane().cdm('\u{035c}'), "W͜e0+x͜e1+y͜e2+z͜e3");
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
            vectors.insert(grade, Self { map });
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
    /// The polarity.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert_eq!(Vee::plane().pol(), Vee::direction().swp());
    /// ```
    #[must_use]
    #[allow(clippy::missing_panics_doc)]
    pub fn pol(self) -> Self {
        self * Self::new([(Symbol::NIL, B::basis().next_back().expect("empty basis"))])
    }
    /// The mixed-grade squared norm (i.e., a Study number).
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// format_eq!(Vee::plane().squared_norm(), "xx+yy+zz");
    /// format_eq!(Vee::point().squared_norm(), "ww");
    /// format_eq!(Vee::line().squared_norm(), "xx+yy+zz+(-2Xx-2Yy-2Zz)I");
    /// format_eq!(Vee::displacement().squared_norm(), "xx+yy+zz");
    /// format_eq!(Vee::moment().squared_norm(), "");
    /// ```
    #[must_use]
    pub fn squared_norm(self) -> Self {
        self.clone() * self.rev()
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
        Self { map }
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
        Self { map }
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
        if self
            .grade()
            .zip(other.grade())
            .map(|(lhs_grade, rhs_grade)| (lhs_grade * rhs_grade) % 2)
            == Some(1)
        {
            self = -self;
        }
        other.clone() * self * other.rev()
    }
}

impl<B: Algebra> Shr for Multivector<B> {
    type Output = Self;

    fn shr(self, other: Self) -> Self::Output {
        (self | other.clone()) * other.rev()
    }
}

impl<B: Algebra> Display for Multivector<B> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, (b, p)) in self.map.iter().enumerate() {
            let is_one = *b == B::default();
            let is_sum = p.map.len() > 1;
            if (f.alternate() && !is_one && is_sum) || (!f.alternate() && i > 0) {
                write!(f, "+")?;
            }
            if !is_one && is_sum {
                write!(f, "(")?;
            }
            if self
                .map
                .values()
                .flat_map(|pol| pol.map.keys())
                .filter_map(|sym| {
                    sym.map
                        .keys()
                        .rfind(|sym| sym.is_pin())
                        .or_else(|| sym.map.keys().last())
                })
                .all(Symbol::is_pin)
            {
                let mut map = BTreeMap::<_, Polynomial>::new();
                for (mut s, c) in p.map.clone() {
                    let key = *s.map.keys().rfind(|sym| sym.is_pin()).unwrap();
                    let pin = s.map.remove_entry(&key).unwrap();
                    assert!(map.entry(pin).or_default().map.insert(s, c).is_none());
                }
                let len = if map.len() % B::N as usize == 0 {
                    B::N as usize
                } else {
                    map.len()
                };
                for (i, ((s, e), p)) in map
                    .iter()
                    .take(len)
                    .cycle()
                    .skip(i)
                    .take(len)
                    .chain(map.iter().skip(len).cycle().skip(i).take(len))
                    .enumerate()
                {
                    if f.alternate() || i > 0 {
                        write!(f, "+")?;
                    }
                    write!(f, "[")?;
                    Display::fmt(p, f)?;
                    write!(f, "]")?;
                    (0..e.get()).map(|_| s).try_for_each(|s| write!(f, "{s}"))?;
                }
            } else {
                Display::fmt(p, f)?;
            }
            if !is_one && is_sum {
                write!(f, ")")?;
            }
            if !is_one {
                Display::fmt(b, f)?;
            }
            if f.alternate() {
                writeln!(f)?;
            }
        }
        Ok(())
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
    /// Extends the symbol space.
    #[must_use]
    #[inline]
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
    #[inline]
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
    #[inline]
    pub fn swp(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, c)| {
            map.insert(s.swp(), c);
            map
        });
        Self { map }
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
        for (s, c) in other.map {
            *self.map.entry(s).or_default() += c;
        }
        self.map.retain(|_s, c| c.p() != 0);
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
        for (s, c) in other.map {
            *self.map.entry(s).or_default() -= c;
        }
        self.map.retain(|_s, c| c.p() != 0);
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
        for (lhs_s, lhs_c) in &self.map {
            for (rhs_s, rhs_c) in &other.map {
                let s = lhs_s.clone() * rhs_s.clone();
                let c = *lhs_c * *rhs_c;
                *map.entry(s).or_default() += c;
            }
        }
        map.retain(|_s, c: &mut Rational| c.p() != 0);
        Self { map }
    }
}

impl MulAssign for Polynomial {
    fn mul_assign(&mut self, other: Self) {
        *self = take(self) * other;
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
    fn mul_assign(&mut self, other: i32) {
        if other == 0 {
            self.map = BTreeMap::default();
        } else {
            self.map.values_mut().for_each(|c| *c *= other);
        }
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
    fn div_assign(&mut self, other: i32) {
        assert_ne!(other, 0, "division by zero");
        self.map.values_mut().for_each(|c| *c /= other);
    }
}

impl Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.map.iter().enumerate().try_for_each(|(i, (s, c))| {
            if !f.sign_plus() && !s.map.is_empty() && c.abs().is_one() {
                if c.is_negative() {
                    write!(f, "-")?;
                } else if f.alternate() || i > 0 {
                    write!(f, "+")?;
                }
            } else {
                write!(f, "{c:+}")?;
            }
            Display::fmt(s, f)
        })
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

    /// Finds the irreducable fraction of numerator `p` and denominator `q`.
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
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.p, f)?;
        if self.q != 1 {
            write!(f, "/{}", self.q)?;
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
    /// Extends the symbol space.
    #[must_use]
    #[inline]
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
    #[inline]
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
    #[inline]
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

impl Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.map
            .iter()
            .flat_map(|(s, e)| (0..e.get()).map(move |_| s))
            .try_for_each(|s| write!(f, "{s}"))
    }
}

/// Symbol as Unicode character with optional *combining diacritical mark*.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Symbol {
    var: char,
    alt: char,
    cdm: char,
}

impl Symbol {
    /// Unicode null (i.e., `'\0'`).
    pub const NIL: char = '\0';

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
        }
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
    /// Creates symbol for variable `var`.
    #[must_use]
    #[inline]
    pub const fn new(var: char) -> Self {
        Self {
            var,
            alt: Self::NIL,
            cdm: Self::NIL,
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

impl From<&str> for Symbol {
    #[inline]
    fn from(sym: &str) -> Self {
        let mut sym = sym.chars();
        let var = sym.next().unwrap_or_default();
        assert_eq!(sym.next(), None, "multi-character symbol");
        Self::new(var)
    }
}

impl From<char> for Symbol {
    #[inline]
    fn from(var: char) -> Self {
        Self::new(var)
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
        }
    }
}

impl Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.var != Self::NIL {
            write!(f, "{}", self.var)?;
            if self.alt != Self::NIL {
                write!(f, "{}", self.alt)?;
            }
            if self.cdm != Self::NIL {
                write!(f, "{}", self.cdm)?;
            }
        }
        Ok(())
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
/// Multivector for Elliptic 4D PGA (experimental).
pub type PgaE4 = Multivector<pga::PgaE4>;
/// Multivector for Elliptic 5D PGA (experimental, no inverse).
pub type PgaE5 = Multivector<pga::PgaE5>;
/// Multivector for Elliptic 6D PGA (experimental, no inverse).
pub type PgaE6 = Multivector<pga::PgaE6>;
/// Multivector for Elliptic 7D PGA (experimental, no inverse).
pub type PgaE7 = Multivector<pga::PgaE7>;

/// Multivector for Hyperbolic 0D PGA.
pub type PgaH0 = Multivector<pga::PgaH0>;
/// Multivector for Hyperbolic 1D PGA.
pub type PgaH1 = Multivector<pga::PgaH1>;
/// Multivector for Hyperbolic 2D PGA.
pub type PgaH2 = Multivector<pga::PgaH2>;
/// Multivector for Hyperbolic 3D PGA.
pub type PgaH3 = Multivector<pga::PgaH3>;
/// Multivector for Hyperbolic 4D PGA (experimental).
pub type PgaH4 = Multivector<pga::PgaH4>;
/// Multivector for Hyperbolic 5D PGA (experimental, no inverse).
pub type PgaH5 = Multivector<pga::PgaH5>;
/// Multivector for Hyperbolic 6D PGA (experimental, no inverse).
pub type PgaH6 = Multivector<pga::PgaH6>;
/// Multivector for Hyperbolic 7D PGA (experimental, no inverse).
pub type PgaH7 = Multivector<pga::PgaH7>;

/// Multivector for Parabolic (Euclidean) 0D PGA.
pub type PgaP0 = Multivector<pga::PgaP0>;
/// Multivector for Parabolic (Euclidean) 1D PGA.
pub type PgaP1 = Multivector<pga::PgaP1>;
/// Multivector for Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Multivector<pga::PgaP2>;
/// Multivector for Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Multivector<pga::PgaP3>;
/// Multivector for Parabolic (Euclidean) 4D PGA (experimental).
pub type PgaP4 = Multivector<pga::PgaP4>;
/// Multivector for Parabolic (Euclidean) 5D PGA (experimental, no inverse).
pub type PgaP5 = Multivector<pga::PgaP5>;
/// Multivector for Parabolic (Euclidean) 6D PGA (experimental, no inverse).
pub type PgaP6 = Multivector<pga::PgaP6>;
/// Multivector for Parabolic (Euclidean) 7D PGA (experimental, no inverse).
pub type PgaP7 = Multivector<pga::PgaP7>;
