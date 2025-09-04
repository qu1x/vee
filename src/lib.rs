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
//! using [`rug`]). Currently, the Pistachio flavor -- Projective Geometric Algebra (PGA) -- is
//! implemented for $`D \equiv N + 1 \le 6 `$ in all three metrics, i.e., elliptic, hyperbolic, and
//! parabolic (Euclidean).[^2] The 5D PGA (i.e., $`N = 5`$) is incomplete as there is no inverse
//! based on Study numbers but it provides dimension-agnostic insights regarding duality and the
//! choice of basis blades.
//!
//! [^1]: S. De Keninck and M. Roelfs, “Normalization, square roots, and the exponential and
//! logarithmic maps in geometric algebras of less than 6D”, [Mathematical Methods in the Applied
//! Sciences 47, 1425–1441](https://doi.org/10.1002/mma.8639).
//! [^2]: M. Roelfs and S. De Keninck, “Graded Symmetry Groups: Plane and Simple”, [Advances in
//! Applied Clifford Algebras 33](https://doi.org/10.1007/s00006-023-01269-9).
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
//! Generates the expression for rotating a plane in [`PgaP3`], i.e., Parabolic (Euclidean) 3D PGA.
//! The [`Multivector::pin()`] method prefixes symbols of [`Multivector::plane()`] with `"~"` to
//! distinguish them from the symbols of [`Multivector::rotator()`].
//!
//! ```
//! use vee::PgaP3 as Vee;
//!
//! assert_eq!(format!("{:#}", Vee::plane().pin() << Vee::rotator()), concat!(
//!   "+(+[+1vv+1xx+1yy+1zz]~W)e0\n",
//!   "+(+[+2vz+2xy]~y+[-2vy+2xz]~z+[+1vv+1xx-1yy-1zz]~x)e1\n",
//!   "+(+[+2vx+2yz]~z+[-2vz+2xy]~x+[+1vv-1xx+1yy-1zz]~y)e2\n",
//!   "+(+[+2vy+2xz]~x+[-2vx+2yz]~y+[+1vv-1xx-1yy+1zz]~z)e3\n",
//! ));
//! ```
//!
//! The symbols are assigned to basis blades such that lowercase symbols are dual to their
//! corresponding uppercase symbols. For blades containing $`\e_0`$, uppercase symbols are used. The
//! [`Multivector::swp()`] method swaps lowercase and uppercase symbols. This is useful for testing
//! duality equivalences.
//!
//! ```
//! use vee::PgaP3 as Vee;
//!
//! assert_eq!(Vee::plane().to_string(), "We0+xe1+ye2+ze3");
//! assert_eq!(Vee::point().to_string(), "we123+Xe032+Ye013+Ze021");
//!
//! assert_ne!(!Vee::plane(), Vee::point());
//! assert_eq!(!Vee::plane(), Vee::point().swp());
//! ```

use change_case::swap_case;
use core::{
    fmt::{self, Debug, Display},
    iter::FromIterator,
    mem::take,
    num::NonZeroI32,
    ops::{
        Add, AddAssign, BitAnd, BitOr, BitXor, Div, DivAssign, Mul, MulAssign, Neg, Not, Rem, Shl,
        Sub, SubAssign,
    },
};
use num_rational::Ratio;
use smartstring::alias::String;
use std::collections::{BTreeMap, BTreeSet};

trait Choose {
    type Output;

    #[must_use]
    fn choose(self, other: Self) -> Self::Output;
}

impl Choose for usize {
    type Output = Self;

    #[inline]
    fn choose(self, other: Self) -> Self::Output {
        let (n, k) = (self, other);
        if n < k {
            0
        } else {
            (0..k).fold(1, |r, i| r * (n - i) / (i + 1))
        }
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
    const N: usize;

    /// The ordered basis (i.e., all basis blades).
    #[must_use]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self>;
    /// The grade.
    #[must_use]
    fn grade(&self) -> usize;
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
        let sig = if self.grade().choose(2) & 1 == 0 {
            1
        } else {
            -1
        };
        (sig, self)
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
    /// use vee::{PgaP3 as Vee, pga::Pga};
    ///
    /// let plane = Vee::new([
    ///     ("W", Pga::new("e0")),
    ///     ("x", Pga::new("e1")),
    ///     ("y", Pga::new("e2")),
    ///     ("z", Pga::new("e3")),
    /// ]);
    ///
    /// assert_eq!(plane, Vee::plane());
    /// assert_eq!(plane.to_string(), "We0+xe1+ye2+ze3");
    /// ```
    #[must_use]
    #[inline]
    pub fn new<E, S>(iter: E) -> Self
    where
        E: IntoIterator<Item = (S, B)>,
        S: Into<String>,
    {
        iter.into_iter().map(|(s, b)| ([[s]], b)).collect()
    }
    /// Adds prefix `"L"` to all symbols pinning this multivector as left-hand side.
    ///
    /// Calls <code>[Self::sym]\(\"L\"\)</code>.
    #[must_use]
    #[inline]
    pub fn lhs(self) -> Self {
        self.sym("L")
    }
    /// Adds prefix `"R"` to all symbols pinning this multivector as right-hand side.
    ///
    /// Calls <code>[Self::sym]\(\"R\"\)</code>.
    #[must_use]
    #[inline]
    pub fn rhs(self) -> Self {
        self.sym("R")
    }
    /// Adds prefix `"~"` to all symbols pinning this multivector as being sandwiched.
    ///
    /// Calls <code>[Self::sym]\(\"~\"\)</code>.
    #[must_use]
    #[inline]
    pub fn pin(self) -> Self {
        self.sym("~")
    }
    /// Adds `prefix` to all symbols.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert_eq!(Vee::plane().to_string(), "We0+xe1+ye2+ze3");
    /// assert_eq!(Vee::plane().sym("S").to_string(), "SWe0+Sxe1+Sye2+Sze3");
    /// ```
    #[must_use]
    pub fn sym(mut self, prefix: &str) -> Self {
        self.map.values_mut().for_each(|p| *p = take(p).sym(prefix));
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
    pub fn grades(&self) -> BTreeSet<usize> {
        self.map.keys().map(Algebra::grade).collect()
    }
    /// Returns the grade or `None` if it is of mixed-grade.
    #[must_use]
    pub fn grade(&self) -> Option<usize> {
        let grades = self.grades();
        (grades.len() == 1)
            .then(|| grades.first().copied())
            .flatten()
    }
    /// Collects the vectors per grade.
    #[must_use]
    pub fn vectors(/*mut*/ self) -> BTreeMap<usize, Self> {
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
    pub fn vector(mut self, grade: usize) -> Self {
        self.map.retain(|b, _p| grade == b.grade());
        self
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
    /// The mixed-grade squared norm (i.e., a Study number).
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert_eq!(Vee::plane().squared_norm().to_string(), "xx+yy+zz");
    /// assert_eq!(Vee::point().squared_norm().to_string(), "ww");
    /// assert_eq!(Vee::line().squared_norm().to_string(), "xx+yy+zz+(-Xx-Yy-Zz)I");
    /// assert_eq!(Vee::displacement().squared_norm().to_string(), "xx+yy+zz");
    /// assert_eq!(Vee::moment().squared_norm().to_string(), "");
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
    S: Into<String>,
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
        assert!(other != 0, "division by zero");
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

impl<B: Algebra> Display for Multivector<B> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (num, (b, p)) in self.map.iter().enumerate() {
            let is_one = *b == B::default();
            let is_sum = p.map.len() > 1;
            if (f.alternate() && !is_one && is_sum) || (!f.alternate() && num > 0) {
                write!(f, "+")?;
            }
            if !is_one && is_sum {
                write!(f, "(")?;
            }
            if self
                .map
                .values()
                .flat_map(|pol| pol.map.keys())
                .filter_map(|sym| sym.map.keys().last())
                .all(|sym| sym.starts_with('~'))
            {
                let mut map = BTreeMap::<_, Polynomial>::new();
                for (mut s, c) in p.map.clone() {
                    let pin = s.map.pop_last().unwrap();
                    assert!(map.entry(pin).or_default().map.insert(s, c).is_none());
                }
                let len = if map.len() % B::N == 0 {
                    B::N
                } else {
                    map.len()
                };
                for (num, ((s, e), p)) in map
                    .iter()
                    .take(len)
                    .cycle()
                    .skip(num)
                    .take(len)
                    .chain(map.iter().skip(len).cycle().skip(num).take(len))
                    .enumerate()
                {
                    if f.alternate() || num > 0 {
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
/// A Laurent polynomial $`P_b`$ is realized as the sum of products of a rational coefficient
/// $`C_m`$ and a primitive Laurent [`Monomial`] $`M_m`$ (i.e., an element of an
/// ordered polynomial basis).
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
    pub map: BTreeMap<Monomial, Ratio<i32>>,
}

impl Polynomial {
    /// Adds `prefix` to all symbols.
    #[must_use]
    #[inline]
    pub fn sym(self, prefix: &str) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (s, c)| {
            map.insert(s.sym(prefix), c);
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
    S: Into<String>,
{
    fn from_iter<P: IntoIterator<Item = M>>(iter: P) -> Self {
        Self {
            map: iter
                .into_iter()
                .map(|sym| (Monomial::from_iter(sym), Ratio::from_integer(1)))
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
            *self.map.entry(s).or_insert_with(|| Ratio::from_integer(0)) += c;
        }
        self.map.retain(|_s, c| *c.numer() != 0);
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
            *self.map.entry(s).or_insert_with(|| Ratio::from_integer(0)) -= c;
        }
        self.map.retain(|_s, c| *c.numer() != 0);
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
                let c = lhs_c * rhs_c;
                *map.entry(s).or_insert_with(|| Ratio::from_integer(0)) += c;
            }
        }
        map.retain(|_s, c| *c.numer() != 0);
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
        assert!(other != 0, "division by zero");
        self.map.values_mut().for_each(|c| *c /= other);
    }
}

impl Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut iter = self.map.iter();
        if !f.alternate() {
            if let Some((s, c)) = iter.next() {
                if !s.map.is_empty() || c.numer().abs() == 1 {
                    if *c.numer() < 1 {
                        write!(f, "-")?;
                    }
                } else {
                    write!(f, "{c:+}")?;
                }
                Display::fmt(s, f)?;
            }
        }
        iter.try_for_each(|(s, c)| {
            if !f.alternate() && (!s.map.is_empty() || c.numer().abs() == 1) {
                if *c.numer() < 1 {
                    write!(f, "-")?;
                } else {
                    write!(f, "+")?;
                }
            } else {
                write!(f, "{c:+}")?;
            }
            Display::fmt(s, f)?;
            Ok(())
        })?;
        Ok(())
    }
}

/// Uniquely reduced form of a symbolic monomial expression.
///
/// A primitive Laurent monomial $`M_m`$ is realized as the product of symbols $`S_s`$ with
/// individual non-zero exponents $`E_s`$ where [`String`] is an element of an ordered set of
/// multi-character symbols.
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
    pub map: BTreeMap<String, NonZeroI32>,
}

impl Monomial {
    /// Adds `prefix` to all symbols.
    #[must_use]
    #[inline]
    pub fn sym(self, prefix: &str) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (mut s, e)| {
            if !s.starts_with(prefix) {
                s.insert_str(0, prefix);
            }
            map.insert(s, e);
            map
        });
        Self { map }
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    #[inline]
    pub fn swp(self) -> Self {
        let map = BTreeMap::new();
        let map = self.map.into_iter().fold(map, |mut map, (mut s, e)| {
            s = swap_case(&s).into();
            map.insert(s, e);
            map
        });
        Self { map }
    }
}

impl<S> FromIterator<S> for Monomial
where
    S: Into<String>,
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
