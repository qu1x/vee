// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{
    Algebra, Factorization, Integer, NegAssign, Polynomial, Rational, Rev, RevAssign, Symbol, Tree,
};
use core::{
    fmt::{self, Alignment, Debug, Display, LowerHex, Octal},
    mem::replace,
    ops::{
        Add, AddAssign, BitAnd, BitOr, BitXor, Div, DivAssign, Mul, MulAssign, Neg, Not, Rem, Shl,
        Shr, Sub, SubAssign,
    },
};
use std::{
    collections::{BTreeMap, BTreeSet, btree_map::Entry},
    iter::{Product, Sum},
};

/// Uniquely reduced form of a symbolic multivector expression.
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
///   * `"{:0}"` for zero newlines (and no alignment environment in case of $`\LaTeX`$),
///   * `"{:.1}"` for floating points,
///   * `"{:<}"` for omitting plus signs,
///   * `"{:>}"` for omitting plus signs and surrounding operators with spaces,
///   * `"{:^}"` for dereferencing input and output fields implying `"{:>#}"`,
///   * `"{:$^}"` for $`\LaTeX`$ where the [`width`](std::fmt#width) parameter as in
///     `r"  \boldsymbol\ell = {:$^2}"` indents successive lines by additional two spaces,
///   * `"{:$>}"` for $`\LaTeX`$ omitting top alignment argument,
///   * `"{:$<}"` for $`\LaTeX`$ omitting environment begin and end.
///
/// Generate code form (i.e., generic statements and Rust) with:
///
///   * `"{:x}"` for factorization of pinned symbols and GCDs,
///   * `"{:-x}"` for factorization of pinned symbols and GCDs inclusive the predominant sign,
///   * `"{:+x}"` for expanded form (i.e., no factorization),
///   * `"{:#x}"` for Rust instead of generic statements,
///   * `"{:^x}"` for Rust dereferencing input and output fields.
///
/// where the [`width`](std::fmt#width) parameter as in `"{:^#4x}"` indents the code by four spaces.
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
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
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
    #[inline]
    fn from(q: Rational) -> Self {
        Self::from((B::scalar(), Polynomial::from(q)))
    }
}

impl<B: Algebra> TryFrom<Symbol> for Multivector<B> {
    type Error = Symbol;

    fn try_from(s: Symbol) -> Result<Self, Self::Error> {
        if s.is_vec() {
            Ok(Self {
                map: BTreeMap::from([(s.try_into()?, Polynomial::one())]),
                onc: false,
            })
        } else {
            Ok(Self::new([(s, B::scalar())]))
        }
    }
}

impl<B: Algebra> Multivector<B> {
    /// Creates a new multivector from `iter` over <code>([Into]<[Symbol]>, B)</code>.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq, pga::PgaP3 as Bee};
    ///
    /// let plane = Vee::new([
    ///     (('W', "e0"), Bee::new("e0").unwrap()),
    ///     (('x', "e1"), Bee::new("e1").unwrap()),
    ///     (('y', "e2"), Bee::new("e2").unwrap()),
    ///     (('z', "e3"), Bee::new("e3").unwrap()),
    /// ]);
    ///
    /// assert_eq!(
    ///     plane,
    ///     Vee::new([Bee::e0(), Bee::e1(), Bee::e2(), Bee::e3()]),
    /// );
    /// assert_eq!(plane, Vee::plane());
    /// format_eq!(plane, ["+We0", "+xe1", "+ye2", "+ze3"]);
    /// ```
    #[must_use]
    #[inline]
    pub fn new<I, S>(iter: I) -> Self
    where
        I: IntoIterator<Item = (S, B)>,
        S: Into<Symbol>,
    {
        iter.into_iter()
            .map(|(s, b)| (b, [([(s, Integer::ONE)], Rational::ONE)]))
            .collect()
    }
    /// Appends Unicode *combining dot above* (i.e., `"◌̇"`) to all symbols.
    ///
    /// This is orthogonal to [`Self::cdm()`] extending the symbol space.
    #[must_use]
    pub fn alt(mut self) -> Self {
        self.map.values_mut().for_each(Polynomial::alt_assign);
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
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// format_eq!(Vee::plane(), ["+We0", "+xe1", "+ye2", "+ze3"]);
    /// format_eq!(Vee::plane().cdm('\u{035c}'), ["+W͜e0", "+x͜e1", "+y͜e2", "+z͜e3"]);
    /// ```
    #[must_use]
    pub fn cdm(mut self, mark: char) -> Self {
        self.map.values_mut().for_each(|p| p.cdm_assign(mark));
        self
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(mut self) -> Self {
        self.map.values_mut().for_each(Polynomial::swp_assign);
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
    /// Whether all grades are odd.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert!(Vee::flector().is_odd());
    /// ```
    ///
    /// Mixed-parity multivectors are neither [`Self::is_odd()`] nor [`Self::is_even()`] whereas the
    /// zero multivector is both.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert!(Vee::flector().is_odd());
    ///
    /// let mixed_parity = Vee::motor() + Vee::flector();
    ///
    /// assert!(!mixed_parity.is_odd());
    /// assert!(!mixed_parity.is_even());
    ///
    /// let zero = Vee::zero();
    ///
    /// assert!(zero.is_odd());
    /// assert!(zero.is_even());
    /// ```
    #[must_use]
    pub fn is_odd(&self) -> bool {
        self.map.keys().map(B::grade).all(|grade| grade & 1 != 0)
    }
    /// Whether all grades are even.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// assert!(Vee::motor().is_even());
    /// ```
    ///
    /// See [`Self::is_odd()`] for mixed-parity and zero multivectors.
    #[must_use]
    pub fn is_even(&self) -> bool {
        self.map.keys().map(B::grade).all(|grade| grade & 1 == 0)
    }
    /// Whether being an entity (i.e., having unique symbols and exactly one per basis blade).
    #[must_use]
    pub fn is_entity(&self) -> bool {
        let mut set = BTreeSet::new();
        for p in self.map.values() {
            if let Some(m) = p.map.keys().next().filter(|_| p.map.len() == 1)
                && let Some(s) = m.map.keys().next().filter(|_| m.map.len() == 1)
                && set.insert(s)
            {
                continue;
            }
            return false;
        }
        true
    }
    /// Collects the vectors per grade.
    #[must_use]
    pub fn vectors(mut self) -> BTreeMap<u32, Self> {
        let mut vectors = BTreeMap::new();
        for grade in self.grades() {
            let map = self
                .map
                .extract_if(.., |b, _p| b.grade() == grade)
                .collect();
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
    pub fn rev(self) -> Self {
        Rev::rev(self)
    }
    /// The zero.
    #[must_use]
    #[inline]
    pub const fn zero() -> Self {
        Self {
            map: BTreeMap::new(),
            onc: false,
        }
    }
    /// The one.
    #[must_use]
    pub fn one() -> Self {
        Self::from(Rational::ONE)
    }
    /// Evaluates each symbol `S` of map `M` as respective rational `Q`.
    #[must_use]
    pub fn eval<M, S, Q>(mut self, map: M) -> Self
    where
        M: IntoIterator<Item = (S, Q)>,
        S: Into<Symbol>,
        Q: TryInto<Rational>,
    {
        // Using inner non-generic function as monomorphization barrier.
        fn eval(vec: Vec<&mut Polynomial>, map: &BTreeMap<Symbol, Option<Rational>>) {
            for old_p in vec {
                let mut new_p = Polynomial::zero();
                for (old_m, old_q) in &old_p.map {
                    let mut new_m = old_m.clone();
                    let mut new_q = Some(*old_q);
                    for (old_s, old_z) in &old_m.map {
                        if let Some(map_q) = map.get(old_s) {
                            new_m.map.remove(old_s);
                            let old_z = old_z
                                .get()
                                .try_into()
                                .expect("attempt to raise rational to power with overflow");
                            new_q = new_q
                                .zip(*map_q)
                                .map(|(new_q, map_q)| new_q * map_q.pow(old_z));
                        }
                    }
                    match new_p.map.entry(new_m) {
                        Entry::Occupied(mut entry) => {
                            match new_q.and_then(|new_q| *entry.get() + new_q) {
                                Some(q) => *entry.get_mut() = q,
                                None => {
                                    entry.remove();
                                }
                            }
                        }
                        Entry::Vacant(entry) => {
                            if let Some(new_q) = new_q {
                                entry.insert(new_q);
                            }
                        }
                    }
                }
                *old_p = new_p;
            }
        }
        eval(
            self.map.values_mut().collect(),
            &map.into_iter()
                .map(|(s, q)| (s.into(), q.try_into().ok()))
                .collect(),
        );
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
    /// The mixed-grade squared norm (i.e., a generalized complex number).
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// format_eq!(Vee::plane().norm_squared(), ["+xx+yy+zz"]);
    /// format_eq!(Vee::point().norm_squared(), ["+ww"]);
    /// format_eq!(Vee::line().norm_squared(), ["+xx+yy+zz", "+2(-Xx-Yy-Zz)I"]);
    /// format_eq!(Vee::displacement().norm_squared(), ["+xx+yy+zz"]);
    /// format_eq!(Vee::moment().norm_squared(), []);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm_squared(self) -> Self {
        self.clone() * self.rev()
    }
    /// Leverages orthonormalization conditions.
    ///
    /// Assumes <code>[Self::norm_squared]\(self\) == [Self::one()]</code>.
    #[must_use]
    #[inline]
    pub const fn unit(mut self) -> Self {
        self.onc = true;
        self
    }
    /// Applies `lhs == rhs` condition to `self`.
    ///
    ///  1. Factors the `gcd` inclusive the predominant sign of `lhs`.
    ///  2. The remaining polynomial is matched with each remaining polynomial of
    ///     <code>self.map.into_values().map([Factorization::from])</code>.
    ///  3. The matched polynomials are replaced with `rhs / gcd` of the respective `lhs` vector.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let l = || Vee::line().lhs();
    /// let r = || Vee::line().rhs();
    /// let mul = l() * r();
    /// let lhs = (l() | r()) * 2 + (l() % r()) * 4 + (l() ^ r()) * 8;
    /// let rhs = Vee::scalar() * 3 + Vee::line() * 5 + Vee::pseudoscalar() * 9;
    ///
    /// format_eq!("{:#x}", mul.cond(&lhs, &rhs), [
    ///     "let e = 3.0 / 2.0 * v;",
    ///     "let e01 = 5.0 / 4.0 * v01;",
    ///     "let e02 = 5.0 / 4.0 * v02;",
    ///     "let e03 = 5.0 / 4.0 * v03;",
    ///     "let e23 = 5.0 / 4.0 * v23;",
    ///     "let e31 = 5.0 / 4.0 * v31;",
    ///     "let e12 = 5.0 / 4.0 * v12;",
    ///     "let e0123 = 9.0 / 8.0 * v0123;",
    /// ]);
    /// ```
    ///
    /// See [`PgaP4::simple_motor()`] where [`Self::cond()`] is used in [`Self::shl`] whenever
    /// [`Self::onc`] is set with [`Self::unit()`] to eliminate orthonormalization conditions.
    ///
    /// [`PgaP4::simple_motor()`]: struct.Multivector.html#method.simple_motor-1
    #[must_use]
    pub fn cond(mut self, lhs: &Self, rhs: &Self) -> Self {
        for (lhs_b, lhs_p) in lhs.map.clone() {
            let rhs_p = rhs.map.get(&lhs_b).cloned();
            let lhs_g = lhs_p.signed_gcd().unwrap_or(Rational::ONE);
            let lhs_p = lhs_p / lhs_g;
            self.map.retain(|_b, p| {
                if p.is_zero() {
                    true
                } else {
                    let mut f = Factorization::from(p.clone());
                    f.map.retain(|_m, (p, _q)| {
                        if p == &lhs_p {
                            rhs_p.as_ref().is_some_and(|rhs_p| {
                                *p = rhs_p.clone() / lhs_g;
                                true
                            })
                        } else {
                            true
                        }
                    });
                    if f.is_zero() {
                        false
                    } else {
                        *p = f.into();
                        true
                    }
                }
            });
        }
        self
    }
    /// Omits zero vectors.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq, pga::PgaP3 as Bee};
    ///
    /// let norm = Vee::norm().eval([(Bee::e0123(), 0)]);
    ///
    /// format_eq!("{:#x}", norm, [
    ///     "let e = v;",
    ///     "let e0123 = 0.0;",
    /// ]);
    ///
    /// let simple_norm = norm.omit();
    ///
    /// format_eq!("{:#x}", simple_norm, [
    ///     "let e = v;",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn omit(mut self) -> Self {
        self.map.retain(|_b, p| !p.is_zero());
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
                    .fold((0, 0), |(mut f_muls, f_adds), (m, (p, q))| {
                        if !m.is_one() {
                            f_muls += 1;
                        }
                        if !q.abs().is_one() {
                            f_muls += 1;
                        }
                        let (p_muls, p_adds) = p.ops();
                        (f_muls + p_muls, f_adds + p_adds)
                    });
            (v_muls + f_muls, v_adds + f_adds)
        })
    }
}

impl<B, P> From<(B, P)> for Multivector<B>
where
    B: Algebra,
    P: Into<Polynomial>,
{
    #[inline]
    fn from((b, p): (B, P)) -> Self {
        Self {
            map: BTreeMap::from([(b, p.into())]),
            onc: false,
        }
    }
}

impl<B, P, M, Q, S, Z> FromIterator<(B, P)> for Multivector<B>
where
    B: Algebra,
    P: IntoIterator<Item = (M, Q)>,
    M: IntoIterator<Item = (S, Z)>,
    Q: TryInto<Rational>,
    S: Into<Symbol>,
    Z: TryInto<Integer>,
{
    #[inline]
    fn from_iter<I: IntoIterator<Item = (B, P)>>(iter: I) -> Self {
        let map = iter
            .into_iter()
            .map(|(b, p)| (b, Polynomial::from_iter(p)))
            .collect();
        Self { map, onc: false }
    }
}

impl<B: Algebra> Sum for Multivector<B> {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Self::add).unwrap_or_else(Self::zero)
    }
}

impl<B: Algebra> Product for Multivector<B> {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Self::mul).unwrap_or_else(Self::one)
    }
}

impl<'a, B: Algebra> Product<&'a Self> for Multivector<B> {
    #[inline]
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Self::mul)
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
            match self.map.entry(b) {
                Entry::Occupied(mut entry) => {
                    if let Some(p) = replace(entry.get_mut(), Polynomial::zero()) + p {
                        *entry.get_mut() = p;
                    } else {
                        entry.remove();
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(p);
                }
            }
        }
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
            match self.map.entry(b) {
                Entry::Occupied(mut entry) => {
                    if let Some(p) = replace(entry.get_mut(), Polynomial::zero()) - p {
                        *entry.get_mut() = p;
                    } else {
                        entry.remove();
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(-p);
                }
            }
        }
    }
}

impl<B: Algebra> Neg for Multivector<B> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self::Output {
        self.neg_assign();
        self
    }
}

impl<B: Algebra> NegAssign for Multivector<B> {
    #[inline]
    fn neg_assign(&mut self) {
        self.map.values_mut().for_each(Polynomial::neg_assign);
    }
}

impl<B: Algebra> Rev for Multivector<B> {
    type Output = Self;

    #[inline]
    fn rev(mut self) -> Self::Output {
        self.rev_assign();
        self
    }
}

impl<B: Algebra> RevAssign for Multivector<B> {
    #[inline]
    fn rev_assign(&mut self) {
        self.map.iter_mut().for_each(|(b, p)| {
            let (s, _b) = b.rev();
            if s < 0 {
                p.neg_assign();
            }
        });
    }
}

impl<B: Algebra> Mul for Multivector<B> {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl<B: Algebra> MulAssign for Multivector<B> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * &other;
    }
}

impl<B: Algebra> Mul<&Self> for Multivector<B> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<B: Algebra> MulAssign<&Self> for Multivector<B> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        *self = &*self * other;
    }
}

impl<B: Algebra> Mul for &Multivector<B> {
    type Output = Multivector<B>;

    fn mul(self, other: Self) -> Self::Output {
        let mut mul = Multivector::zero();
        for (&lhs_b, lhs_p) in &self.map {
            for (&rhs_b, rhs_p) in &other.map {
                let (s, b) = lhs_b * rhs_b;
                let q = Rational::new_integer(s.into());
                let Some(p) = q.map(|q| lhs_p.clone() * rhs_p.clone() * q) else {
                    continue;
                };
                match mul.map.entry(b) {
                    Entry::Occupied(mut entry) => {
                        if let Some(p) = replace(entry.get_mut(), Polynomial::zero()) + p {
                            *entry.get_mut() = p;
                        } else {
                            entry.remove();
                        }
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(p);
                    }
                }
            }
        }
        mul
    }
}

impl<B: Algebra> Mul<Multivector<B>> for &Multivector<B> {
    type Output = Multivector<B>;

    #[inline]
    fn mul(self, other: Multivector<B>) -> Self::Output {
        self * &other
    }
}

impl<B: Algebra, Q: TryInto<Rational>> Mul<Q> for Multivector<B> {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: Q) -> Self::Output {
        self *= other;
        self
    }
}

impl<B: Algebra, Q: TryInto<Rational>> MulAssign<Q> for Multivector<B> {
    #[inline]
    fn mul_assign(&mut self, other: Q) {
        if let Ok(other) = other.try_into() {
            self.map.values_mut().for_each(|p| *p *= other);
        } else {
            self.map.clear();
        }
    }
}

impl<B: Algebra, Q: TryInto<Rational>> Div<Q> for Multivector<B> {
    type Output = Self;

    #[inline]
    fn div(mut self, other: Q) -> Self::Output {
        self /= other;
        self
    }
}

impl<B: Algebra, Q: TryInto<Rational>> DivAssign<Q> for Multivector<B> {
    #[inline]
    fn div_assign(&mut self, other: Q) {
        let other = other
            .try_into()
            .unwrap_or_else(|_| panic!("attempt to divide multivector by zero"));
        self.map.values_mut().for_each(|p| *p /= other);
    }
}

impl<B: Algebra> BitOr for Multivector<B> {
    type Output = Self;

    fn bitor(self, other: Self) -> Self::Output {
        let mut v = Self::zero();
        for (lhs_grade, lhs_vector) in self.vectors() {
            for (rhs_grade, rhs_vector) in other.clone().vectors() {
                v += (lhs_vector.clone() * rhs_vector).vector(rhs_grade.abs_diff(lhs_grade));
            }
        }
        v
    }
}

impl<B: Algebra> BitXor for Multivector<B> {
    type Output = Self;

    fn bitxor(self, other: Self) -> Self::Output {
        let mut v = Self::zero();
        for (lhs_grade, lhs_vector) in self.vectors() {
            for (rhs_grade, rhs_vector) in other.clone().vectors() {
                v += (lhs_vector.clone() * rhs_vector).vector(lhs_grade + rhs_grade);
            }
        }
        v
    }
}

impl<B: Algebra> Not for Multivector<B> {
    type Output = Self;

    fn not(self) -> Self::Output {
        let map = self
            .map
            .into_iter()
            .map(|(b, p)| {
                let (s, b) = !b;
                (b, p * Rational::new_integer(s.into()).unwrap())
            })
            .collect();
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
        (self.clone() * other.clone() - other * self) / Rational::TWO
    }
}

impl<B: Algebra> Shl for Multivector<B> {
    type Output = Self;

    fn shl(mut self, other: Self) -> Self::Output {
        let onc = other.onc.then(|| {
            let lhs = other.clone().norm_squared();
            (self.clone() | (Self::one() - lhs.clone()), lhs, Self::one())
        });
        if self.is_odd() && other.is_odd() {
            self.neg_assign();
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
            shr.cond(&other.norm_squared(), &Self::one())
        } else {
            shr
        }
    }
}

impl<B: Algebra> Display for Multivector<B> {
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fn traverse<'a>(
            fmt: &mut fmt::Formatter,
            tree: &Tree,
            depth: usize,
            grasp: bool,
            group: bool,
            close: bool,
            mut defer: &'a str,
        ) -> Result<&'a str, fmt::Error> {
            let math = fmt.fill() == '$';
            let wide = matches!(fmt.align(), Some(Alignment::Center | Alignment::Right)) || math;
            let code = fmt.align() == Some(Alignment::Center) || fmt.alternate() || math;
            let rust = wide && !math && fmt.precision().is_some();
            match tree {
                Tree::Add(siblings) => {
                    if siblings.is_empty() && depth != 0 {
                        write!(fmt, "0")?;
                        if fmt.precision().is_some() {
                            write!(fmt, ".0")?;
                        }
                    } else {
                        let grasp = grasp || !(defer.is_empty() || defer == "+" || defer == " + ");
                        if grasp {
                            write!(fmt, "{defer}")?;
                            defer = "";
                            write!(fmt, "(")?;
                        }
                        for (index, sibling) in siblings.iter().enumerate() {
                            if fmt.align().is_none() || index != 0 {
                                defer = if wide { " + " } else { "+" };
                            }
                            if let Some(sym) = sibling.as_sym()
                                && sym.is_scalar()
                            {
                                defer = if fmt.precision().is_some() {
                                    if defer == "+" { "+1.0" } else { "1.0" }
                                } else {
                                    if defer == "+" { "+1" } else { "1" }
                                };
                            }
                            let close = index + 1 != siblings.len();
                            defer = traverse(fmt, sibling, depth + 1, grasp, false, close, defer)?;
                            write!(fmt, "{defer}")?;
                        }
                        if grasp {
                            write!(fmt, ")")?;
                        }
                    }
                    defer = "";
                }
                Tree::Mul(siblings) => {
                    let is_num = siblings
                        .first()
                        .and_then(Tree::as_num)
                        .is_some_and(|num| num.abs().is_one());
                    let is_one = siblings
                        .last()
                        .and_then(Tree::as_sym)
                        .is_some_and(|sym| (sym.is_scalar()) && depth <= 1 && siblings.len() == 2);
                    for (index, sibling) in siblings.iter().enumerate() {
                        let group = group || index > 0;
                        defer = if index == 0 {
                            defer
                        } else if code && defer.is_empty() && !is_one {
                            if math {
                                " "
                            } else if wide {
                                " * "
                            } else {
                                "*"
                            }
                        } else if is_num && is_one {
                            if fmt.precision().is_some() {
                                "1.0"
                            } else {
                                "1"
                            }
                        } else {
                            ""
                        };
                        let grasp = !(index == 0 && is_one);
                        defer = traverse(fmt, sibling, depth + 1, grasp, group, close, defer)?;
                    }
                    defer = "";
                }
                Tree::Num(num) => {
                    if num.is_negative() {
                        if group {
                            if code {
                                if wide {
                                    write!(fmt, " * ")?;
                                } else {
                                    write!(fmt, "*")?;
                                }
                            }
                            if !rust {
                                write!(fmt, "(")?;
                            }
                            write!(fmt, "-")?;
                        } else if wide && !defer.is_empty() {
                            write!(fmt, " - ")?;
                        } else {
                            write!(fmt, "-")?;
                        }

                        defer = "";
                    }
                    if !group && num.abs().is_one() {
                        if !num.is_negative() && fmt.align().is_none() {
                            write!(fmt, "+")?;
                        }
                        defer = if fmt.precision().is_some() {
                            "1.0"
                        } else {
                            "1"
                        };
                    } else {
                        write!(fmt, "{defer}")?;
                        Display::fmt(&num.abs(), fmt)?;
                        defer = "";
                    }
                    if num.is_negative() && group && !rust {
                        write!(fmt, ")")?;
                    }
                }
                Tree::Sym(sym) => {
                    write!(fmt, "{defer}")?;
                    if sym.is_vec() && math && !fmt.sign_aware_zero_pad() {
                        if sym.is_scalar() {
                            write!(fmt, " ")?;
                        }
                        write!(fmt, "& ")?;
                    }
                    defer = "";
                    if !sym.is_scalar() || depth == 0 {
                        Display::fmt(sym, fmt)?;
                    }
                    if !fmt.sign_aware_zero_pad() && sym.is_vec() {
                        if close && math {
                            if !sym.is_scalar() {
                                write!(fmt, " ")?;
                            }
                            write!(fmt, "\\\\\n ")?;
                            if let Some(width) = fmt.width() {
                                write!(fmt, "{:width$}", "")?;
                            }
                        } else {
                            writeln!(fmt)?;
                        }
                    }
                }
            }
            if depth == 0 {
                if defer == "1" || defer == "1.0" {
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
        let math = fmt.fill() == '$';
        let wide = matches!(fmt.align(), Some(Alignment::Center | Alignment::Right));
        if math && !fmt.sign_aware_zero_pad() && wide {
            let align = if fmt.align() == Some(Alignment::Center) {
                "[t]"
            } else {
                ""
            };
            writeln!(fmt, "\\begin{{aligned}}{align}")?;
        }
        if math && !fmt.sign_aware_zero_pad() {
            write!(fmt, "  ")?;
            if let Some(width) = fmt.width()
                && wide
            {
                write!(fmt, "{:width$}", "")?;
            }
        }
        traverse(fmt, &tree, 0, false, false, false, defer)?;
        if math && !fmt.sign_aware_zero_pad() && wide {
            if let Some(width) = fmt.width() {
                write!(fmt, "{:width$}", "")?;
            }
            writeln!(fmt, "\\end{{aligned}}")?;
        }
        Ok(())
    }
}

impl<B: Algebra> LowerHex for Multivector<B> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        for (b, p) in &self.map {
            let map = BTreeMap::from([(B::scalar(), p.clone())]);
            let v = Self { map, onc: self.onc };
            if let Some(width) = fmt.width() {
                write!(fmt, "{:width$}", "")?;
            }
            let deref = fmt.align() == Some(Alignment::Center);
            let rust = fmt.alternate() || deref;
            if fmt.sign_plus() {
                if rust {
                    if deref {
                        writeln!(fmt, "{b:^} = {v:^+0.1};")?;
                    } else {
                        writeln!(fmt, "let {b:#} = {v:>+#0.1};")?;
                    }
                } else {
                    writeln!(fmt, "{b:#}={v:<+#0}")?;
                }
            } else if fmt.sign_minus() {
                if rust {
                    if deref {
                        writeln!(fmt, "{b:^} = {v:^-0.1};")?;
                    } else {
                        writeln!(fmt, "let {b:#} = {v:>-#0.1};")?;
                    }
                } else {
                    writeln!(fmt, "{b:#}={v:<-#0}")?;
                }
            } else if rust {
                if deref {
                    writeln!(fmt, "{b:^} = {v:^0.1};")?;
                } else {
                    writeln!(fmt, "let {b:#} = {v:>#0.1};")?;
                }
            } else {
                writeln!(fmt, "{b:#}={v:<#0}")?;
            }
        }
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
