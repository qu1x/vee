// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Factorization, Integer, Monomial, NegAssign, Rational, Symbol};
use core::ops::{Add, Div, DivAssign, Mul, MulAssign, Neg, Sub};
use std::{
    collections::{BTreeMap, btree_map::Entry},
    iter::{Product, Sum},
};

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
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Polynomial {
    /// Symbolic storage.
    pub map: BTreeMap<Monomial, Rational>,
}

impl Polynomial {
    /// Creates a new polynomial from `iter` over <code>[Into]<[Symbol]></code>.
    #[must_use]
    #[inline]
    pub fn new<I>(iter: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<Symbol>,
    {
        iter.into_iter()
            .map(|s| ([(s, Integer::ONE)], Rational::ONE))
            .collect()
    }
    /// The zero.
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert_eq!(Polynomial::zero() * Polynomial::zero(), Polynomial::zero());
    /// assert_eq!(
    ///     Polynomial::zero() + Polynomial::zero(),
    ///     Some(Polynomial::zero())
    /// );
    /// assert_eq!(
    ///     Polynomial::zero() - Polynomial::zero(),
    ///     Some(Polynomial::zero())
    /// );
    /// ```
    #[must_use]
    #[inline]
    pub const fn zero() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }
    /// The one.
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert_eq!(Polynomial::one() * Polynomial::one(), Polynomial::one());
    ///
    /// assert_eq!(Polynomial::one() * Polynomial::zero(), Polynomial::zero());
    /// assert_eq!(
    ///     Polynomial::one() + Polynomial::zero(),
    ///     Some(Polynomial::one())
    /// );
    ///
    /// assert_eq!(
    ///     Polynomial::zero() + Polynomial::one(),
    ///     Some(Polynomial::one())
    /// );
    /// assert_eq!(
    ///     Polynomial::one() + Polynomial::zero(),
    ///     Some(Polynomial::one())
    /// );
    ///
    /// assert_eq!(
    ///     Polynomial::zero() - Polynomial::one(),
    ///     Some(-Polynomial::one())
    /// );
    /// assert_eq!(
    ///     Polynomial::one() - Polynomial::zero(),
    ///     Some(Polynomial::one())
    /// );
    /// ```
    #[must_use]
    #[inline]
    pub fn one() -> Self {
        Self::from((Monomial::one(), Rational::ONE))
    }
    /// Whether this polynomial is zero.
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert!(Polynomial::zero().is_zero());
    /// ```
    #[must_use]
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty()
    }
    /// Whether this polynomial is one.
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert!(Polynomial::one().is_one());
    /// ```
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.map
            .first_key_value()
            .is_some_and(|(m, q)| self.map.len() == 1 && m.is_one() && q.is_one())
    }
    /// Extends the symbol space.
    #[must_use]
    pub fn alt(&self) -> Self {
        Self {
            map: self.map.iter().map(|(m, &q)| (m.alt(), q)).collect(),
        }
    }
    #[inline]
    pub(crate) fn alt_assign(&mut self) {
        *self = self.alt();
    }
    /// Appends combining diacritical `mark` to all symbols.
    #[must_use]
    pub fn cdm(&self, mark: char) -> Self {
        Self {
            map: self.map.iter().map(|(m, &q)| (m.cdm(mark), q)).collect(),
        }
    }
    #[inline]
    pub(crate) fn cdm_assign(&mut self, mark: char) {
        *self = self.cdm(mark);
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(&self) -> Self {
        Self {
            map: self.map.iter().map(|(m, &q)| (m.swp(), q)).collect(),
        }
    }
    #[inline]
    pub(crate) fn swp_assign(&mut self) {
        *self = self.swp();
    }
    /// Returns the number of `(multiplications, additions)`.
    #[must_use]
    pub fn ops(&self) -> (usize, usize) {
        (
            self.map
                .iter()
                .map(|(m, q)| {
                    (m.map
                        .values()
                        .filter_map(|z| usize::try_from(z.get()).ok())
                        .sum::<usize>()
                        + usize::from(!q.abs().is_one()))
                    .saturating_sub(1)
                })
                .sum::<usize>(),
            self.map.len().saturating_sub(1),
        )
    }
    /// GCD of coefficients.
    #[must_use]
    pub fn gcd(&self) -> Option<Rational> {
        Rational::gcd_reduce(self.map.values().copied())
    }
    /// LCM of coefficients.
    #[must_use]
    pub fn lcm(&self) -> Option<Rational> {
        Rational::lcm_reduce(self.map.values().copied())
    }
    /// GCD and predominant sign of coefficients.
    #[must_use]
    pub fn signed_gcd(&self) -> Option<Rational> {
        Rational::signed_gcd_reduce(self.map.values().copied())
    }
    /// LCD and predominant sign of coefficients.
    #[must_use]
    pub fn signed_lcd(&self) -> Option<Rational> {
        Rational::signed_lcm_reduce(self.map.values().copied())
    }
}

impl From<Factorization> for Polynomial {
    fn from(f: Factorization) -> Self {
        let map = f
            .map
            .into_iter()
            .flat_map(|(f_m, (f_p, f_q))| {
                f_p.map
                    .into_iter()
                    .map(move |(m, q)| (f_m.clone() * m, q * f_q * f.gcd))
            })
            .collect();
        Self { map }
    }
}

impl From<Rational> for Polynomial {
    #[inline]
    fn from(q: Rational) -> Self {
        Self {
            map: BTreeMap::from([(Monomial::one(), q)]),
        }
    }
}

impl<M, Q> From<(M, Q)> for Polynomial
where
    M: Into<Monomial>,
    Q: TryInto<Rational>,
{
    #[inline]
    fn from((m, q): (M, Q)) -> Self {
        q.try_into().ok().map_or_else(Self::zero, |q| Self {
            map: BTreeMap::from([(m.into(), q)]),
        })
    }
}

impl<M, Q, S, Z> FromIterator<(M, Q)> for Polynomial
where
    M: IntoIterator<Item = (S, Z)>,
    Q: TryInto<Rational>,
    S: Into<Symbol>,
    Z: TryInto<Integer>,
{
    #[inline]
    fn from_iter<I: IntoIterator<Item = (M, Q)>>(iter: I) -> Self {
        let map = iter
            .into_iter()
            .filter_map(|(m, q)| q.try_into().ok().map(|q| (Monomial::from_iter(m), q)))
            .collect();
        Self { map }
    }
}

impl Sum for Polynomial {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |a, b| a.add(b).unwrap_or_else(Self::zero))
    }
}

impl Sum<(Monomial, Rational)> for Polynomial {
    #[inline]
    fn sum<I: Iterator<Item = (Monomial, Rational)>>(iter: I) -> Self {
        iter.fold(Self::zero(), |mut p, (m, q)| {
            p.add_entry(m, q);
            p
        })
    }
}

impl Product for Polynomial {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Self::mul).unwrap_or_else(Self::one)
    }
}

impl<'a> Product<&'a Self> for Polynomial {
    #[inline]
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Self::mul)
    }
}

impl Add for Polynomial {
    type Output = Option<Self>;

    fn add(mut self, other: Self) -> Self::Output {
        if self.is_zero() {
            Some(other)
        } else {
            for (m, q) in other.map {
                self.add_entry(m, q);
            }
            (!self.is_zero()).then_some(self)
        }
    }
}

impl Polynomial {
    #[inline]
    fn add_entry(&mut self, m: Monomial, q: Rational) {
        match self.map.entry(m) {
            Entry::Occupied(mut entry) => match *entry.get() + q {
                Some(sum) => *entry.get_mut() = sum,
                None => {
                    entry.remove();
                }
            },
            Entry::Vacant(entry) => {
                entry.insert(q);
            }
        }
    }
}

impl Sub for Polynomial {
    type Output = Option<Self>;

    fn sub(mut self, other: Self) -> Self::Output {
        if self.is_zero() {
            Some(-other)
        } else {
            for (m, q) in other.map {
                self.add_entry(m, -q);
            }
            (!self.is_zero()).then_some(self)
        }
    }
}

impl Neg for Polynomial {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self::Output {
        self.neg_assign();
        self
    }
}

impl NegAssign for Polynomial {
    #[inline]
    fn neg_assign(&mut self) {
        self.map.values_mut().for_each(Rational::neg_assign);
    }
}

impl Mul for Polynomial {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}

impl MulAssign for Polynomial {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * &other;
    }
}

impl Mul<&Self> for Polynomial {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl MulAssign<&Self> for Polynomial {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        *self = &*self * other;
    }
}

impl Mul for &Polynomial {
    type Output = Polynomial;

    fn mul(self, other: Self) -> Self::Output {
        self.map
            .iter()
            .flat_map(|(lhs_m, lhs_q)| {
                other
                    .map
                    .iter()
                    .map(move |(rhs_m, rhs_q)| (lhs_m.clone() * rhs_m, *lhs_q * *rhs_q))
            })
            .sum()
    }
}

impl Mul<Polynomial> for &Polynomial {
    type Output = Polynomial;

    #[inline]
    fn mul(self, other: Polynomial) -> Self::Output {
        self * &other
    }
}

impl<Q: TryInto<Rational>> Mul<Q> for Polynomial {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: Q) -> Self::Output {
        self *= other;
        self
    }
}

impl<Q: TryInto<Rational>> MulAssign<Q> for Polynomial {
    #[inline]
    fn mul_assign(&mut self, other: Q) {
        if let Ok(other) = other.try_into() {
            self.map.values_mut().for_each(|q| *q *= other);
        } else {
            self.map.clear();
        }
    }
}

impl<Q: TryInto<Rational>> Div<Q> for Polynomial {
    type Output = Self;

    #[inline]
    fn div(mut self, other: Q) -> Self::Output {
        self /= other;
        self
    }
}

impl<Q: TryInto<Rational>> DivAssign<Q> for Polynomial {
    #[inline]
    fn div_assign(&mut self, other: Q) {
        let other = other
            .try_into()
            .unwrap_or_else(|_| panic!("attempt to divide polynomial by zero"));
        self.map.values_mut().for_each(|q| *q /= other);
    }
}
