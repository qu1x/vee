// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Integer, Inv, InvAssign, NegAssign, Symbol};
use core::ops::{Div, DivAssign, Mul, MulAssign};
use std::{
    collections::{BTreeMap, btree_map::Entry},
    iter::Product,
};

/// Uniquely reduced form of a symbolic monomial expression.
///
/// A primitive Laurent monomial $`M_m`$ is realized as the product of <code>[Symbol]s</code>
/// $`S_s`$ with individual non-zero exponents $`E_s`$.
///
/// ```math
/// M_m \equiv \prod_s S_s^{E_s}
/// ```
///
/// All operators (e.g., [`Mul`], [`Div`]) implemented for [`Monomial`] reduce an arbitrary
/// expression into this unique form.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Monomial {
    /// Symbolic storage.
    pub map: BTreeMap<Symbol, Integer>,
}

impl Monomial {
    /// Creates a new monomial from `iter` over <code>[Into]<[Symbol]></code>.
    #[must_use]
    #[inline]
    pub fn new<I>(iter: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<Symbol>,
    {
        iter.into_iter().map(|s| (s, Integer::ONE)).collect()
    }
    /// The one.
    ///
    /// ```
    /// use vee::Monomial;
    ///
    /// assert_eq!(Monomial::one() * Monomial::one(), Monomial::one());
    /// ```
    #[must_use]
    #[inline]
    pub const fn one() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }
    /// Whether this monomial is one.
    #[must_use]
    #[inline]
    pub fn is_one(&self) -> bool {
        self.map.is_empty()
    }
    /// Extends the symbol space.
    #[must_use]
    pub fn alt(&self) -> Self {
        Self {
            map: self.map.iter().map(|(s, &z)| (s.alt(), z)).collect(),
        }
    }
    /// Appends combining diacritical `mark` to all symbols.
    #[must_use]
    pub fn cdm(&self, mark: char) -> Self {
        Self {
            map: self.map.iter().map(|(s, &z)| (s.cdm(mark), z)).collect(),
        }
    }
    /// Swaps lowercase and uppercase symbols.
    #[must_use]
    pub fn swp(&self) -> Self {
        Self {
            map: self.map.iter().map(|(s, &z)| (!s, z)).collect(),
        }
    }
}

impl<S, Z> From<(S, Z)> for Monomial
where
    S: Into<Symbol>,
    Z: TryInto<Integer>,
{
    #[inline]
    fn from((s, z): (S, Z)) -> Self {
        z.try_into().ok().map_or_else(Self::one, |z| Self {
            map: BTreeMap::from([(s.into(), z)]),
        })
    }
}

impl<S, Z> FromIterator<(S, Z)> for Monomial
where
    S: Into<Symbol>,
    Z: TryInto<Integer>,
{
    #[inline]
    fn from_iter<I: IntoIterator<Item = (S, Z)>>(iter: I) -> Self {
        let map = iter
            .into_iter()
            .filter_map(|(s, z)| z.try_into().ok().map(|z| (s.into(), z)))
            .collect();
        Self { map }
    }
}

impl Product for Monomial {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(Self::mul).unwrap_or_else(Self::one)
    }
}

impl<'a> Product<&'a Self> for Monomial {
    #[inline]
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Self::mul)
    }
}

impl Product<(Symbol, Integer)> for Monomial {
    #[inline]
    fn product<I: Iterator<Item = (Symbol, Integer)>>(iter: I) -> Self {
        iter.fold(Self::one(), |mut m, (s, z)| {
            m.mul_entry(s, z);
            m
        })
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

impl Mul<&Self> for Monomial {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: &Self) -> Self::Output {
        self *= other;
        self
    }
}

impl Mul<Monomial> for &Monomial {
    type Output = Monomial;

    #[inline]
    fn mul(self, mut other: Monomial) -> Self::Output {
        other *= self;
        other
    }
}

impl MulAssign for Monomial {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl MulAssign<&Self> for Monomial {
    fn mul_assign(&mut self, other: &Self) {
        for (&s, &z) in &other.map {
            self.mul_entry(s, z);
        }
    }
}

impl Monomial {
    #[inline]
    fn mul_entry(&mut self, s: Symbol, z: Integer) {
        match self.map.entry(s) {
            #[allow(clippy::suspicious_op_assign_impl)]
            Entry::Occupied(mut entry) => {
                if let Some(z) = *entry.get() + z {
                    *entry.get_mut() = z;
                } else {
                    entry.remove();
                }
            }
            Entry::Vacant(entry) => {
                entry.insert(z);
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

impl Div<&Self> for Monomial {
    type Output = Self;

    #[inline]
    fn div(mut self, other: &Self) -> Self::Output {
        self /= other;
        self
    }
}

impl Div<Monomial> for &Monomial {
    type Output = Monomial;

    #[inline]
    fn div(self, mut other: Monomial) -> Self::Output {
        other /= self;
        other
    }
}

impl DivAssign for Monomial {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        *self /= &other;
    }
}

impl DivAssign<&Self> for Monomial {
    fn div_assign(&mut self, other: &Self) {
        for (&s, &z) in &other.map {
            match self.map.entry(s) {
                #[allow(clippy::suspicious_op_assign_impl)]
                Entry::Occupied(mut entry) => {
                    if let Some(z) = *entry.get() - z {
                        *entry.get_mut() = z;
                    } else {
                        entry.remove();
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(-z);
                }
            }
        }
    }
}

impl Inv for Monomial {
    type Output = Self;

    #[inline]
    fn inv(mut self) -> Self::Output {
        self.inv_assign();
        self
    }
}

impl InvAssign for Monomial {
    #[inline]
    fn inv_assign(&mut self) {
        self.map.values_mut().for_each(Integer::neg_assign);
    }
}
