// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::Symbol;
use core::{
    num::NonZeroI32,
    ops::{Div, DivAssign, Mul, MulAssign},
};
use std::collections::BTreeMap;

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
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Monomial {
    /// Symbolic storage.
    pub map: BTreeMap<Symbol, NonZeroI32>,
}

impl Monomial {
    /// The one.
    ///
    /// ```
    /// use vee::Monomial;
    ///
    /// assert_eq!(Monomial::one() * Monomial::one(), Monomial::one());
    /// ```
    #[must_use]
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub const fn one() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }
    /// The inverse.
    #[must_use]
    pub fn inv(mut self) -> Self {
        self.map.values_mut().for_each(|e| *e = -*e);
        self
    }
    /// Whether this monomial is one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.map.is_empty()
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
    #[inline]
    #[allow(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}
