// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Factor, Factorization, Monomial, Rational, Symbol};
use core::{
    mem::take,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use std::collections::BTreeMap;

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
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert_eq!(Polynomial::one() * Polynomial::one(), Polynomial::one());
    ///
    /// assert_eq!(
    ///     Polynomial::one() * Polynomial::default(),
    ///     Polynomial::default()
    /// );
    /// assert_eq!(
    ///     Polynomial::default() + Polynomial::default(),
    ///     Polynomial::default()
    /// );
    /// assert_eq!(Polynomial::one() + Polynomial::default(), Polynomial::one());
    /// ```
    #[must_use]
    pub fn one() -> Self {
        Self {
            map: BTreeMap::from([(Monomial::one(), Rational::ONE)]),
        }
    }
    /// Whether this polynomial is zero.
    ///
    /// ```
    /// use vee::Polynomial;
    ///
    /// assert!(Polynomial::default().is_zero());
    /// ```
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty() || self.map.iter().all(|(_m, c)| c.is_zero())
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
        Rational::gcd_reduce(self.map.values().copied())
    }
    /// LCM of coefficients.
    #[must_use]
    pub fn lcm(&self) -> Rational {
        Rational::lcm_reduce(self.map.values().copied())
    }
    /// GCD and predominant sign of coefficients.
    #[must_use]
    pub fn signed_gcd(&self) -> Rational {
        Rational::signed_gcd_reduce(self.map.values().copied())
    }
    /// LCD and predominant sign of coefficients.
    #[must_use]
    pub fn signed_lcd(&self) -> Rational {
        Rational::signed_lcm_reduce(self.map.values().copied())
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
