// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Factor, Monomial, Polynomial, Rational};
use std::collections::BTreeMap;

/// Uniquely reduced but **volatile** form of symbolic polynomial factorization.
///
/// Factors pinned symbols as monomials and factors the greatest common divisors (GCDs) of the
/// remaining polynomials and the GCD among them. Optionally, the GCDs are
/// [`signed`](`Factor::signed_reduce()`) comprising the factored predominant sign.
///
/// Initially, a factorization is uniquely reduced but in contrast to [`Polynomial`], the invariants
/// are no longer enforced by storage, making this form volatile. As the members are public, the
/// invariants are unguarded. For instance, manipulating a polynomial in [`Self::map`], potentially
/// invalidates the [`Self::gcd`].
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
                let mut g = Monomial {
                    map: m.map.extract_if(.., |s, _e| s.is_pin()).collect(),
                };
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
        f.map.values_mut().for_each(|(p, r)| {
            *r = if signed || p.map.len() == 1 {
                p.signed_gcd()
            } else {
                p.gcd()
            };
            *p /= *r;
        });
        f.gcd = if signed {
            Rational::signed_gcd_reduce(f.map.values().map(|(_p, r)| *r))
        } else {
            Rational::gcd_reduce(f.map.values().map(|(_p, r)| *r))
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
                .all(|(_m, (p, c))| p.is_zero() || c.is_zero())
    }
    /// Whether this factorization is one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        let mut sum = 0;
        self.gcd.is_one()
            && self
                .map
                .iter()
                .filter(|(_m, (p, c))| !(p.is_zero() || c.is_zero()))
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
