// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Monomial, Polynomial, Rational};
use core::ops::DivAssign;
use std::collections::{BTreeMap, btree_map::Entry};

/// Uniquely reduced but **volatile** form of symbolic polynomial factorization.
///
/// Factors pinned symbols as monomials and factors the greatest common divisors (GCDs) of the
/// remaining polynomials and the GCD among them. Optionally, the GCDs are
/// [`signed`](`Rational::signed_gcd_reduce()`) comprising the factored predominant sign.
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
        let mut f = p.map.into_iter().fold(Self::zero(), |mut f, (mut m, q)| {
            let f_m = Monomial {
                map: m.map.extract_if(.., |s, _z| s.is_pin()).collect(),
            };
            match f.map.entry(f_m) {
                Entry::Occupied(mut entry) => {
                    let (f_p, f_q) = entry.get_mut();
                    match f_p.map.entry(m) {
                        Entry::Occupied(mut entry) => match *entry.get() + q {
                            Some(q) => *entry.get_mut() = q,
                            None => {
                                entry.remove();
                            }
                        },
                        Entry::Vacant(entry) => {
                            entry.insert(q);
                        }
                    }
                    *f_q = Rational::ONE;
                }
                Entry::Vacant(entry) => {
                    let map = BTreeMap::from([(m, q)]);
                    entry.insert((Polynomial { map }, Rational::ONE));
                }
            }
            f
        });
        f.map.values_mut().for_each(|(p, q)| {
            *q = if signed || p.map.len() == 1 {
                p.signed_gcd()
            } else {
                p.gcd()
            }
            .unwrap_or(Rational::ONE);
            p.div_assign(*q);
        });
        f.gcd = if signed {
            Rational::signed_gcd_reduce(f.map.values().map(|(_p, q)| *q))
        } else {
            Rational::gcd_reduce(f.map.values().map(|(_p, q)| *q))
        }
        .unwrap_or(Rational::ONE);
        f.map.values_mut().for_each(|(_p, q)| q.div_assign(f.gcd));
        f
    }
    /// The zero.
    #[must_use]
    #[inline]
    pub const fn zero() -> Self {
        Self {
            map: BTreeMap::new(),
            gcd: Rational::ONE,
        }
    }
    /// The one.
    #[must_use]
    #[inline]
    pub fn one() -> Self {
        Self {
            map: BTreeMap::from([(Monomial::one(), (Polynomial::one(), Rational::ONE))]),
            gcd: Rational::ONE,
        }
    }
    /// Whether this factorization is zero.
    ///
    /// ```
    /// use vee::Factorization;
    ///
    /// assert!(Factorization::zero().is_zero());
    /// ```
    #[must_use]
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.map.is_empty()
    }
    /// Whether this factorization is one.
    ///
    /// ```
    /// use vee::Factorization;
    ///
    /// assert!(Factorization::one().is_one());
    /// ```
    #[must_use]
    pub fn is_one(&self) -> bool {
        !self.is_zero()
            && self.gcd.is_one()
            && self
                .map
                .iter()
                .all(|(m, (p, q))| m.is_one() && p.is_one() && q.is_one())
    }
}

impl From<Polynomial> for Factorization {
    #[inline]
    fn from(p: Polynomial) -> Self {
        Self::new(p, true)
    }
}
