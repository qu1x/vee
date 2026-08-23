// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![deny(clippy::arithmetic_side_effects)]

use core::ops::Neg;

pub trait Factor
where
    Self: Copy + Neg<Output = Self>,
{
    #[must_use]
    fn is_negative(self) -> bool;
    #[must_use]
    #[inline]
    fn is_positive(self) -> bool {
        !self.is_negative()
    }
    #[must_use]
    #[inline]
    fn abs(self) -> Self {
        if self.is_negative() { self.neg() } else { self }
    }
    #[must_use]
    fn gcd(self, other: Self) -> Self;
    #[must_use]
    fn lcm(self, other: Self) -> Self;
    #[must_use]
    #[inline]
    fn gcd_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        let mut iter = iter.into_iter().filter_map(|q| q.try_into().ok());
        iter.next().map(|q| iter.fold(q.abs(), Self::gcd))
    }
    #[must_use]
    #[inline]
    fn lcm_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        let mut iter = iter.into_iter().map(|q| q.try_into().ok());
        iter.next()
            .map(|q| iter.try_fold(q?.abs(), |l, q| q.map(|q| Self::lcm(l, q))))?
    }
    #[must_use]
    #[inline]
    fn signed_gcd_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        iter.into_iter()
            .filter_map(|q| q.try_into().ok())
            .map(|q| {
                (
                    q.abs(),
                    usize::from(q.is_negative()),
                    usize::from(q.is_positive()),
                )
            })
            .reduce(|(acc, neg, pos), (q, n, p)| {
                (Self::gcd(acc, q), neg.strict_add(n), pos.strict_add(p))
            })
            .map(|(acc, neg, pos)| if neg > pos { acc.neg() } else { acc })
    }
    #[must_use]
    #[inline]
    fn signed_lcm_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        let mut iter = iter.into_iter().map(|q| {
            q.try_into().ok().map(|q| {
                (
                    q.abs(),
                    usize::from(q.is_negative()),
                    usize::from(q.is_positive()),
                )
            })
        });
        iter.next()
            .map(|a| {
                iter.try_fold(a?, |(acc, neg, pos), q| {
                    q.map(|(q, n, p)| (Self::lcm(acc, q), neg.strict_add(n), pos.strict_add(p)))
                })
            })?
            .map(|(acc, neg, pos)| if neg > pos { acc.neg() } else { acc })
    }
}
