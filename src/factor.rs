// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::Rational;
use core::{
    cmp::min,
    mem::swap,
    ops::{Div, Mul, Neg},
};

/// Finds the greatest common divisor (GCD) and the least common multiple (LCM).
pub trait Factor
where
    Self: Copy + Mul<Output = Self> + Div<Output = Self> + PartialEq + Eq + Default,
{
    /// The zero constant.
    const ZERO: Self;
    /// The one constant.
    const ONE: Self;

    /// The three-way parity.[^1]
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(Rational::new(2, 1).parity(), Some(false)); // even
    /// assert_eq!(Rational::new(1, 1).parity(), Some(true)); // odd
    /// assert_eq!(Rational::new(1, 2).parity(), None); // none
    /// ```
    ///
    /// [^1]: P. Lynch and M. Mackey, “Parity and Partition of the Rational Numbers”, [The College
    /// Mathematics Journal, 55(5), 387–399](https://doi.org/10.1080/07468342.2024.2311632).
    #[must_use]
    fn parity(self) -> Option<bool>;
    /// Whether the parity is odd, i.e., <code>[Self::parity()] == [Some]\([true]\)</code>.
    #[must_use]
    #[inline]
    fn is_odd(self) -> bool {
        self.parity() == Some(true)
    }
    /// Whether the parity is even, i.e., <code>[Self::parity()] == [Some]\([false]\)</code>.
    #[must_use]
    #[inline]
    fn is_even(self) -> bool {
        self.parity() == Some(false)
    }

    /// The absolute.
    #[must_use]
    fn abs(self) -> Self;

    /// Finds the greatest common divisor (GCD) of `self` and `other`.
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(2i32, (-2).gcd(4));
    /// assert_eq!(Rational::from(2), Rational::from(-2).gcd(4.into()));
    /// assert_eq!(Rational::from(2), Rational::from(-2).gcd((-4).into()));
    /// assert_eq!(Rational::from(2), Rational::from(2).gcd((-4).into()));
    /// ```
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
        let gcd = self.gcd(other);
        let lcm = (gcd != Self::ZERO)
            .then(|| self * (other / gcd))
            .map_or(Self::ZERO, Self::abs);
        (gcd, lcm)
    }
    /// Finds the GCD of iterator over [`Self`].
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(Rational::ZERO, Rational::gcd_reduce([]));
    /// assert_eq!(
    ///     Rational::from(2),
    ///     Rational::gcd_reduce([Rational::from(-2)])
    /// );
    /// assert_eq!(
    ///     Rational::new(3, 1),
    ///     Rational::gcd_reduce([
    ///         Rational::new(-3, 1),
    ///         Rational::new(-6, 1),
    ///         Rational::new(9, 1),
    ///     ])
    /// );
    /// ```
    #[inline]
    #[must_use]
    fn gcd_reduce(r: impl IntoIterator<Item = Self>) -> Self {
        r.into_iter()
            .reduce(Self::gcd)
            .map_or(Self::ZERO, Self::abs)
    }
    /// Finds the LCM of iterator over [`Self`].
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(Rational::ONE, Rational::lcm_reduce([]));
    /// assert_eq!(
    ///     Rational::from(2),
    ///     Rational::lcm_reduce([Rational::from(-2)])
    /// );
    /// assert_eq!(
    ///     Rational::new(18, 1),
    ///     Rational::lcm_reduce([
    ///         Rational::new(-3, 1),
    ///         Rational::new(-6, 1),
    ///         Rational::new(9, 1),
    ///     ])
    /// );
    /// ```
    #[inline]
    #[must_use]
    fn lcm_reduce(r: impl IntoIterator<Item = Self>) -> Self {
        r.into_iter().reduce(Self::lcm).map_or(Self::ONE, Self::abs)
    }
    /// Finds the [`Self::gcd()`] or [`Self::lcm()`] (specified with function `f` and respective
    /// identity `i`, i.e., [`Self::ZERO`] or [`Self::ONE`]) of iterator `r` over [`Self`].
    ///
    /// See [`Self::signed_gcd_reduce()`] and [`Self::signed_lcm_reduce()`] for convenience.
    #[must_use]
    fn signed_reduce(f: fn(Self, Self) -> Self, i: Self, r: impl IntoIterator<Item = Self>) -> Self
    where
        Self: Neg<Output = Self> + PartialOrd + Ord,
    {
        let (acc, neg, len) = r.into_iter().fold((i, 0, 0), |(acc, neg, len), r| {
            (f(acc, r), neg + usize::from(r < Self::ZERO), len + 1)
        });
        if neg > len / 2 { -acc } else { acc }
    }
    /// Finds the [`Self::gcd()`] and the predominant sign of iterator over [`Self`].
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(Rational::ZERO, Rational::signed_gcd_reduce([]));
    /// assert_eq!(
    ///     Rational::from(-2),
    ///     Rational::signed_gcd_reduce([Rational::from(-2)])
    /// );
    /// assert_eq!(
    ///     Rational::new(-3, 1),
    ///     Rational::signed_gcd_reduce([
    ///         Rational::new(-3, 1),
    ///         Rational::new(-6, 1),
    ///         Rational::new(9, 1),
    ///     ])
    /// );
    /// assert_eq!(
    ///     Rational::new(3, 1),
    ///     Rational::signed_gcd_reduce(
    ///         ([
    ///             Rational::new(-3, 1),
    ///             Rational::new(-6, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(9, 1),
    ///         ])
    ///     )
    /// );
    /// assert_eq!(
    ///     Rational::new(-3, 1),
    ///     Rational::signed_gcd_reduce(
    ///         ([
    ///             Rational::new(-3, 1),
    ///             Rational::new(-6, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(-9, 1),
    ///         ])
    ///     )
    /// );
    /// ```
    #[inline]
    #[must_use]
    fn signed_gcd_reduce(r: impl IntoIterator<Item = Self>) -> Self
    where
        Self: Neg<Output = Self> + PartialOrd + Ord,
    {
        Self::signed_reduce(Self::gcd, Self::ZERO, r)
    }
    /// Finds the [`Self::lcm()`] and the predominant sign of iterator over [`Self`].
    ///
    /// ```
    /// use vee::{Factor, Rational};
    ///
    /// assert_eq!(Rational::ONE, Rational::signed_lcm_reduce([]));
    /// assert_eq!(
    ///     Rational::from(-2),
    ///     Rational::signed_lcm_reduce([Rational::from(-2)])
    /// );
    /// assert_eq!(
    ///     Rational::new(-18, 1),
    ///     Rational::signed_lcm_reduce([
    ///         Rational::new(-3, 1),
    ///         Rational::new(-6, 1),
    ///         Rational::new(9, 1),
    ///     ])
    /// );
    /// assert_eq!(
    ///     Rational::new(18, 1),
    ///     Rational::signed_lcm_reduce(
    ///         ([
    ///             Rational::new(-3, 1),
    ///             Rational::new(-6, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(9, 1),
    ///         ])
    ///     )
    /// );
    /// assert_eq!(
    ///     Rational::new(-18, 1),
    ///     Rational::signed_lcm_reduce(
    ///         ([
    ///             Rational::new(-3, 1),
    ///             Rational::new(-6, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(9, 1),
    ///             Rational::new(-9, 1),
    ///         ])
    ///     )
    /// );
    /// ```
    #[inline]
    #[must_use]
    fn signed_lcm_reduce(r: impl IntoIterator<Item = Self>) -> Self
    where
        Self: Neg<Output = Self> + PartialOrd + Ord,
    {
        Self::signed_reduce(Self::lcm, Self::ONE, r)
    }
}

impl Factor for Rational {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    #[inline]
    fn parity(self) -> Option<bool> {
        self.q.is_odd().then_some(self.p.is_odd())
    }
    #[inline]
    fn abs(self) -> Self {
        Self::abs(&self)
    }
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
    const ONE: Self = 1;

    #[inline]
    fn parity(self) -> Option<bool> {
        Some(self & 1 != 0)
    }
    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }
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
    const ONE: Self = 1;

    #[inline]
    fn parity(self) -> Option<bool> {
        Some(self & 1 != 0)
    }
    #[inline]
    fn abs(self) -> Self {
        self
    }
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
