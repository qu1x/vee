// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![deny(clippy::arithmetic_side_effects)]

use super::{Factor, Inv, NegAssign};
use core::{
    mem::swap,
    ops::{Add, Div, Mul, MulAssign, Neg, Sub},
};

/// 64-bit non-zero integer with solely strict arithmetic.
///
/// The [`None`] variant of <code>[Option]<[Integer]></code> represents zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(transparent)]
pub struct Integer(pub(crate) i64);

impl Integer {
    /// The one.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::ONE + Integer::MIN, Some(-Integer::MAX));
    /// ```
    pub const ONE: Self = Self(1);
    /// The two.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::ONE + Integer::ONE, Some(Integer::TWO));
    /// ```
    pub const TWO: Self = Self(2);
    /// The minimum.
    pub const MIN: Self = Self(i64::MIN);
    /// The maximum.
    pub const MAX: Self = Self(i64::MAX);

    /// The number of bits.
    pub const BITS: u32 = i64::BITS;

    /// The non-zero integer or [`None`] if zero.
    #[must_use]
    #[inline]
    pub const fn new(m: i64) -> Option<Self> {
        if m != 0 { Some(Self(m)) } else { None }
    }
    /// The inner integer.
    #[must_use]
    #[inline]
    pub const fn get(self) -> i64 {
        self.0
    }
    /// Whether the integer is [`Self::ONE`].
    #[must_use]
    #[inline]
    pub const fn is_one(self) -> bool {
        self.eq(Self::ONE)
    }
    /// Whether the parity is even.
    #[must_use]
    #[inline]
    pub const fn is_even(self) -> bool {
        self.0 & 1 == 0
    }
    /// Whether the parity is odd.
    #[must_use]
    #[inline]
    pub const fn is_odd(self) -> bool {
        self.0 & 1 != 0
    }
    /// Whether the integer is positive.
    #[must_use]
    #[inline]
    pub const fn is_positive(self) -> bool {
        self.0 > 0
    }
    /// Whether the integer is negative.
    #[must_use]
    #[inline]
    pub const fn is_negative(self) -> bool {
        self.0 < 0
    }
    /// The absolute.
    #[must_use]
    #[inline]
    pub const fn abs(self) -> Self {
        Self(self.0.strict_abs())
    }
    /// The negation.
    #[must_use]
    #[inline]
    pub const fn neg(self) -> Self {
        Self(self.0.strict_neg())
    }
    /// The inverse.
    #[must_use]
    #[inline]
    pub const fn inv(self) -> Option<Self> {
        Self::ONE.div(self)
    }
    /// The power.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// assert_eq!(z(3).pow(2), z(9));
    /// ```
    #[must_use]
    #[inline]
    pub const fn pow(self, exp: u32) -> Self {
        Self(self.0.strict_pow(exp))
    }
    /// The equality.
    #[must_use]
    #[inline]
    pub const fn eq(self, other: Self) -> bool {
        self.0 == other.0
    }
    /// The multiplication.
    #[must_use]
    #[inline]
    pub const fn mul(self, other: Self) -> Self {
        Self(self.0.strict_mul(other.0))
    }
    /// The division.
    #[must_use]
    #[inline]
    pub const fn div(self, other: Self) -> Option<Self> {
        Self::new(self.0.strict_div(other.0))
    }
    /// The addition.
    #[must_use]
    #[inline]
    pub const fn add(self, other: Self) -> Option<Self> {
        Self::new(self.0.strict_add(other.0))
    }
    /// The subtraction.
    #[must_use]
    #[inline]
    pub const fn sub(self, other: Self) -> Option<Self> {
        Self::new(self.0.strict_sub(other.0))
    }
    /// The greatest common divisor (GCD).
    ///
    /// See [`Self::gcd_lcm()`] for calculating [`Self::gcd()`] and [`Self::lcm()`] by reusing the
    /// former for the latter.
    ///
    /// # Examples
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// assert_eq!(z(12).gcd(z(18)), z(6));
    /// assert_eq!(z(17).gcd(z(5)), z(1)); // coprime
    /// assert_eq!(z(-12).gcd(z(18)), z(6)); // sign is ignored
    /// assert_eq!(z(-12).gcd(z(-18)), z(6)); // sign is ignored
    /// assert_eq!(z(7).gcd(z(7)), z(7)); // equal inputs
    /// assert_eq!(z(1).gcd(z(9_223_372_036_854_775_807)), z(1)); // GCD with 1
    /// assert_eq!(z(4_294_967_311).gcd(z(4_294_967_357)), z(1)); // works here but overflows LCM
    /// assert_eq!(z(4_294_967_291).gcd(z(4_294_967_279)), z(1)); // works here but overflows LCM
    /// ```
    ///
    /// # Panics
    ///
    /// Panics on the one input pair whose GCD exceeds [`i64::MAX`].
    ///
    /// ```should_panic
    /// use vee::Integer;
    ///
    /// let _ = Integer::MIN.gcd(Integer::MIN); // `Integer::MIN` overflows `Integer::MAX`
    /// ```
    #[must_use]
    #[inline]
    #[allow(clippy::arithmetic_side_effects, clippy::many_single_char_names)]
    pub const fn gcd(self, other: Self) -> Self {
        let mut u = self.0.unsigned_abs();
        let mut v = other.0.unsigned_abs();
        let i = u.trailing_zeros();
        let j = v.trailing_zeros();
        u >>= i;
        v >>= j;
        while u != v {
            if u < v {
                swap(&mut u, &mut v);
            }
            u -= v;
            u >>= u.trailing_zeros();
        }
        let g = u << if i < j { i } else { j };
        assert!(
            g <= i64::MAX.cast_unsigned(),
            "attempt to calculate the GCD with overflow"
        );
        Self(g.cast_signed())
    }
    /// The least common multiple (LCM).
    ///
    /// Equals <code>`Self::gcd_lcm()`.1</code>.
    ///
    /// # Examples
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// assert_eq!(z(4).lcm(z(6)), z(12));
    /// assert_eq!(z(17).lcm(z(5)), z(85)); // LCM is product for coprimes
    /// assert_eq!(z(-4).lcm(z(6)), z(12)); // sign is ignored
    /// assert_eq!(z(-4).lcm(z(-6)), z(12)); // sign is ignored
    /// assert_eq!(z(5).lcm(z(5)), z(5)); // equal inputs
    /// ```
    ///
    /// # Panics
    ///
    /// Panics for LCM exceeding [`i64::MAX`].
    ///
    /// ```should_panic
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// // Both primes are well within `i64::MAX` but LCM overflows `u64::MAX`.
    /// let _ = z(4_294_967_311).lcm(z(4_294_967_357));
    /// ```
    ///
    /// ```should_panic
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// // Both primes are well within `i64::MAX` but LCM overflows `i64::MAX`.
    /// let _ = z(4_294_967_291).lcm(z(4_294_967_279));
    /// ```
    #[must_use]
    #[inline]
    pub const fn lcm(self, other: Self) -> Self {
        self.gcd_lcm(other).1
    }
    /// The [`Self::gcd()`] and [`Self::lcm()`] calculated by reusing the former for the latter.
    ///
    /// # Examples
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// let z = |m| Integer::new(m).unwrap();
    ///
    /// assert_eq!(z(12).gcd_lcm(z(18)), (z(6), z(36)));
    /// assert_eq!(z(-12).gcd_lcm(z(18)), (z(6), z(36))); // both always positive
    /// assert_eq!(z(9).gcd_lcm(z(9)), (z(9), z(9))); // equal inputs
    /// ```
    ///
    /// # Panics
    ///
    /// The GCD can succeed while the LCM still overflows, see when [`Self::lcm()`] panics. The same
    /// case panics here too, since both results are calculated at once.
    #[must_use]
    #[inline]
    pub const fn gcd_lcm(self, other: Self) -> (Self, Self) {
        let g = self.gcd(other);
        #[allow(clippy::arithmetic_side_effects)]
        match (self.0.unsigned_abs() / g.0.cast_unsigned()).checked_mul(other.0.unsigned_abs()) {
            Some(l) if l <= i64::MAX.cast_unsigned() => (g, Self(l.cast_signed())),
            _ => panic!("attempt to calculate the LCM with overflow"),
        }
    }
    /// The [`Self::gcd()`] of `iter` over <code>[TryInto]<[Self]></code>.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::gcd_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Integer::gcd_reduce([3, 0]), Integer::new(3));
    /// assert_eq!(Integer::gcd_reduce([0, 3]), Integer::new(3));
    /// assert_eq!(Integer::gcd_reduce([-2]), Integer::new(2));
    /// assert_eq!(Integer::gcd_reduce([-3, -6, 9]), Integer::new(3));
    /// ```
    #[must_use]
    #[inline]
    pub fn gcd_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        Factor::gcd_reduce(iter)
    }
    /// The [`Self::lcm()`] of `iter` over <code>[TryInto]<[Self]></code>.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::lcm_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Integer::lcm_reduce([3, 0]), None);
    /// assert_eq!(Integer::lcm_reduce([0, 3]), None);
    /// assert_eq!(Integer::lcm_reduce([-2]), Integer::new(2));
    /// assert_eq!(Integer::lcm_reduce([-3, -6, 9]), Integer::new(18));
    /// ```
    #[must_use]
    #[inline]
    pub fn lcm_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        Factor::lcm_reduce(iter)
    }
    /// The [`Self::gcd()`] inclusive the predominant sign of `iter` over
    /// <code>[TryInto]<[Self]></code>.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::signed_gcd_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Integer::signed_gcd_reduce([-3, 0]), Integer::new(-3));
    /// assert_eq!(Integer::signed_gcd_reduce([0, -3]), Integer::new(-3));
    /// assert_eq!(Integer::signed_gcd_reduce([-2]), Integer::new(-2));
    /// assert_eq!(Integer::signed_gcd_reduce([-3, -6, 9]), Integer::new(-3));
    /// assert_eq!(
    ///     Integer::signed_gcd_reduce(([-3, -6, 9, 9])),
    ///     Integer::new(3),
    /// );
    /// assert_eq!(
    ///     Integer::signed_gcd_reduce(([-3, -6, 9, 9, -9])),
    ///     Integer::new(-3),
    /// );
    /// ```
    #[must_use]
    #[inline]
    pub fn signed_gcd_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        Factor::signed_gcd_reduce(iter)
    }
    /// The [`Self::lcm()`] inclusive the predominant sign of `iter` over
    /// <code>[TryInto]<[Self]></code>.
    ///
    /// ```
    /// use vee::Integer;
    ///
    /// assert_eq!(Integer::signed_lcm_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Integer::signed_lcm_reduce([-3, 0]), None);
    /// assert_eq!(Integer::signed_lcm_reduce([0, -3]), None);
    /// assert_eq!(Integer::signed_lcm_reduce([-2]), Integer::new(-2));
    /// assert_eq!(Integer::signed_lcm_reduce([-3, -6, 9]), Integer::new(-18));
    /// assert_eq!(
    ///     Integer::signed_lcm_reduce(([-3, -6, 9, 9])),
    ///     Integer::new(18),
    /// );
    /// assert_eq!(
    ///     Integer::signed_lcm_reduce(([-3, -6, 9, 9, -9])),
    ///     Integer::new(-18),
    /// );
    /// ```
    #[must_use]
    #[inline]
    pub fn signed_lcm_reduce<I>(iter: I) -> Option<Self>
    where
        I: IntoIterator,
        I::Item: TryInto<Self>,
    {
        Factor::signed_lcm_reduce(iter)
    }
}

impl TryFrom<i64> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: i64) -> Result<Self, Self::Error> {
        Self::new(i).ok_or(())
    }
}

impl TryFrom<u64> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: u64) -> Result<Self, Self::Error> {
        let i = i
            .try_into()
            .expect("attempt to create integer with overflow");
        Self::new(i).ok_or(())
    }
}

impl TryFrom<i32> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: i32) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl TryFrom<u32> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: u32) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl TryFrom<i16> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: i16) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl TryFrom<u16> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: u16) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl TryFrom<i8> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: i8) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl TryFrom<u8> for Integer {
    type Error = ();

    #[inline]
    fn try_from(i: u8) -> Result<Self, Self::Error> {
        Self::new(i.into()).ok_or(())
    }
}

impl Factor for Integer {
    #[inline]
    fn is_negative(self) -> bool {
        self.is_negative()
    }
    #[inline]
    fn is_positive(self) -> bool {
        self.is_positive()
    }
    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }
    #[inline]
    fn gcd(self, other: Self) -> Self {
        self.gcd(other)
    }
    #[inline]
    fn lcm(self, other: Self) -> Self {
        self.lcm(other)
    }
}

impl Add for Integer {
    type Output = Option<Self>;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        self.add(other)
    }
}

impl Sub for Integer {
    type Output = Option<Self>;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        self.sub(other)
    }
}

impl Neg for Integer {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        self.neg()
    }
}

impl NegAssign for Integer {
    #[inline]
    fn neg_assign(&mut self) {
        *self = self.neg();
    }
}

impl Mul for Integer {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        self.mul(other)
    }
}

impl MulAssign for Integer {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = self.mul(other);
    }
}

impl Div for Integer {
    type Output = Option<Self>;

    #[inline]
    fn div(self, other: Self) -> Self::Output {
        self.div(other)
    }
}

impl Inv for Integer {
    type Output = Option<Self>;

    #[inline]
    fn inv(self) -> Self::Output {
        self.inv()
    }
}
