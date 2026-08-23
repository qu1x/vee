// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![deny(clippy::arithmetic_side_effects)]

use super::{Factor, Integer, Inv, InvAssign, NegAssign};
use core::{
    cmp::Ordering,
    fmt::{self, Alignment, Display},
    ops::{Add, Div, DivAssign, Mul, MulAssign, Neg, Sub},
};

/// 64-bit non-zero rational in canonical form with solely strict arithmetic.
///
/// Comprises two non-zero 64-bit <code>[Integer]s</code>, the numerator and the denominator.
///
/// The [`None`] variant of <code>[Option]<[Rational]></code> represents zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rational(Integer, Integer);

impl Rational {
    /// The one.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::ONE + Rational::MIN, Some(-Rational::MAX));
    /// ```
    pub const ONE: Self = Self::from_integer(Integer::ONE);
    /// The two.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::ONE + Rational::ONE, Some(Rational::TWO));
    /// ```
    pub const TWO: Self = Self::from_integer(Integer::TWO);
    /// The minimum.
    pub const MIN: Self = Self::from_integer(Integer::MIN);
    /// The maximum.
    pub const MAX: Self = Self::from_integer(Integer::MAX);

    /// The number of [`Integer`] bits.
    pub const BITS: u32 = Integer::BITS;

    /// The non-zero rational `q = m / n` in canonical form or [`None`] if zero.
    ///
    /// # Panics
    ///
    /// Panics if denominator `n` is zero regardless of numerator `m`.
    ///
    /// ```should_panic
    /// use vee::Rational;
    ///
    /// let _ = Rational::new(0, 0);
    /// ```
    #[must_use]
    #[inline]
    pub const fn new(m: i64, n: i64) -> Option<Self> {
        let n = Integer::new(n).expect("attempt to reduce rational with division by zero");
        if let Some(m) = Integer::new(m) {
            Some(Self::reduce(m, n))
        } else {
            None
        }
    }
    /// The non-zero integer in canonical form.
    #[must_use]
    #[inline]
    pub const fn new_integer(m: i64) -> Option<Self> {
        if let Some(m) = Integer::new(m) {
            Some(Self::from_integer(m))
        } else {
            None
        }
    }
    /// The non-zero integer in canonical form.
    #[must_use]
    #[inline]
    pub const fn from_integer(m: Integer) -> Self {
        Self(m, Integer::ONE)
    }
    /// The non-zero rational in canonical form.
    #[must_use]
    #[inline]
    #[allow(clippy::arithmetic_side_effects)]
    pub const fn reduce(mut m: Integer, mut n: Integer) -> Self {
        if !n.is_one() {
            let mut g = m.gcd(n);
            if n.is_negative() {
                g = g.neg();
            }
            m.0 /= g.0;
            n.0 /= g.0;
        }
        Self(m, n)
    }
    /// The numerator.
    #[must_use]
    #[inline]
    pub const fn numerator(self) -> Integer {
        self.0
    }
    /// The denominator.
    #[must_use]
    #[inline]
    pub const fn denominator(self) -> Integer {
        self.1
    }
    /// Whether the rational is [`Self::ONE`].
    #[must_use]
    #[inline]
    pub const fn is_one(self) -> bool {
        self.eq(Self::ONE)
    }
    /// The three-way parity.[^1]
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// let q = |m, n| Rational::new(m, n).unwrap();
    ///
    /// assert_eq!(q(2, 1).parity(), Some(false)); // even
    /// assert_eq!(q(1, 1).parity(), Some(true)); // odd
    /// assert_eq!(q(1, 2).parity(), None); // none
    /// ```
    ///
    /// [^1]: P. Lynch and M. Mackey, “Parity and Partition of the Rational Numbers”, [The College
    /// Mathematics Journal, 55(5), 387–399](https://doi.org/10.1080/07468342.2024.2311632).
    #[must_use]
    #[inline]
    pub const fn parity(self) -> Option<bool> {
        if self.1.is_odd() {
            Some(self.0.is_odd())
        } else {
            None
        }
    }
    /// Whether the parity is odd, i.e., <code>[Self::parity()] == [Some]\([true]\)</code>.
    #[must_use]
    #[inline]
    pub const fn is_odd(self) -> bool {
        matches!(self.parity(), Some(true))
    }
    /// Whether the parity is even, i.e., <code>[Self::parity()] == [Some]\([false]\)</code>.
    #[must_use]
    #[inline]
    pub const fn is_even(self) -> bool {
        matches!(self.parity(), Some(false))
    }
    /// The integer if [`Self::is_integer()`].
    #[must_use]
    #[inline]
    pub const fn integer(self) -> Option<Integer> {
        if self.is_integer() {
            Some(self.0)
        } else {
            None
        }
    }
    /// Whether the rational is an integer.
    #[must_use]
    #[inline]
    pub const fn is_integer(self) -> bool {
        self.1.eq(Integer::ONE)
    }
    /// Whether the rational is a rational.
    #[must_use]
    #[inline]
    pub const fn is_rational(self) -> bool {
        !self.is_integer()
    }
    /// Whether the rational is positive.
    #[must_use]
    #[inline]
    pub const fn is_positive(self) -> bool {
        self.0.is_positive()
    }
    /// Whether the rational is negative.
    #[must_use]
    #[inline]
    pub const fn is_negative(self) -> bool {
        self.0.is_negative()
    }
    /// The absolute.
    #[must_use]
    #[inline]
    pub const fn abs(self) -> Self {
        Self(self.0.abs(), self.1)
    }
    /// The negation.
    #[must_use]
    #[inline]
    pub const fn neg(self) -> Self {
        Self(self.0.neg(), self.1)
    }
    /// The inverse.
    #[must_use]
    #[inline]
    pub const fn inv(self) -> Self {
        if self.0.is_negative() {
            Self(self.1.neg(), self.0.neg())
        } else {
            Self(self.1, self.0)
        }
    }
    /// The power.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// let q = |m, n| Rational::new(m, n).unwrap();
    ///
    /// assert_eq!(q(-3, 4).pow(2), q(9, 16));
    /// assert_eq!(q(-3, 4).pow(3), q(-27, 64));
    /// assert_eq!(q(3, 4).pow(0), q(1, 1));
    /// assert_eq!(q(3, 4).pow(-2), q(16, 9));
    /// ```
    #[must_use]
    #[inline]
    pub const fn pow(self, exp: i32) -> Self {
        if exp == 0 {
            Self::ONE
        } else {
            let abs = exp.unsigned_abs();
            let pow = Self(self.0.pow(abs), self.1.pow(abs));
            if exp.is_negative() { pow.inv() } else { pow }
        }
    }
    /// The equality.
    #[must_use]
    #[inline]
    pub const fn eq(self, other: Self) -> bool {
        self.0.eq(other.0) && self.1.eq(other.1)
    }
    /// The multiplication.
    #[must_use]
    #[inline]
    pub const fn mul(self, other: Self) -> Self {
        Self::reduce(self.0.mul(other.0), self.1.mul(other.1))
    }
    /// The division.
    #[must_use]
    #[inline]
    pub const fn div(self, other: Self) -> Self {
        Self::reduce(self.0.mul(other.1), self.1.mul(other.0))
    }
    /// The addition.
    #[must_use]
    #[inline]
    pub const fn add(self, other: Self) -> Option<Self> {
        if let Some(m) = self.0.mul(other.1).add(self.1.mul(other.0)) {
            Some(Self::reduce(m, self.1.mul(other.1)))
        } else {
            None
        }
    }
    /// The subtraction.
    #[must_use]
    #[inline]
    pub const fn sub(self, other: Self) -> Option<Self> {
        if let Some(m) = self.0.mul(other.1).sub(self.1.mul(other.0)) {
            Some(Self::reduce(m, self.1.mul(other.1)))
        } else {
            None
        }
    }
    /// The greatest common divisor (GCD).
    ///
    /// See [`Self::gcd_lcm()`] for calculating [`Self::gcd()`] and [`Self::lcm()`] by reusing the
    /// former for the latter.
    #[must_use]
    #[inline]
    pub const fn gcd(self, other: Self) -> Self {
        Self(self.0.gcd(other.0), self.1.lcm(other.1))
    }
    /// The least common multiple (LCM).
    ///
    /// See [`Self::gcd_lcm()`] for calculating [`Self::gcd()`] and [`Self::lcm()`] by reusing the
    /// former for the latter.
    #[must_use]
    #[inline]
    pub const fn lcm(self, other: Self) -> Self {
        Self(self.0.lcm(other.0), self.1.gcd(other.1))
    }
    /// The [`Self::gcd()`] and [`Self::lcm()`] calculated by reusing the former for the latter.
    #[must_use]
    #[inline]
    pub const fn gcd_lcm(self, other: Self) -> (Self, Self) {
        let (g_0, l_0) = self.0.gcd_lcm(other.0);
        let (g_1, l_1) = self.1.gcd_lcm(other.1);
        (Self(g_0, l_1), Self(l_0, g_1))
    }
    /// The [`Self::gcd()`] of `iter` over <code>[TryInto]<[Self]></code>.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::gcd_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Rational::gcd_reduce([3, 0]), Rational::new_integer(3));
    /// assert_eq!(Rational::gcd_reduce([0, 3]), Rational::new_integer(3));
    /// assert_eq!(Rational::gcd_reduce([-2]), Rational::new_integer(2));
    /// assert_eq!(
    ///     Rational::gcd_reduce([(-3, 1), (-6, 1), (9, 1)]),
    ///     Rational::new(3, 1)
    /// );
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
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::lcm_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Rational::lcm_reduce([3, 0]), None);
    /// assert_eq!(Rational::lcm_reduce([0, 3]), None);
    /// assert_eq!(Rational::lcm_reduce([-2]), Rational::new_integer(2));
    /// assert_eq!(
    ///     Rational::lcm_reduce([(-3, 1), (-6, 1), (9, 1)]),
    ///     Rational::new(18, 1)
    /// );
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
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::signed_gcd_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(
    ///     Rational::signed_gcd_reduce([-3, 0]),
    ///     Rational::new_integer(-3)
    /// );
    /// assert_eq!(
    ///     Rational::signed_gcd_reduce([0, -3]),
    ///     Rational::new_integer(-3)
    /// );
    /// assert_eq!(Rational::signed_gcd_reduce([-2]), Rational::new_integer(-2));
    /// assert_eq!(
    ///     Rational::signed_gcd_reduce([(-3, 1), (-6, 1), (9, 1)]),
    ///     Rational::new(-3, 1),
    /// );
    /// assert_eq!(
    ///     Rational::signed_gcd_reduce(([(-3, 1), (-6, 1), (9, 1), (9, 1)])),
    ///     Rational::new(3, 1),
    /// );
    /// assert_eq!(
    ///     Rational::signed_gcd_reduce(([(-3, 1), (-6, 1), (9, 1), (9, 1), (-9, 1)])),
    ///     Rational::new(-3, 1),
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
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::signed_lcm_reduce(<[i64; 0]>::default()), None);
    /// assert_eq!(Rational::signed_lcm_reduce([-3, 0]), None);
    /// assert_eq!(Rational::signed_lcm_reduce([0, -3]), None);
    /// assert_eq!(Rational::signed_lcm_reduce([-2]), Rational::new_integer(-2));
    /// assert_eq!(
    ///     Rational::signed_lcm_reduce([(-3, 1), (-6, 1), (9, 1)]),
    ///     Rational::new(-18, 1),
    /// );
    /// assert_eq!(
    ///     Rational::signed_lcm_reduce(([(-3, 1), (-6, 1), (9, 1), (9, 1)])),
    ///     Rational::new(18, 1),
    /// );
    /// assert_eq!(
    ///     Rational::signed_lcm_reduce(([(-3, 1), (-6, 1), (9, 1), (9, 1), (-9, 1)])),
    ///     Rational::new(-18, 1),
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

impl From<Integer> for Rational {
    #[inline]
    fn from(m: Integer) -> Self {
        Self::from_integer(m)
    }
}

impl From<(Integer, Integer)> for Rational {
    #[inline]
    fn from((m, n): (Integer, Integer)) -> Self {
        Self::reduce(m, n)
    }
}

impl TryFrom<i64> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: i64) -> Result<Self, Self::Error> {
        Self::new_integer(i).ok_or(())
    }
}

impl TryFrom<(i64, i64)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (i64, i64)) -> Result<Self, Self::Error> {
        Self::new(m, n).ok_or(())
    }
}

impl TryFrom<u64> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: u64) -> Result<Self, Self::Error> {
        let i = i
            .try_into()
            .expect("attempt to create rational with overflow");
        Self::new_integer(i).ok_or(())
    }
}

impl TryFrom<(u64, u64)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (u64, u64)) -> Result<Self, Self::Error> {
        let (m, n) = m
            .try_into()
            .and_then(|m| n.try_into().map(|n| (m, n)))
            .expect("attempt to create rational with overflow");
        Self::new(m, n).ok_or(())
    }
}

impl TryFrom<i32> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: i32) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(i32, i32)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (i32, i32)) -> Result<Self, Self::Error> {
        Self::new(m.into(), n.into()).ok_or(())
    }
}

impl TryFrom<u32> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: u32) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(u32, u32)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (u32, u32)) -> Result<Self, Self::Error> {
        Self::new(m.into(), n.into()).ok_or(())
    }
}

impl TryFrom<i16> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: i16) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(i16, i16)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (i16, i16)) -> Result<Self, Self::Error> {
        Self::new(m.into(), n.into()).ok_or(())
    }
}

impl TryFrom<u16> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: u16) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(u16, u16)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (u16, u16)) -> Result<Self, Self::Error> {
        let m = m.into();
        Self::new(m, n.into()).ok_or(())
    }
}

impl TryFrom<i8> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: i8) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(i8, i8)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (i8, i8)) -> Result<Self, Self::Error> {
        Self::new(m.into(), n.into()).ok_or(())
    }
}

impl TryFrom<u8> for Rational {
    type Error = ();

    #[inline]
    fn try_from(i: u8) -> Result<Self, Self::Error> {
        Self::new_integer(i.into()).ok_or(())
    }
}

impl TryFrom<(u8, u8)> for Rational {
    type Error = ();

    #[inline]
    fn try_from((m, n): (u8, u8)) -> Result<Self, Self::Error> {
        Self::new(m.into(), n.into()).ok_or(())
    }
}

impl Factor for Rational {
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

impl Add for Rational {
    type Output = Option<Self>;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        self.add(other)
    }
}

impl Sub for Rational {
    type Output = Option<Self>;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        self.sub(other)
    }
}

impl Neg for Rational {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        self.neg()
    }
}

impl NegAssign for Rational {
    #[inline]
    fn neg_assign(&mut self) {
        self.0.neg_assign();
    }
}

impl Mul for Rational {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        self.mul(other)
    }
}

impl MulAssign for Rational {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = self.mul(other);
    }
}

impl Div for Rational {
    type Output = Self;

    #[inline]
    fn div(self, other: Self) -> Self::Output {
        self.div(other)
    }
}

impl DivAssign for Rational {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        *self = self.div(other);
    }
}

impl Inv for Rational {
    type Output = Self;

    #[inline]
    fn inv(self) -> Self::Output {
        self.inv()
    }
}

impl InvAssign for Rational {
    #[inline]
    fn inv_assign(&mut self) {
        *self = self.inv();
    }
}

impl PartialOrd for Rational {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.mul(other.1).cmp(&other.0.mul(self.1))
    }
}

impl Display for Rational {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let suffix = if fmt.precision().is_some() { ".0" } else { "" };
        write!(fmt, "{}{suffix}", self.numerator().get())?;
        if self.is_rational() {
            let math = fmt.fill() == '$';
            let wide = matches!(fmt.align(), Some(Alignment::Center | Alignment::Right)) || math;
            if wide {
                write!(fmt, " / {}", self.denominator().get())?;
            } else {
                write!(fmt, "/{}", self.denominator().get())?;
            }
            write!(fmt, "{suffix}")?;
        }
        Ok(())
    }
}
