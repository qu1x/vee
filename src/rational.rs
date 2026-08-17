// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::Factor;
use core::{
    cmp::Ordering,
    fmt::{self, Alignment, Display},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// Rational number in canonical form.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rational {
    pub(crate) p: i32,
    pub(crate) q: i32,
}

impl Rational {
    /// The zero constant.
    pub const ZERO: Self = Self { p: 0, q: 1 };
    /// The one constant.
    pub const ONE: Self = Self { p: 1, q: 1 };

    /// Finds the irreducible fraction of numerator `p` and denominator `q`.
    ///
    /// # Panics
    ///
    /// Panics if denominator is zero.
    #[must_use]
    pub fn new(mut p: i32, mut q: i32) -> Self {
        assert_ne!(q, 0, "division by zero");
        if p == 0 {
            Self::ZERO
        } else if p == q {
            Self::ONE
        } else {
            let mut g = p.gcd(q);
            if q < 0 {
                g = -g;
            }
            p /= g;
            q /= g;
            Self { p, q }
        }
    }
    /// The numerator.
    #[must_use]
    #[inline]
    pub const fn p(&self) -> i32 {
        self.p
    }
    /// The denominator.
    #[must_use]
    #[inline]
    pub const fn q(&self) -> i32 {
        self.q
    }
    /// Inverts the fraction.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// let pos = Rational::new(3, 4);
    /// let neg = -pos;
    ///
    /// assert_eq!(pos * pos.inv(), Rational::ONE);
    /// assert_eq!(neg * neg.inv(), Rational::ONE);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if numerator is zero or [`i32::MIN`].
    ///
    /// ```should_panic
    /// use vee::Rational;
    ///
    /// Rational::ZERO.inv();
    /// ```
    ///
    /// ```should_panic
    /// use vee::Rational;
    ///
    /// Rational::from(i32::MIN).inv();
    /// ```
    #[must_use]
    pub fn inv(&self) -> Self {
        assert_ne!(self.p, 0, "division by zero");
        Self {
            p: self.q * self.p.signum(),
            q: self.p.strict_abs(),
        }
    }
    /// The absolute.
    ///
    /// # Panics
    ///
    /// Panics if numerator is [`i32::MIN`].
    ///
    /// ```should_panic
    /// use vee::Rational;
    ///
    /// Rational::from(i32::MIN).abs();
    /// ```
    #[must_use]
    #[inline]
    pub const fn abs(&self) -> Self {
        Self {
            p: self.p.strict_abs(),
            q: self.q,
        }
    }
    /// Whether this rational number is rational.
    #[must_use]
    #[inline]
    pub const fn is_rational(&self) -> bool {
        !self.is_integer()
    }
    /// Whether this rational number is integer.
    #[must_use]
    #[inline]
    pub const fn is_integer(&self) -> bool {
        self.q == 1
    }
    /// Whether this rational number is negative.
    #[must_use]
    #[inline]
    pub const fn is_negative(&self) -> bool {
        self.p.is_negative()
    }
    /// Whether this rational number is positive.
    #[must_use]
    #[inline]
    pub const fn is_positive(&self) -> bool {
        self.p.is_positive()
    }
    /// Whether this rational number is one.
    #[must_use]
    #[inline]
    pub const fn is_one(&self) -> bool {
        self.p == 1 && self.q == 1
    }
    /// Whether this rational number is zero.
    #[must_use]
    #[inline]
    pub const fn is_zero(&self) -> bool {
        self.p == 0
    }
    /// Power.
    ///
    /// ```
    /// use vee::Rational;
    ///
    /// assert_eq!(Rational::new(3, 4).pow(2), Rational::new(9, 16));
    /// assert_eq!(Rational::new(3, 4).pow(-2), Rational::new(16, 9));
    /// ```
    #[must_use]
    pub fn pow(&self, exp: i32) -> Self {
        let abs = exp.unsigned_abs();
        let pow = Self::new(self.p().pow(abs), self.q().pow(abs));
        if exp.is_negative() { pow.inv() } else { pow }
    }
}

impl Default for Rational {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl From<i32> for Rational {
    #[inline]
    fn from(p: i32) -> Self {
        Self { p, q: 1 }
    }
}

impl From<(i32, i32)> for Rational {
    #[inline]
    fn from((p, q): (i32, i32)) -> Self {
        Self::new(p, q)
    }
}

impl Add for Rational {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(self.p * other.q + self.q * other.p, self.q * other.q)
    }
}

impl AddAssign for Rational {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Add<i32> for Rational {
    type Output = Self;

    #[inline]
    fn add(self, other: i32) -> Self {
        Self::new(self.p + self.q * other, self.q)
    }
}

impl AddAssign<i32> for Rational {
    #[inline]
    fn add_assign(&mut self, other: i32) {
        *self = *self + other;
    }
}

impl Add<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn add(self, other: Rational) -> Rational {
        Rational::new(self * other.q + other.p, other.q)
    }
}

impl Sub for Rational {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(self.p * other.q - self.q * other.p, self.q * other.q)
    }
}

impl SubAssign for Rational {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Rational {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            p: self.p.strict_neg(),
            q: self.q,
        }
    }
}

impl Sub<i32> for Rational {
    type Output = Self;

    #[inline]
    fn sub(self, other: i32) -> Self {
        Self::new(self.p - self.q * other, self.q)
    }
}

impl SubAssign<i32> for Rational {
    #[inline]
    fn sub_assign(&mut self, other: i32) {
        *self = *self - other;
    }
}

impl Sub<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn sub(self, other: Rational) -> Rational {
        Rational::new(self * other.q - other.p, other.q)
    }
}

impl Mul for Rational {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self {
        Self::new(self.p * other.p, self.q * other.q)
    }
}

impl MulAssign for Rational {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Mul<i32> for Rational {
    type Output = Self;

    #[inline]
    fn mul(self, other: i32) -> Self {
        Self::new(self.p * other, self.q)
    }
}

impl MulAssign<i32> for Rational {
    #[inline]
    fn mul_assign(&mut self, other: i32) {
        *self = *self * other;
    }
}

impl Mul<Rational> for i32 {
    type Output = Rational;

    #[inline]
    fn mul(self, other: Rational) -> Rational {
        Rational::new(self * other.p, other.q)
    }
}

impl Div for Rational {
    type Output = Self;

    #[inline]
    fn div(self, other: Self) -> Self {
        Self::new(self.p * other.q, self.q * other.p)
    }
}

impl DivAssign for Rational {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}

impl Div<i32> for Rational {
    type Output = Self;

    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: i32) -> Self {
        Self::new(self.p, self.q * other)
    }
}

impl DivAssign<i32> for Rational {
    #[inline]
    fn div_assign(&mut self, other: i32) {
        *self = *self / other;
    }
}

impl Div<Rational> for i32 {
    type Output = Rational;

    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: Rational) -> Rational {
        Rational::new(self * other.q, other.p)
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
        (self.p * other.q).cmp(&(other.p * self.q))
    }
}

impl Display for Rational {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let suffix = if fmt.precision().is_some() { ".0" } else { "" };
        write!(fmt, "{}{suffix}", self.p)?;
        if self.is_rational() {
            let math = fmt.fill() == '$';
            let wide = matches!(fmt.align(), Some(Alignment::Center | Alignment::Right)) || math;
            if wide {
                write!(fmt, " / {}", self.q)?;
            } else {
                write!(fmt, "/{}", self.q)?;
            }
            write!(fmt, "{suffix}")?;
        }
        Ok(())
    }
}
