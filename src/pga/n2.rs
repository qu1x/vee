// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB2, LUT2, 2, [
    ('v', e),
    ('W', e0),
    ('x', e1),
    ('y', e2),
    ('Y', e01),
    ('X', e20),
    ('w', e12),
    ('V', e012),
]);

/// The named entities of the PGA with embedded dimension $`N = 2`$.
impl<const M: i8> Multivector<Pga<M, 2>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{012}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of moment $`\ell_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn moment() -> Self {
        Self::e0()
    }
    /// The multivector of displacement $`\ell_0 \equiv x\e_1 + y\e_2`$.
    #[must_use]
    #[inline]
    pub fn displacement() -> Self {
        Self::e1() + Self::e2()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+X͔Y͕-X͕Y͔)e0",
    ///     "+(+Y͔w͕-Y͕w͔)e1",
    ///     "+(-X͔w͕+X͕w͔)e2",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::moment() + Self::displacement()
    }
    /// The multivector of direction $`P_\infty \equiv X\e_{20} + Y\e_{01}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e20() + Self::e01()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{12}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e12()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let point = Vee::line().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-W͔y͕+W͕y͔)e20",
    ///     "+(+W͔x͕-W͕x͔)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + P_0`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "+x͔x͕+y͔y͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-w͔w͕",
    ///     "+(-Y͔w͕+Y͕w͔)e20",
    ///     "+(+X͔w͕-X͕w͔)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::direction()
    }
    /// The multivector of motor $`m \equiv s + P`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "+x͔x͕+y͔y͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-W͔y͕+W͕y͔)e20",
    ///     "+(+W͔x͕-W͕x͔)e01",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::motor();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+vv+ww)w͓e12",
    ///     "+(+(+vv-ww)X͓+2vwY͓+2(+Xw-Yv)w͓)e20",
    ///     "+(-2vwX͓+(+vv-ww)Y͓+2(+Xv+Yw)w͓)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of rotoreflector $`f_r \equiv \ell_0`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let rotoreflector = Vee::displacement().lhs() * Vee::rotator().rhs();
    ///
    /// assert_eq!(
    ///     rotoreflector.basis_blades(),
    ///     Vee::rotoreflector().basis_blades()
    /// );
    /// format_eq!(rotoreflector, [
    ///     "+(+v͕x͔-w͕y͔)e1",
    ///     "+(+v͕y͔+w͕x͔)e2",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(-xx-yy)w͓e12",
    ///     "+(+(+xx-yy)X͓+2xyY͓)e20",
    ///     "+(+2xyX͓+(-xx+yy)Y͓)e01",
    /// ]);
    /// ```
    #[must_use]
    pub fn rotoreflector() -> Self {
        Self::displacement()
    }
    /// The multivector of flector $`f \equiv \ell + S`$.
    ///
    /// ```
    /// use vee::{PgaP2 as Vee, format_eq};
    ///
    /// let flector = Vee::line().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+W͔v͕+X͕y͔-Y͕x͔)e0",
    ///     "+(+v͕x͔-w͕y͔)e1",
    ///     "+(+v͕y͔+w͕x͔)e2",
    ///     "+(+W͔w͕+X͕x͔+Y͕y͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::flector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(-xx-yy)w͓e12",
    ///     "+(+(+xx-yy)X͓+2xyY͓+2(+Vy+Wx)w͓)e20",
    ///     "+(+2xyX͓+(-xx+yy)Y͓+2(-Vx+Wy)w͓)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::line() + Self::pseudoscalar()
    }
}
