// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB3, LUT3, 3, [
    ('v', e),
    ('W', e0),
    ('x', e1),
    ('y', e2),
    ('z', e3),
    ('X', e01),
    ('Y', e02),
    ('Z', e03),
    ('z', e12),
    ('y', e31),
    ('x', e23),
    ('Z', e021),
    ('Y', e013),
    ('X', e032),
    ('w', e123),
    ('V', e0123),
]);

/// The named entities of the PGA with embedded dimension $`N = 3`$.
impl<const M: i8> Multivector<Pga<M, 3>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{0123}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0123()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let norm_squared = Vee::line().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+xx+yy+zz",
    ///     "+2(-Xx-Yy-Zz)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of bias $`p_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`p_0 \equiv x\e_1 + y\e_2 + z\e_3`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let plane = Vee::line().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    /// format_eq!(plane, [
    ///     "+(-X͔X͕-Y͔Y͕-Z͔Z͕)e0",
    ///     "+(+X͔w͕-Y͕z͔+Z͕y͔)e1",
    ///     "+(+X͕z͔+Y͔w͕-Z͕x͔)e2",
    ///     "+(-X͕y͔+Y͕x͔+Z͔w͕)e3",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of displacement $`\ell_0 \equiv x\e_{23} + y\e_{31} + z\e_{12}`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// // A line through the origin as the join of a point and the origin.
    /// let displacement = Vee::point().lhs() & Vee::weight().rhs();
    ///
    /// assert_eq!(
    ///     displacement.basis_blades(),
    ///     Vee::displacement().basis_blades()
    /// );
    /// format_eq!(displacement, [
    ///     "-X͔w͕e23",
    ///     "-Y͔w͕e31",
    ///     "-Z͔w͕e12",
    /// ]);
    ///
    /// // A line through the origin as the meet of two planes through the origin.
    /// let displacement = Vee::normal().lhs() ^ Vee::normal().rhs();
    ///
    /// assert_eq!(
    ///     displacement.basis_blades(),
    ///     Vee::displacement().basis_blades()
    /// );
    /// format_eq!(displacement, [
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn displacement() -> Self {
        Self::e23() + Self::e31() + Self::e12()
    }
    /// The multivector of moment $`\ell_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03}`$.
    #[must_use]
    #[inline]
    pub fn moment() -> Self {
        Self::e01() + Self::e02() + Self::e03()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// // A line as the join of two points.
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+Y͔Z͕-Y͕Z͔)e01",
    ///     "+(-X͔Z͕+X͕Z͔)e02",
    ///     "+(+X͔Y͕-X͕Y͔)e03",
    ///     "+(-X͔w͕+X͕w͔)e23",
    ///     "+(-Y͔w͕+Y͕w͔)e31",
    ///     "+(-Z͔w͕+Z͕w͔)e12",
    /// ]);
    ///
    /// // A line as the meet of two planes.
    /// let line = Vee::plane().lhs() ^ Vee::plane().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::displacement() + Self::moment()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{123}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e123()
    }
    /// The multivector of direction $`P_\infty \equiv X\e_{032} + Y\e_{013} + Z\e_{021}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e032() + Self::e013() + Self::e021()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// // A point as the meet of a plane and a line.
    /// let point = Vee::plane().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
    ///     "+(-W͔x͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-W͔y͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-W͔z͕-X͕y͔+Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + \ell_0`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "-x͔x͕-y͔y͕-z͔z͕",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::displacement()
    }
    /// The multivector of translator $`t \equiv s + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-w͔w͕",
    ///     "+(+X͔w͕-X͕w͔)e01",
    ///     "+(+Y͔w͕-Y͕w͔)e02",
    ///     "+(+Z͔w͕-Z͕w͔)e03",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + \ell`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let simple_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(
    ///     simple_motor.basis_blades(),
    ///     Vee::simple_motor().basis_blades()
    /// );
    /// format_eq!(simple_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::line()
    }
    /// The multivector of motor $`m \equiv s + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "-x͔x͕-y͔y͕-z͔z͕",
    ///     "+(-Y͔z͕+Y͕z͔+Z͔y͕-Z͕y͔)e01",
    ///     "+(+X͔z͕-X͕z͔-Z͔x͕+Z͕x͔)e02",
    ///     "+(-X͔y͕+X͕y͔+Y͔x͕-Y͕x͔)e03",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(+X͔x͕+X͕x͔+Y͔y͕+Y͕y͔+Z͔z͕+Z͕z͔)I",
    /// ]);
    ///
    /// let norm_squared = motor.norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+x͔x͔x͕x͕+x͔x͔y͕y͕+x͔x͔z͕z͕+x͕x͕y͔y͔+x͕x͕z͔z͔+y͔y͔y͕y͕+y͔y͔z͕z͕+y͕y͕z͔z͔+z͔z͔z͕z͕",
    ///     "+2(-X͔x͔x͕x͕-X͔x͔y͕y͕-X͔x͔z͕z͕\
    ///         -X͕x͔x͔x͕-X͕x͕y͔y͔-X͕x͕z͔z͔\
    ///         -Y͔x͕x͕y͔-Y͔y͔y͕y͕-Y͔y͔z͕z͕\
    ///         -Y͕x͔x͔y͕-Y͕y͔y͔y͕-Y͕y͕z͔z͔\
    ///         -Z͔x͕x͕z͔-Z͔y͕y͕z͔-Z͔z͔z͕z͕\
    ///         -Z͕x͔x͔z͕-Z͕y͔y͔z͕-Z͕z͔z͔z͕)I",
    /// ]);
    ///
    /// let motor = Vee::rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕z͔-Z͕y͔)e01",
    ///     "+(-X͕z͔+Y͕v͔+Z͕x͔)e02",
    ///     "+(+X͕y͔-Y͕x͔+Z͕v͔)e03",
    ///     "+v͕x͔e23",
    ///     "+v͕y͔e31",
    ///     "+v͕z͔e12",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔)I",
    /// ]);
    ///
    /// let norm_squared = motor.norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().vector(0).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+v͔v͔v͕v͕+v͕v͕x͔x͔+v͕v͕y͔y͔+v͕v͕z͔z͔",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::line()
    }
    /// The multivector of rotoreflector $`f_r \equiv p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::rotator().rhs();
    ///
    /// assert_eq!(
    ///     rotoreflector.basis_blades(),
    ///     Vee::rotoreflector().basis_blades()
    /// );
    /// format_eq!(rotoreflector, [
    ///     "+(+v͕x͔-y͔z͕+y͕z͔)e1",
    ///     "+(+v͕y͔+x͔z͕-x͕z͔)e2",
    ///     "+(+v͕z͔-x͔y͕+x͕y͔)e3",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv p + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     transflector.basis_blades(),
    ///     Vee::transflector().basis_blades()
    /// );
    /// format_eq!(transflector, [
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+(-Y͕z͔+Z͕y͔)e032",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(-X͕y͔+Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::plane() + Self::direction()
    }
    /// The multivector of flector $`f \equiv p + P`$.
    ///
    /// ```
    /// use vee::{PgaP3 as Vee, format_eq};
    ///
    /// let flector = Vee::plane().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔)e0",
    ///     "+(+v͕x͔-y͔z͕+y͕z͔)e1",
    ///     "+(+v͕y͔+x͔z͕-x͕z͔)e2",
    ///     "+(+v͕z͔-x͔y͕+x͕y͔)e3",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
    ///     "+(+V͕x͔-W͔x͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(+V͕y͔-W͔y͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+V͕z͔-W͔z͕-X͕y͔+Y͕x͔)e021",
    /// ]);
    ///
    /// let plane = Vee::plane().pin() << Vee::flector();
    ///
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    /// format_eq!(plane, [
    ///     "+(+(+ww+xx+yy+zz)W͓+2(-Wx+Xw+Yz-Zy)x͓+2(-Wy-Xz+Yw+Zx)y͓+2(-Wz+Xy-Yx+Zw)z͓)e0",
    ///     "+(+(-ww-xx+yy+zz)x͓+2(+wz-xy)y͓+2(-wy-xz)z͓)e1",
    ///     "+(+2(-wz-xy)x͓+(-ww+xx-yy+zz)y͓+2(+wx-yz)z͓)e2",
    ///     "+(+2(+wy-xz)x͓+2(-wx-yz)y͓+(-ww+xx+yy-zz)z͓)e3",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::plane() + Self::point()
    }
}
