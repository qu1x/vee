// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB4, LUT4, 4, [
    ('v', e),
    ('W', e0),
    ('x', e1),
    ('y', e2),
    ('z', e3),
    ('ð', e4),
    ('X', e01),
    ('Y', e02),
    ('Z', e03),
    ('Ð', e04),
    ('a', e12),
    ('b', e13),
    ('c', e14),
    ('d', e23),
    ('e', e24),
    ('f', e34),
    ('F', e012),
    ('E', e031),
    ('D', e014),
    ('C', e023),
    ('B', e042),
    ('A', e034),
    ('ð', e132),
    ('z', e124),
    ('y', e143),
    ('x', e234),
    ('Ð', e0123),
    ('Z', e0142),
    ('Y', e0134),
    ('X', e0243),
    ('w', e1234),
    ('V', e01234),
]);

/// The named entities of the PGA with embedded dimension $`N = 4`$.
impl<const M: i8> Multivector<Pga<M, 4>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{01234}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01234()
    }
    /// The multivector of norm $`n \equiv s + P`$.
    ///
    /// Quadvector $`P`$ does square to a scalar, therefore $`n`$ is a generalized complex number.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let quadvector_norm_squared = Vee::point().norm_squared();
    ///
    /// assert_eq!(
    ///     quadvector_norm_squared.basis_blades(),
    ///     Vee::scalar().basis_blades()
    /// );
    /// format_eq!(quadvector_norm_squared, ["+ww"]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of bias $`v_\infty \equiv w\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of plane displacement
    /// $`p_0 \equiv a\e_{12} + b\e_{13} + c\e_{14} + d\e_{23} + e\e_{24} + f\e_{34}`$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        (Self::e12() + Self::e13() + Self::e14() + Self::e23() + Self::e24() + Self::e34()).alt()
    }
    /// The multivector of plane moment
    /// $`p_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{04}`$.
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        (Self::e01() + Self::e02() + Self::e03() + Self::e04()).alt()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_displacement() + Self::plane_moment()
    }
    /// The multivector of line displacement
    /// $`\ell_0 \equiv x\e_{132} + y\e_{124} + z\e_{143} + ð\e_{234}`$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        (Self::e132() + Self::e124() + Self::e143() + Self::e234()).alt()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{012} + B\e_{031} + C\e_{014} + D\e_{023} + E\e_{042} + F\e_{034}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        (Self::e012() + Self::e031() + Self::e014() + Self::e023() + Self::e042() + Self::e034())
            .alt()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_displacement() + Self::line_moment()
    }
    /// The multivector of weight $`P_0 \equiv W\e_{1234}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e1234()
    }
    /// The multivector of direction
    /// $`P_\infty \equiv X\e_{0123} + Y\e_{0142} + Z\e_{0134} + Ð\e_{0243}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e0123() + Self::e0142() + Self::e0134() + Self::e0243()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of single rotator $`r_1 \equiv s + \ell_0`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(
    ///     single_rotator.basis_blades(),
    ///     Vee::single_rotator().basis_blades()
    /// );
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(
    ///     single_rotator.basis_blades(),
    ///     Vee::single_rotator().basis_blades()
    /// );
    /// format_eq!(single_rotator, [
    ///     "-ẋ͔ẋ͕-ẏ͔ẏ͕-ż͔ż͕-ð͔̇ð͕̇",
    ///     "+(-ẋ͔ẏ͕+ẋ͕ẏ͔)e12",
    ///     "+(-ẋ͔ż͕+ẋ͕ż͔)e13",
    ///     "+(-ẋ͔ð͕̇+ẋ͕ð͔̇)e14",
    ///     "+(-ẏ͔ż͕+ẏ͕ż͔)e23",
    ///     "+(-ẏ͔ð͕̇+ẏ͕ð͔̇)e24",
    ///     "+(-ż͔ð͕̇+ż͕ð͔̇)e34",
    /// ]);
    ///
    /// let norm_squared = Vee::single_rotator().norm_squared();
    /// assert_eq!(
    ///     norm_squared.basis_blades(),
    ///     (Vee::scalar() + Vee::weight()).basis_blades()
    /// );
    /// format_eq!(norm_squared, [
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv",
    ///     "+2(-ȧḟ+ḃė-ċḋ)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::plane_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + \ell_0 + P_0`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(
    ///     double_rotator.basis_blades(),
    ///     Vee::double_rotator().basis_blades()
    /// );
    /// format_eq!(double_rotator, [
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕+v͔v͕",
    ///     "+(+ȧ͔v͕+ȧ͕v͔-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔+ḃ͔v͕+ḃ͕v͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔+ċ͔v͕+ċ͕v͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔+ḋ͔v͕+ḋ͕v͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔+ė͔v͕+ė͕v͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔+ḟ͔v͕+ḟ͕v͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(
    ///     double_rotator.basis_blades(),
    ///     Vee::double_rotator().basis_blades()
    /// );
    /// format_eq!(double_rotator, [
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕",
    ///     "+(-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    /// ]);
    ///
    /// let norm_squared = Vee::double_rotator().norm_squared();
    /// assert_eq!(
    ///     norm_squared.basis_blades(),
    ///     (Vee::scalar() + Vee::weight()).basis_blades()
    /// );
    /// format_eq!(norm_squared, [
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv+ww",
    ///     "+2(-ȧḟ+ḃė-ċḋ+vw)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + p_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "+w͔w͕",
    ///     "+(-X͔w͕+X͕w͔)e01",
    ///     "+(-Y͔w͕+Y͕w͔)e02",
    ///     "+(-Z͔w͕+Z͕w͔)e03",
    ///     "+(+w͔Ð͕-w͕Ð͔)e04",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::plane_moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let simple_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(simple_motor.basis_blades(), Vee::simple_motor().basis_blades());
    /// format_eq!(simple_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    /// ]);
    ///
    /// let norm_squared = Vee::simple_motor().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     // Scalar condition.
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv",
    ///     // Point condition.
    ///     "+2(-ȧḟ+ḃė-ċḋ)e1234", // Weight condition.
    ///     "+2(+Ẏḟ-Żė+ḋÐ̇)e0243", // Direction condition.
    ///     "+2(-Ẋḟ+Żċ-ḃÐ̇)e0134", // Direction condition.
    ///     "+2(+Ẋė-Ẏċ+ȧÐ̇)e0142", // Direction condition.
    ///     "+2(-Ẋḋ+Ẏḃ-Żȧ)e0123", // Direction condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(-ȧḟ+ḃė-ċḋ)w͓", // Vanishes with weight condition.
    ///     "+(+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv)w͓e1234",
    ///     "+(+(-ȧȧ-ḃḃ-ċċ+ḋḋ+ėė+ḟḟ+vv)X͓+2(+ȧv-ḃḋ-ċė)Y͓+2(+ȧḋ+ḃv-ċḟ)Z͓+2(-Ẋv-Ẏȧ-Żḃ-ċÐ̇)w͓\
    ///        +2(+ȧė+ḃḟ+ċv)Ð͓)e0243",
    ///     "+(+2(-ȧv-ḃḋ-ċė)X͓+(-ȧȧ+ḃḃ+ċċ-ḋḋ-ėė+ḟḟ+vv)Y͓+2(-ȧḃ+ḋv-ėḟ)Z͓+2(+Ẋȧ-Ẏv-Żḋ-ėÐ̇)w͓\
    ///        +2(-ȧċ+ḋḟ+ėv)Ð͓)e0134",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)X͓+2(-ȧḃ-ḋv-ėḟ)Y͓+(+ȧȧ-ḃḃ+ċċ-ḋḋ+ėė-ḟḟ+vv)Z͓+2(+Ẋḃ+Ẏḋ-Żv-ḟÐ̇)w͓\
    ///        +2(-ḃċ-ḋė+ḟv)Ð͓)e0142",
    ///     "+(+2(+ȧė+ḃḟ-ċv)X͓+2(-ȧċ+ḋḟ-ėv)Y͓+2(-ḃċ-ḋė-ḟv)Z͓+2(+Ẋċ+Ẏė+Żḟ-vÐ̇)w͓\
    ///      +(+ȧȧ+ḃḃ-ċċ+ḋḋ-ėė-ḟḟ+vv)Ð͓)e0123",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor().unit();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+w͓e1234",
    ///     "+(+(+1-2ȧȧ-2ḃḃ-2ċċ)X͓+2(+ȧv-ḃḋ-ċė)Y͓+2(+ȧḋ+ḃv-ċḟ)Z͓+2(-Ẋv-Ẏȧ-Żḃ-ċÐ̇)w͓+2(+ȧė+ḃḟ+ċv)Ð͓)e0243",
    ///     "+(+2(-ȧv-ḃḋ-ċė)X͓+(+1-2ȧȧ-2ḋḋ-2ėė)Y͓+2(-ȧḃ+ḋv-ėḟ)Z͓+2(+Ẋȧ-Ẏv-Żḋ-ėÐ̇)w͓+2(-ȧċ+ḋḟ+ėv)Ð͓)e0134",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)X͓+2(-ȧḃ-ḋv-ėḟ)Y͓+(+1-2ḃḃ-2ḋḋ-2ḟḟ)Z͓+2(+Ẋḃ+Ẏḋ-Żv-ḟÐ̇)w͓+2(-ḃċ-ḋė+ḟv)Ð͓)e0142",
    ///     "+(+2(+ȧė+ḃḟ-ċv)X͓+2(-ȧċ+ḋḟ-ėv)Y͓+2(-ḃċ-ḋė-ḟv)Z͓+2(+Ẋċ+Ẏė+Żḟ-vÐ̇)w͓+(+1-2ċċ-2ėė-2ḟḟ)Ð͓)e0123",
    /// ]);
    ///
    /// let line = Vee::line().pin() << Vee::simple_motor();
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// let plane = Vee::plane().pin() << Vee::simple_motor();
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    ///
    /// let volume = Vee::volume().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(volume.basis_blades(), (Vee::pseudoscalar() + Vee::volume()).basis_blades());
    /// format_eq!(volume, [
    ///     "+(+(+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv)W͓+2(+Ẋv-Ẏȧ-Żḃ-ċÐ̇)x͓+2(+Ẋȧ+Ẏv-Żḋ-ėÐ̇)y͓+2(+Ẋḃ+Ẏḋ+Żv-ḟÐ̇)z͓\
    ///      +2(+Ẋċ+Ẏė+Żḟ+vÐ̇)ð͓)e0",
    ///     "+(+(-ȧȧ-ḃḃ-ċċ+ḋḋ+ėė+ḟḟ+vv)x͓+2(+ȧv-ḃḋ-ċė)y͓+2(+ȧḋ+ḃv-ċḟ)z͓+2(+ȧė+ḃḟ+ċv)ð͓)e1",
    ///     "+(+2(-ȧv-ḃḋ-ċė)x͓+(-ȧȧ+ḃḃ+ċċ-ḋḋ-ėė+ḟḟ+vv)y͓+2(-ȧḃ+ḋv-ėḟ)z͓+2(-ȧċ+ḋḟ+ėv)ð͓)e2",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)x͓+2(-ȧḃ-ḋv-ėḟ)y͓+(+ȧȧ-ḃḃ+ċċ-ḋḋ+ėė-ḟḟ+vv)z͓+2(-ḃċ-ḋė+ḟv)ð͓)e3",
    ///     "+(+2(+ȧė+ḃḟ-ċv)x͓+2(-ȧċ+ḋḟ-ėv)y͓+2(-ḃċ-ḋė-ḟv)z͓+(+ȧȧ+ḃḃ-ċċ+ḋḋ-ėė-ḟḟ+vv)ð͓)e4",
    ///     // Vanishes with point condition.
    ///     "+2(+(-ȧḟ+ḃė-ċḋ)W͓+(+Ẏḟ-Żė+ḋÐ̇)x͓+(-Ẋḟ+Żċ-ḃÐ̇)y͓+(+Ẋė-Ẏċ+ȧÐ̇)z͓+(-Ẋḋ+Ẏḃ-Żȧ)ð͓)I",
    /// ]);
    ///
    /// let volume = Vee::volume().pin() << Vee::simple_motor().unit();
    ///
    /// assert_eq!(volume.basis_blades(), Vee::volume().basis_blades());
    /// format_eq!(volume, [
    ///     "+(+W͓+2(+Ẋv-Ẏȧ-Żḃ-ċÐ̇)x͓+2(+Ẋȧ+Ẏv-Żḋ-ėÐ̇)y͓+2(+Ẋḃ+Ẏḋ+Żv-ḟÐ̇)z͓+2(+Ẋċ+Ẏė+Żḟ+vÐ̇)ð͓)e0",
    ///     "+(+(+1-2ȧȧ-2ḃḃ-2ċċ)x͓+2(+ȧv-ḃḋ-ċė)y͓+2(+ȧḋ+ḃv-ċḟ)z͓+2(+ȧė+ḃḟ+ċv)ð͓)e1",
    ///     "+(+2(-ȧv-ḃḋ-ċė)x͓+(+1-2ȧȧ-2ḋḋ-2ėė)y͓+2(-ȧḃ+ḋv-ėḟ)z͓+2(-ȧċ+ḋḟ+ėv)ð͓)e2",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)x͓+2(-ȧḃ-ḋv-ėḟ)y͓+(+1-2ḃḃ-2ḋḋ-2ḟḟ)z͓+2(-ḃċ-ḋė+ḟv)ð͓)e3",
    ///     "+(+2(+ȧė+ḃḟ-ċv)x͓+2(-ȧċ+ḋḟ-ėv)y͓+2(-ḃċ-ḋė-ḟv)z͓+(+1-2ċċ-2ėė-2ḟḟ)ð͓)e4",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::plane()
    }
    /// The multivector of single motor $`m_1 \equiv s + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     single_motor.basis_blades(),
    ///     Vee::single_motor().basis_blades()
    /// );
    /// format_eq!(single_motor, [
    ///     "+v͔v͕",
    ///     "+(+Ẋ͕v͔+Ẏ͕ȧ͔+Ż͕ḃ͔+ċ͔Ð͕̇)e01",
    ///     "+(-Ẋ͕ȧ͔+Ẏ͕v͔+Ż͕ḋ͔+ė͔Ð͕̇)e02",
    ///     "+(-Ẋ͕ḃ͔-Ẏ͕ḋ͔+Ż͕v͔+ḟ͔Ð͕̇)e03",
    ///     "+(-Ẋ͕ċ͔-Ẏ͕ė͔-Ż͕ḟ͔+v͔Ð͕̇)e04",
    ///     "+ȧ͔v͕e12",
    ///     "+ḃ͔v͕e13",
    ///     "+ċ͔v͕e14",
    ///     "+ḋ͔v͕e23",
    ///     "+ė͔v͕e24",
    ///     "+ḟ͔v͕e34",
    ///     "+(-Ẏ͕ḟ͔+Ż͕ė͔-ḋ͔Ð͕̇)e0243",
    ///     "+(+Ẋ͕ḟ͔-Ż͕ċ͔+ḃ͔Ð͕̇)e0134",
    ///     "+(-Ẋ͕ė͔+Ẏ͕ċ͔-ȧ͔Ð͕̇)e0142",
    ///     "+(+Ẋ͕ḋ͔-Ẏ͕ḃ͔+Ż͕ȧ͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(
    ///     single_motor.basis_blades(),
    ///     Vee::single_motor().basis_blades()
    /// );
    /// format_eq!(single_motor, [
    ///     "-ẋ͔ẋ͕-ẏ͔ẏ͕-ż͔ż͕-ð͔̇ð͕̇",
    ///     "+(+Ȧ͔ẏ͕-Ȧ͕ẏ͔+Ḃ͔ż͕-Ḃ͕ż͔+Ċ͔ð͕̇-Ċ͕ð͔̇)e01",
    ///     "+(-Ȧ͔ẋ͕+Ȧ͕ẋ͔+Ḋ͔ż͕-Ḋ͕ż͔+Ė͔ð͕̇-Ė͕ð͔̇)e02",
    ///     "+(-Ḃ͔ẋ͕+Ḃ͕ẋ͔-Ḋ͔ẏ͕+Ḋ͕ẏ͔+Ḟ͔ð͕̇-Ḟ͕ð͔̇)e03",
    ///     "+(-Ċ͔ẋ͕+Ċ͕ẋ͔-Ė͔ẏ͕+Ė͕ẏ͔-Ḟ͔ż͕+Ḟ͕ż͔)e04",
    ///     "+(-ẋ͔ẏ͕+ẋ͕ẏ͔)e12",
    ///     "+(-ẋ͔ż͕+ẋ͕ż͔)e13",
    ///     "+(-ẋ͔ð͕̇+ẋ͕ð͔̇)e14",
    ///     "+(-ẏ͔ż͕+ẏ͕ż͔)e23",
    ///     "+(-ẏ͔ð͕̇+ẏ͕ð͔̇)e24",
    ///     "+(-ż͔ð͕̇+ż͕ð͔̇)e34",
    ///     "+(-Ḋ͔ð͕̇-Ḋ͕ð͔̇+Ė͔ż͕+Ė͕ż͔-Ḟ͔ẏ͕-Ḟ͕ẏ͔)e0243",
    ///     "+(+Ḃ͔ð͕̇+Ḃ͕ð͔̇-Ċ͔ż͕-Ċ͕ż͔+Ḟ͔ẋ͕+Ḟ͕ẋ͔)e0134",
    ///     "+(-Ȧ͔ð͕̇-Ȧ͕ð͔̇+Ċ͔ẏ͕+Ċ͕ẏ͔-Ė͔ẋ͕-Ė͕ẋ͔)e0142",
    ///     "+(+Ȧ͔ż͕+Ȧ͕ż͔-Ḃ͔ẏ͕-Ḃ͕ẏ͔+Ḋ͔ẋ͕+Ḋ͕ẋ͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::plane() + Self::direction()
    }
    /// The multivector of double motor $`m_2 \equiv s + p + P`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     double_motor.basis_blades(),
    ///     Vee::double_motor().basis_blades()
    /// );
    /// format_eq!(double_motor, [
    ///     "+v͔v͕",
    ///     "+(+Ẋ͕v͔+Ẏ͕ȧ͔+Ż͕ḃ͔+ċ͔Ð͕̇)e01",
    ///     "+(-Ẋ͕ȧ͔+Ẏ͕v͔+Ż͕ḋ͔+ė͔Ð͕̇)e02",
    ///     "+(-Ẋ͕ḃ͔-Ẏ͕ḋ͔+Ż͕v͔+ḟ͔Ð͕̇)e03",
    ///     "+(-Ẋ͕ċ͔-Ẏ͕ė͔-Ż͕ḟ͔+v͔Ð͕̇)e04",
    ///     "+ȧ͔v͕e12",
    ///     "+ḃ͔v͕e13",
    ///     "+ċ͔v͕e14",
    ///     "+ḋ͔v͕e23",
    ///     "+ė͔v͕e24",
    ///     "+ḟ͔v͕e34",
    ///     "+v͕w͔e1234",
    ///     "+(+Ẋ͕w͔-Ẏ͕ḟ͔+Ż͕ė͔-ḋ͔Ð͕̇)e0243",
    ///     "+(+Ẋ͕ḟ͔+Ẏ͕w͔-Ż͕ċ͔+ḃ͔Ð͕̇)e0134",
    ///     "+(-Ẋ͕ė͔+Ẏ͕ċ͔+Ż͕w͔-ȧ͔Ð͕̇)e0142",
    ///     "+(+Ẋ͕ḋ͔-Ẏ͕ḃ͔+Ż͕ȧ͔+w͔Ð͕̇)e0123",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(
    ///     double_motor.basis_blades(),
    ///     Vee::double_motor().basis_blades()
    /// );
    /// format_eq!(double_motor, [
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕",
    ///     "+(-Ẏ͔ȧ͕+Ẏ͕ȧ͔-Ż͔ḃ͕+Ż͕ḃ͔+ċ͔Ð͕̇-ċ͕Ð͔̇)e01",
    ///     "+(+Ẋ͔ȧ͕-Ẋ͕ȧ͔-Ż͔ḋ͕+Ż͕ḋ͔+ė͔Ð͕̇-ė͕Ð͔̇)e02",
    ///     "+(+Ẋ͔ḃ͕-Ẋ͕ḃ͔+Ẏ͔ḋ͕-Ẏ͕ḋ͔+ḟ͔Ð͕̇-ḟ͕Ð͔̇)e03",
    ///     "+(+Ẋ͔ċ͕-Ẋ͕ċ͔+Ẏ͔ė͕-Ẏ͕ė͔+Ż͔ḟ͕-Ż͕ḟ͔)e04",
    ///     "+(-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    ///     "+(-Ẏ͔ḟ͕-Ẏ͕ḟ͔+Ż͔ė͕+Ż͕ė͔-ḋ͔Ð͕̇-ḋ͕Ð͔̇)e0243",
    ///     "+(+Ẋ͔ḟ͕+Ẋ͕ḟ͔-Ż͔ċ͕-Ż͕ċ͔+ḃ͔Ð͕̇+ḃ͕Ð͔̇)e0134",
    ///     "+(-Ẋ͔ė͕-Ẋ͕ė͔+Ẏ͔ċ͕+Ẏ͕ċ͔-ȧ͔Ð͕̇-ȧ͕Ð͔̇)e0142",
    ///     "+(+Ẋ͔ḋ͕+Ẋ͕ḋ͔-Ẏ͔ḃ͕-Ẏ͕ḃ͔+Ż͔ȧ͕+Ż͕ȧ͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::plane() + Self::point()
    }
    /// The multivector of rotoreflector $`f_r \equiv v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// format_eq!(rotoreflector, [
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    /// ]);
    ///
    /// let norm_squared = Vee::rotoreflector().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+xx+ẋẋ+yy+ẏẏ+zz+żż+ðð+ð̇ð̇",
    ///     "+2(-xẋ-yẏ-zż-ðð̇)e1234", // Weight condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(+xẋ+yẏ+zż+ðð̇)w͓", // Vanishes with weight condition.
    ///     "+(-xx-ẋẋ-yy-ẏẏ-zz-żż-ðð-ð̇ð̇)w͓e1234",
    ///     "+(+(+xx-ẋẋ-yy+ẏẏ-zz+żż-ðð+ð̇ð̇)X͓+2(+xy-ẋẏ+zð̇-żð)Y͓+2(+xz-ẋż-yð̇+ẏð)Z͓\
    ///      +2(+xð-ẋð̇+yż-ẏz)Ð͓)e0243",
    ///     "+(+2(+xy-ẋẏ-zð̇+żð)X͓+(-xx+ẋẋ+yy-ẏẏ-zz+żż-ðð+ð̇ð̇)Y͓+2(+xð̇-ẋð+yz-ẏż)Z͓\
    ///      +2(-xż+ẋz+yð-ẏð̇)Ð͓)e0134",
    ///     "+(+2(+xz-ẋż+yð̇-ẏð)X͓+2(-xð̇+ẋð+yz-ẏż)Y͓+(-xx+ẋẋ-yy+ẏẏ+zz-żż-ðð+ð̇ð̇)Z͓\
    ///      +2(+xẏ-ẋy+zð-żð̇)Ð͓)e0142",
    ///     "+(+2(+xð-ẋð̇-yż+ẏz)X͓+2(+xż-ẋz+yð-ẏð̇)Y͓+2(-xẏ+ẋy+zð-żð̇)Z͓\
    ///      +(-xx+ẋẋ-yy+ẏẏ-zz+żż+ðð-ð̇ð̇)Ð͓)e0123",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector().unit();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+1-2xx-2ẋẋ-2yy-2ẏẏ-2zz-2żż-2ðð-2ð̇ð̇)w͓e1234",
    ///     "+(+(+1-2ẋẋ-2yy-2zz-2ðð)X͓+2(+xy-ẋẏ+zð̇-żð)Y͓+2(+xz-ẋż-yð̇+ẏð)Z͓+2(+xð-ẋð̇+yż-ẏz)Ð͓)e0243",
    ///     "+(+2(+xy-ẋẏ-zð̇+żð)X͓+(+1-2xx-2ẏẏ-2zz-2ðð)Y͓+2(+xð̇-ẋð+yz-ẏż)Z͓+2(-xż+ẋz+yð-ẏð̇)Ð͓)e0134",
    ///     "+(+2(+xz-ẋż+yð̇-ẏð)X͓+2(-xð̇+ẋð+yz-ẏż)Y͓+(+1-2xx-2yy-2żż-2ðð)Z͓+2(+xẏ-ẋy+zð-żð̇)Ð͓)e0142",
    ///     "+(+2(+xð-ẋð̇-yż+ẏz)X͓+2(+xż-ẋz+yð-ẏð̇)Y͓+2(-xẏ+ẋy+zð-żð̇)Z͓+(+1-2xx-2yy-2zz-2ð̇ð̇)Ð͓)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::line_displacement()
    }
    /// The multivector of transflector $`f_t \equiv v + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     transflector.basis_blades(),
    ///     Vee::transflector().basis_blades()
    /// );
    /// format_eq!(transflector, [
    ///     "+(-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+(+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-Ẏ͕ð͔+y͔Ð͕̇)e042",
    ///     "+(+Ẏ͕z͔-Ż͕y͔)e023",
    ///     "+(+Ẋ͕ð͔-x͔Ð͕̇)e014",
    ///     "+(-Ẋ͕z͔+Ż͕x͔)e031",
    ///     "+(+Ẋ͕y͔-Ẏ͕x͔)e012",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume() + Self::line_moment()
    }
    /// The multivector of flector $`f_s \equiv v + \ell`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let flector = Vee::volume().lhs() * Vee::simple_motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::simple_flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+W͔v͕-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    ///     "+(+W͔ḟ͕+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-W͔ė͕-Ẏ͕ð͔+y͔Ð͕̇)e042",
    ///     "+(+W͔ḋ͕+Ẏ͕z͔-Ż͕y͔)e023",
    ///     "+(+W͔ċ͕+Ẋ͕ð͔-x͔Ð͕̇)e014",
    ///     "+(-W͔ḃ͕-Ẋ͕z͔+Ż͕x͔)e031",
    ///     "+(+W͔ȧ͕+Ẋ͕y͔-Ẏ͕x͔)e012",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_flector() -> Self {
        Self::volume() + Self::line()
    }
    /// The multivector of flector $`f \equiv v + \ell + S`$.
    ///
    /// ```
    /// use vee::{PgaP4 as Vee, format_eq};
    ///
    /// let flector = Vee::volume().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+W͔v͕-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    ///     "+(+W͔ḟ͕+X͕y͔-Y͕x͔+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-W͔ė͕+X͕z͔-Ẏ͕ð͔-Z͕x͔+y͔Ð͕̇)e042",
    ///     "+(+W͔ḋ͕+X͕ð͔+Ẏ͕z͔-Ż͕y͔-x͔Ð͕)e023",
    ///     "+(+W͔ċ͕+Ẋ͕ð͔+Y͕z͔-Z͕y͔-x͔Ð͕̇)e014",
    ///     "+(-W͔ḃ͕-Ẋ͕z͔+Y͕ð͔+Ż͕x͔-y͔Ð͕)e031",
    ///     "+(+W͔ȧ͕+Ẋ͕y͔-Ẏ͕x͔+Z͕ð͔-z͔Ð͕)e012",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::volume() + Self::line() + Self::pseudoscalar()
    }
}
