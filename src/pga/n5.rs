// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB5, LUT5, 5, [
    ('v', e),
    ('W', e0),
    ('x', e1),
    ('y', e2),
    ('z', e3),
    ('ð', e4),
    ('ø', e5),
    ('X', e01),
    ('Y', e02),
    ('Z', e03),
    ('Ð', e04),
    ('Ø', e05),
    ('a', e12),
    ('b', e13),
    ('c', e14),
    ('d', e15),
    ('e', e23),
    ('f', e24),
    ('g', e25),
    ('h', e34),
    ('i', e35),
    ('j', e45),
    ('A', e012),
    ('B', e013),
    ('C', e014),
    ('D', e015),
    ('E', e023),
    ('F', e024),
    ('G', e025),
    ('H', e034),
    ('I', e035),
    ('J', e045),
    ('j', e123),
    ('i', e142),
    ('h', e125),
    ('g', e134),
    ('f', e153),
    ('e', e145),
    ('d', e243),
    ('c', e235),
    ('b', e254),
    ('a', e345),
    ('J', e0123),
    ('I', e0142),
    ('H', e0125),
    ('G', e0134),
    ('F', e0153),
    ('E', e0145),
    ('D', e0243),
    ('C', e0235),
    ('B', e0254),
    ('A', e0345),
    ('ø', e1234),
    ('ð', e1253),
    ('z', e1245),
    ('y', e1354),
    ('x', e2345),
    ('Ø', e01243),
    ('Ð', e01235),
    ('Z', e01254),
    ('Y', e01345),
    ('X', e02354),
    ('w', e12345),
    ('V', e012345),
]);

/// The named entities of the PGA with embedded dimension $`N = 5`$.
impl<const M: i8> Multivector<Pga<M, 5>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{012345}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012345()
    }
    /// The multivector of norm $`n \equiv s + \ell`$.
    ///
    /// Quadvector $`\ell`$ does not square to a scalar, therefore $`n`$ is **not** a generalized
    /// complex number.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let quadvector_norm_squared = Vee::line().norm_squared();
    ///
    /// assert_eq!(
    ///     quadvector_norm_squared.basis_blades(),
    ///     (Vee::scalar() + Vee::line_moment()).basis_blades()
    /// );
    /// format_eq!(quadvector_norm_squared, [
    ///     "+xx+yy+zz+ðð+øø",
    ///     "+2(-Hø+Ið-Jz)e0345",
    ///     "+2(+Fø-Gð+Jy)e0254",
    ///     "+2(-Eø+Gz-Iy)e0235",
    ///     "+2(+Eð-Fz+Hy)e0243",
    ///     "+2(-Cø+Dð-Jx)e0145",
    ///     "+2(+Bø-Dz+Ix)e0153",
    ///     "+2(-Bð+Cz-Hx)e0134",
    ///     "+2(-Aø+Dy-Gx)e0125",
    ///     "+2(+Að-Cy+Fx)e0142",
    ///     "+2(-Az+By-Ex)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::line()
    }
    /// The multivector of bias $`v^4_\infty \equiv w\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v^4_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4 + ø\e_5`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4() + Self::e5()
    }
    /// The multivector of $`4`$-volume $`v^4 \equiv v^4_0 + v^4_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume4() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of volume displacement.
    ///
    /// ```math
    /// v_0 \equiv a\e_{12} + b\e_{13} + c\e_{14} + d\e_{15} + e\e_{23}
    ///          + f\e_{24} + g\e_{25} + h\e_{34} + i\e_{35} + j\e_{45}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e12()
            + Self::e13()
            + Self::e14()
            + Self::e15()
            + Self::e23()
            + Self::e24()
            + Self::e25()
            + Self::e34()
            + Self::e35()
            + Self::e45()
    }
    /// The multivector of volume moment.
    ///
    /// ```math
    /// v_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{04} + Ø\e_{05}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e04() + Self::e05()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::volume_displacement() + Self::volume_moment()
    }
    /// The multivector of plane displacement.
    ///
    /// ```math
    /// p_0 \equiv a\e_{123} + b\e_{142} + c\e_{125} + d\e_{134} + e\e_{153}
    ///          + f\e_{145} + g\e_{243} + h\e_{235} + i\e_{254} + j\e_{345}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e123()
            + Self::e142()
            + Self::e125()
            + Self::e134()
            + Self::e153()
            + Self::e145()
            + Self::e243()
            + Self::e235()
            + Self::e254()
            + Self::e345()
    }
    /// The multivector of plane moment.
    ///
    /// ```math
    /// p_\infty \equiv A\e_{012} + B\e_{013} + C\e_{014} + D\e_{015} + E\e_{023}
    ///               + F\e_{024} + G\e_{025} + H\e_{034} + I\e_{035} + J\e_{045}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e012()
            + Self::e013()
            + Self::e014()
            + Self::e015()
            + Self::e023()
            + Self::e024()
            + Self::e025()
            + Self::e034()
            + Self::e035()
            + Self::e045()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_displacement() + Self::plane_moment()
    }
    /// The multivector of line displacement.
    ///
    /// ```math
    /// \ell_0 \equiv x\e_{1234} + y\e_{1253} + z\e_{1245} + ð\e_{1354} + ø\e_{2345}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e1234() + Self::e1253() + Self::e1245() + Self::e1354() + Self::e2345()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{0123} + B\e_{0142} + C\e_{0125} + D\e_{0134} + E\e_{0153}
    ///                  + F\e_{0145} + G\e_{0243} + H\e_{0235} + I\e_{0254} + J\e_{0345}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e0123()
            + Self::e0142()
            + Self::e0125()
            + Self::e0134()
            + Self::e0153()
            + Self::e0145()
            + Self::e0243()
            + Self::e0235()
            + Self::e0254()
            + Self::e0345()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_displacement() + Self::line_moment()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{12345}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e12345()
    }
    /// The multivector of direction.
    ///
    /// ```math
    /// P_\infty \equiv X\e_{01243} + Y\e_{01235} + Z\e_{01254} + Ð\e_{01345} + Ø\e_{02354}
    /// ```
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e01243() + Self::e01235() + Self::e01254() + Self::e01345() + Self::e02354()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of single rotator $`r_1 \equiv s + v_0`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(
    ///     single_rotator.basis_blades(),
    ///     Vee::single_rotator().basis_blades()
    /// );
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(
    ///     single_rotator.basis_blades(),
    ///     Vee::single_rotator().basis_blades()
    /// );
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::volume_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(
    ///     double_rotator.basis_blades(),
    ///     Vee::double_rotator().basis_blades()
    /// );
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕+v͔v͕",
    ///     "+(+a͔v͕+a͕v͔-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔+b͔v͕+b͕v͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔+c͔v͕+c͕v͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔+d͔v͕+d͕v͔)e15",
    ///     "+(-a͔b͕+a͕b͔+e͔v͕+e͕v͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔+f͔v͕+f͕v͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔+g͔v͕+g͕v͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔+h͔v͕+h͕v͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔+i͔v͕+i͕v͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔+j͔v͕+j͕v͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(
    ///     double_rotator.basis_blades(),
    ///     Vee::double_rotator().basis_blades()
    /// );
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(
    ///     double_rotator.basis_blades(),
    ///     Vee::double_rotator().basis_blades()
    /// );
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::volume_displacement() + Self::line_displacement()
    }
    /// The multivector of translator $`t \equiv s + v_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
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
    ///     "+(+w͔Ø͕-w͕Ø͔)e05",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::volume_moment()
    }
    /// The multivector of simple single motor $`m_{s1} \equiv s + v`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let simple_single_motor = Vee::volume4().lhs() * Vee::volume4().rhs();
    ///
    /// assert_eq!(
    ///     simple_single_motor.basis_blades(),
    ///     Vee::simple_single_motor().basis_blades()
    /// );
    /// format_eq!(simple_single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+W͔ø͕-W͕ø͔)e05",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_motor() -> Self {
        Self::scalar() + Self::volume()
    }
    /// The multivector of single motor $`m_1 \equiv s + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     single_motor.basis_blades(),
    ///     Vee::single_motor().basis_blades()
    /// );
    /// format_eq!(single_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕a͔+Z͕b͔+c͔Ð͕+d͔Ø͕)e01",
    ///     "+(-X͕a͔+Y͕v͔+Z͕e͔+f͔Ð͕+g͔Ø͕)e02",
    ///     "+(-X͕b͔-Y͕e͔+Z͕v͔+h͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕c͔-Y͕f͔-Z͕h͔+j͔Ø͕+v͔Ð͕)e04",
    ///     "+(-X͕d͔-Y͕g͔-Z͕i͔-j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e12",
    ///     "+b͔v͕e13",
    ///     "+c͔v͕e14",
    ///     "+d͔v͕e15",
    ///     "+e͔v͕e23",
    ///     "+f͔v͕e24",
    ///     "+g͔v͕e25",
    ///     "+h͔v͕e34",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+(+Z͕j͔+h͔Ø͕-i͔Ð͕)e0345",
    ///     "+(-Y͕j͔-f͔Ø͕+g͔Ð͕)e0254",
    ///     "+(+Y͕i͔-Z͕g͔+e͔Ø͕)e0235",
    ///     "+(-Y͕h͔+Z͕f͔-e͔Ð͕)e0243",
    ///     "+(+X͕j͔+c͔Ø͕-d͔Ð͕)e0145",
    ///     "+(-X͕i͔+Z͕d͔-b͔Ø͕)e0153",
    ///     "+(+X͕h͔-Z͕c͔+b͔Ð͕)e0134",
    ///     "+(+X͕g͔-Y͕d͔+a͔Ø͕)e0125",
    ///     "+(-X͕f͔+Y͕c͔-a͔Ð͕)e0142",
    ///     "+(+X͕e͔-Y͕b͔+Z͕a͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(
    ///     single_motor.basis_blades(),
    ///     Vee::single_motor().basis_blades()
    /// );
    /// format_eq!(single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(-A͔y͕+A͕y͔-B͔z͕+B͕z͔-C͔ð͕+C͕ð͔-D͔ø͕+D͕ø͔)e01",
    ///     "+(+A͔x͕-A͕x͔-E͔z͕+E͕z͔-F͔ð͕+F͕ð͔-G͔ø͕+G͕ø͔)e02",
    ///     "+(+B͔x͕-B͕x͔+E͔y͕-E͕y͔-H͔ð͕+H͕ð͔-I͔ø͕+I͕ø͔)e03",
    ///     "+(+C͔x͕-C͕x͔+F͔y͕-F͕y͔+H͔z͕-H͕z͔-J͔ø͕+J͕ø͔)e04",
    ///     "+(+D͔x͕-D͕x͔+G͔y͕-G͕y͔+I͔z͕-I͕z͔+J͔ð͕-J͕ð͔)e05",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(-H͔ø͕-H͕ø͔+I͔ð͕+I͕ð͔-J͔z͕-J͕z͔)e0345",
    ///     "+(+F͔ø͕+F͕ø͔-G͔ð͕-G͕ð͔+J͔y͕+J͕y͔)e0254",
    ///     "+(-E͔ø͕-E͕ø͔+G͔z͕+G͕z͔-I͔y͕-I͕y͔)e0235",
    ///     "+(+E͔ð͕+E͕ð͔-F͔z͕-F͕z͔+H͔y͕+H͕y͔)e0243",
    ///     "+(-C͔ø͕-C͕ø͔+D͔ð͕+D͕ð͔-J͔x͕-J͕x͔)e0145",
    ///     "+(+B͔ø͕+B͕ø͔-D͔z͕-D͕z͔+I͔x͕+I͕x͔)e0153",
    ///     "+(-B͔ð͕-B͕ð͔+C͔z͕+C͕z͔-H͔x͕-H͕x͔)e0134",
    ///     "+(-A͔ø͕-A͕ø͔+D͔y͕+D͕y͔-G͔x͕-G͕x͔)e0125",
    ///     "+(+A͔ð͕+A͕ð͔-C͔y͕-C͕y͔+F͔x͕+F͕x͔)e0142",
    ///     "+(-A͔z͕-A͕z͔+B͔y͕+B͕y͔-E͔x͕-E͕x͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line_moment()
    }
    /// The multivector of simple double motor $`m_{s2} \equiv s + v + \ell`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let simple_double_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(
    ///     simple_double_motor.basis_blades(),
    ///     Vee::simple_double_motor().basis_blades()
    /// );
    /// format_eq!(simple_double_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-Y͔a͕+Y͕a͔-Z͔b͕+Z͕b͔+c͔Ð͕-c͕Ð͔+d͔Ø͕-d͕Ø͔)e01",
    ///     "+(+X͔a͕-X͕a͔-Z͔e͕+Z͕e͔+f͔Ð͕-f͕Ð͔+g͔Ø͕-g͕Ø͔)e02",
    ///     "+(+X͔b͕-X͕b͔+Y͔e͕-Y͕e͔+h͔Ð͕-h͕Ð͔+i͔Ø͕-i͕Ø͔)e03",
    ///     "+(+X͔c͕-X͕c͔+Y͔f͕-Y͕f͔+Z͔h͕-Z͕h͔+j͔Ø͕-j͕Ø͔)e04",
    ///     "+(+X͔d͕-X͕d͔+Y͔g͕-Y͕g͔+Z͔i͕-Z͕i͔-j͔Ð͕+j͕Ð͔)e05",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    ///     "+(+Z͔j͕+Z͕j͔+h͔Ø͕+h͕Ø͔-i͔Ð͕-i͕Ð͔)e0345",
    ///     "+(-Y͔j͕-Y͕j͔-f͔Ø͕-f͕Ø͔+g͔Ð͕+g͕Ð͔)e0254",
    ///     "+(+Y͔i͕+Y͕i͔-Z͔g͕-Z͕g͔+e͔Ø͕+e͕Ø͔)e0235",
    ///     "+(-Y͔h͕-Y͕h͔+Z͔f͕+Z͕f͔-e͔Ð͕-e͕Ð͔)e0243",
    ///     "+(+X͔j͕+X͕j͔+c͔Ø͕+c͕Ø͔-d͔Ð͕-d͕Ð͔)e0145",
    ///     "+(-X͔i͕-X͕i͔+Z͔d͕+Z͕d͔-b͔Ø͕-b͕Ø͔)e0153",
    ///     "+(+X͔h͕+X͕h͔-Z͔c͕-Z͕c͔+b͔Ð͕+b͕Ð͔)e0134",
    ///     "+(+X͔g͕+X͕g͔-Y͔d͕-Y͕d͔+a͔Ø͕+a͕Ø͔)e0125",
    ///     "+(-X͔f͕-X͕f͔+Y͔c͕+Y͕c͔-a͔Ð͕-a͕Ð͔)e0142",
    ///     "+(+X͔e͕+X͕e͔-Y͔b͕-Y͕b͔+Z͔a͕+Z͕a͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line()
    }
    /// The multivector of double motor $`m_2 \equiv s + v + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕a͔+Z͕b͔+c͔Ð͕+d͔Ø͕)e01",
    ///     "+(-X͕a͔+Y͕v͔+Z͕e͔+f͔Ð͕+g͔Ø͕)e02",
    ///     "+(-X͕b͔-Y͕e͔+Z͕v͔+h͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕c͔-Y͕f͔-Z͕h͔+j͔Ø͕+v͔Ð͕)e04",
    ///     "+(-X͕d͔-Y͕g͔-Z͕i͔-j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e12",
    ///     "+b͔v͕e13",
    ///     "+c͔v͕e14",
    ///     "+d͔v͕e15",
    ///     "+e͔v͕e23",
    ///     "+f͔v͕e24",
    ///     "+g͔v͕e25",
    ///     "+h͔v͕e34",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+v͕x͔e2345",
    ///     "+v͕y͔e1354",
    ///     "+v͕z͔e1245",
    ///     "+v͕ð͔e1253",
    ///     "+v͕ø͔e1234",
    ///     "+(+X͕y͔-Y͕x͔+Z͕j͔+h͔Ø͕-i͔Ð͕)e0345",
    ///     "+(+X͕z͔-Y͕j͔-Z͕x͔-f͔Ø͕+g͔Ð͕)e0254",
    ///     "+(+X͕ð͔+Y͕i͔-Z͕g͔+e͔Ø͕-x͔Ð͕)e0235",
    ///     "+(+X͕ø͔-Y͕h͔+Z͕f͔-e͔Ð͕-x͔Ø͕)e0243",
    ///     "+(+X͕j͔+Y͕z͔-Z͕y͔+c͔Ø͕-d͔Ð͕)e0145",
    ///     "+(-X͕i͔+Y͕ð͔+Z͕d͔-b͔Ø͕-y͔Ð͕)e0153",
    ///     "+(+X͕h͔+Y͕ø͔-Z͕c͔+b͔Ð͕-y͔Ø͕)e0134",
    ///     "+(+X͕g͔-Y͕d͔+Z͕ð͔+a͔Ø͕-z͔Ð͕)e0125",
    ///     "+(-X͕f͔+Y͕c͔+Z͕ø͔-a͔Ð͕-z͔Ø͕)e0142",
    ///     "+(+X͕e͔-Y͕b͔+Z͕a͔+Ð͕ø͔-Ø͕ð͔)e0123",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔+Ø͕ø͔)I",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-E͔j͕+E͕j͔+F͔i͕-F͕i͔-G͔h͕+G͕h͔-H͔g͕+H͕g͔+I͔f͕-I͕f͔-J͔e͕+J͕e͔)e01",
    ///     "+(+B͔j͕-B͕j͔-C͔i͕+C͕i͔+D͔h͕-D͕h͔+H͔d͕-H͕d͔-I͔c͕+I͕c͔+J͔b͕-J͕b͔)e02",
    ///     "+(-A͔j͕+A͕j͔+C͔g͕-C͕g͔-D͔f͕+D͕f͔-F͔d͕+F͕d͔+G͔c͕-G͕c͔-J͔a͕+J͕a͔)e03",
    ///     "+(+A͔i͕-A͕i͔-B͔g͕+B͕g͔+D͔e͕-D͕e͔+E͔d͕-E͕d͔-G͔b͕+G͕b͔+I͔a͕-I͕a͔)e04",
    ///     "+(-A͔h͕+A͕h͔+B͔f͕-B͕f͔-C͔e͕+C͕e͔-E͔c͕+E͕c͔+F͔b͕-F͕b͔-H͔a͕+H͕a͔)e05",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    ///     "+(-B͔e͕-B͕e͔-C͔f͕-C͕f͔-D͔g͕-D͕g͔+E͔b͕+E͕b͔+F͔c͕+F͕c͔+G͔d͕+G͕d͔)e0345",
    ///     "+(+A͔e͕+A͕e͔-C͔h͕-C͕h͔-D͔i͕-D͕i͔-E͔a͕-E͕a͔+H͔c͕+H͕c͔+I͔d͕+I͕d͔)e0254",
    ///     "+(+A͔f͕+A͕f͔+B͔h͕+B͕h͔-D͔j͕-D͕j͔-F͔a͕-F͕a͔-H͔b͕-H͕b͔+J͔d͕+J͕d͔)e0235",
    ///     "+(+A͔g͕+A͕g͔+B͔i͕+B͕i͔+C͔j͕+C͕j͔-G͔a͕-G͕a͔-I͔b͕-I͕b͔-J͔c͕-J͕c͔)e0243",
    ///     "+(-A͔b͕-A͕b͔+B͔a͕+B͕a͔-F͔h͕-F͕h͔-G͔i͕-G͕i͔+H͔f͕+H͕f͔+I͔g͕+I͕g͔)e0145",
    ///     "+(-A͔c͕-A͕c͔+C͔a͕+C͕a͔+E͔h͕+E͕h͔-G͔j͕-G͕j͔-H͔e͕-H͕e͔+J͔g͕+J͕g͔)e0153",
    ///     "+(-A͔d͕-A͕d͔+D͔a͕+D͕a͔+E͔i͕+E͕i͔+F͔j͕+F͕j͔-I͔e͕-I͕e͔-J͔f͕-J͕f͔)e0134",
    ///     "+(-B͔c͕-B͕c͔+C͔b͕+C͕b͔-E͔f͕-E͕f͔+F͔e͕+F͕e͔-I͔j͕-I͕j͔+J͔i͕+J͕i͔)e0125",
    ///     "+(-B͔d͕-B͕d͔+D͔b͕+D͕b͔-E͔g͕-E͕g͔+G͔e͕+G͕e͔+H͔j͕+H͕j͔-J͔h͕-J͕h͔)e0142",
    ///     "+(-C͔d͕-C͕d͔+D͔c͕+D͕c͔-F͔g͕-F͕g͔+G͔f͕+G͕f͔-H͔i͕-H͕i͔+I͔h͕+I͕h͔)e0123",
    ///     "+(+A͔a͕-A͕a͔+B͔b͕-B͕b͔+C͔c͕-C͕c͔+D͔d͕-D͕d͔+E͔e͕-E͕e͔+F͔f͕-F͕f͔+G͔g͕-G͕g͔+H͔h͕-H͕h͔+I͔i͕-I͕i͔+J͔j͕-J͕j͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::double_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::point() + Vee::volume4()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(+(+ay+bz+cð+dø+ej-fi+gh-vx)X͓+(-ax-bj+ci-dh+ez+fð+gø-vy)Y͓\
    ///         +(+aj-bx-cg+df-ey+hð+iø-vz)Z͓+(-Aa-Bb-Cc-Dd-Ee-Ff-Gg-Hh-Ii-Jj+Vv+Xx+Yy+Zz+Ðð+Øø)w͓\
    ///         +(-ai+bg-cx-de-fy-hz+jø-vð)Ð͓+(+ah-bf+ce-dx-gy-iz-jð-vø)Ø͓)e0",
    ///     "+2(+ay+bz+cð+dø-ej+fi-gh+vx)w͓e1",
    ///     "+2(-ax+bj-ci+dh+ez+fð+gø+vy)w͓e2",
    ///     "+2(-aj-bx+cg-df-ey+hð+iø+vz)w͓e3",
    ///     "+2(+ai-bg-cx+de-fy-hz+jø+vð)w͓e4",
    ///     "+2(-ah+bf-ce-dx-gy-iz-jð+vø)w͓e5",
    ///     "+(+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+vv+xx+yy+zz+ðð+øø)w͓e12345",
    ///     "+(+(-aa-bb-cc-dd+ee+ff+gg+hh+ii+jj+vv+xx-yy-zz-ðð-øø)X͓+2(+av-be-cf-dg+hø-ið+jz+xy)Y͓\
    ///        +2(+ae+bv-ch-di-fø+gð-jy+xz)Z͓+2(+Ay+Bz+Cð+Dø-Ej+Fi-Gh-Hg+If-Je-Vx-Xv-Ya-Zb-cÐ-dØ)w͓\
    ///        +2(+af+bh+cv-dj+eø-gz+iy+xð)Ð͓+2(+ag+bi+cj+dv-eð+fz-hy+xø)Ø͓)e02354",
    ///     "+(+2(-av-be-cf-dg-hø+ið-jz+xy)X͓+(-aa+bb+cc+dd-ee-ff-gg+hh+ii+jj+vv-xx+yy-zz-ðð-øø)Y͓\
    ///        +2(-ab+cø-dð+ev-fh-gi+jx+yz)Z͓+2(-Ax+Bj-Ci+Dh+Ez+Fð+Gø+Hd-Ic+Jb-Vy+Xa-Yv-Ze-fÐ-gØ)w͓\
    ///        +2(-ac-bø+dz+eh+fv-gj-ix+yð)Ð͓+2(-ad+bð-cz+ei+fj+gv+hx+yø)Ø͓)e01345",
    ///     "+(+2(+ae-bv-ch-di+fø-gð+jy+xz)X͓+2(-ab-cø+dð-ev-fh-gi-jx+yz)Y͓\
    ///        +(+aa-bb+cc+dd-ee+ff+gg-hh-ii+jj+vv-xx-yy+zz-ðð-øø)Z͓\
    ///        +2(-Aj-Bx+Cg-Df-Ey-Fd+Gc+Hð+Iø-Ja-Vz+Xb+Ye-Zv-hÐ-iØ)w͓+2(+aø-bc-dy-ef+gx+hv-ij+zð)Ð͓\
    ///        +2(-að-bd+cy-eg-fx+hj+iv+zø)Ø͓)e01254",
    ///     "+(+2(+af+bh-cv-dj-eø+gz-iy+xð)X͓+2(-ac+bø-dz+eh-fv-gj+ix+yð)Y͓\
    ///        +2(-aø-bc+dy-ef-gx-hv-ij+zð)Z͓+2(+Ai-Bg-Cx+De+Ed-Fy-Gb-Hz+Ia+Jø-Vð+Xc+Yf+Zh-jØ-vÐ)w͓\
    ///        +(+aa+bb-cc+dd+ee-ff+gg-hh+ii-jj+vv-xx-yy-zz+ðð-øø)Ð͓\
    ///        +2(+az-by-cd+ex-fg-hi+jv+ðø)Ø͓)e01235",
    ///     "+(+2(+ag+bi+cj-dv+eð-fz+hy+xø)X͓+2(-ad-bð+cz+ei+fj-gv-hx+yø)Y͓\
    ///        +2(+að-bd-cy-eg+fx+hj-iv+zø)Z͓+2(-Ah+Bf-Ce-Dx-Ec+Fb-Gy-Ha-Iz-Jð-Vø+Xd+Yg+Zi+jÐ-vØ)w͓\
    ///        +2(-az+by-cd-ex-fg-hi-jv+ðø)Ð͓\
    ///        +(+aa+bb+cc-dd+ee+ff-gg+hh-ii-jj+vv-xx-yy-zz-ðð+øø)Ø͓)e01243",
    /// ]);
    ///
    /// let norm_squared = Vee::double_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+vv+xx+yy+zz+ðð+øø",
    ///     "+2(+ay+bz+cð+dø-ej+fi-gh+vx)e2345",
    ///     "+2(-ax+bj-ci+dh+ez+fð+gø+vy)e1354",
    ///     "+2(-aj-bx+cg-df-ey+hð+iø+vz)e1245",
    ///     "+2(+ai-bg-cx+de-fy-hz+jø+vð)e1253",
    ///     "+2(-ah+bf-ce-dx-gy-iz-jð+vø)e1234",
    ///     "+2(+Av+Be+Cf+Dg-Eb-Fc-Gd-Hø+Ið-Jz+Va-Xy+Yx-Zj-hØ+iÐ)e0345",
    ///     "+2(-Ae+Bv+Ch+Di+Ea+Fø-Gð-Hc-Id+Jy+Vb-Xz+Yj+Zx+fØ-gÐ)e0254",
    ///     "+2(-Af-Bh+Cv+Dj-Eø+Fa+Gz+Hb-Iy-Jd+Vc-Xð-Yi+Zg-eØ+xÐ)e0235",
    ///     "+2(-Ag-Bi-Cj+Dv+Eð-Fz+Ga+Hy+Ib+Jc+Vd-Xø+Yh-Zf+eÐ+xØ)e0243",
    ///     "+2(+Ab-Ba-Cø+Dð+Ev+Fh+Gi-Hf-Ig-Jx+Ve-Xj-Yz+Zy-cØ+dÐ)e0145",
    ///     "+2(+Ac+Bø-Ca-Dz-Eh+Fv+Gj+He+Ix-Jg+Vf+Xi-Yð-Zd+bØ+yÐ)e0153",
    ///     "+2(+Ad-Bð+Cz-Da-Ei-Fj+Gv-Hx+Ie+Jf+Vg-Xh-Yø+Zc-bÐ+yØ)e0134",
    ///     "+2(-Aø+Bc-Cb+Dy+Ef-Fe-Gx+Hv+Ij-Ji+Vh-Xg+Yd-Zð-aØ+zÐ)e0125",
    ///     "+2(+Að+Bd-Cy-Db+Eg+Fx-Ge-Hj+Iv+Jh+Vi+Xf-Yc-Zø+aÐ+zØ)e0142",
    ///     "+2(-Az+By+Cd-Dc-Ex+Fg-Gf+Hi-Ih+Jv+Vj-Xe+Yb-Za-Ðø+Øð)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
    /// The multivector of single rotoreflector $`f_{r1} \equiv v^4_0 + p_0`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let single_rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(
    ///     single_rotoreflector.basis_blades(),
    ///     Vee::single_rotoreflector().basis_blades()
    /// );
    /// format_eq!(single_rotoreflector, [
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement()
    }
    /// The multivector of double rotoreflector $`f_{r2} \equiv v^4_0 + p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let double_rotoreflector = Vee::normal().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(
    ///     double_rotoreflector.basis_blades(),
    ///     Vee::double_rotoreflector().basis_blades()
    /// );
    /// format_eq!(double_rotoreflector, [
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv v^4 + p_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(
    ///     transflector.basis_blades(),
    ///     Vee::transflector().basis_blades()
    /// );
    /// format_eq!(transflector, [
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+v͕ø͔e5",
    ///     "+(+X͕y͔-Y͕x͔)e012",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+Y͕z͔-Z͕y͔)e023",
    ///     "+(+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ð͕ø͔-Ø͕ð͔)e045",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume4() + Self::plane_moment()
    }
    /// The multivector of simple single flector $`f_{s1} \equiv v^4 + p`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let simple_single_flector = Vee::volume4().lhs() * Vee::simple_single_motor().rhs();
    ///
    /// assert_eq!(
    ///     simple_single_flector.basis_blades(),
    ///     Vee::simple_single_flector().basis_blades()
    /// );
    /// format_eq!(simple_single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_flector() -> Self {
        Self::volume4() + Self::plane()
    }
    /// The multivector of single flector $`f_1 \equiv v^4 + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let single_flector = Vee::volume4().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(
    ///     single_flector.basis_blades(),
    ///     Vee::single_flector().basis_blades()
    /// );
    /// format_eq!(single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::volume4() + Self::plane() + Self::direction()
    }
    /// The multivector of simple double flector $`f_{s2} \equiv v^4 + p + P`$.
    ///
    /// ```
    /// use vee::{PgaP5 as Vee, format_eq};
    ///
    /// let simple_double_flector = Vee::volume4().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(
    ///     simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades()
    /// );
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔-W͔x͕)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔-W͔y͕)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔-W͔z͕)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔-W͔ð͕)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔-W͔ø͕)e01243",
    /// ]);
    ///
    /// let simple_double_flector = Vee::volume4().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(
    ///     simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades()
    /// );
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔+V͕x͔-W͔x͕)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔+V͕y͔-W͔y͕)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔+V͕z͔-W͔z͕)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔+V͕ð͔-W͔ð͕)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔+V͕ø͔-W͔ø͕)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_flector() -> Self {
        Self::volume4() + Self::plane() + Self::point()
    }
}
