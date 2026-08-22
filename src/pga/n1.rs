// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB1, LUT1, 1, [
    ('v', e),
    ('W', e0),
    ('w', e1),
    ('V', e01),
]);

/// The named entities of the PGA with embedded dimension $`N = 1`$.
impl<const M: i8> Multivector<Pga<M, 1>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{01}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of weight $`P_0 \equiv w\e_1`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e1()
    }
    /// The multivector of direction $`P_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e0()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of translator $`t \equiv s + S`$.
    ///
    /// ```
    /// use vee::{PgaP1 as Vee, format_eq};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "+w͔w͕",
    ///     "+(+W͔w͕-W͕w͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::translator();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+vvW͓+2Vvw͓)e0",
    ///     "+vvw͓e1",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::norm()
    }
}
