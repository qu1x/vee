// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{BasisBlade, Multivector, Pga, Symbol};

#[rustfmt::skip]
basis!(TAB0, LUT0, 0, [
    ('v', e),
    ('V', e0),
]);

/// The named entities of the PGA with embedded dimension $`N = 0`$.
impl<const M: i8> Multivector<Pga<M, 0>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_0`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
}
