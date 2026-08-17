// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Choose, Factor, Symbol};
use core::{
    fmt::{Debug, Display},
    ops::{Mul, Not},
};

/// A geometric algebra defined by a flavor's basis (i.e., all its basis blades).
///
/// Implementations for this trait (e.g., [`Pga`]) define the basis of a particular flavor from
/// which the Cayley table can be constructed. The generators of a basis (i.e., the basis blades
/// of <code>[Self::grade()] == 1</code>) define the whole algebra whereas the required [`Mul`]
/// operator implements the signature-aware anti-commutative multiplication for the chosen storage.
///
/// [`Pga`]: `super::pga::Pga`
pub trait Algebra
where
    Self: Copy
        + Clone
        + Eq
        + PartialEq
        + Ord
        + PartialOrd
        + Default
        + Into<Symbol>
        + TryFrom<Symbol, Error = Symbol>
        + Debug
        + Display
        + Mul<Output = (i8, Self)>
        + Not<Output = (i8, Self)>,
{
    /// The embedded dimension.
    ///
    /// Not to be confused with the embedding dimension, e.g., `N + 1` is the embedding dimension
    /// for one-up flavors of embedded dimension `N`.
    const N: u32;

    /// The ordered basis (i.e., all basis blades).
    #[must_use]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self>;
    /// The scalar.
    #[must_use]
    fn scalar() -> Self;
    /// The pseudoscalar.
    #[must_use]
    fn pseudoscalar() -> Self;
    /// The grade.
    #[must_use]
    fn grade(&self) -> u32;
    /// The number $`n`$ of basis blades with the same grade $`g`$ of this basis blade.
    ///
    /// ```math
    /// n = { N + 1 \choose g }
    /// ```
    #[must_use]
    fn blade_len(&self) -> usize;
    /// The reverse.
    #[must_use]
    #[inline]
    fn rev(self) -> (i8, Self) {
        (1 - i8::from(self.grade().choose(2).is_odd()) * 2, self)
    }
}
