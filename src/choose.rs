// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Finds the binomial coefficient.
pub trait Choose {
    /// The output type.
    type Output;

    /// Finds the binomial coefficient as in `self` over `other`.
    #[must_use]
    fn choose(self, other: Self) -> Self::Output;
}

impl Choose for u32 {
    type Output = Self;

    fn choose(self, other: Self) -> Self::Output {
        let (n, k) = (self, other);
        if n < k {
            0
        } else {
            (0..k).fold(1, |r, i| r * (n - i) / (i + 1))
        }
    }
}
