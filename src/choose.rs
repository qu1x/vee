// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// The binomial coefficient `n` over `k` up to `n = 34` for all `k`.
///
/// ```
/// use vee::choose;
///
/// assert_eq!(choose(34, 16), 2_203_961_430);
/// assert_eq!(choose(34, 17), 2_333_606_220);
/// assert_eq!(choose(34, 18), 2_203_961_430);
/// ```
///
/// # Panics
///
/// Panics on overflow.
///
/// ```should_panic
/// use vee::choose;
///
/// let _ = choose(35, 17);
/// ```
#[must_use]
pub const fn choose(n: u32, k: u32) -> u32 {
    if k > n {
        0
    } else {
        let (n, k) = (n as u64, if k > n - k { n - k } else { k } as u64);
        let (mut r, mut i) = (1, 0);
        while i < k {
            r = r * (n - i) / (i + 1);
            i += 1;
        }
        assert!(
            r <= u32::MAX as u64,
            "attempt to calculate the binomial coefficient with overflow"
        );
        #[allow(clippy::cast_possible_truncation)]
        {
            r as u32
        }
    }
}
