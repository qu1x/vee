// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(missing_docs)]

use vee::choose;

fn checked_triangle(n: u32) -> Option<Vec<Vec<u32>>> {
    let n = usize::try_from(n).ok()?.checked_add(1)?;
    let mut rows = Vec::with_capacity(n);
    rows.push(vec![1u32]);
    for n in 0..n - 1 {
        let last = &rows[n];
        let k = last.len() + 1;
        let mut next = Vec::with_capacity(k);
        next.push(1);
        for k in 1..k - 1 {
            next.push(last[k - 1].checked_add(last[k])?);
        }
        next.push(1);
        rows.push(next);
    }
    Some(rows)
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn triangle() {
    let triangle = checked_triangle(34).unwrap();
    for n in 0..triangle.len() as u32 {
        for k in 0..=n {
            assert_eq!(
                choose(n, k),
                triangle[n as usize][k as usize],
                "mismatch at (n, k) = ({n}, {k})"
            );
        }
    }
}
