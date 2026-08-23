// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(missing_docs)]

use std::{
    fs::{read_to_string, write},
    path::Path,
};
use vee::pga::Pga;

macro_rules! mul {
    ($($f:ident: $t:expr, ($m:expr, $n:expr);)*) => {
        $(
            #[test]
            fn $f() {
                mul($t, &Pga::<$m, $n>::table().unwrap());
            }
        )*
    };
}

fn mul(alias: &str, table: &str) {
    let path = Path::new("tests")
        .join("mul")
        .join(alias)
        .with_extension("ct");
    if let Ok(text) = read_to_string(&path) {
        assert_eq!(table, text);
    } else {
        write(&path, table).unwrap();
    }
}

const E: i8 = 1;
const H: i8 = -1;
const P: i8 = 0;

#[cfg(feature = "rudimentary")]
mul! {
    pga_e0: "PgaE0", (E, 0);
    pga_h0: "PgaH0", (H, 0);
    pga_p0: "PgaP0", (P, 0);

    pga_e1: "PgaE1", (E, 1);
    pga_h1: "PgaH1", (H, 1);
    pga_p1: "PgaP1", (P, 1);
}

mul! {
    pga_e2: "PgaE2", (E, 2);
    pga_h2: "PgaH2", (H, 2);
    pga_p2: "PgaP2", (P, 2);

    pga_e3: "PgaE3", (E, 3);
    pga_h3: "PgaH3", (H, 3);
    pga_p3: "PgaP3", (P, 3);

    pga_e4: "PgaE4", (E, 4);
    pga_h4: "PgaH4", (H, 4);
    pga_p4: "PgaP4", (P, 4);
}

#[cfg(feature = "exploratory")]
mul! {
    pga_e5: "PgaE5", (E, 5);
    pga_h5: "PgaH5", (H, 5);
    pga_p5: "PgaP5", (P, 5);

    pga_e6: "PgaE6", (E, 6);
    pga_h6: "PgaH6", (H, 6);
    pga_p6: "PgaP6", (P, 6);

    pga_e7: "PgaE7", (E, 7);
    pga_h7: "PgaH7", (H, 7);
    pga_p7: "PgaP7", (P, 7);
}
