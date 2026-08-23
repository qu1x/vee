// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(missing_docs)]

use vee::{Algebra, Multivector, PgaP2, PgaP3, PgaP4};
#[cfg(feature = "rudimentary")]
use vee::{PgaP0, PgaP1};
#[cfg(feature = "exploratory")]
use vee::{PgaP5, PgaP6, PgaP7};

fn mul<B: Algebra>(vec: &[Multivector<B>], ops: &[(usize, usize)]) {
    let mut idx = 0;
    for row in 0..vec.len() {
        for col in 0..=row {
            assert_eq!(
                (vec[row].clone().lhs() * vec[col].clone().rhs()).ops(),
                ops[idx],
                "(vec[{row}] * vec[{col}]).ops() == ops[{idx}]"
            );
            idx += 1;
        }
    }
}

#[cfg(feature = "rudimentary")]
#[test]
fn pga_p0() {
    #[rustfmt::skip]
    let vec = [
        PgaP0::scalar(),
        PgaP0::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        (1,0),
        (1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[cfg(feature = "rudimentary")]
#[test]
fn pga_p1() {
    #[rustfmt::skip]
    let vec = [
        PgaP1::scalar(),
        PgaP1::point(),
        PgaP1::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        (1,0),
        (2,0),(3,1),
        (1,0),(1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[test]
fn pga_p2() {
    let vec = [
        PgaP2::scalar(),
        PgaP2::line(),
        PgaP2::point(),
        PgaP2::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        (1,0),
        (3,0),(8,4),
        (3,0),(7,3),(5,2),
        (1,0),(2,0),(1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[test]
fn pga_p3() {
    let vec = [
        PgaP3::scalar(),
        PgaP3::plane(),
        PgaP3::line(),
        PgaP3::point(),
        PgaP3::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        (1,0),
        (4,0),(15, 8),
        (6,0),(21,13),(27,19),
        (4,0),(13, 6),(15, 8),(7,3),
        (1,0),( 3, 0),( 3, 0),(1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[test]
fn pga_p4() {
    let vec = [
        PgaP4::scalar(),
        PgaP4::volume(),
        PgaP4::plane(),
        PgaP4::line(),
        PgaP4::point(),
        PgaP4::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        ( 1,0),
        ( 5,0),(24,13),
        (10,0),(46,31),(84,68),
        (10,0),(44,29),(76,60),(64,49),
        ( 5,0),(21,10),(34,20),(26,15),(9,4),
        ( 1,0),( 4, 0),( 6, 0),( 4, 0),(1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p5() {
    let vec = [
        PgaP5::scalar(),
        PgaP5::volume4(),
        PgaP5::volume(),
        PgaP5::plane(),
        PgaP5::line(),
        PgaP5::point(),
        PgaP5::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        ( 1,0),
        ( 6,0),( 35,19),
        (15,0),( 85,59),(200,169),
        (20,0),(110,80),(250,218),(300,268),
        (15,0),( 80,54),(175,144),(200,169),(125,99),
        ( 6,0),( 31,15),( 65, 40),( 70, 45),( 40,24),(11,5),
        ( 1,0),(  5, 0),( 10,  0),( 10,  0),(  5, 0),( 1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p6() {
    let vec = [
        PgaP6::scalar(),
        PgaP6::volume5(),
        PgaP6::volume4(),
        PgaP6::volume(),
        PgaP6::plane(),
        PgaP6::line(),
        PgaP6::point(),
        PgaP6::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        ( 1,0),
        ( 7,0),( 48, 26),
        (21,0),(141, 99),(405,348),
        (35,0),(230,174),(645,582),(1000,936),
        (35,0),(225,169),(615,552),( 925,861),(825,762),
        (21,0),(132, 90),(351,294),( 510,448),(435,378),(216,174),
        ( 7,0),( 43, 21),(111, 70),( 155,105),(125, 84),( 57, 35),(13,6),
        ( 1,0),(  6,  0),( 15,  0),(  20,  0),( 15,  0),(  6,  0),( 1,0),(0,0),
    ];
    mul(&vec, &ops);
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p7() {
    let vec = [
        PgaP7::scalar(),
        PgaP7::volume6(),
        PgaP7::volume5(),
        PgaP7::volume4(),
        PgaP7::volume(),
        PgaP7::plane(),
        PgaP7::line(),
        PgaP7::point(),
        PgaP7::pseudoscalar(),
    ];
    #[rustfmt::skip]
    let ops = [
        (1, 0),
        (8, 0),( 63, 34),
        (28,0),(217,153),( 735, 636),
        (56,0),(427,329),(1421,1301),(2695,2568),
        (70,0),(525,413),(1715,1589),(3185,3057),(3675,3547),
        (56,0),(413,315),(1323,1203),(2401,2274),(2695,2568),(1911,1791),
        (28,0),(203,139),( 637, 538),(1127,1008),(1225,1106),( 833, 734),(343,279),
        (8, 0),( 57, 28),( 175, 112),( 301, 210),( 315, 224),( 203, 140),( 77, 48),(15,7),
        (1, 0),(  7,  0),(  21,   0),(  35,   0),(  35,   0),(  21,   0),(  7,  0),( 1,0),(0,0),
    ];
    mul(&vec, &ops);
}
