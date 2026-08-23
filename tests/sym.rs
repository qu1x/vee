// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(missing_docs)]

#[cfg(feature = "rudimentary")]
use vee::{PgaP0, PgaP1};
use vee::{PgaP2, PgaP3, PgaP4};
#[cfg(feature = "exploratory")]
use vee::{PgaP5, PgaP6, PgaP7};

#[cfg(feature = "rudimentary")]
#[test]
#[allow(clippy::iter_on_single_items)]
fn pga_p0() {
    #[rustfmt::skip]
    assert!(
        [
            PgaP0::norm(),
        ]
        .iter()
        .all(PgaP0::is_entity)
    );
}

#[cfg(feature = "rudimentary")]
#[test]
fn pga_p1() {
    #[rustfmt::skip]
    assert!(
        [
            PgaP1::point(),
            PgaP1::norm(),
            PgaP1::translator(),
        ]
        .iter()
        .all(PgaP1::is_entity)
    );
}

#[test]
fn pga_p2() {
    assert!(
        [
            PgaP2::line(),
            PgaP2::point(),
            PgaP2::norm(),
            PgaP2::rotator(),
            PgaP2::translator(),
            PgaP2::motor(),
            PgaP2::rotoreflector(),
            PgaP2::flector(),
        ]
        .iter()
        .all(PgaP2::is_entity)
    );
}

#[test]
fn pga_p3() {
    assert!(
        [
            PgaP3::plane(),
            PgaP3::line(),
            PgaP3::point(),
            PgaP3::norm(),
            PgaP3::rotator(),
            PgaP3::translator(),
            PgaP3::simple_motor(),
            PgaP3::motor(),
            PgaP3::rotoreflector(),
            PgaP3::transflector(),
            PgaP3::flector(),
        ]
        .iter()
        .all(PgaP3::is_entity)
    );
}

#[test]
fn pga_p4() {
    assert!(
        [
            PgaP4::volume(),
            PgaP4::plane(),
            PgaP4::line(),
            PgaP4::point(),
            PgaP4::norm(),
            PgaP4::single_rotator(),
            PgaP4::double_rotator(),
            PgaP4::translator(),
            PgaP4::simple_motor(),
            PgaP4::single_motor(),
            PgaP4::double_motor(),
            PgaP4::rotoreflector(),
            PgaP4::transflector(),
            PgaP4::simple_flector(),
            PgaP4::flector(),
        ]
        .iter()
        .all(PgaP4::is_entity)
    );
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p5() {
    assert!(
        [
            PgaP5::volume4(),
            PgaP5::volume(),
            PgaP5::plane(),
            PgaP5::line(),
            PgaP5::point(),
            PgaP5::norm(),
            PgaP5::single_rotator(),
            PgaP5::double_rotator(),
            PgaP5::translator(),
            PgaP5::simple_single_motor(),
            PgaP5::single_motor(),
            PgaP5::double_motor(),
            PgaP5::single_rotoreflector(),
            PgaP5::double_rotoreflector(),
            PgaP5::transflector(),
            PgaP5::simple_single_flector(),
            PgaP5::single_flector(),
            PgaP5::simple_double_flector(),
        ]
        .iter()
        .all(PgaP5::is_entity)
    );
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p6() {
    assert!(
        [
            PgaP6::volume5(),
            PgaP6::volume4(),
            PgaP6::volume(),
            PgaP6::plane(),
            PgaP6::line(),
            PgaP6::point(),
            PgaP6::norm(),
            PgaP6::single_rotator(),
            PgaP6::double_rotator(),
            PgaP6::triple_rotator(),
            PgaP6::translator(),
            PgaP6::single_motor(),
            PgaP6::simple_double_motor(),
            PgaP6::double_motor(),
            PgaP6::simple_triple_motor(),
            PgaP6::triple_motor(),
            PgaP6::single_rotoreflector(),
            PgaP6::double_rotoreflector(),
            PgaP6::transflector(),
            PgaP6::single_flector(),
            PgaP6::simple_double_flector(),
            PgaP6::double_flector(),
            PgaP6::triple_flector(),
        ]
        .iter()
        .all(PgaP6::is_entity)
    );
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p7() {
    assert!(
        [
            PgaP7::volume6(),
            PgaP7::volume5(),
            PgaP7::volume4(),
            PgaP7::volume(),
            PgaP7::plane(),
            PgaP7::line(),
            PgaP7::point(),
            PgaP7::norm(),
            PgaP7::single_rotator(),
            PgaP7::double_rotator(),
            PgaP7::triple_rotator(),
            PgaP7::translator(),
            PgaP7::simple_single_motor(),
            PgaP7::single_motor(),
            PgaP7::simple_double_motor(),
            PgaP7::double_motor(),
            PgaP7::simple_triple_motor(),
            PgaP7::triple_motor(),
            PgaP7::single_rotoreflector(),
            PgaP7::double_rotoreflector(),
            PgaP7::triple_rotoreflector(),
            PgaP7::transflector(),
            PgaP7::simple_single_flector(),
            PgaP7::single_flector(),
            PgaP7::simple_double_flector(),
            PgaP7::double_flector(),
            PgaP7::simple_triple_flector(),
        ]
        .iter()
        .all(PgaP7::is_entity)
    );
}
