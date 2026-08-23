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
fn pga_p0() {
    assert_eq!(!PgaP0::norm(), PgaP0::norm().swp());
}

#[cfg(feature = "rudimentary")]
#[test]
fn pga_p1() {
    assert_eq!(
        !PgaP1::point(),
        (PgaP1::weight() - PgaP1::direction()).swp()
    );
    assert_eq!(!!PgaP1::point(), -PgaP1::point());
}

#[test]
fn pga_p2() {
    assert_eq!(!PgaP2::line(), PgaP2::point().swp());
    assert_eq!(!PgaP2::point(), PgaP2::line().swp());
}

#[test]
fn pga_p3() {
    assert_eq!(!PgaP3::plane(), PgaP3::point().swp());
    assert_eq!(!PgaP3::line(), PgaP3::line().swp());
    assert_eq!(!PgaP3::point(), -PgaP3::plane().swp());
}

#[test]
fn pga_p4() {
    assert_eq!(!PgaP4::volume(), PgaP4::point().swp());
    assert_eq!(!PgaP4::plane(), PgaP4::line().swp());
    assert_eq!(!PgaP4::line(), PgaP4::plane().swp());
    assert_eq!(!PgaP4::point(), PgaP4::volume().swp());
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p5() {
    assert_eq!(!PgaP5::volume4(), PgaP5::point().swp());
    assert_eq!(!PgaP5::volume(), PgaP5::line().swp());
    assert_eq!(
        !PgaP5::plane(),
        (PgaP5::plane_displacement() - PgaP5::plane_moment()).swp()
    );
    assert_eq!(!!PgaP5::plane(), -PgaP5::plane());
    assert_eq!(!PgaP5::line(), PgaP5::volume().swp());
    assert_eq!(!PgaP5::point(), -PgaP5::volume4().swp());
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p6() {
    assert_eq!(!PgaP6::volume5().alt(), PgaP6::point().swp());
    assert_eq!(!PgaP6::volume4().alt(), PgaP6::line().swp());
    assert_eq!(!PgaP6::volume().alt(), PgaP6::plane().swp());
    assert_eq!(!PgaP6::plane(), PgaP6::volume().alt().swp());
    assert_eq!(!PgaP6::line(), PgaP6::volume4().alt().swp());
    assert_eq!(!PgaP6::point(), PgaP6::volume5().alt().swp());
}

#[cfg(feature = "exploratory")]
#[test]
fn pga_p7() {
    assert_eq!(!PgaP7::volume6(), PgaP7::point().swp());
    assert_eq!(!PgaP7::volume5(), PgaP7::line().swp());
    assert_eq!(!PgaP7::volume4(), PgaP7::plane().swp());
    assert_eq!(!PgaP7::volume(), PgaP7::volume().swp());
    assert_eq!(!PgaP7::plane(), -PgaP7::volume4().swp());
    assert_eq!(!PgaP7::line(), PgaP7::volume5().swp());
    assert_eq!(!PgaP7::point(), -PgaP7::volume6().swp());
}
