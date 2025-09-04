// Copyright © 2025 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Pistachio Flavor -- Projective Geometric Algebra (PGA)

#![allow(clippy::zero_prefixed_literal)]

use super::{Algebra, Choose, Multivector};
use core::{
    cmp::Ordering,
    fmt::{self, Debug, Display, Write},
    ops::{Mul, Not},
};
use smartstring::alias::String;

/// Basis blade of Elliptic 0D PGA.
pub type PgaE0 = Pga<1, 0>;
/// Basis blade of Elliptic 1D PGA.
pub type PgaE1 = Pga<1, 1>;
/// Basis blade of Elliptic 2D PGA.
pub type PgaE2 = Pga<1, 2>;
/// Basis blade of Elliptic 3D PGA.
pub type PgaE3 = Pga<1, 3>;
/// Basis blade of Elliptic 4D PGA (experimental).
pub type PgaE4 = Pga<1, 4>;
/// Basis blade of Elliptic 5D PGA (experimental, no inverse).
pub type PgaE5 = Pga<1, 5>;

/// Basis blade of Hyperbolic 0D PGA.
pub type PgaH0 = Pga<-1, 0>;
/// Basis blade of Hyperbolic 1D PGA.
pub type PgaH1 = Pga<-1, 1>;
/// Basis blade of Hyperbolic 2D PGA.
pub type PgaH2 = Pga<-1, 2>;
/// Basis blade of Hyperbolic 3D PGA.
pub type PgaH3 = Pga<-1, 3>;
/// Basis blade of Hyperbolic 4D PGA (experimental).
pub type PgaH4 = Pga<-1, 4>;
/// Basis blade of Hyperbolic 5D PGA (experimental, no inverse).
pub type PgaH5 = Pga<-1, 5>;

/// Basis blade of Parabolic (Euclidean) 0D PGA.
pub type PgaP0 = Pga<0, 0>;
/// Basis blade of Parabolic (Euclidean) 1D PGA.
pub type PgaP1 = Pga<0, 1>;
/// Basis blade of Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Pga<0, 2>;
/// Basis blade of Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Pga<0, 3>;
/// Basis blade of Parabolic (Euclidean) 4D PGA (experimental).
pub type PgaP4 = Pga<0, 4>;
/// Basis blade of Parabolic (Euclidean) 5D PGA (experimental, no inverse).
pub type PgaP5 = Pga<0, 5>;

/// Basis blade of PGA with metric $`M\in\{\pm 1,0\}`$ and embedded dimension $`N\in[0, 5]`$.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub struct Pga<const M: i8, const N: usize> {
    pub(crate) idx: u8,
}

/// Flavor-specific methods.
impl<const M: i8, const N: usize> Pga<M, N> {
    /// Creates basis blade from $`\e`$-notation.
    ///
    /// ```
    /// use vee::pga::PgaP3;
    ///
    /// let e = PgaP3::new("e");
    /// let e0 = PgaP3::new("e0");
    /// let e1 = PgaP3::new("e1");
    /// let e2 = PgaP3::new("e2");
    /// let e12 = PgaP3::new("e12");
    ///
    /// assert_eq!(e1 * e2, (1, e12));
    /// assert_eq!(e2 * e1, (-1, e12));
    /// assert_eq!(e0 * e0, (0, e));
    /// ```
    #[must_use]
    pub const fn new(sym: &'static str) -> Self {
        Self {
            #[allow(clippy::cast_possible_truncation)]
            idx: BasisBlade::new(N as u8, sym).idx,
        }
    }
    /// Constructs Cayley table.
    #[must_use]
    pub fn table() -> String {
        let mut tab = String::new();
        for a in Self::basis() {
            for b in Self::basis() {
                let (s, c) = a * b;
                let e = if s == 0 {
                    "0".into()
                } else {
                    format!("{}{}", if s > 0 { " " } else { "-" }, c)
                };
                write!(&mut tab, "{e:>0$}", N + 3).unwrap();
            }
            writeln!(&mut tab).unwrap();
        }
        tab
    }
}

impl<const M: i8, const N: usize> Algebra for Pga<M, N> {
    const N: usize = N;

    #[inline]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self> {
        TAB[N].iter().map(|b| Self { idx: b.idx })
    }
    #[inline]
    fn blade_len(&self) -> usize {
        (N + 1).choose(self.grade())
    }
    #[inline]
    fn grade(&self) -> usize {
        self.idx.count_ones() as usize
    }
}

impl<const M: i8, const N: usize> Mul for Pga<M, N> {
    type Output = (i8, Self);

    fn mul(self, other: Self) -> Self::Output {
        let [lhs, rhs] = [self, other].map(|b| b.idx);
        let mul = Self { idx: lhs ^ rhs };
        let cnt = ((1..=N).fold(0, |p, n| p ^ (lhs >> n)) & rhs).count_ones()
            + [self, other, mul]
                .map(|b| u32::from(LUT[N][b.idx as usize].cnt))
                .into_iter()
                .sum::<u32>();
        let sig = if cnt & 1 == 0 { 1 } else { -1 };
        let sig = if lhs & rhs & 1 == 0 { sig } else { sig * M };
        (sig, mul)
    }
}

impl<const M: i8, const N: usize> Not for Pga<M, N> {
    type Output = (i8, Self);

    #[inline]
    fn not(self) -> Self::Output {
        let not = Self {
            idx: !self.idx & ((1 << (N + 1)) - 1),
        };
        let (sig, _pss) = self * not;
        (sig, not)
    }
}

impl<const M: i8, const N: usize> Ord for Pga<M, N> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        let [lhs, rhs] = [self, other].map(|b| b.idx as usize);
        self.grade()
            .cmp(&other.grade())
            .then(LUT[N][lhs].idx.cmp(&LUT[N][rhs].idx))
    }
}

impl<const M: i8, const N: usize> PartialOrd for Pga<M, N> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const M: i8, const N: usize> Display for Pga<M, N> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(LUT[N][self.idx as usize].sym, f)
    }
}

#[derive(Debug, Copy, Clone)]
struct BasisBlade {
    sym: &'static str,
    cnt: u8,
    idx: u8,
}

impl BasisBlade {
    const fn new(dim: u8, mut sym: &'static str) -> Self {
        let s = sym.as_bytes();
        assert!(s[0] == b'e');
        let mut i = 1;
        let mut idx = 0;
        let mut cnt = 0;
        let mut lst = 0;
        while i < sym.len() {
            assert!(b'0' <= s[i] && s[i] <= b'9');
            idx |= 1 << (s[i] - b'0');
            if lst > b'0' && lst > s[i] {
                cnt += 1;
            }
            lst = s[i];
            i += 1;
        }
        match sym.len() - 1 {
            0 => sym = "1",
            l if l == dim as usize + 1 => sym = "I",
            _ => {}
        }
        Self { sym, cnt, idx }
    }
    const fn tab<const LEN: usize>(sym: [&'static str; LEN]) -> [Self; LEN] {
        let mut tab = [Self {
            sym: "",
            cnt: 0,
            idx: 0,
        }; LEN];
        let mut pss = 0;
        let mut i = 0;
        while i < LEN {
            #[allow(clippy::cast_possible_truncation)]
            let len = sym[i].len() as u8;
            if len > pss {
                pss = len;
            }
            i += 1;
        }
        let mut i = 0;
        while i < LEN {
            tab[i] = Self::new(pss - 2, sym[i]);
            i += 1;
        }
        tab
    }
    const fn lut<const LEN: usize>(tab: [Self; LEN]) -> [Self; LEN] {
        let mut lut = tab;
        let mut idx = 0;
        while idx < LEN {
            lut[tab[idx].idx as usize] = Self {
                sym: tab[idx].sym,
                cnt: tab[idx].cnt,
                #[allow(clippy::cast_possible_truncation)]
                idx: if idx >= LEN / 2 {
                    LEN - idx + LEN / 2
                } else {
                    idx
                } as u8,
            };
            idx += 1;
        }
        lut
    }
}

const TAB: [&[BasisBlade]; 6] = [&TAB0, &TAB1, &TAB2, &TAB3, &TAB4, &TAB5];
const LUT: [&[BasisBlade]; 6] = [&LUT0, &LUT1, &LUT2, &LUT3, &LUT4, &LUT5];

macro_rules! basis {
    ($t:ident, $u:ident, $n:tt, [$(($s:tt, $b:tt),)*]) => {
        /// The basis blades of the PGA with embedded dimension $`N = 0`$.
        ///
        /// ```gdef
        /// \gdef\idx#1{\expandafter\sub#1\relax}
        /// \gdef\sub#1#2\relax{#2}
        /// \gdef\fmt#1{\e_{\idx{#1}}}
        /// \gdef\e{\boldsymbol e}
        /// ```
        impl<const M: i8> Multivector<Pga<M, $n>> {
            $(
                #[doc = concat!(
                    "The multivector of basis blade $`",
                    stringify!($s),
                    "\\fmt{",
                    stringify!($b),
                    "}`$.",
                )]
                #[must_use]
                #[inline]
                pub fn $b() -> Self {
                    Self::new([(stringify!($s), const { Pga::new(stringify!($b)) })])
                }
            )*
        }
        const $t: [BasisBlade; count!($(($s, $b)),*)] = BasisBlade::tab([$(stringify!($b),)*]);
        const $u: [BasisBlade; count!($(($s, $b)),*)] = BasisBlade::lut($t);
    };
}

macro_rules! count {
    () => (0);
    ($head:expr $(, $tail:expr)*) => (1 + count!($($tail),*));
}

#[rustfmt::skip]
basis!(TAB0, LUT0, 0, [
    (v, e),
    (V, e0),
]);
#[rustfmt::skip]
basis!(TAB1, LUT1, 1, [
    (v, e),
    (W, e0),
    (w, e1),
    (V, e01),
]);
#[rustfmt::skip]
basis!(TAB2, LUT2, 2, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (Y, e01),
    (X, e20),
    (w, e12),
    (V, e012),
]);
#[rustfmt::skip]
basis!(TAB3, LUT3, 3, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (z, e3),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (z, e12),
    (y, e31),
    (x, e23),
    (Z, e021),
    (Y, e013),
    (X, e032),
    (w, e123),
    (V, e0123),
]);
#[rustfmt::skip]
basis!(TAB4, LUT4, 4, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (z, e3),
    (þ, e4),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Þ, e40),
    (a, e23),
    (b, e31),
    (c, e12),
    (d, e41),
    (e, e42),
    (f, e43),
    (F, e021),
    (E, e013),
    (D, e032),
    (C, e034),
    (B, e024),
    (A, e014),
    (þ, e123),
    (z, e124),
    (y, e314),
    (x, e234),
    (Þ, e0123),
    (Z, e0214),
    (Y, e0134),
    (X, e0324),
    (w, e1234),
    (V, e01234),
]);
#[rustfmt::skip]
basis!(TAB5, LUT5, 5, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (z, e3),
    (þ, e4),
    (ð, e5),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Þ, e40),
    (Ð, e05),
    (a, e23),
    (b, e31),
    (c, e12),
    (d, e41),
    (e, e42),
    (f, e43),
    (g, e15),
    (h, e25),
    (i, e35),
    (j, e45),
    (A, e015),
    (B, e052),
    (C, e035),
    (D, e054),
    (E, e014),
    (F, e042),
    (G, e034),
    (H, e032),
    (I, e013),
    (J, e021),
    (j, e345),
    (i, e245),
    (h, e145),
    (g, e152),
    (f, e315),
    (e, e253),
    (d, e123),
    (c, e124),
    (b, e134), 
    (a, e234),
    (J, e0123),
    (I, e0214),
    (H, e0134),
    (G, e0324),
    (F, e0215),
    (E, e0135),
    (D, e0325),
    (C, e0345),
    (B, e0245),
    (A, e0145),
    (ð, e1234),
    (þ, e1235),
    (z, e1245),
    (y, e3145),
    (x, e2345),
    (Ð, e01243),
    (Þ, e01235),
    (Z, e02145),
    (Y, e01345),
    (X, e03245),
    (w, e12345),
    (V, e012345),
]);

/// The named entities of the PGA with embedded dimension $`N = 0`$.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 0>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_0`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0()
    }
    /// The multivector of norm $`n \equiv n_0 + n_\infty`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 1`$.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 1>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_{01}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01()
    }
    /// The multivector of norm $`n \equiv n_0 + n_\infty`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of weight $`P_0 \equiv w\e_1`.
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
    /// The multivector of translator $`t \equiv v + V\e_{01}`$.
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::norm()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 2`$.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 2>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_{012}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012()
    }
    /// The multivector of norm $`n \equiv n_0 + n_\infty`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of moment $`\ell_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn moment() -> Self {
        Self::e0()
    }
    /// The multivector of displacement $`\ell_0 \equiv x\e_1 + y\e_2`$.
    #[must_use]
    #[inline]
    pub fn displacement() -> Self {
        Self::e1() + Self::e2()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::moment() + Self::displacement()
    }
    /// The multivector of direction $`P_\infty \equiv X\e_{20}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e20() + Self::e01()
    }
    /// The multivector of weight $`P_0 \equiv W\e_{12}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e12()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv n_0 + P_0`$.
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::weight()
    }
    /// The multivector of translator $`t \equiv n_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::direction()
    }
    /// The multivector of motor $`m \equiv n_0 + P`$.
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of rotoflector $`f_r \equiv \ell_0 + P_0`$.
    #[must_use]
    #[inline]
    pub fn rotoflector() -> Self {
        Self::displacement() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv \ell + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::line() + Self::direction()
    }
    /// The multivector of flector $`f \equiv \ell + P`$.
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::line() + Self::point()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 3`$.
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 3>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_{0123}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0123()
    }
    /// The multivector of norm $`n \equiv n_0 + n_\infty`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of bias $`p_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`p_0 \equiv x\e_1 + y\e_2 + z\e_3`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of displacement $`\ell \equiv x\e_{23} + y\e_{31} + z\e_{12}`$.
    #[must_use]
    #[inline]
    pub fn displacement() -> Self {
        Self::e23() + Self::e31() + Self::e12()
    }
    /// The multivector of moment $`\ell \equiv X\e_{01} + Y\e_{02} + Z\e_{03}`$.
    #[must_use]
    #[inline]
    pub fn moment() -> Self {
        Self::e01() + Self::e02() + Self::e03()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::displacement() + Self::moment()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{123}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e123()
    }
    /// The multivector of direction $`P_\infty \equiv X\e_{032} + Y\e_{013} + Z\e_{021}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e032() + Self::e013() + Self::e021()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv n_0 + \ell_0`$.
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::displacement()
    }
    /// The multivector of translator $`t \equiv n_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::moment()
    }
    /// The multivector of motor $`m \equiv n + \ell`$.
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::line()
    }
    /// The multivector of rotoflector $`f_r \equiv p_0 + P_0`$.
    #[must_use]
    #[inline]
    pub fn rotoflector() -> Self {
        Self::normal() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv p + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::plane() + Self::direction()
    }
    /// The multivector of flector $`f \equiv p + P`$.
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::plane() + Self::point()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 4`$ (experimental).
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 4>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_{01234}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01234()
    }
    /// The multivector of norm $`n \equiv n_0 + P`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of bias $`v_\infty \equiv w\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v_0 \equiv x\e_1 + y\e_2 + z\e_3 + þ\e_4`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of plane displacement
    /// $`p_0 \equiv a\e_{23} + b\e_{31} + c\e_{12} + d\e_{41} + e\e_{42} + f\e_{43}`$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e23() + Self::e31() + Self::e12() + Self::e41() + Self::e42() + Self::e43()
    }
    /// The multivector of plane moment
    /// $`p_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Þ\e_{40}`$.
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e40()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_displacement() + Self::plane_moment()
    }
    /// The multivector of line displacement
    /// $`\ell_\infty \equiv x\e_{234} + y\e_{314} + z\e_{124} + þ\e_{123}`$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e234() + Self::e314() + Self::e124() + Self::e123()
    }
    /// The multivector of line moment
    /// $`p_0 \equiv A\e_{014} + B\e_{024} + C\e_{034} + D\e_{035} + E\e_{013} + F\e_{021}`$.
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e014() + Self::e024() + Self::e034() + Self::e032() + Self::e013() + Self::e021()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_displacement() + Self::line_moment()
    }
    /// The multivector of bias $`P_0 \equiv W\e_{1234}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e1234()
    }
    /// The multivector of direction
    /// $`p_\infty \equiv X\e_{0324} + Y\e_{0134} + Z\e_{0214} + Þ\e_{0123}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e0324() + Self::e0134() + Self::e0214() + Self::e0123()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv n_0 + \ell_0`$.
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::plane_displacement()
    }
    /// The multivector of translator $`t \equiv n_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::plane_moment()
    }
    /// The multivector of translator $`m \equiv n + p`$.
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::plane()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 5`$ (experimental).
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 5>> {
    /// The multivector of scalar $`n_0 \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`n_\infty \equiv V\I`$ where $`\I \equiv \e_{012345}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012345()
    }
    /// The multivector of norm $`n \equiv n_0 + n_\infty`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of bias $`h_\infty \equiv w\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`h_0 \equiv x\e_1 + y\e_2 + z\e_3 + þ\e_4 + ð\e_5`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4() + Self::e5()
    }
    /// The multivector of hypervolume $`h \equiv h_0 + h_\infty`$.
    #[must_use]
    #[inline]
    pub fn hypervolume() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of volume displacement.
    ///
    /// ```math
    /// v_0 \equiv a\e_{23} + b\e_{31} + c\e_{12} + d\e_{41} + e\e_{42}
    ///          + f\e_{43} + g\e_{15} + h\e_{25} + i\e_{35} + j\e_{45}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e23()
            + Self::e31()
            + Self::e12()
            + Self::e41()
            + Self::e42()
            + Self::e43()
            + Self::e15()
            + Self::e25()
            + Self::e35()
            + Self::e45()
    }
    /// The multivector of volume moment.
    ///
    /// ```math
    /// v_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Þ\e_{40} + Ð\e_{05}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e40() + Self::e05()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::volume_displacement() + Self::volume_moment()
    }
    /// The multivector of plane displacement.
    ///
    /// ```math
    /// p_0 \equiv a\e_{234} + b\e_{134} + c\e_{124} + d\e_{123} + e\e_{253}
    ///          + f\e_{315} + g\e_{152} + h\e_{145} + i\e_{245} + j\e_{345}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e234()
            + Self::e134()
            + Self::e124()
            + Self::e123()
            + Self::e253()
            + Self::e315()
            + Self::e152()
            + Self::e145()
            + Self::e245()
            + Self::e345()
    }
    /// The multivector of plane moment.
    ///
    /// ```math
    /// p_\infty \equiv A\e_{015} + B\e_{052} + C\e_{035} + D\e_{054} + E\e_{014}
    ///               + F\e_{042} + G\e_{034} + H\e_{032} + I\e_{013} + J\e_{021}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e015()
            + Self::e052()
            + Self::e035()
            + Self::e054()
            + Self::e014()
            + Self::e042()
            + Self::e034()
            + Self::e032()
            + Self::e013()
            + Self::e021()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_displacement() + Self::plane_moment()
    }
    /// The multivector of line displacement.
    ///
    /// ```math
    /// \ell_0 \equiv x\e_{2345} + y\e_{3145} + z\e_{1245} + þ\e_{1235} + ð\e_{1234}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e2345() + Self::e3145() + Self::e1245() + Self::e1235() + Self::e1234()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{0145} + B\e_{0245} + C\e_{0345} + D\e_{0325} + E\e_{0135}
    ///                  + F\e_{0215} + G\e_{0324} + H\e_{0134} + I\e_{0214} + J\e_{0123}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e0145()
            + Self::e0245()
            + Self::e0345()
            + Self::e0325()
            + Self::e0135()
            + Self::e0215()
            + Self::e0324()
            + Self::e0134()
            + Self::e0214()
            + Self::e0123()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_displacement() + Self::line_moment()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{12345}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e12345()
    }
    /// The multivector of direction.
    ///
    /// ```math
    /// P_\infty \equiv X\e_{03245} + Y\e_{01345} + Z\e_{02145} + Þ\e_{01235} + Ð\e_{01243}
    /// ```
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e03245() + Self::e01345() + Self::e02145() + Self::e01235() + Self::e01243()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv n_0 + p_0`$.
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::plane_displacement()
    }
    /// The multivector of translator $`t \equiv n_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::plane_moment()
    }
    /// The multivector of motor $`m \equiv n + p`$.
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::plane()
    }
}

// #[test]
// fn dim() {
//     assert!(TAB.windows(2).all(|tab| {
//         let a = tab[0];
//         let b = tab[1];
//         a.iter().all(|a| b.iter().any(|b| b.sym == a.sym))
//     }));
// }

#[test]
fn not() {
    use super::{PgaP0, PgaP1, PgaP2, PgaP3, PgaP4, PgaP5};

    assert_eq!(!PgaP0::norm(), PgaP0::norm().swp());

    assert_eq!(
        !PgaP1::point(),
        (PgaP1::weight() - PgaP1::direction()).swp()
    );

    assert_eq!(!PgaP2::line(), PgaP2::point().swp());
    assert_eq!(!PgaP2::point(), PgaP2::line().swp());

    assert_eq!(!PgaP3::plane(), PgaP3::point().swp());
    assert_eq!(!PgaP3::line(), PgaP3::line().swp());
    assert_eq!(!PgaP3::point(), -PgaP3::plane().swp());

    assert_eq!(!PgaP4::volume(), PgaP4::point().swp());
    assert_eq!(!PgaP4::plane(), PgaP4::line().swp());
    assert_eq!(!PgaP4::line(), PgaP4::plane().swp());
    assert_eq!(!PgaP4::point(), PgaP4::volume().swp());

    assert_eq!(!PgaP5::hypervolume(), PgaP5::point().swp());
    assert_eq!(!PgaP5::volume(), PgaP5::line().swp());
    assert_eq!(
        !PgaP5::plane(),
        (PgaP5::plane_moment() - PgaP5::plane_displacement()).swp()
    );
    assert_eq!(!!PgaP5::plane(), -PgaP5::plane());
    assert_eq!(!PgaP5::line(), PgaP5::volume().swp());
    assert_eq!(!PgaP5::point(), -PgaP5::hypervolume().swp());
}

#[test]
fn mul() {
    use std::{
        fs::{read_to_string, write},
        path::Path,
    };
    let tables = [
        ("PgaE0", PgaE0::table()),
        ("PgaE1", PgaE1::table()),
        ("PgaE2", PgaE2::table()),
        ("PgaE3", PgaE3::table()),
        ("PgaE4", PgaE4::table()),
        ("PgaE5", PgaE5::table()),
        ("PgaH0", PgaH0::table()),
        ("PgaH1", PgaH1::table()),
        ("PgaH2", PgaH2::table()),
        ("PgaH3", PgaH3::table()),
        ("PgaH4", PgaH4::table()),
        ("PgaH5", PgaH5::table()),
        ("PgaP0", PgaP0::table()),
        ("PgaP1", PgaP1::table()),
        ("PgaP2", PgaP2::table()),
        ("PgaP3", PgaP3::table()),
        ("PgaP4", PgaP4::table()),
        ("PgaP5", PgaP5::table()),
    ];
    for (pga, table) in tables {
        let path = Path::new("tests").join(pga).with_extension("ct");
        if let Ok(text) = read_to_string(&path) {
            assert_eq!(table, text);
        } else {
            write(&path, table).unwrap();
        }
    }
}
