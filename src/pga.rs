// Copyright © 2025 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Pistachio Flavor -- Projective Geometric Algebra (PGA)

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
/// Basis blade of Elliptic 6D PGA (experimental, no inverse).
pub type PgaE6 = Pga<1, 6>;
/// Basis blade of Elliptic 7D PGA (experimental, no inverse).
pub type PgaE7 = Pga<1, 7>;

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
/// Basis blade of Hyperbolic 6D PGA (experimental, no inverse).
pub type PgaH6 = Pga<-1, 6>;
/// Basis blade of Hyperbolic 7D PGA (experimental, no inverse).
pub type PgaH7 = Pga<-1, 7>;

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
/// Basis blade of Parabolic (Euclidean) 6D PGA (experimental, no inverse).
pub type PgaP6 = Pga<0, 6>;
/// Basis blade of Parabolic (Euclidean) 7D PGA (experimental, no inverse).
pub type PgaP7 = Pga<0, 7>;

/// Basis blade of PGA with metric $`M\in\{\pm 1,0\}`$ and embedded dimension $`N\in[0, 7]`$.
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
pub struct Pga<const M: i8, const N: u32> {
    pub(crate) idx: u8,
}

/// Flavor-specific methods.
impl<const M: i8, const N: u32> Pga<M, N> {
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
                write!(&mut tab, "{e:>0$}", N as usize + 3).unwrap();
            }
            writeln!(&mut tab).unwrap();
        }
        tab
    }
}

impl<const M: i8, const N: u32> Algebra for Pga<M, N> {
    const N: u32 = N;

    #[inline]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self> {
        TAB[N as usize].iter().map(|b| Self { idx: b.idx })
    }
    #[inline]
    fn blade_len(&self) -> usize {
        (N + 1).choose(self.grade()) as usize
    }
    #[inline]
    fn grade(&self) -> u32 {
        self.idx.count_ones()
    }
}

impl<const M: i8, const N: u32> Mul for Pga<M, N> {
    type Output = (i8, Self);

    fn mul(self, other: Self) -> Self::Output {
        let [lhs, rhs] = [self, other].map(|b| b.idx);
        let mul = Self { idx: lhs ^ rhs };
        let cnt = ((1..=N).fold(0, |p, n| p ^ (lhs >> n)) & rhs).count_ones()
            + [self, other, mul]
                .map(|b| u32::from(LUT[N as usize][b.idx as usize].cnt))
                .into_iter()
                .sum::<u32>();
        let sig = if cnt & 1 == 0 { 1 } else { -1 };
        let sig = if lhs & rhs & 1 == 0 { sig } else { sig * M };
        (sig, mul)
    }
}

impl<const M: i8, const N: u32> Not for Pga<M, N> {
    type Output = (i8, Self);

    #[inline]
    fn not(self) -> Self::Output {
        let pss = const { u8::MAX >> (u8::BITS - (N + 1)) };
        let not = Self {
            idx: !self.idx & pss,
        };
        let (sig, _pss) = self * not;
        (sig, not)
    }
}

impl<const M: i8, const N: u32> Ord for Pga<M, N> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        let [lhs, rhs] = [self, other].map(|b| b.idx as usize);
        self.grade()
            .cmp(&other.grade())
            .then(LUT[N as usize][lhs].idx.cmp(&LUT[N as usize][rhs].idx))
    }
}

impl<const M: i8, const N: u32> PartialOrd for Pga<M, N> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const M: i8, const N: u32> Display for Pga<M, N> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(LUT[N as usize][self.idx as usize].sym, f)
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

const TAB: [&[BasisBlade]; 8] = [&TAB0, &TAB1, &TAB2, &TAB3, &TAB4, &TAB5, &TAB6, &TAB7];
const LUT: [&[BasisBlade]; 8] = [&LUT0, &LUT1, &LUT2, &LUT3, &LUT4, &LUT5, &LUT6, &LUT7];

macro_rules! basis {
    ($t:ident, $u:ident, $n:tt, $l:tt, [$(($s:tt, $b:tt),)*]) => {
        #[doc = concat!("The basis blades of the PGA with embedded dimension $`N = ", $n, "`$.")]
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
        const $t: [BasisBlade; $l] = BasisBlade::tab([$(stringify!($b),)*]);
        const $u: [BasisBlade; $l] = BasisBlade::lut($t);
    };
}

#[rustfmt::skip]
basis!(TAB0, LUT0, 0, 2, [
    (v, e),
    (V, e0),
]);
#[rustfmt::skip]
basis!(TAB1, LUT1, 1, 4, [
    (v, e),
    (W, e0),
    (w, e1),
    (V, e01),
]);
#[rustfmt::skip]
basis!(TAB2, LUT2, 2, 8, [
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
basis!(TAB3, LUT3, 3, 16, [
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
basis!(TAB4, LUT4, 4, 32, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (z, e3),
    (ð, e4),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e40),
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
    (ð, e123),
    (z, e124),
    (y, e314),
    (x, e234),
    (Ð, e0123),
    (Z, e0214),
    (Y, e0134),
    (X, e0324),
    (w, e1234),
    (V, e01234),
]);
#[rustfmt::skip]
basis!(TAB5, LUT5, 5, 64, [
    (v, e),
    (W, e0),
    (x, e1),
    (y, e2),
    (z, e3),
    (ð, e4),
    (ø, e5),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e40),
    (Ø, e05),
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
    (ø, e1234),
    (ð, e1235),
    (z, e1245),
    (y, e3145),
    (x, e2345),
    (Ø, e01243),
    (Ð, e01235),
    (Z, e02145),
    (Y, e01345),
    (X, e03245),
    (w, e12345),
    (V, e012345),
]);
#[rustfmt::skip]
basis!(TAB6, LUT6, 6, 128, [
    (v, e),
    (W, e0),
    // 6
    (x, e1),
    (y, e2),
    (z, e3),
    (ð, e4),
    (ø, e5),
    (þ, e6),
    // 6
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e40),
    (Ø, e05),
    (Þ, e60),
    // 15
    (α, e23),
    (β, e31),
    (γ, e12),
    (δ, e41),
    (ε, e42),
    (ζ, e43),
    (η, e15),
    (θ, e25),
    (ι, e35),
    (κ, e45),
    (λ, e16),
    (μ, e62),
    (ν, e36),
    (ξ, e64),
    (ο, e56),
    // 15
    (Α, e015),
    (Β, e052),
    (Γ, e035),
    (Δ, e054),
    (Ε, e014),
    (Ζ, e042),
    (Η, e034),
    (Θ, e032),
    (Ι, e013),
    (Κ, e021),
    (Λ, e016),
    (Μ, e062),
    (Ν, e036),
    (Ξ, e064),
    (Ο, e056),
    // 20
    (a, e345),
    (b, e245),
    (c, e145),
    (d, e152),
    (e, e315),
    (f, e253),
    (g, e123),
    (h, e124),
    (i, e134),
    (j, e234),
    (k, e126),
    (l, e163),
    (m, e146),
    (n, e165),
    (o, e236),
    (p, e264),
    (q, e256),
    (r, e346),
    (s, e365),
    (t, e456),
    // 20
    (T, e0123),
    (S, e0124),
    (R, e0125),
    (Q, e0134),
    (P, e0135),
    (O, e0145),
    (N, e0234),
    (M, e0235),
    (L, e0245),
    (K, e0345),
    (J, e0156),
    (I, e0265),
    (H, e0356),
    (G, e0465),
    (F, e0146),
    (E, e0264),
    (D, e0346),
    (C, e0263),
    (B, e0136),
    (A, e0162),
    // 15
    (ο, e1234),
    (ξ, e1235),
    (ν, e1245),
    (μ, e1345),
    (λ, e2345),
    (κ, e3465),
    (ι, e2465),
    (θ, e1465),
    (η, e1256),
    (ζ, e1356),
    (ε, e2356),
    (δ, e1263),
    (γ, e1264),
    (β, e1364),
    (α, e2364),
    // 15
    (Ο, e01234),
    (Ξ, e01235),
    (Ν, e01245),
    (Μ, e01345),
    (Λ, e02345),
    (Κ, e01236),
    (Ι, e01264),
    (Θ, e01346),
    (Η, e02364),
    (Ζ, e01265),
    (Ε, e01356),
    (Δ, e02365),
    (Γ, e03456),
    (Β, e02456),
    (Α, e01456),
    // 6
    (þ, e12345),
    (ø, e12346),
    (ð, e12356),
    (z, e12456),
    (y, e13465),
    (x, e23456),
    // 6
    (Þ, e012345),
    (Ø, e012436),
    (Ð, e012356),
    (Z, e021456),
    (Y, e013456),
    (X, e032456),

    (w, e123456),
    (V, e0123456),
]);
#[rustfmt::skip]
basis!(TAB7, LUT7, 7, 256, [
    // 1
    (v, e),
    // 1
    (W, e0),
    // 7
    (x, e1),
    (y, e2),
    (z, e3),
    (ð, e4),
    (ø, e5),
    (þ, e6),
    (œ, e7),
    // 7
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e04),
    (Ø, e05),
    (Þ, e06),
    (Œ, e07),
    // 21
    (α, e23),
    (β, e31),
    (γ, e12),
    (δ, e41),
    (ε, e42),
    (ζ, e43),
    (η, e15),
    (θ, e25),
    (ι, e35),
    (κ, e45),
    (λ, e16),
    (μ, e62),
    (ν, e36),
    (ξ, e64),
    (ο, e56),
    (π, e17),
    (ρ, e27),
    (σ, e37),
    (τ, e47),
    (υ, e57),
    (φ, e67),
    // 21
    (Α, e012),
    (Β, e013),
    (Γ, e014),
    (Δ, e015),
    (Ε, e016),
    (Ζ, e017),
    (Η, e023),
    (Θ, e024),
    (Ι, e025),
    (Κ, e026),
    (Λ, e027),
    (Μ, e034),
    (Ν, e035),
    (Ξ, e036),
    (Ο, e037),
    (Π, e045),
    (Ρ, e046),
    (Σ, e047),
    (Τ, e056),
    (Υ, e057),
    (Φ, e067),
    // 35
    (a, e123),
    (b, e124),
    (c, e125),
    (d, e126),
    (e, e127),
    (f, e134),
    (g, e135),
    (h, e136),
    (i, e137),
    (j, e145),
    (k, e146),
    (l, e147),
    (m, e156),
    (n, e157),
    (o, e167),
    (p, e234),
    (q, e235),
    (r, e236),
    (s, e237),
    (t, e245),
    (u, e246),
    (á, e247),
    (ä, e256),
    (å, e257),
    (æ, e267),
    (ç, e345),
    (é, e346),
    (ë, e347),
    (í, e356),
    (ï, e357),
    (ñ, e367),
    (ó, e456),
    (ö, e457),
    (ú, e467),
    (ü, e567),
    // 35
    (A, e0123),
    (B, e0124),
    (C, e0125),
    (D, e0126),
    (E, e0127),
    (F, e0134),
    (G, e0135),
    (H, e0136),
    (I, e0137),
    (J, e0145),
    (K, e0146),
    (L, e0147),
    (M, e0156),
    (N, e0157),
    (O, e0167),
    (P, e0234),
    (Q, e0235),
    (R, e0236),
    (S, e0237),
    (T, e0245),
    (U, e0246),
    (Á, e0247),
    (Ä, e0256),
    (Å, e0257),
    (Æ, e0267),
    (Ç, e0345),
    (É, e0346),
    (Ë, e0347),
    (Í, e0356),
    (Ï, e0357),
    (Ñ, e0367),
    (Ó, e0456),
    (Ö, e0457),
    (Ú, e0467),
    (Ü, e0567),
    // 35
    (ü, e1234),
    (ú, e1253),
    (ö, e1236),
    (ó, e1273),
    (ñ, e1245),
    (ï, e1264),
    (í, e1247),
    (ë, e1256),
    (é, e1275),
    (ç, e1267),
    (æ, e1354),
    (å, e1346),
    (ä, e1374),
    (á, e1365),
    (u, e1357),
    (t, e1376),
    (s, e1456),
    (r, e1475),
    (q, e1467),
    (p, e1576),
    (o, e2345),
    (n, e2364),
    (m, e2347),
    (l, e2356),
    (k, e2375),
    (j, e2367),
    (i, e2465),
    (h, e2457),
    (g, e2476),
    (f, e2567),
    (e, e3456),
    (d, e3475),
    (c, e3467),
    (b, e3576),
    (a, e4567),
    // 35
    (Ü, e01243),
    (Ú, e01235),
    (Ö, e01263),
    (Ó, e01237),
    (Ñ, e01254),
    (Ï, e01246),
    (Í, e01274),
    (Ë, e01265),
    (É, e01257),
    (Ç, e01276),
    (Æ, e01345),
    (Å, e01364),
    (Ä, e01347),
    (Á, e01356),
    (U, e01375),
    (T, e01367),
    (S, e01465),
    (R, e01457),
    (Q, e01476),
    (P, e01567),
    (O, e02354),
    (N, e02346),
    (M, e02374),
    (L, e02365),
    (K, e02357),
    (J, e02376),
    (I, e02456),
    (H, e02475),
    (G, e02467),
    (F, e02576),
    (E, e03465),
    (D, e03457),
    (C, e03476),
    (B, e03567),
    (A, e04576),
    // 21
    (φ, e12345),
    (υ, e12364),
    (τ, e12347),
    (σ, e12356),
    (ρ, e12375),
    (π, e12367),
    (ο, e12465),
    (ξ, e12457),
    (ν, e12476),
    (μ, e12567),
    (λ, e13456),
    (κ, e13475),
    (ι, e13467),
    (θ, e13576),
    (η, e14567),
    (ζ, e23465),
    (ε, e23457),
    (δ, e23476),
    (γ, e23567),
    (β, e24576),
    (α, e34567),
    // 21
    (Φ, e012345),
    (Υ, e012364),
    (Τ, e012356),
    (Σ, e012465),
    (Ρ, e013456),
    (Π, e023465),
    (Ο, e012347),
    (Ξ, e012357),
    (Ν, e012457),
    (Μ, e013457),
    (Λ, e023457),
    (Κ, e012367),
    (Ι, e012476),
    (Θ, e013467),
    (Η, e023647),
    (Ζ, e012576),
    (Ε, e013567),
    (Δ, e023657),
    (Γ, e034567),
    (Β, e024567),
    (Α, e014567),
    // 7
    (œ, e123456),
    (þ, e123475),
    (ø, e123467),
    (ð, e123576),
    (z, e124567),
    (y, e134657),
    (x, e234567),
    // 7
    (Œ, e0124356),
    (Þ, e0123457),
    (Ø, e0124367),
    (Ð, e0123567),
    (Z, e0214567),
    (Y, e0134567),
    (X, e0324567),
    // 1
    (w, e1234567),
    // 1
    (V, e01234567),
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
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{01}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::pseudoscalar()
    }
    /// The multivector of weight $`P_0 \equiv w\e_1`$.
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
    /// The multivector of translator $`t \equiv s + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP1 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "+1w͔w͕",
    ///     "+(+1W͔w͕-1W͕w͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::translator();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[+1vv]W͓+[+2Vv]w͓)e0",
    ///     "+[+1vv]w͓e1",
    /// ]);
    /// ```
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
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{012}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012()
    }
    /// The multivector of norm $`n \equiv s + S`$.
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
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+1X͔Y͕-1X͕Y͔)e0",
    ///     "+(+1Y͔w͕-1Y͕w͔)e1",
    ///     "+(-1X͔w͕+1X͕w͔)e2",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::moment() + Self::displacement()
    }
    /// The multivector of direction $`P_\infty \equiv X\e_{20} + Y\e_{01}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e20() + Self::e01()
    }
    /// The multivector of weight $`P_0 \equiv w\e_{12}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e12()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let point = Vee::line().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1W͔y͕+1W͕y͔)e20",
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "+1x͔x͕+1y͔y͕",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-1w͔w͕",
    ///     "+(-1Y͔w͕+1Y͕w͔)e20",
    ///     "+(+1X͔w͕-1X͕w͔)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::direction()
    }
    /// The multivector of motor $`m \equiv s + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "+1x͔x͕+1y͔y͕",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1W͔y͕+1W͕y͔)e20",
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::motor();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[+1vv+1ww]w͓)e12",
    ///     "+(+[+2vw]Y͓+[+2Xw-2Yv]w͓+[+1vv-1ww]X͓)e20",
    ///     "+(+[+2Xv+2Yw]w͓+[-2vw]X͓+[+1vv-1ww]Y͓)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of flector $`f \equiv \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let flector = Vee::line().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+1W͔v͕+1X͕y͔-1Y͕x͔)e0",
    ///     "+(+1v͕x͔-1w͕y͔)e1",
    ///     "+(+1v͕y͔+1w͕x͔)e2",
    ///     "+(+1W͔w͕+1X͕x͔+1Y͕y͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::flector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[-1xx-1yy]w͓)e12",
    ///     "+(+[+2xy]Y͓+[+2Vy+2Wx]w͓+[+1xx-1yy]X͓)e20",
    ///     "+(+[-2Vx+2Wy]w͓+[+2xy]X͓+[-1xx+1yy]Y͓)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::line() + Self::pseudoscalar()
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
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{0123}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0123()
    }
    /// The multivector of norm $`n \equiv s + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let squared_norm = Vee::line().squared_norm();
    ///
    /// assert_eq!(squared_norm.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+1xx+1yy+1zz",
    ///     "+(-2Xx-2Yy-2Zz)I",
    /// ]);
    /// ```
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
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let plane = Vee::line().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    /// format_eq!(plane, [
    ///     "+(-1X͔X͕-1Y͔Y͕-1Z͔Z͕)e0",
    ///     "+(+1X͔w͕-1Y͕z͔+1Z͕y͔)e1",
    ///     "+(+1X͕z͔+1Y͔w͕-1Z͕x͔)e2",
    ///     "+(-1X͕y͔+1Y͕x͔+1Z͔w͕)e3",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of displacement $`\ell_0 \equiv x\e_{23} + y\e_{31} + z\e_{12}`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// // A line through the origin as the join of a point and the origin.
    /// let displacement = Vee::point().lhs() & Vee::weight().rhs();
    ///
    /// assert_eq!(displacement.basis_blades(), Vee::displacement().basis_blades());
    /// format_eq!(displacement, [
    ///     "-1X͔w͕e23",
    ///     "-1Y͔w͕e31",
    ///     "-1Z͔w͕e12",
    /// ]);
    ///
    /// // A line through the origin as the meet of two planes through the origin.
    /// let displacement = Vee::normal().lhs() ^ Vee::normal().rhs();
    ///
    /// assert_eq!(displacement.basis_blades(), Vee::displacement().basis_blades());
    /// format_eq!(displacement, [
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn displacement() -> Self {
        Self::e23() + Self::e31() + Self::e12()
    }
    /// The multivector of moment $`\ell_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03}`$.
    #[must_use]
    #[inline]
    pub fn moment() -> Self {
        Self::e01() + Self::e02() + Self::e03()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// // A line as the join of two points.
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+1Y͔Z͕-1Y͕Z͔)e01",
    ///     "+(-1X͔Z͕+1X͕Z͔)e02",
    ///     "+(+1X͔Y͕-1X͕Y͔)e03",
    ///     "+(-1X͔w͕+1X͕w͔)e23",
    ///     "+(-1Y͔w͕+1Y͕w͔)e31",
    ///     "+(-1Z͔w͕+1Z͕w͔)e12",
    /// ]);
    ///
    /// // A line as the meet of two planes.
    /// let line = Vee::plane().lhs() ^ Vee::plane().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    ///     "+(+1W͔y͕-1W͕y͔)e02",
    ///     "+(+1W͔z͕-1W͕z͔)e03",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    /// ]);
    /// ```
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
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// // A point as the meet of a plane and a line.
    /// let point = Vee::plane().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕)e123",
    ///     "+(-1W͔x͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(-1W͔y͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1W͔z͕-1X͕y͔+1Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    /// ]);
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "-1x͔x͕-1y͔y͕-1z͔z͕",
    ///     "+(-1y͔z͕+1y͕z͔)e23",
    ///     "+(+1x͔z͕-1x͕z͔)e31",
    ///     "+(-1x͔y͕+1x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::displacement()
    }
    /// The multivector of translator $`t \equiv s + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-1w͔w͕",
    ///     "+(+1X͔w͕-1X͕w͔)e01",
    ///     "+(+1Y͔w͕-1Y͕w͔)e02",
    ///     "+(+1Z͔w͕-1Z͕w͔)e03",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + \ell`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let simple_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(simple_motor.basis_blades(), Vee::simple_motor().basis_blades());
    /// format_eq!(format!("{simple_motor:#}"), [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕",
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    ///     "+(+1W͔y͕-1W͕y͔)e02",
    ///     "+(+1W͔z͕-1W͕z͔)e03",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::line()
    }
    /// The multivector of motor $`m \equiv s + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "-1x͔x͕-1y͔y͕-1z͔z͕",
    ///     "+(-1Y͔z͕+1Y͕z͔+1Z͔y͕-1Z͕y͔)e01",
    ///     "+(+1X͔z͕-1X͕z͔-1Z͔x͕+1Z͕x͔)e02",
    ///     "+(-1X͔y͕+1X͕y͔+1Y͔x͕-1Y͕x͔)e03",
    ///     "+(-1y͔z͕+1y͕z͔)e23",
    ///     "+(+1x͔z͕-1x͕z͔)e31",
    ///     "+(-1x͔y͕+1x͕y͔)e12",
    ///     "+(+1X͔x͕+1X͕x͔+1Y͔y͕+1Y͕y͔+1Z͔z͕+1Z͕z͔)I",
    /// ]);
    ///
    /// let motor = Vee::rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "+1v͔v͕",
    ///     "+(+1X͕v͔+1Y͕z͔-1Z͕y͔)e01",
    ///     "+(-1X͕z͔+1Y͕v͔+1Z͕x͔)e02",
    ///     "+(+1X͕y͔-1Y͕x͔+1Z͕v͔)e03",
    ///     "+1v͕x͔e23",
    ///     "+1v͕y͔e31",
    ///     "+1v͕z͔e12",
    ///     "+(+1X͕x͔+1Y͕y͔+1Z͕z͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::line()
    }
    /// The multivector of rotoreflector $`f_r \equiv p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// format_eq!(rotoreflector, [
    ///     "+(+1v͕x͔-1y͔z͕+1y͕z͔)e1",
    ///     "+(+1v͕y͔+1x͔z͕-1x͕z͔)e2",
    ///     "+(+1v͕z͔-1x͔y͕+1x͕y͔)e3",
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(-1X͕x͔-1Y͕y͔-1Z͕z͔)e0",
    ///     "+1v͕x͔e1",
    ///     "+1v͕y͔e2",
    ///     "+1v͕z͔e3",
    ///     "+(-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1X͕y͔+1Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::plane() + Self::direction()
    }
    /// The multivector of flector $`f \equiv p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP3 as Vee};
    ///
    /// let flector = Vee::plane().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(format!("{flector:#}"), [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔)e0",
    ///     "+(+1v͕x͔-1y͔z͕+1y͕z͔)e1",
    ///     "+(+1v͕y͔+1x͔z͕-1x͕z͔)e2",
    ///     "+(+1v͕z͔-1x͔y͕+1x͕y͔)e3",
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕)e123",
    ///     "+(+1V͕x͔-1W͔x͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(+1V͕y͔-1W͔y͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(+1V͕z͔-1W͔z͕-1X͕y͔+1Y͕x͔)e021",
    /// ]);
    /// ```
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
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{01234}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01234()
    }
    /// The multivector of norm $`n \equiv s + P`$.
    ///
    /// Quadvector $`P`$ does square to a scalar, therefore $`n`$ is a Study number.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let quadvector_squared_norm = Vee::point().squared_norm();
    ///
    /// format_eq!(quadvector_squared_norm, "ww");
    /// ```
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
    /// The multivector of normal $`v_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4`$.
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
    /// $`p_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{40}`$.
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
    /// $`\ell_0 \equiv x\e_{234} + y\e_{314} + z\e_{124} + ð\e_{123}`$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e234() + Self::e314() + Self::e124() + Self::e123()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{014} + B\e_{024} + C\e_{034} + D\e_{032} + E\e_{013} + F\e_{021}
    /// ```
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
    /// The multivector of weight $`P_0 \equiv W\e_{1234}`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e1234()
    }
    /// The multivector of direction
    /// $`P_\infty \equiv X\e_{0324} + Y\e_{0134} + Z\e_{0214} + Ð\e_{0123}`$.
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
    /// The multivector of single rotator $`r_1 \equiv s + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "-1x͔x͕-1y͔y͕-1z͔z͕-1ð͔ð͕",
    ///     "+(-1y͔z͕+1y͕z͔)e23",
    ///     "+(+1x͔z͕-1x͕z͔)e31",
    ///     "+(-1x͔y͕+1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    /// ]);
    ///
    /// let squared_norm = Vee::single_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv",
    ///     "+(+2ad+2be+2cf)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::plane_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + \ell_0 + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕+1v͔v͕",
    ///     "+(+1a͔v͕+1a͕v͔-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1b͔v͕+1b͕v͔+1d͔f͕-1d͕f͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔+1c͔v͕+1c͕v͔-1d͔e͕+1d͕e͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔+1d͔v͕+1d͕v͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔+1e͔v͕+1e͕v͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔+1f͔v͕+1f͕v͔)e43",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕",
    ///     "+(-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1d͔f͕-1d͕f͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔-1d͔e͕+1d͕e͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔)e43",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    /// ]);
    ///
    /// let squared_norm = Vee::double_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv+1ww",
    ///     "+(+2ad+2be+2cf+2vw)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "+1w͔w͕",
    ///     "+(-1X͔w͕+1X͕w͔)e01",
    ///     "+(-1Y͔w͕+1Y͕w͔)e02",
    ///     "+(-1Z͔w͕+1Z͕w͔)e03",
    ///     "+(-1w͔Ð͕+1w͕Ð͔)e40",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::plane_moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let simple_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(simple_motor.basis_blades(), Vee::simple_motor().basis_blades());
    /// format_eq!(simple_motor, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕",
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    ///     "+(+1W͔y͕-1W͕y͔)e02",
    ///     "+(+1W͔z͕-1W͕z͔)e03",
    ///     "+(-1W͔ð͕+1W͕ð͔)e40",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    /// ]);
    ///
    /// let squared_norm = Vee::simple_motor().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(squared_norm, [
    ///     // Scalar condition.
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv",
    ///     // Point condition.
    ///     "+(+2ad+2be+2cf)e1234", // Weight condition.
    ///     "+(-2Yf+2Ze-2aÐ)e0324", // Direction condition.
    ///     "+(+2Xf-2Zd-2bÐ)e0134", // Direction condition.
    ///     "+(-2Xe+2Yd-2cÐ)e0214", // Direction condition.
    ///     "+(-2Xa-2Yb-2Zc)e0123", // Direction condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+[+2ad+2be+2cf]w͓", // Vanishes with weight condition.
    ///     "+(+[+1aa+1bb+1cc+1dd+1ee+1ff+1vv]w͓)e1234",
    ///     "+(+[+2ac-2bv-2df]Z͓+[-2Xv-2Yc+2Zb-2dÐ]w͓+[+2bf-2ce-2dv]Ð͓\
    ///        +[+1aa-1bb-1cc-1dd+1ee+1ff+1vv]X͓+[+2ab+2cv-2de]Y͓)e0324",
    ///     "+(+[+2Xc-2Yv-2Za-2eÐ]w͓+[-2af+2cd-2ev]Ð͓+[+2ab-2cv-2de]X͓\
    ///        +[-1aa+1bb-1cc+1dd-1ee+1ff+1vv]Y͓+[+2av+2bc-2ef]Z͓)e0134",
    ///     "+(+[+2ae-2bd-2fv]Ð͓+[+2ac+2bv-2df]X͓+[-2av+2bc-2ef]Y͓\
    ///        +[-1aa-1bb+1cc+1dd+1ee-1ff+1vv]Z͓+[-2Xb+2Ya-2Zv-2fÐ]w͓)e0214",
    ///     "+(+[+2bf-2ce+2dv]X͓+[-2af+2cd+2ev]Y͓+[+2ae-2bd+2fv]Z͓\
    ///        +[-2Xd-2Ye-2Zf+2vÐ]w͓+[+1aa+1bb+1cc-1dd-1ee-1ff+1vv]Ð͓)e0123",
    /// ]);
    ///
    /// let line = Vee::line().pin() << Vee::simple_motor();
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// let plane = Vee::plane().pin() << Vee::simple_motor();
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    ///
    /// let volume = Vee::volume().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(volume.basis_blades(), (Vee::pseudoscalar() + Vee::volume()).basis_blades());
    /// format_eq!(volume, [
    ///     "+(+[+1aa+1bb+1cc+1dd+1ee+1ff+1vv]W͓\
    ///        +[+2Xv-2Yc+2Zb-2dÐ]x͓+[+2Xc+2Yv-2Za-2eÐ]y͓+[-2Xb+2Ya+2Zv-2fÐ]z͓+[-2Xd-2Ye-2Zf-2vÐ]ð͓)e0",
    ///     "+(+[+2ab+2cv-2de]y͓+[+2ac-2bv-2df]z͓+[+2bf-2ce-2dv]ð͓+[+1aa-1bb-1cc-1dd+1ee+1ff+1vv]x͓)e1",
    ///     "+(+[+2av+2bc-2ef]z͓+[-2af+2cd-2ev]ð͓+[+2ab-2cv-2de]x͓+[-1aa+1bb-1cc+1dd-1ee+1ff+1vv]y͓)e2",
    ///     "+(+[+2ae-2bd-2fv]ð͓+[+2ac+2bv-2df]x͓+[-2av+2bc-2ef]y͓+[-1aa-1bb+1cc+1dd+1ee-1ff+1vv]z͓)e3",
    ///     "+(+[+2bf-2ce+2dv]x͓+[-2af+2cd+2ev]y͓+[+2ae-2bd+2fv]z͓+[+1aa+1bb+1cc-1dd-1ee-1ff+1vv]ð͓)e4",
    ///     // Vanishes with point condition.
    ///     "+(+[+2ad+2be+2cf]W͓+[-2Yf+2Ze-2aÐ]x͓+[+2Xf-2Zd-2bÐ]y͓+[-2Xe+2Yd-2cÐ]z͓+[-2Xa-2Yb-2Zc]ð͓)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::plane()
    }
    /// The multivector of single motor $`m_1 \equiv s + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+1v͔v͕",
    ///     "+(+1X͕v͔+1Y͕c͔-1Z͕b͔+1d͔Ð͕)e01",
    ///     "+(-1X͕c͔+1Y͕v͔+1Z͕a͔+1e͔Ð͕)e02",
    ///     "+(+1X͕b͔-1Y͕a͔+1Z͕v͔+1f͔Ð͕)e03",
    ///     "+(-1X͕d͔-1Y͕e͔-1Z͕f͔+1v͔Ð͕)e40",
    ///     "+1a͔v͕e23",
    ///     "+1b͔v͕e31",
    ///     "+1c͔v͕e12",
    ///     "+1d͔v͕e41",
    ///     "+1e͔v͕e42",
    ///     "+1f͔v͕e43",
    ///     "+(+1Y͕f͔-1Z͕e͔+1a͔Ð͕)e0324",
    ///     "+(-1X͕f͔+1Z͕d͔+1b͔Ð͕)e0134",
    ///     "+(+1X͕e͔-1Y͕d͔+1c͔Ð͕)e0214",
    ///     "+(+1X͕a͔+1Y͕b͔+1Z͕c͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "-1x͔x͕-1y͔y͕-1z͔z͕-1ð͔ð͕",
    ///     "+(-1B͔z͕+1B͕z͔+1C͔y͕-1C͕y͔+1D͔ð͕-1D͕ð͔)e01",
    ///     "+(+1A͔z͕-1A͕z͔-1C͔x͕+1C͕x͔+1E͔ð͕-1E͕ð͔)e02",
    ///     "+(-1A͔y͕+1A͕y͔+1B͔x͕-1B͕x͔+1F͔ð͕-1F͕ð͔)e03",
    ///     "+(-1D͔x͕+1D͕x͔-1E͔y͕+1E͕y͔-1F͔z͕+1F͕z͔)e40",
    ///     "+(-1y͔z͕+1y͕z͔)e23",
    ///     "+(+1x͔z͕-1x͕z͔)e31",
    ///     "+(-1x͔y͕+1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    ///     "+(+1A͔ð͕+1A͕ð͔-1E͔z͕-1E͕z͔+1F͔y͕+1F͕y͔)e0324",
    ///     "+(+1B͔ð͕+1B͕ð͔+1D͔z͕+1D͕z͔-1F͔x͕-1F͕x͔)e0134",
    ///     "+(+1C͔ð͕+1C͕ð͔-1D͔y͕-1D͕y͔+1E͔x͕+1E͕x͔)e0214",
    ///     "+(+1A͔x͕+1A͕x͔+1B͔y͕+1B͕y͔+1C͔z͕+1C͕z͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::plane() + Self::direction()
    }
    /// The multivector of double motor $`m_2 \equiv s + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "+1v͔v͕",
    ///     "+(+1X͕v͔+1Y͕c͔-1Z͕b͔+1d͔Ð͕)e01",
    ///     "+(-1X͕c͔+1Y͕v͔+1Z͕a͔+1e͔Ð͕)e02",
    ///     "+(+1X͕b͔-1Y͕a͔+1Z͕v͔+1f͔Ð͕)e03",
    ///     "+(-1X͕d͔-1Y͕e͔-1Z͕f͔+1v͔Ð͕)e40",
    ///     "+1a͔v͕e23",
    ///     "+1b͔v͕e31",
    ///     "+1c͔v͕e12",
    ///     "+1d͔v͕e41",
    ///     "+1e͔v͕e42",
    ///     "+1f͔v͕e43",
    ///     "+1v͕w͔e1234",
    ///     "+(+1X͕w͔+1Y͕f͔-1Z͕e͔+1a͔Ð͕)e0324",
    ///     "+(-1X͕f͔+1Y͕w͔+1Z͕d͔+1b͔Ð͕)e0134",
    ///     "+(+1X͕e͔-1Y͕d͔+1Z͕w͔+1c͔Ð͕)e0214",
    ///     "+(+1X͕a͔+1Y͕b͔+1Z͕c͔-1w͔Ð͕)e0123",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕",
    ///     "+(-1Y͔c͕+1Y͕c͔+1Z͔b͕-1Z͕b͔+1d͔Ð͕-1d͕Ð͔)e01",
    ///     "+(+1X͔c͕-1X͕c͔-1Z͔a͕+1Z͕a͔+1e͔Ð͕-1e͕Ð͔)e02",
    ///     "+(-1X͔b͕+1X͕b͔+1Y͔a͕-1Y͕a͔+1f͔Ð͕-1f͕Ð͔)e03",
    ///     "+(+1X͔d͕-1X͕d͔+1Y͔e͕-1Y͕e͔+1Z͔f͕-1Z͕f͔)e40",
    ///     "+(-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1d͔f͕-1d͕f͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔-1d͔e͕+1d͕e͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔)e43",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    ///     "+(+1Y͔f͕+1Y͕f͔-1Z͔e͕-1Z͕e͔+1a͔Ð͕+1a͕Ð͔)e0324",
    ///     "+(-1X͔f͕-1X͕f͔+1Z͔d͕+1Z͕d͔+1b͔Ð͕+1b͕Ð͔)e0134",
    ///     "+(+1X͔e͕+1X͕e͔-1Y͔d͕-1Y͕d͔+1c͔Ð͕+1c͕Ð͔)e0214",
    ///     "+(+1X͔a͕+1X͕a͔+1Y͔b͕+1Y͕b͔+1Z͔c͕+1Z͕c͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::plane() + Self::point()
    }
    /// The multivector of rotoreflector $`f_r \equiv v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// format_eq!(rotoreflector, [
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔+1v͕ð͔)e4",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔)e234",
    ///     "+(+1b͕ð͔-1d͕z͔+1f͕x͔)e314",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔)e123",
    /// ]);
    ///
    /// let squared_norm = Vee::rotoreflector().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+2xx+2yy+2zz+2ðð",
    ///     "+(-2xx-2yy-2zz+2ðð)e1234", // Weight condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+[+2xx+2yy+2zz-2ðð]w͓", // Vanishes with weight condition.
    ///     "+(+[-2xx-2yy-2zz-2ðð]w͓)e1234",
    ///     "+(+[+4xð]Ð͓+[-4zð]Y͓+[+4yð]Z͓)e0324",
    ///     "+(+[+4zð]X͓+[-4xð]Z͓+[+4yð]Ð͓)e0134",
    ///     "+(+[+4xð]Y͓+[+4zð]Ð͓+[-4yð]X͓)e0214",
    ///     "+(+[+4zð]Z͓+[+4xð]X͓+[+4yð]Y͓)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::line_displacement()
    }
    /// The multivector of transflector $`f_t \equiv v + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔)e0",
    ///     "+1v͕x͔e1",
    ///     "+1v͕y͔e2",
    ///     "+1v͕z͔e3",
    ///     "+1v͕ð͔e4",
    ///     "+(+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(+1Y͕ð͔+1y͔Ð͕)e024",
    ///     "+(+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1X͕y͔+1Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume() + Self::line_moment()
    }
    /// The multivector of flector $`f \equiv v + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let flector = Vee::volume().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔)e0",
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔+1v͕ð͔)e4",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔)e234",
    ///     "+(+1b͕ð͔-1d͕z͔+1f͕x͔)e314",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔)e123",
    ///     "+(-1W͔d͕+1X͕ð͔+1Y͕z͔-1Z͕y͔+1x͔Ð͕)e014",
    ///     "+(-1W͔e͕-1X͕z͔+1Y͕ð͔+1Z͕x͔+1y͔Ð͕)e024",
    ///     "+(-1W͔f͕+1X͕y͔-1Y͕x͔+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1W͔a͕-1X͕ð͔-1Y͕z͔+1Z͕y͔+1x͔Ð͕)e032",
    ///     "+(-1W͔b͕+1X͕z͔-1Y͕ð͔-1Z͕x͔+1y͔Ð͕)e013",
    ///     "+(-1W͔c͕-1X͕y͔+1Y͕x͔-1Z͕ð͔+1z͔Ð͕)e021",
    ///     "+(+1X͕x͔+1Y͕y͔+1Z͕z͔+1Ð͕ð͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::volume() + Self::line() + Self::pseudoscalar()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 5`$ (experimental, no inverse).
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
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{012345}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e012345()
    }
    /// The multivector of norm $`n \equiv s + \ell`$.
    ///
    /// Quadvector $`\ell`$ does not square to a scalar, therefore $`n`$ is **not** a Study number.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let quadvector_squared_norm = Vee::line().squared_norm();
    ///
    /// assert_eq!(quadvector_squared_norm.basis_blades(),
    ///     (Vee::scalar() + Vee::line_moment()).basis_blades());
    /// format_eq!(quadvector_squared_norm, [
    ///     "+1xx+1yy+1zz+1ðð+1øø",
    ///     "+(+2Dø-2Gð-2Jx)e0145",
    ///     "+(+2Eø-2Hð-2Jy)e0245",
    ///     "+(+2Fø-2Ið-2Jz)e0345",
    ///     "+(+2Aø-2Hz+2Iy)e0325",
    ///     "+(+2Bø+2Gz-2Ix)e0135",
    ///     "+(+2Cø-2Gy+2Hx)e0215",
    ///     "+(-2Að+2Ez-2Fy)e0324",
    ///     "+(-2Bð-2Dz+2Fx)e0134",
    ///     "+(-2Cð+2Dy-2Ex)e0214",
    ///     "+(-2Ax-2By-2Cz)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::line()
    }
    /// The multivector of bias $`v^4_\infty \equiv w\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v^4_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4 + þ\e_5`$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4() + Self::e5()
    }
    /// The multivector of $`4`$-volume $`v^4 \equiv v^4_0 + v^4_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume4() -> Self {
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
    /// v_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{40} + Þ\e_{05}
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
    /// \ell_0 \equiv x\e_{2345} + y\e_{3145} + z\e_{1245} + ð\e_{1235} + þ\e_{1234}
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
    /// P_\infty \equiv X\e_{03245} + Y\e_{01345} + Z\e_{02145} + Ð\e_{01235} + Þ\e_{01243}
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
    /// The multivector of single rotator $`r_1 \equiv s + v_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕+1ø͔ø͕",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    ///     "+(+1x͔ø͕-1x͕ø͔)e15",
    ///     "+(+1y͔ø͕-1y͕ø͔)e25",
    ///     "+(+1z͔ø͕-1z͕ø͔)e35",
    ///     "+(+1ð͔ø͕-1ð͕ø͔)e45",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕+1ø͔ø͕",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(+1x͔ð͕-1x͕ð͔)e41",
    ///     "+(+1y͔ð͕-1y͕ð͔)e42",
    ///     "+(+1z͔ð͕-1z͕ð͔)e43",
    ///     "+(+1x͔ø͕-1x͕ø͔)e15",
    ///     "+(+1y͔ø͕-1y͕ø͔)e25",
    ///     "+(+1z͔ø͕-1z͕ø͔)e35",
    ///     "+(-1ð͔ø͕+1ð͕ø͔)e45",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::volume_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕-1g͔g͕-1h͔h͕-1i͔i͕-1j͔j͕+1v͔v͕",
    ///     "+(+1a͔v͕+1a͕v͔-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔-1h͔i͕+1h͕i͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1b͔v͕+1b͕v͔+1d͔f͕-1d͕f͔+1g͔i͕-1g͕i͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔+1c͔v͕+1c͕v͔-1d͔e͕+1d͕e͔-1g͔h͕+1g͕h͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔+1d͔v͕+1d͕v͔+1g͔j͕-1g͕j͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔+1e͔v͕+1e͕v͔+1h͔j͕-1h͕j͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔+1f͔v͕+1f͕v͔+1i͔j͕-1i͕j͔)e43",
    ///     "+(-1b͔i͕+1b͕i͔+1c͔h͕-1c͕h͔-1d͔j͕+1d͕j͔+1g͔v͕+1g͕v͔)e15",
    ///     "+(+1a͔i͕-1a͕i͔-1c͔g͕+1c͕g͔-1e͔j͕+1e͕j͔+1h͔v͕+1h͕v͔)e25",
    ///     "+(-1a͔h͕+1a͕h͔+1b͔g͕-1b͕g͔-1f͔j͕+1f͕j͔+1i͔v͕+1i͕v͔)e35",
    ///     "+(+1d͔g͕-1d͕g͔+1e͔h͕-1e͕h͔+1f͔i͕-1f͕i͔+1j͔v͕+1j͕v͔)e45",
    ///     "+(+1a͔j͕+1a͕j͔+1e͔i͕+1e͕i͔-1f͔h͕-1f͕h͔)e2345",
    ///     "+(+1b͔j͕+1b͕j͔-1d͔i͕-1d͕i͔+1f͔g͕+1f͕g͔)e3145",
    ///     "+(+1c͔j͕+1c͕j͔+1d͔h͕+1d͕h͔-1e͔g͕-1e͕g͔)e1245",
    ///     "+(+1a͔g͕+1a͕g͔+1b͔h͕+1b͕h͔+1c͔i͕+1c͕i͔)e1235",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕-1g͔g͕-1h͔h͕-1i͔i͕-1j͔j͕",
    ///     "+(-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔-1h͔i͕+1h͕i͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1d͔f͕-1d͕f͔+1g͔i͕-1g͕i͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔-1d͔e͕+1d͕e͔-1g͔h͕+1g͕h͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔+1g͔j͕-1g͕j͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔+1h͔j͕-1h͕j͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔+1i͔j͕-1i͕j͔)e43",
    ///     "+(-1b͔i͕+1b͕i͔+1c͔h͕-1c͕h͔-1d͔j͕+1d͕j͔)e15",
    ///     "+(+1a͔i͕-1a͕i͔-1c͔g͕+1c͕g͔-1e͔j͕+1e͕j͔)e25",
    ///     "+(-1a͔h͕+1a͕h͔+1b͔g͕-1b͕g͔-1f͔j͕+1f͕j͔)e35",
    ///     "+(+1d͔g͕-1d͕g͔+1e͔h͕-1e͕h͔+1f͔i͕-1f͕i͔)e45",
    ///     "+(+1a͔j͕+1a͕j͔+1e͔i͕+1e͕i͔-1f͔h͕-1f͕h͔)e2345",
    ///     "+(+1b͔j͕+1b͕j͔-1d͔i͕-1d͕i͔+1f͔g͕+1f͕g͔)e3145",
    ///     "+(+1c͔j͕+1c͕j͔+1d͔h͕+1d͕h͔-1e͔g͕-1e͕g͔)e1245",
    ///     "+(+1a͔g͕+1a͕g͔+1b͔h͕+1b͕h͔+1c͔i͕+1c͕i͔)e1235",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕-1g͔g͕-1h͔h͕-1i͔i͕-1j͔j͕",
    ///     "+(+1b͔c͕-1b͕c͔+1f͔g͕-1f͕g͔-1i͔j͕+1i͕j͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1e͔g͕-1e͕g͔+1h͔j͕-1h͕j͔)e31",
    ///     "+(+1a͔b͕-1a͕b͔+1e͔f͕-1e͕f͔-1h͔i͕+1h͕i͔)e12",
    ///     "+(-1a͔d͕+1a͕d͔+1f͔j͕-1f͕j͔+1g͔i͕-1g͕i͔)e41",
    ///     "+(+1b͔d͕-1b͕d͔+1e͔j͕-1e͕j͔-1g͔h͕+1g͕h͔)e42",
    ///     "+(-1c͔d͕+1c͕d͔-1e͔i͕+1e͕i͔-1f͔h͕+1f͕h͔)e43",
    ///     "+(-1b͔j͕+1b͕j͔-1c͔i͕+1c͕i͔+1d͔e͕-1d͕e͔)e15",
    ///     "+(-1a͔j͕+1a͕j͔+1c͔h͕-1c͕h͔-1d͔f͕+1d͕f͔)e25",
    ///     "+(+1a͔i͕-1a͕i͔+1b͔h͕-1b͕h͔+1d͔g͕-1d͕g͔)e35",
    ///     "+(+1a͔e͕-1a͕e͔+1b͔f͕-1b͕f͔+1c͔g͕-1c͕g͔)e45",
    ///     "+(-1b͔g͕-1b͕g͔+1c͔f͕+1c͕f͔+1d͔h͕+1d͕h͔)e2345",
    ///     "+(-1a͔g͕-1a͕g͔+1c͔e͕+1c͕e͔+1d͔i͕+1d͕i͔)e3145",
    ///     "+(-1a͔f͕-1a͕f͔+1b͔e͕+1b͕e͔+1d͔j͕+1d͕j͔)e1245",
    ///     "+(-1a͔h͕-1a͕h͔+1b͔i͕+1b͕i͔-1c͔j͕-1c͕j͔)e1235",
    ///     "+(-1e͔h͕-1e͕h͔+1f͔i͕+1f͕i͔-1g͔j͕-1g͕j͔)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::volume_displacement() + Self::line_displacement()
    }
    /// The multivector of translator $`t \equiv s + v_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "+1w͔w͕",
    ///     "+(-1X͔w͕+1X͕w͔)e01",
    ///     "+(-1Y͔w͕+1Y͕w͔)e02",
    ///     "+(-1Z͔w͕+1Z͕w͔)e03",
    ///     "+(-1w͔Ð͕+1w͕Ð͔)e40",
    ///     "+(+1w͔Ø͕-1w͕Ø͔)e05",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::volume_moment()
    }
    /// The multivector of simple single motor $`m_{s1} \equiv s + v`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let simple_single_motor = Vee::volume4().lhs() * Vee::volume4().rhs();
    ///
    /// assert_eq!(simple_single_motor.basis_blades(), Vee::simple_single_motor().basis_blades());
    /// format_eq!(simple_single_motor, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕+1ø͔ø͕",
    ///     "+(+1W͔x͕-1W͕x͔)e01",
    ///     "+(+1W͔y͕-1W͕y͔)e02",
    ///     "+(+1W͔z͕-1W͕z͔)e03",
    ///     "+(-1W͔ð͕+1W͕ð͔)e40",
    ///     "+(+1W͔ø͕-1W͕ø͔)e05",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(-1x͔ð͕+1x͕ð͔)e41",
    ///     "+(-1y͔ð͕+1y͕ð͔)e42",
    ///     "+(-1z͔ð͕+1z͕ð͔)e43",
    ///     "+(+1x͔ø͕-1x͕ø͔)e15",
    ///     "+(+1y͔ø͕-1y͕ø͔)e25",
    ///     "+(+1z͔ø͕-1z͕ø͔)e35",
    ///     "+(+1ð͔ø͕-1ð͕ø͔)e45",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_motor() -> Self {
        Self::scalar() + Self::volume()
    }
    /// The multivector of single motor $`m_1 \equiv s + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+1v͔v͕",
    ///     "+(+1X͕v͔+1Y͕c͔-1Z͕b͔+1d͔Ð͕+1g͔Ø͕)e01",
    ///     "+(-1X͕c͔+1Y͕v͔+1Z͕a͔+1e͔Ð͕+1h͔Ø͕)e02",
    ///     "+(+1X͕b͔-1Y͕a͔+1Z͕v͔+1f͔Ð͕+1i͔Ø͕)e03",
    ///     "+(-1X͕d͔-1Y͕e͔-1Z͕f͔-1j͔Ø͕+1v͔Ð͕)e40",
    ///     "+(-1X͕g͔-1Y͕h͔-1Z͕i͔+1j͔Ð͕+1v͔Ø͕)e05",
    ///     "+1a͔v͕e23",
    ///     "+1b͔v͕e31",
    ///     "+1c͔v͕e12",
    ///     "+1d͔v͕e41",
    ///     "+1e͔v͕e42",
    ///     "+1f͔v͕e43",
    ///     "+1g͔v͕e15",
    ///     "+1h͔v͕e25",
    ///     "+1i͔v͕e35",
    ///     "+1j͔v͕e45",
    ///     "+(+1X͕j͔-1d͔Ø͕+1g͔Ð͕)e0145",
    ///     "+(+1Y͕j͔-1e͔Ø͕+1h͔Ð͕)e0245",
    ///     "+(+1Z͕j͔-1f͔Ø͕+1i͔Ð͕)e0345",
    ///     "+(-1Y͕i͔+1Z͕h͔-1a͔Ø͕)e0325",
    ///     "+(+1X͕i͔-1Z͕g͔-1b͔Ø͕)e0135",
    ///     "+(-1X͕h͔+1Y͕g͔-1c͔Ø͕)e0215",
    ///     "+(+1Y͕f͔-1Z͕e͔+1a͔Ð͕)e0324",
    ///     "+(-1X͕f͔+1Z͕d͔+1b͔Ð͕)e0134",
    ///     "+(+1X͕e͔-1Y͕d͔+1c͔Ð͕)e0214",
    ///     "+(+1X͕a͔+1Y͕b͔+1Z͕c͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+1x͔x͕+1y͔y͕+1z͔z͕+1ð͔ð͕+1ø͔ø͕",
    ///     "+(+1B͔z͕-1B͕z͔-1C͔y͕+1C͕y͔-1D͔ð͕+1D͕ð͔-1G͔ø͕+1G͕ø͔)e01",
    ///     "+(-1A͔z͕+1A͕z͔+1C͔x͕-1C͕x͔-1E͔ð͕+1E͕ð͔-1H͔ø͕+1H͕ø͔)e02",
    ///     "+(+1A͔y͕-1A͕y͔-1B͔x͕+1B͕x͔-1F͔ð͕+1F͕ð͔-1I͔ø͕+1I͕ø͔)e03",
    ///     "+(+1D͔x͕-1D͕x͔+1E͔y͕-1E͕y͔+1F͔z͕-1F͕z͔+1J͔ø͕-1J͕ø͔)e40",
    ///     "+(+1G͔x͕-1G͕x͔+1H͔y͕-1H͕y͔+1I͔z͕-1I͕z͔-1J͔ð͕+1J͕ð͔)e05",
    ///     "+(+1y͔z͕-1y͕z͔)e23",
    ///     "+(-1x͔z͕+1x͕z͔)e31",
    ///     "+(+1x͔y͕-1x͕y͔)e12",
    ///     "+(+1x͔ð͕-1x͕ð͔)e41",
    ///     "+(+1y͔ð͕-1y͕ð͔)e42",
    ///     "+(+1z͔ð͕-1z͕ð͔)e43",
    ///     "+(+1x͔ø͕-1x͕ø͔)e15",
    ///     "+(+1y͔ø͕-1y͕ø͔)e25",
    ///     "+(+1z͔ø͕-1z͕ø͔)e35",
    ///     "+(-1ð͔ø͕+1ð͕ø͔)e45",
    ///     "+(+1D͔ø͕+1D͕ø͔-1G͔ð͕-1G͕ð͔-1J͔x͕-1J͕x͔)e0145",
    ///     "+(+1E͔ø͕+1E͕ø͔-1H͔ð͕-1H͕ð͔-1J͔y͕-1J͕y͔)e0245",
    ///     "+(+1F͔ø͕+1F͕ø͔-1I͔ð͕-1I͕ð͔-1J͔z͕-1J͕z͔)e0345",
    ///     "+(+1A͔ø͕+1A͕ø͔-1H͔z͕-1H͕z͔+1I͔y͕+1I͕y͔)e0325",
    ///     "+(+1B͔ø͕+1B͕ø͔+1G͔z͕+1G͕z͔-1I͔x͕-1I͕x͔)e0135",
    ///     "+(+1C͔ø͕+1C͕ø͔-1G͔y͕-1G͕y͔+1H͔x͕+1H͕x͔)e0215",
    ///     "+(-1A͔ð͕-1A͕ð͔+1E͔z͕+1E͕z͔-1F͔y͕-1F͕y͔)e0324",
    ///     "+(-1B͔ð͕-1B͕ð͔-1D͔z͕-1D͕z͔+1F͔x͕+1F͕x͔)e0134",
    ///     "+(-1C͔ð͕-1C͕ð͔+1D͔y͕+1D͕y͔-1E͔x͕-1E͕x͔)e0214",
    ///     "+(-1A͔x͕-1A͕x͔-1B͔y͕-1B͕y͔-1C͔z͕-1C͕z͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line_moment()
    }
    /// The multivector of simple double motor $`m_{s2} \equiv s + v + \ell`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let simple_double_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(simple_double_motor.basis_blades(), Vee::simple_double_motor().basis_blades());
    /// format_eq!(simple_double_motor, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕-1g͔g͕-1h͔h͕-1i͔i͕-1j͔j͕",
    ///     "+(-1Y͔c͕+1Y͕c͔+1Z͔b͕-1Z͕b͔+1d͔Ð͕-1d͕Ð͔+1g͔Ø͕-1g͕Ø͔)e01",
    ///     "+(+1X͔c͕-1X͕c͔-1Z͔a͕+1Z͕a͔+1e͔Ð͕-1e͕Ð͔+1h͔Ø͕-1h͕Ø͔)e02",
    ///     "+(-1X͔b͕+1X͕b͔+1Y͔a͕-1Y͕a͔+1f͔Ð͕-1f͕Ð͔+1i͔Ø͕-1i͕Ø͔)e03",
    ///     "+(+1X͔d͕-1X͕d͔+1Y͔e͕-1Y͕e͔+1Z͔f͕-1Z͕f͔-1j͔Ø͕+1j͕Ø͔)e40",
    ///     "+(+1X͔g͕-1X͕g͔+1Y͔h͕-1Y͕h͔+1Z͔i͕-1Z͕i͔+1j͔Ð͕-1j͕Ð͔)e05",
    ///     "+(-1b͔c͕+1b͕c͔-1e͔f͕+1e͕f͔-1h͔i͕+1h͕i͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1d͔f͕-1d͕f͔+1g͔i͕-1g͕i͔)e31",
    ///     "+(-1a͔b͕+1a͕b͔-1d͔e͕+1d͕e͔-1g͔h͕+1g͕h͔)e12",
    ///     "+(-1b͔f͕+1b͕f͔+1c͔e͕-1c͕e͔+1g͔j͕-1g͕j͔)e41",
    ///     "+(+1a͔f͕-1a͕f͔-1c͔d͕+1c͕d͔+1h͔j͕-1h͕j͔)e42",
    ///     "+(-1a͔e͕+1a͕e͔+1b͔d͕-1b͕d͔+1i͔j͕-1i͕j͔)e43",
    ///     "+(-1b͔i͕+1b͕i͔+1c͔h͕-1c͕h͔-1d͔j͕+1d͕j͔)e15",
    ///     "+(+1a͔i͕-1a͕i͔-1c͔g͕+1c͕g͔-1e͔j͕+1e͕j͔)e25",
    ///     "+(-1a͔h͕+1a͕h͔+1b͔g͕-1b͕g͔-1f͔j͕+1f͕j͔)e35",
    ///     "+(+1d͔g͕-1d͕g͔+1e͔h͕-1e͕h͔+1f͔i͕-1f͕i͔)e45",
    ///     "+(+1a͔j͕+1a͕j͔+1e͔i͕+1e͕i͔-1f͔h͕-1f͕h͔)e2345",
    ///     "+(+1b͔j͕+1b͕j͔-1d͔i͕-1d͕i͔+1f͔g͕+1f͕g͔)e3145",
    ///     "+(+1c͔j͕+1c͕j͔+1d͔h͕+1d͕h͔-1e͔g͕-1e͕g͔)e1245",
    ///     "+(+1a͔g͕+1a͕g͔+1b͔h͕+1b͕h͔+1c͔i͕+1c͕i͔)e1235",
    ///     "+(-1a͔d͕-1a͕d͔-1b͔e͕-1b͕e͔-1c͔f͕-1c͕f͔)e1234",
    ///     "+(+1X͔j͕+1X͕j͔-1d͔Ø͕-1d͕Ø͔+1g͔Ð͕+1g͕Ð͔)e0145",
    ///     "+(+1Y͔j͕+1Y͕j͔-1e͔Ø͕-1e͕Ø͔+1h͔Ð͕+1h͕Ð͔)e0245",
    ///     "+(+1Z͔j͕+1Z͕j͔-1f͔Ø͕-1f͕Ø͔+1i͔Ð͕+1i͕Ð͔)e0345",
    ///     "+(-1Y͔i͕-1Y͕i͔+1Z͔h͕+1Z͕h͔-1a͔Ø͕-1a͕Ø͔)e0325",
    ///     "+(+1X͔i͕+1X͕i͔-1Z͔g͕-1Z͕g͔-1b͔Ø͕-1b͕Ø͔)e0135",
    ///     "+(-1X͔h͕-1X͕h͔+1Y͔g͕+1Y͕g͔-1c͔Ø͕-1c͕Ø͔)e0215",
    ///     "+(+1Y͔f͕+1Y͕f͔-1Z͔e͕-1Z͕e͔+1a͔Ð͕+1a͕Ð͔)e0324",
    ///     "+(-1X͔f͕-1X͕f͔+1Z͔d͕+1Z͕d͔+1b͔Ð͕+1b͕Ð͔)e0134",
    ///     "+(+1X͔e͕+1X͕e͔-1Y͔d͕-1Y͕d͔+1c͔Ð͕+1c͕Ð͔)e0214",
    ///     "+(+1X͔a͕+1X͕a͔+1Y͔b͕+1Y͕b͔+1Z͔c͕+1Z͕c͔)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line()
    }
    /// The multivector of double motor $`m_2 \equiv s + v + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "+1v͔v͕",
    ///     "+(+1X͕v͔+1Y͕c͔-1Z͕b͔+1d͔Ð͕+1g͔Ø͕)e01",
    ///     "+(-1X͕c͔+1Y͕v͔+1Z͕a͔+1e͔Ð͕+1h͔Ø͕)e02",
    ///     "+(+1X͕b͔-1Y͕a͔+1Z͕v͔+1f͔Ð͕+1i͔Ø͕)e03",
    ///     "+(-1X͕d͔-1Y͕e͔-1Z͕f͔-1j͔Ø͕+1v͔Ð͕)e40",
    ///     "+(-1X͕g͔-1Y͕h͔-1Z͕i͔+1j͔Ð͕+1v͔Ø͕)e05",
    ///     "+1a͔v͕e23",
    ///     "+1b͔v͕e31",
    ///     "+1c͔v͕e12",
    ///     "+1d͔v͕e41",
    ///     "+1e͔v͕e42",
    ///     "+1f͔v͕e43",
    ///     "+1g͔v͕e15",
    ///     "+1h͔v͕e25",
    ///     "+1i͔v͕e35",
    ///     "+1j͔v͕e45",
    ///     "+1v͕x͔e2345",
    ///     "+1v͕y͔e3145",
    ///     "+1v͕z͔e1245",
    ///     "+1v͕ð͔e1235",
    ///     "+1v͕ø͔e1234",
    ///     "+(+1X͕j͔+1Y͕z͔-1Z͕y͔-1d͔Ø͕+1g͔Ð͕)e0145",
    ///     "+(-1X͕z͔+1Y͕j͔+1Z͕x͔-1e͔Ø͕+1h͔Ð͕)e0245",
    ///     "+(+1X͕y͔-1Y͕x͔+1Z͕j͔-1f͔Ø͕+1i͔Ð͕)e0345",
    ///     "+(+1X͕ð͔-1Y͕i͔+1Z͕h͔-1a͔Ø͕-1x͔Ð͕)e0325",
    ///     "+(+1X͕i͔+1Y͕ð͔-1Z͕g͔-1b͔Ø͕-1y͔Ð͕)e0135",
    ///     "+(-1X͕h͔+1Y͕g͔+1Z͕ð͔-1c͔Ø͕-1z͔Ð͕)e0215",
    ///     "+(+1X͕ø͔+1Y͕f͔-1Z͕e͔+1a͔Ð͕-1x͔Ø͕)e0324",
    ///     "+(-1X͕f͔+1Y͕ø͔+1Z͕d͔+1b͔Ð͕-1y͔Ø͕)e0134",
    ///     "+(+1X͕e͔-1Y͕d͔+1Z͕ø͔+1c͔Ð͕-1z͔Ø͕)e0214",
    ///     "+(+1X͕a͔+1Y͕b͔+1Z͕c͔-1Ð͕ø͔+1Ø͕ð͔)e0123",
    ///     "+(+1X͕x͔+1Y͕y͔+1Z͕z͔+1Ð͕ð͔+1Ø͕ø͔)I",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-1a͔a͕-1b͔b͕-1c͔c͕-1d͔d͕-1e͔e͕-1f͔f͕-1g͔g͕-1h͔h͕-1i͔i͕-1j͔j͕",
    ///     "+(-1B͔g͕+1B͕g͔+1C͔f͕-1C͕f͔+1D͔h͕-1D͕h͔+1F͔c͕-1F͕c͔-1G͔b͕+1G͕b͔+1H͔d͕-1H͕d͔)e01",
    ///     "+(-1A͔g͕+1A͕g͔+1C͔e͕-1C͕e͔+1D͔i͕-1D͕i͔+1E͔c͕-1E͕c͔-1G͔a͕+1G͕a͔+1I͔d͕-1I͕d͔)e02",
    ///     "+(-1A͔f͕+1A͕f͔+1B͔e͕-1B͕e͔+1D͔j͕-1D͕j͔+1E͔b͕-1E͕b͔-1F͔a͕+1F͕a͔+1J͔d͕-1J͕d͔)e03",
    ///     "+(-1A͔h͕+1A͕h͔+1B͔i͕-1B͕i͔-1C͔j͕+1C͕j͔-1H͔a͕+1H͕a͔+1I͔b͕-1I͕b͔-1J͔c͕+1J͕c͔)e40",
    ///     "+(-1E͔h͕+1E͕h͔+1F͔i͕-1F͕i͔-1G͔j͕+1G͕j͔-1H͔e͕+1H͕e͔+1I͔f͕-1I͕f͔-1J͔g͕+1J͕g͔)e05",
    ///     "+(+1b͔c͕-1b͕c͔+1f͔g͕-1f͕g͔-1i͔j͕+1i͕j͔)e23",
    ///     "+(+1a͔c͕-1a͕c͔+1e͔g͕-1e͕g͔+1h͔j͕-1h͕j͔)e31",
    ///     "+(+1a͔b͕-1a͕b͔+1e͔f͕-1e͕f͔-1h͔i͕+1h͕i͔)e12",
    ///     "+(-1a͔d͕+1a͕d͔+1f͔j͕-1f͕j͔+1g͔i͕-1g͕i͔)e41",
    ///     "+(+1b͔d͕-1b͕d͔+1e͔j͕-1e͕j͔-1g͔h͕+1g͕h͔)e42",
    ///     "+(-1c͔d͕+1c͕d͔-1e͔i͕+1e͕i͔-1f͔h͕+1f͕h͔)e43",
    ///     "+(-1b͔j͕+1b͕j͔-1c͔i͕+1c͕i͔+1d͔e͕-1d͕e͔)e15",
    ///     "+(-1a͔j͕+1a͕j͔+1c͔h͕-1c͕h͔-1d͔f͕+1d͕f͔)e25",
    ///     "+(+1a͔i͕-1a͕i͔+1b͔h͕-1b͕h͔+1d͔g͕-1d͕g͔)e35",
    ///     "+(+1a͔e͕-1a͕e͔+1b͔f͕-1b͕f͔+1c͔g͕-1c͕g͔)e45",
    ///     "+(-1b͔g͕-1b͕g͔+1c͔f͕+1c͕f͔+1d͔h͕+1d͕h͔)e2345",
    ///     "+(-1a͔g͕-1a͕g͔+1c͔e͕+1c͕e͔+1d͔i͕+1d͕i͔)e3145",
    ///     "+(-1a͔f͕-1a͕f͔+1b͔e͕+1b͕e͔+1d͔j͕+1d͕j͔)e1245",
    ///     "+(-1a͔h͕-1a͕h͔+1b͔i͕+1b͕i͔-1c͔j͕-1c͕j͔)e1235",
    ///     "+(-1e͔h͕-1e͕h͔+1f͔i͕+1f͕i͔-1g͔j͕-1g͕j͔)e1234",
    ///     "+(-1B͔c͕-1B͕c͔+1C͔b͕+1C͕b͔-1F͔g͕-1F͕g͔+1G͔f͕+1G͕f͔+1I͔j͕+1I͕j͔-1J͔i͕-1J͕i͔)e0145",
    ///     "+(-1A͔c͕-1A͕c͔+1C͔a͕+1C͕a͔-1E͔g͕-1E͕g͔+1G͔e͕+1G͕e͔-1H͔j͕-1H͕j͔+1J͔h͕+1J͕h͔)e0245",
    ///     "+(-1A͔b͕-1A͕b͔+1B͔a͕+1B͕a͔-1E͔f͕-1E͕f͔+1F͔e͕+1F͕e͔+1H͔i͕+1H͕i͔-1I͔h͕-1I͕h͔)e0345",
    ///     "+(+1A͔d͕+1A͕d͔-1D͔a͕-1D͕a͔-1F͔j͕-1F͕j͔-1G͔i͕-1G͕i͔+1I͔g͕+1I͕g͔+1J͔f͕+1J͕f͔)e0325",
    ///     "+(-1B͔d͕-1B͕d͔+1D͔b͕+1D͕b͔-1E͔j͕-1E͕j͔+1G͔h͕+1G͕h͔-1H͔g͕-1H͕g͔+1J͔e͕+1J͕e͔)e0135",
    ///     "+(+1C͔d͕+1C͕d͔-1D͔c͕-1D͕c͔+1E͔i͕+1E͕i͔+1F͔h͕+1F͕h͔-1H͔f͕-1H͕f͔-1I͔e͕-1I͕e͔)e0215",
    ///     "+(+1B͔j͕+1B͕j͔+1C͔i͕+1C͕i͔-1D͔e͕-1D͕e͔+1E͔d͕+1E͕d͔-1I͔c͕-1I͕c͔-1J͔b͕-1J͕b͔)e0324",
    ///     "+(+1A͔j͕+1A͕j͔-1C͔h͕-1C͕h͔+1D͔f͕+1D͕f͔-1F͔d͕-1F͕d͔+1H͔c͕+1H͕c͔-1J͔a͕-1J͕a͔)e0134",
    ///     "+(-1A͔i͕-1A͕i͔-1B͔h͕-1B͕h͔-1D͔g͕-1D͕g͔+1G͔d͕+1G͕d͔+1H͔b͕+1H͕b͔+1I͔a͕+1I͕a͔)e0214",
    ///     "+(-1A͔e͕-1A͕e͔-1B͔f͕-1B͕f͔-1C͔g͕-1C͕g͔+1E͔a͕+1E͕a͔+1F͔b͕+1F͕b͔+1G͔c͕+1G͕c͔)e0123",
    ///     "+(-1A͔a͕+1A͕a͔-1B͔b͕+1B͕b͔-1C͔c͕+1C͕c͔-1D͔d͕+1D͕d͔-1E͔e͕+1E͕e͔-1F͔f͕+1F͕f͔-1G͔g͕+1G͕g͔-1H͔h͕+1H͕h͔-1I͔i͕+1I͕i͔-1J͔j͕+1J͕j͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
    /// The multivector of single rotoreflector $`f_{r1} \equiv h_0 + p_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(single_rotoreflector.basis_blades(), Vee::single_rotoreflector().basis_blades());
    /// format_eq!(single_rotoreflector, [
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔)e345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement()
    }
    /// The multivector of double rotoreflector $`f_{r2} \equiv h_0 + p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_rotoreflector = Vee::normal().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(double_rotoreflector.basis_blades(), Vee::double_rotoreflector().basis_blades());
    /// format_eq!(double_rotoreflector, [
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔+1x͔ø͕-1x͕ø͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔-1y͔ø͕+1y͕ø͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔+1z͔ø͕-1z͕ø͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔-1ð͔ø͕-1ð͕ø͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔-1x͔ð͕-1x͕ð͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔+1y͔ð͕+1y͕ð͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔-1z͔ð͕-1z͕ð͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔-1y͔z͕+1y͕z͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔+1x͔z͕-1x͕z͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔-1x͔y͕+1x͕y͔)e345",
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕-1ð͔ð͕+1ø͔ø͕)e12345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv h + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔-1Ø͕ø͔)e0",
    ///     "+1v͕x͔e1",
    ///     "+1v͕y͔e2",
    ///     "+1v͕z͔e3",
    ///     "+1v͕ð͔e4",
    ///     "+1v͕ø͔e5",
    ///     "+(+1X͕ø͔-1x͔Ø͕)e015",
    ///     "+(-1Y͕ø͔+1y͔Ø͕)e052",
    ///     "+(+1Z͕ø͔-1z͔Ø͕)e035",
    ///     "+(+1Ð͕ø͔+1Ø͕ð͔)e054",
    ///     "+(+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(-1Y͕ð͔-1y͔Ð͕)e042",
    ///     "+(+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1X͕y͔+1Y͕x͔)e021",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume4() + Self::plane_moment()
    }
    /// The multivector of simple single flector $`f_{s1} \equiv h + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let simple_single_flector = Vee::volume4().lhs() * Vee::simple_single_motor().rhs();
    ///
    /// assert_eq!(simple_single_flector.basis_blades(),
    ///     Vee::simple_single_flector().basis_blades());
    /// format_eq!(simple_single_flector, [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔-1Ø͕ø͔)e0",
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1W͔g͕+1X͕ø͔-1x͔Ø͕)e015",
    ///     "+(-1W͔h͕-1Y͕ø͔+1y͔Ø͕)e052",
    ///     "+(+1W͔i͕+1Z͕ø͔-1z͔Ø͕)e035",
    ///     "+(-1W͔j͕+1Ð͕ø͔+1Ø͕ð͔)e054",
    ///     "+(-1W͔d͕+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(+1W͔e͕-1Y͕ð͔-1y͔Ð͕)e042",
    ///     "+(-1W͔f͕+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1W͔a͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(-1W͔b͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1W͔c͕-1X͕y͔+1Y͕x͔)e021",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔)e345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_flector() -> Self {
        Self::volume4() + Self::plane()
    }
    /// The multivector of single flector $`f_1 \equiv h + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_flector = Vee::volume4().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(single_flector.basis_blades(), Vee::single_flector().basis_blades());
    /// format_eq!(single_flector, [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔-1Ø͕ø͔)e0",
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1A͕ð͔+1E͕z͔-1F͕y͔+1W͔g͕+1X͕ø͔-1x͔Ø͕)e015",
    ///     "+(-1B͕ð͔+1D͕z͔-1F͕x͔-1W͔h͕-1Y͕ø͔+1y͔Ø͕)e052",
    ///     "+(+1C͕ð͔+1D͕y͔-1E͕x͔+1W͔i͕+1Z͕ø͔-1z͔Ø͕)e035",
    ///     "+(+1A͕x͔+1B͕y͔+1C͕z͔-1W͔j͕+1Ð͕ø͔+1Ø͕ð͔)e054",
    ///     "+(-1A͕ø͔+1H͕z͔-1I͕y͔-1W͔d͕+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(+1B͕ø͔+1G͕z͔-1I͕x͔+1W͔e͕-1Y͕ð͔-1y͔Ð͕)e042",
    ///     "+(-1C͕ø͔+1G͕y͔-1H͕x͔-1W͔f͕+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1D͕ø͔-1G͕ð͔+1J͕x͔-1W͔a͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(-1E͕ø͔-1H͕ð͔+1J͕y͔-1W͔b͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1F͕ø͔-1I͕ð͔+1J͕z͔-1W͔c͕-1X͕y͔+1Y͕x͔)e021",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔)e345",
    ///     "+(-1B͕z͔+1C͕y͔-1D͕ð͔+1G͕ø͔)e03245",
    ///     "+(+1A͕z͔-1C͕x͔-1E͕ð͔+1H͕ø͔)e01345",
    ///     "+(-1A͕y͔+1B͕x͔-1F͕ð͔+1I͕ø͔)e02145",
    ///     "+(+1D͕x͔+1E͕y͔+1F͕z͔+1J͕ø͔)e01235",
    ///     "+(-1G͕x͔-1H͕y͔-1I͕z͔-1J͕ð͔)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::volume4() + Self::plane() + Self::direction()
    }
    /// The multivector of double flector $`f_2 \equiv h + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_flector = Vee::volume4().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// format_eq!(double_flector, [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔-1Ø͕ø͔)e0",
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1A͕ð͔+1E͕z͔-1F͕y͔+1W͔g͕+1X͕ø͔-1x͔Ø͕)e015",
    ///     "+(-1B͕ð͔+1D͕z͔-1F͕x͔-1W͔h͕-1Y͕ø͔+1y͔Ø͕)e052",
    ///     "+(+1C͕ð͔+1D͕y͔-1E͕x͔+1W͔i͕+1Z͕ø͔-1z͔Ø͕)e035",
    ///     "+(+1A͕x͔+1B͕y͔+1C͕z͔-1W͔j͕+1Ð͕ø͔+1Ø͕ð͔)e054",
    ///     "+(-1A͕ø͔+1H͕z͔-1I͕y͔-1W͔d͕+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(+1B͕ø͔+1G͕z͔-1I͕x͔+1W͔e͕-1Y͕ð͔-1y͔Ð͕)e042",
    ///     "+(-1C͕ø͔+1G͕y͔-1H͕x͔-1W͔f͕+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1D͕ø͔-1G͕ð͔+1J͕x͔-1W͔a͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(-1E͕ø͔-1H͕ð͔+1J͕y͔-1W͔b͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1F͕ø͔-1I͕ð͔+1J͕z͔-1W͔c͕-1X͕y͔+1Y͕x͔)e021",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔+1x͔ø͕-1x͕ø͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔-1y͔ø͕+1y͕ø͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔+1z͔ø͕-1z͕ø͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔-1ð͔ø͕-1ð͕ø͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔-1x͔ð͕-1x͕ð͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔+1y͔ð͕+1y͕ð͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔-1z͔ð͕-1z͕ð͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔-1y͔z͕+1y͕z͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔+1x͔z͕-1x͕z͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔-1x͔y͕+1x͕y͔)e345",
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕-1ð͔ð͕+1ø͔ø͕)e12345",
    ///     "+(-1B͕z͔+1C͕y͔-1D͕ð͔+1G͕ø͔-1W͔x͕)e03245",
    ///     "+(+1A͕z͔-1C͕x͔-1E͕ð͔+1H͕ø͔-1W͔y͕)e01345",
    ///     "+(-1A͕y͔+1B͕x͔-1F͕ð͔+1I͕ø͔-1W͔z͕)e02145",
    ///     "+(+1D͕x͔+1E͕y͔+1F͕z͔+1J͕ø͔+1W͔ð͕)e01235",
    ///     "+(-1G͕x͔-1H͕y͔-1I͕z͔-1J͕ð͔-1W͔ø͕)e01243",
    /// ]);
    ///
    /// let double_flector = Vee::volume4().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// format_eq!(double_flector, [
    ///     "+(+1W͔v͕-1X͕x͔-1Y͕y͔-1Z͕z͔+1Ð͕ð͔-1Ø͕ø͔)e0",
    ///     "+(+1b͕z͔-1c͕y͔+1d͕ð͔-1g͕ø͔+1v͕x͔)e1",
    ///     "+(-1a͕z͔+1c͕x͔+1e͕ð͔-1h͕ø͔+1v͕y͔)e2",
    ///     "+(+1a͕y͔-1b͕x͔+1f͕ð͔-1i͕ø͔+1v͕z͔)e3",
    ///     "+(-1d͕x͔-1e͕y͔-1f͕z͔-1j͕ø͔+1v͕ð͔)e4",
    ///     "+(+1g͕x͔+1h͕y͔+1i͕z͔+1j͕ð͔+1v͕ø͔)e5",
    ///     "+(+1A͕ð͔+1E͕z͔-1F͕y͔+1W͔g͕+1X͕ø͔-1x͔Ø͕)e015",
    ///     "+(-1B͕ð͔+1D͕z͔-1F͕x͔-1W͔h͕-1Y͕ø͔+1y͔Ø͕)e052",
    ///     "+(+1C͕ð͔+1D͕y͔-1E͕x͔+1W͔i͕+1Z͕ø͔-1z͔Ø͕)e035",
    ///     "+(+1A͕x͔+1B͕y͔+1C͕z͔-1W͔j͕+1Ð͕ø͔+1Ø͕ð͔)e054",
    ///     "+(-1A͕ø͔+1H͕z͔-1I͕y͔-1W͔d͕+1X͕ð͔+1x͔Ð͕)e014",
    ///     "+(+1B͕ø͔+1G͕z͔-1I͕x͔+1W͔e͕-1Y͕ð͔-1y͔Ð͕)e042",
    ///     "+(-1C͕ø͔+1G͕y͔-1H͕x͔-1W͔f͕+1Z͕ð͔+1z͔Ð͕)e034",
    ///     "+(-1D͕ø͔-1G͕ð͔+1J͕x͔-1W͔a͕-1Y͕z͔+1Z͕y͔)e032",
    ///     "+(-1E͕ø͔-1H͕ð͔+1J͕y͔-1W͔b͕+1X͕z͔-1Z͕x͔)e013",
    ///     "+(-1F͕ø͔-1I͕ð͔+1J͕z͔-1W͔c͕-1X͕y͔+1Y͕x͔)e021",
    ///     "+(+1a͕ð͔+1e͕z͔-1f͕y͔+1x͔ø͕-1x͕ø͔)e234",
    ///     "+(-1b͕ð͔+1d͕z͔-1f͕x͔-1y͔ø͕+1y͕ø͔)e134",
    ///     "+(+1c͕ð͔+1d͕y͔-1e͕x͔+1z͔ø͕-1z͕ø͔)e124",
    ///     "+(+1a͕x͔+1b͕y͔+1c͕z͔-1ð͔ø͕-1ð͕ø͔)e123",
    ///     "+(-1a͕ø͔+1h͕z͔-1i͕y͔-1x͔ð͕-1x͕ð͔)e253",
    ///     "+(+1b͕ø͔+1g͕z͔-1i͕x͔+1y͔ð͕+1y͕ð͔)e315",
    ///     "+(-1c͕ø͔+1g͕y͔-1h͕x͔-1z͔ð͕-1z͕ð͔)e152",
    ///     "+(-1d͕ø͔-1g͕ð͔+1j͕x͔-1y͔z͕+1y͕z͔)e145",
    ///     "+(-1e͕ø͔-1h͕ð͔+1j͕y͔+1x͔z͕-1x͕z͔)e245",
    ///     "+(-1f͕ø͔-1i͕ð͔+1j͕z͔-1x͔y͕+1x͕y͔)e345",
    ///     "+(+1x͔x͕+1y͔y͕+1z͔z͕-1ð͔ð͕+1ø͔ø͕)e12345",
    ///     "+(-1B͕z͔+1C͕y͔-1D͕ð͔+1G͕ø͔+1V͕x͔-1W͔x͕)e03245",
    ///     "+(+1A͕z͔-1C͕x͔-1E͕ð͔+1H͕ø͔+1V͕y͔-1W͔y͕)e01345",
    ///     "+(-1A͕y͔+1B͕x͔-1F͕ð͔+1I͕ø͔+1V͕z͔-1W͔z͕)e02145",
    ///     "+(+1D͕x͔+1E͕y͔+1F͕z͔+1J͕ø͔+1V͕ð͔+1W͔ð͕)e01235",
    ///     "+(-1G͕x͔-1H͕y͔-1I͕z͔-1J͕ð͔+1V͕ø͔-1W͔ø͕)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_flector() -> Self {
        Self::volume4() + Self::plane() + Self::point()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 6`$ (experimental, no inverse).
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 6>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{0123456}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e0123456()
    }
    /// The multivector of bias $`v^5_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v^5_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4() + Self::e5() + Self::e6()
    }
    /// The multivector of $`5`$-volume $`v^5 \equiv v^5_0 + v^5_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume5() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of $`4`$-volume moment $`v^4_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume4_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e40() + Self::e05() + Self::e60()
    }
    /// The multivector of $`4`$-volume displacement $`v^4_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume4_displacement() -> Self {
        Self::e12()
            + Self::e31()
            + Self::e41()
            + Self::e15()
            + Self::e16()
            + Self::e23()
            + Self::e42()
            + Self::e25()
            + Self::e62()
            + Self::e43()
            + Self::e35()
            + Self::e36()
            + Self::e45()
            + Self::e64()
            + Self::e56()
    }
    /// The multivector of $`4`$-volume $`v^4 \equiv v^4_0 + v^4_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume4() -> Self {
        Self::volume4_moment() + Self::volume4_displacement()
    }
    /// The multivector of volume moment $`v_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume_moment() -> Self {
        Self::e021()
            + Self::e013()
            + Self::e014()
            + Self::e015()
            + Self::e016()
            + Self::e032()
            + Self::e042()
            + Self::e052()
            + Self::e062()
            + Self::e034()
            + Self::e035()
            + Self::e036()
            + Self::e054()
            + Self::e064()
            + Self::e056()
    }
    /// The multivector of volume displacement $`v_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e123()
            + Self::e124()
            + Self::e152()
            + Self::e126()
            + Self::e134()
            + Self::e315()
            + Self::e163()
            + Self::e145()
            + Self::e146()
            + Self::e165()
            + Self::e234()
            + Self::e253()
            + Self::e236()
            + Self::e245()
            + Self::e264()
            + Self::e256()
            + Self::e345()
            + Self::e346()
            + Self::e365()
            + Self::e456()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::volume_moment() + Self::volume_displacement()
    }
    /// The multivector of plane moment $`p_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e0123()
            + Self::e0124()
            + Self::e0125()
            + Self::e0162()
            + Self::e0134()
            + Self::e0135()
            + Self::e0136()
            + Self::e0145()
            + Self::e0146()
            + Self::e0156()
            + Self::e0234()
            + Self::e0235()
            + Self::e0263()
            + Self::e0245()
            + Self::e0264()
            + Self::e0265()
            + Self::e0345()
            + Self::e0346()
            + Self::e0356()
            + Self::e0465()
    }
    /// The multivector of plane displacement $`p_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e1234()
            + Self::e1235()
            + Self::e1263()
            + Self::e1245()
            + Self::e1264()
            + Self::e1256()
            + Self::e1345()
            + Self::e1364()
            + Self::e1356()
            + Self::e1465()
            + Self::e2345()
            + Self::e2364()
            + Self::e2356()
            + Self::e2465()
            + Self::e3465()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_moment() + Self::plane_displacement()
    }
    /// The multivector of line moment $`\ell_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e01234()
            + Self::e01235()
            + Self::e01236()
            + Self::e01245()
            + Self::e01264()
            + Self::e01265()
            + Self::e01345()
            + Self::e01346()
            + Self::e01356()
            + Self::e01456()
            + Self::e02345()
            + Self::e02364()
            + Self::e02365()
            + Self::e02456()
            + Self::e03456()
    }
    /// The multivector of line displacement $`ell_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e12345()
            + Self::e12346()
            + Self::e12356()
            + Self::e12456()
            + Self::e13465()
            + Self::e23456()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_moment() + Self::line_displacement()
    }
    /// The multivector of direction $`P_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e012345()
            + Self::e012436()
            + Self::e012356()
            + Self::e021456()
            + Self::e013456()
            + Self::e032456()
    }
    /// The multivector of weight $`P_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e123456()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::direction() + Self::weight()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 7`$ (experimental, no inverse).
///
/// ```gdef
/// \gdef\e{
///   \boldsymbol e
/// }
/// \gdef\I{
///   \boldsymbol I
/// }
/// ```
impl<const M: i8> Multivector<Pga<M, 7>> {
    /// The multivector of scalar $`s \equiv v\e`$ where $`\e \equiv 1`$.
    #[must_use]
    #[inline]
    pub fn scalar() -> Self {
        Self::e()
    }
    /// The multivector of pseudoscalar $`S \equiv V\I`$ where $`\I \equiv \e_{0123456}`$.
    #[must_use]
    #[inline]
    pub fn pseudoscalar() -> Self {
        Self::e01234567()
    }
    /// The multivector of bias $`v^6_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v^6_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn normal() -> Self {
        Self::e1() + Self::e2() + Self::e3() + Self::e4() + Self::e5() + Self::e6() + Self::e7()
    }
    /// The multivector of $`6`$-volume $`v^6 \equiv v^6_0 + v^6_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume6() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of $`5`$-volume moment $`v^5_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume5_moment() -> Self {
        Self::e01()
            + Self::e02()
            + Self::e03()
            + Self::e04()
            + Self::e05()
            + Self::e06()
            + Self::e07()
    }
    /// The multivector of $`5`$-volume displacement $`v^5_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume5_displacement() -> Self {
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
            + Self::e16()
            + Self::e62()
            + Self::e36()
            + Self::e64()
            + Self::e56()
            + Self::e17()
            + Self::e27()
            + Self::e37()
            + Self::e47()
            + Self::e57()
            + Self::e67()
    }
    /// The multivector of $`5`$-volume $`v^5 \equiv v^5_0 + v^5_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume5() -> Self {
        Self::volume5_moment() + Self::volume5_displacement()
    }
    /// The multivector of $`4`$-volume moment $`v^4_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume4_moment() -> Self {
        Self::e012()
            + Self::e013()
            + Self::e014()
            + Self::e015()
            + Self::e016()
            + Self::e017()
            + Self::e023()
            + Self::e024()
            + Self::e025()
            + Self::e026()
            + Self::e027()
            + Self::e034()
            + Self::e035()
            + Self::e036()
            + Self::e037()
            + Self::e045()
            + Self::e046()
            + Self::e047()
            + Self::e056()
            + Self::e057()
            + Self::e067()
    }
    /// The multivector of $`4`$-volume displacement $`v^4_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume4_displacement() -> Self {
        Self::e123()
            + Self::e124()
            + Self::e125()
            + Self::e126()
            + Self::e127()
            + Self::e134()
            + Self::e135()
            + Self::e136()
            + Self::e137()
            + Self::e145()
            + Self::e146()
            + Self::e147()
            + Self::e156()
            + Self::e157()
            + Self::e167()
            + Self::e234()
            + Self::e235()
            + Self::e236()
            + Self::e237()
            + Self::e245()
            + Self::e246()
            + Self::e247()
            + Self::e256()
            + Self::e257()
            + Self::e267()
            + Self::e345()
            + Self::e346()
            + Self::e347()
            + Self::e356()
            + Self::e357()
            + Self::e367()
            + Self::e456()
            + Self::e457()
            + Self::e467()
            + Self::e567()
    }
    /// The multivector of $`4`$-volume $`v^4 \equiv v^4_0 + v^4_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume4() -> Self {
        Self::volume4_moment() + Self::volume4_displacement()
    }
    /// The multivector of volume moment $`v_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume_moment() -> Self {
        Self::e0123()
            + Self::e0124()
            + Self::e0125()
            + Self::e0126()
            + Self::e0127()
            + Self::e0134()
            + Self::e0135()
            + Self::e0136()
            + Self::e0137()
            + Self::e0145()
            + Self::e0146()
            + Self::e0147()
            + Self::e0156()
            + Self::e0157()
            + Self::e0167()
            + Self::e0234()
            + Self::e0235()
            + Self::e0236()
            + Self::e0237()
            + Self::e0245()
            + Self::e0246()
            + Self::e0247()
            + Self::e0256()
            + Self::e0257()
            + Self::e0267()
            + Self::e0345()
            + Self::e0346()
            + Self::e0347()
            + Self::e0356()
            + Self::e0357()
            + Self::e0367()
            + Self::e0456()
            + Self::e0457()
            + Self::e0467()
            + Self::e0567()
    }
    /// The multivector of volume displacement $`v_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e1234()
            + Self::e1253()
            + Self::e1236()
            + Self::e1273()
            + Self::e1245()
            + Self::e1264()
            + Self::e1247()
            + Self::e1256()
            + Self::e1275()
            + Self::e1267()
            + Self::e1354()
            + Self::e1346()
            + Self::e1374()
            + Self::e1365()
            + Self::e1357()
            + Self::e1376()
            + Self::e1456()
            + Self::e1475()
            + Self::e1467()
            + Self::e1576()
            + Self::e2345()
            + Self::e2364()
            + Self::e2347()
            + Self::e2356()
            + Self::e2375()
            + Self::e2367()
            + Self::e2465()
            + Self::e2457()
            + Self::e2476()
            + Self::e2567()
            + Self::e3456()
            + Self::e3475()
            + Self::e3467()
            + Self::e3576()
            + Self::e4567()
    }
    /// The multivector of volume $`v \equiv v_0 + v_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume() -> Self {
        Self::volume_moment() + Self::volume_displacement()
    }
    /// The multivector of plane moment $`p_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e01243()
            + Self::e01235()
            + Self::e01263()
            + Self::e01237()
            + Self::e01254()
            + Self::e01246()
            + Self::e01274()
            + Self::e01265()
            + Self::e01257()
            + Self::e01276()
            + Self::e01345()
            + Self::e01364()
            + Self::e01347()
            + Self::e01356()
            + Self::e01375()
            + Self::e01367()
            + Self::e01465()
            + Self::e01457()
            + Self::e01476()
            + Self::e01567()
            + Self::e02354()
            + Self::e02346()
            + Self::e02374()
            + Self::e02365()
            + Self::e02357()
            + Self::e02376()
            + Self::e02456()
            + Self::e02475()
            + Self::e02467()
            + Self::e02576()
            + Self::e03465()
            + Self::e03457()
            + Self::e03476()
            + Self::e03567()
            + Self::e04576()
    }
    /// The multivector of plane displacement $`p_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e12345()
            + Self::e12364()
            + Self::e12347()
            + Self::e12356()
            + Self::e12375()
            + Self::e12367()
            + Self::e12465()
            + Self::e12457()
            + Self::e12476()
            + Self::e12567()
            + Self::e13456()
            + Self::e13475()
            + Self::e13467()
            + Self::e13576()
            + Self::e14567()
            + Self::e23465()
            + Self::e23457()
            + Self::e23476()
            + Self::e23567()
            + Self::e24576()
            + Self::e34567()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_moment() + Self::plane_displacement()
    }
    /// The multivector of line moment $`\ell_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e012345()
            + Self::e012364()
            + Self::e012347()
            + Self::e012356()
            + Self::e012357()
            + Self::e012367()
            + Self::e012465()
            + Self::e012457()
            + Self::e012476()
            + Self::e012576()
            + Self::e013456()
            + Self::e013457()
            + Self::e013467()
            + Self::e013567()
            + Self::e014567()
            + Self::e023465()
            + Self::e023457()
            + Self::e023647()
            + Self::e023657()
            + Self::e024567()
            + Self::e034567()
    }
    /// The multivector of line displacement $`ell_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e123456()
            + Self::e123475()
            + Self::e123467()
            + Self::e123576()
            + Self::e124567()
            + Self::e134657()
            + Self::e234567()
    }
    /// The multivector of line $`\ell \equiv \ell_0 + \ell_\infty`$.
    #[must_use]
    #[inline]
    pub fn line() -> Self {
        Self::line_moment() + Self::line_displacement()
    }
    // 7
    /// The multivector of direction $`P_\infty \equiv `$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e0124356()
            + Self::e0123457()
            + Self::e0124367()
            + Self::e0123567()
            + Self::e0214567()
            + Self::e0134567()
            + Self::e0324567()
    }
    /// The multivector of weight $`P_0 \equiv `$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e1234567()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::direction() + Self::weight()
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
    use super::{PgaP0, PgaP1, PgaP2, PgaP3, PgaP4, PgaP5, PgaP6, PgaP7};

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

    assert_eq!(!PgaP5::volume4(), PgaP5::point().swp());
    assert_eq!(!PgaP5::volume(), PgaP5::line().swp());
    assert_eq!(
        !PgaP5::plane(),
        (PgaP5::plane_moment() - PgaP5::plane_displacement()).swp()
    );
    assert_eq!(!!PgaP5::plane(), -PgaP5::plane());
    assert_eq!(!PgaP5::line(), PgaP5::volume().swp());
    assert_eq!(!PgaP5::point(), -PgaP5::volume4().swp());

    assert_eq!(!PgaP6::volume5(), PgaP6::point().swp());
    assert_eq!(!PgaP6::volume4(), PgaP6::line().swp());
    assert_eq!(!PgaP6::volume(), PgaP6::plane().swp());
    assert_eq!(!PgaP6::plane(), PgaP6::volume().swp());
    assert_eq!(!PgaP6::line(), PgaP6::volume4().swp());
    assert_eq!(!PgaP6::point(), PgaP6::volume5().swp());

    assert_eq!(!PgaP7::volume6(), PgaP7::point().swp());
    assert_eq!(!PgaP7::volume5(), PgaP7::line().swp());
    assert_eq!(!PgaP7::volume4(), PgaP7::plane().swp());
    assert_eq!(!PgaP7::volume(), PgaP7::volume().swp());
    assert_eq!(!PgaP7::plane(), -PgaP7::volume4().swp());
    assert_eq!(!PgaP7::line(), PgaP7::volume5().swp());
    assert_eq!(!PgaP7::point(), -PgaP7::volume6().swp());
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
        ("PgaE6", PgaE6::table()),
        ("PgaE7", PgaE7::table()),
        ("PgaH0", PgaH0::table()),
        ("PgaH1", PgaH1::table()),
        ("PgaH2", PgaH2::table()),
        ("PgaH3", PgaH3::table()),
        ("PgaH4", PgaH4::table()),
        ("PgaH5", PgaH5::table()),
        ("PgaH6", PgaH6::table()),
        ("PgaH7", PgaH7::table()),
        ("PgaP0", PgaP0::table()),
        ("PgaP1", PgaP1::table()),
        ("PgaP2", PgaP2::table()),
        ("PgaP3", PgaP3::table()),
        ("PgaP4", PgaP4::table()),
        ("PgaP5", PgaP5::table()),
        ("PgaP6", PgaP6::table()),
        ("PgaP7", PgaP7::table()),
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
