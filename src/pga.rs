// Copyright © 2025 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Planed-Based Pistachio Flavor -- Projective Geometric Algebra (PGA)

use super::{Algebra, Choose, Multivector};
use core::{
    cmp::Ordering,
    fmt::{self, Debug, Display, Error, Write},
    ops::{Mul, Not},
};

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
    ///
    /// # Errors
    ///
    /// Fails in case of formatting errors.
    pub fn table() -> Result<String, Error> {
        let basis_len = Self::basis().len();
        let blade_len = N as usize + 3;
        let table_len = blade_len * basis_len.pow(2) + basis_len;
        let mut table = String::with_capacity(table_len);
        let mut blade = String::with_capacity(blade_len);
        for row in Self::basis() {
            for col in Self::basis() {
                let (sig, mul) = row * col;
                blade.clear();
                if sig == 0 {
                    write!(&mut blade, "0")?;
                } else {
                    let sig = if sig > 0 { " " } else { "-" };
                    write!(&mut blade, "{sig}{mul}")?;
                }
                write!(&mut table, "{blade:>blade_len$}")?;
            }
            writeln!(&mut table)?;
        }
        debug_assert_eq!(table.len(), table_len);
        Ok(table)
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
    // 1
    (v, e),
    // 1
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
    // 1
    (w, e123456),
    // 1
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
    ///     "+w͔w͕",
    ///     "+(+W͔w͕-W͕w͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::translator();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[+vv]W͓+[+2Vv]w͓)e0",
    ///     "+[+vv]w͓e1",
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
    ///     "+(+X͔Y͕-X͕Y͔)e0",
    ///     "+(+Y͔w͕-Y͕w͔)e1",
    ///     "+(-X͔w͕+X͕w͔)e2",
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
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-W͔y͕+W͕y͔)e20",
    ///     "+(+W͔x͕-W͕x͔)e01",
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
    ///     "+x͔x͕+y͔y͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
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
    ///     "-w͔w͕",
    ///     "+(-Y͔w͕+Y͕w͔)e20",
    ///     "+(+X͔w͕-X͕w͔)e01",
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
    ///     "+x͔x͕+y͔y͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-W͔y͕+W͕y͔)e20",
    ///     "+(+W͔x͕-W͕x͔)e01",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::motor();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[+vv+ww]w͓)e12",
    ///     "+(+[+2vw]Y͓+[+2Xw-2Yv]w͓+[+vv-ww]X͓)e20",
    ///     "+(+[+2Xv+2Yw]w͓+[-2vw]X͓+[+vv-ww]Y͓)e01",
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
    ///     "+(+W͔v͕+X͕y͔-Y͕x͔)e0",
    ///     "+(+v͕x͔-w͕y͔)e1",
    ///     "+(+v͕y͔+w͕x͔)e2",
    ///     "+(+W͔w͕+X͕x͔+Y͕y͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::flector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+[-xx-yy]w͓)e12",
    ///     "+(+[+2xy]Y͓+[+2Vy+2Wx]w͓+[+xx-yy]X͓)e20",
    ///     "+(+[-2Vx+2Wy]w͓+[+2xy]X͓+[-xx+yy]Y͓)e01",
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
    ///     "+xx+yy+zz",
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
    ///     "+(-X͔X͕-Y͔Y͕-Z͔Z͕)e0",
    ///     "+(+X͔w͕-Y͕z͔+Z͕y͔)e1",
    ///     "+(+X͕z͔+Y͔w͕-Z͕x͔)e2",
    ///     "+(-X͕y͔+Y͕x͔+Z͔w͕)e3",
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
    ///     "-X͔w͕e23",
    ///     "-Y͔w͕e31",
    ///     "-Z͔w͕e12",
    /// ]);
    ///
    /// // A line through the origin as the meet of two planes through the origin.
    /// let displacement = Vee::normal().lhs() ^ Vee::normal().rhs();
    ///
    /// assert_eq!(displacement.basis_blades(), Vee::displacement().basis_blades());
    /// format_eq!(displacement, [
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
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
    ///     "+(+Y͔Z͕-Y͕Z͔)e01",
    ///     "+(-X͔Z͕+X͕Z͔)e02",
    ///     "+(+X͔Y͕-X͕Y͔)e03",
    ///     "+(-X͔w͕+X͕w͔)e23",
    ///     "+(-Y͔w͕+Y͕w͔)e31",
    ///     "+(-Z͔w͕+Z͕w͔)e12",
    /// ]);
    ///
    /// // A line as the meet of two planes.
    /// let line = Vee::plane().lhs() ^ Vee::plane().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// format_eq!(line, [
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
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
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
    ///     "+(-W͔x͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-W͔y͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-W͔z͕-X͕y͔+Y͕x͔)e021",
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
    ///     "+x͔x͕+y͔y͕+z͔z͕",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    /// ]);
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// format_eq!(rotator, [
    ///     "-x͔x͕-y͔y͕-z͔z͕",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
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
    ///     "-w͔w͕",
    ///     "+(+X͔w͕-X͕w͔)e01",
    ///     "+(+Y͔w͕-Y͕w͔)e02",
    ///     "+(+Z͔w͕-Z͕w͔)e03",
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
    /// format_eq!(simple_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
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
    ///     "-x͔x͕-y͔y͕-z͔z͕",
    ///     "+(-Y͔z͕+Y͕z͔+Z͔y͕-Z͕y͔)e01",
    ///     "+(+X͔z͕-X͕z͔-Z͔x͕+Z͕x͔)e02",
    ///     "+(-X͔y͕+X͕y͔+Y͔x͕-Y͕x͔)e03",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(+X͔x͕+X͕x͔+Y͔y͕+Y͕y͔+Z͔z͕+Z͕z͔)I",
    /// ]);
    ///
    /// let motor = Vee::rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// format_eq!(motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕z͔-Z͕y͔)e01",
    ///     "+(-X͕z͔+Y͕v͔+Z͕x͔)e02",
    ///     "+(+X͕y͔-Y͕x͔+Z͕v͔)e03",
    ///     "+v͕x͔e23",
    ///     "+v͕y͔e31",
    ///     "+v͕z͔e12",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔)I",
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
    ///     "+(+v͕x͔-y͔z͕+y͕z͔)e1",
    ///     "+(+v͕y͔+x͔z͕-x͕z͔)e2",
    ///     "+(+v͕z͔-x͔y͕+x͕y͔)e3",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
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
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+(-Y͕z͔+Z͕y͔)e032",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(-X͕y͔+Y͕x͔)e021",
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
    /// format_eq!(flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔)e0",
    ///     "+(+v͕x͔-y͔z͕+y͕z͔)e1",
    ///     "+(+v͕y͔+x͔z͕-x͕z͔)e2",
    ///     "+(+v͕z͔-x͔y͕+x͕y͔)e3",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕)e123",
    ///     "+(+V͕x͔-W͔x͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(+V͕y͔-W͔y͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+V͕z͔-W͔z͕-X͕y͔+Y͕x͔)e021",
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
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "-x͔x͕-y͔y͕-z͔z͕-ð͔ð͕",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    /// ]);
    ///
    /// let squared_norm = Vee::single_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+aa+bb+cc+dd+ee+ff+vv",
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
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕+v͔v͕",
    ///     "+(+a͔v͕+a͕v͔-b͔c͕+b͕c͔-e͔f͕+e͕f͔)e23",
    ///     "+(+a͔c͕-a͕c͔+b͔v͕+b͕v͔+d͔f͕-d͕f͔)e31",
    ///     "+(-a͔b͕+a͕b͔+c͔v͕+c͕v͔-d͔e͕+d͕e͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔+d͔v͕+d͕v͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔+e͔v͕+e͕v͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔+f͔v͕+f͕v͔)e43",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔)e23",
    ///     "+(+a͔c͕-a͕c͔+d͔f͕-d͕f͔)e31",
    ///     "+(-a͔b͕+a͕b͔-d͔e͕+d͕e͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔)e43",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    /// ]);
    ///
    /// let squared_norm = Vee::double_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(squared_norm, [
    ///     "+aa+bb+cc+dd+ee+ff+vv+ww",
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
    ///     "+w͔w͕",
    ///     "+(-X͔w͕+X͕w͔)e01",
    ///     "+(-Y͔w͕+Y͕w͔)e02",
    ///     "+(-Z͔w͕+Z͕w͔)e03",
    ///     "+(-w͔Ð͕+w͕Ð͔)e40",
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
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(-W͔ð͕+W͕ð͔)e40",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    /// ]);
    ///
    /// let squared_norm = Vee::simple_motor().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(squared_norm, [
    ///     // Scalar condition.
    ///     "+aa+bb+cc+dd+ee+ff+vv",
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
    ///     "+(+[+aa+bb+cc+dd+ee+ff+vv]w͓)e1234",
    ///     "+(+[+2ac-2bv-2df]Z͓+[-2Xv-2Yc+2Zb-2dÐ]w͓+[+2bf-2ce-2dv]Ð͓\
    ///        +[+aa-bb-cc-dd+ee+ff+vv]X͓+[+2ab+2cv-2de]Y͓)e0324",
    ///     "+(+[+2Xc-2Yv-2Za-2eÐ]w͓+[-2af+2cd-2ev]Ð͓+[+2ab-2cv-2de]X͓\
    ///        +[-aa+bb-cc+dd-ee+ff+vv]Y͓+[+2av+2bc-2ef]Z͓)e0134",
    ///     "+(+[+2ae-2bd-2fv]Ð͓+[+2ac+2bv-2df]X͓+[-2av+2bc-2ef]Y͓\
    ///        +[-aa-bb+cc+dd+ee-ff+vv]Z͓+[-2Xb+2Ya-2Zv-2fÐ]w͓)e0214",
    ///     "+(+[+2bf-2ce+2dv]X͓+[-2af+2cd+2ev]Y͓+[+2ae-2bd+2fv]Z͓\
    ///        +[-2Xd-2Ye-2Zf+2vÐ]w͓+[+aa+bb+cc-dd-ee-ff+vv]Ð͓)e0123",
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
    ///     "+(+[+aa+bb+cc+dd+ee+ff+vv]W͓\
    ///        +[+2Xv-2Yc+2Zb-2dÐ]x͓+[+2Xc+2Yv-2Za-2eÐ]y͓+[-2Xb+2Ya+2Zv-2fÐ]z͓+[-2Xd-2Ye-2Zf-2vÐ]ð͓)e0",
    ///     "+(+[+2ab+2cv-2de]y͓+[+2ac-2bv-2df]z͓+[+2bf-2ce-2dv]ð͓+[+aa-bb-cc-dd+ee+ff+vv]x͓)e1",
    ///     "+(+[+2av+2bc-2ef]z͓+[-2af+2cd-2ev]ð͓+[+2ab-2cv-2de]x͓+[-aa+bb-cc+dd-ee+ff+vv]y͓)e2",
    ///     "+(+[+2ae-2bd-2fv]ð͓+[+2ac+2bv-2df]x͓+[-2av+2bc-2ef]y͓+[-aa-bb+cc+dd+ee-ff+vv]z͓)e3",
    ///     "+(+[+2bf-2ce+2dv]x͓+[-2af+2cd+2ev]y͓+[+2ae-2bd+2fv]z͓+[+aa+bb+cc-dd-ee-ff+vv]ð͓)e4",
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
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕c͔-Z͕b͔+d͔Ð͕)e01",
    ///     "+(-X͕c͔+Y͕v͔+Z͕a͔+e͔Ð͕)e02",
    ///     "+(+X͕b͔-Y͕a͔+Z͕v͔+f͔Ð͕)e03",
    ///     "+(-X͕d͔-Y͕e͔-Z͕f͔+v͔Ð͕)e40",
    ///     "+a͔v͕e23",
    ///     "+b͔v͕e31",
    ///     "+c͔v͕e12",
    ///     "+d͔v͕e41",
    ///     "+e͔v͕e42",
    ///     "+f͔v͕e43",
    ///     "+(+Y͕f͔-Z͕e͔+a͔Ð͕)e0324",
    ///     "+(-X͕f͔+Z͕d͔+b͔Ð͕)e0134",
    ///     "+(+X͕e͔-Y͕d͔+c͔Ð͕)e0214",
    ///     "+(+X͕a͔+Y͕b͔+Z͕c͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "-x͔x͕-y͔y͕-z͔z͕-ð͔ð͕",
    ///     "+(-B͔z͕+B͕z͔+C͔y͕-C͕y͔+D͔ð͕-D͕ð͔)e01",
    ///     "+(+A͔z͕-A͕z͔-C͔x͕+C͕x͔+E͔ð͕-E͕ð͔)e02",
    ///     "+(-A͔y͕+A͕y͔+B͔x͕-B͕x͔+F͔ð͕-F͕ð͔)e03",
    ///     "+(-D͔x͕+D͕x͔-E͔y͕+E͕y͔-F͔z͕+F͕z͔)e40",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(+x͔z͕-x͕z͔)e31",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    ///     "+(+A͔ð͕+A͕ð͔-E͔z͕-E͕z͔+F͔y͕+F͕y͔)e0324",
    ///     "+(+B͔ð͕+B͕ð͔+D͔z͕+D͕z͔-F͔x͕-F͕x͔)e0134",
    ///     "+(+C͔ð͕+C͕ð͔-D͔y͕-D͕y͔+E͔x͕+E͕x͔)e0214",
    ///     "+(+A͔x͕+A͕x͔+B͔y͕+B͕y͔+C͔z͕+C͕z͔)e0123",
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
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕c͔-Z͕b͔+d͔Ð͕)e01",
    ///     "+(-X͕c͔+Y͕v͔+Z͕a͔+e͔Ð͕)e02",
    ///     "+(+X͕b͔-Y͕a͔+Z͕v͔+f͔Ð͕)e03",
    ///     "+(-X͕d͔-Y͕e͔-Z͕f͔+v͔Ð͕)e40",
    ///     "+a͔v͕e23",
    ///     "+b͔v͕e31",
    ///     "+c͔v͕e12",
    ///     "+d͔v͕e41",
    ///     "+e͔v͕e42",
    ///     "+f͔v͕e43",
    ///     "+v͕w͔e1234",
    ///     "+(+X͕w͔+Y͕f͔-Z͕e͔+a͔Ð͕)e0324",
    ///     "+(-X͕f͔+Y͕w͔+Z͕d͔+b͔Ð͕)e0134",
    ///     "+(+X͕e͔-Y͕d͔+Z͕w͔+c͔Ð͕)e0214",
    ///     "+(+X͕a͔+Y͕b͔+Z͕c͔-w͔Ð͕)e0123",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕",
    ///     "+(-Y͔c͕+Y͕c͔+Z͔b͕-Z͕b͔+d͔Ð͕-d͕Ð͔)e01",
    ///     "+(+X͔c͕-X͕c͔-Z͔a͕+Z͕a͔+e͔Ð͕-e͕Ð͔)e02",
    ///     "+(-X͔b͕+X͕b͔+Y͔a͕-Y͕a͔+f͔Ð͕-f͕Ð͔)e03",
    ///     "+(+X͔d͕-X͕d͔+Y͔e͕-Y͕e͔+Z͔f͕-Z͕f͔)e40",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔)e23",
    ///     "+(+a͔c͕-a͕c͔+d͔f͕-d͕f͔)e31",
    ///     "+(-a͔b͕+a͕b͔-d͔e͕+d͕e͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔)e43",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    ///     "+(+Y͔f͕+Y͕f͔-Z͔e͕-Z͕e͔+a͔Ð͕+a͕Ð͔)e0324",
    ///     "+(-X͔f͕-X͕f͔+Z͔d͕+Z͕d͔+b͔Ð͕+b͕Ð͔)e0134",
    ///     "+(+X͔e͕+X͕e͔-Y͔d͕-Y͕d͔+c͔Ð͕+c͕Ð͔)e0214",
    ///     "+(+X͔a͕+X͕a͔+Y͔b͕+Y͕b͔+Z͔c͕+Z͕c͔)e0123",
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
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔+v͕ð͔)e4",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔)e234",
    ///     "+(+b͕ð͔-d͕z͔+f͕x͔)e314",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔)e123",
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
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+(+X͕ð͔+x͔Ð͕)e014",
    ///     "+(+Y͕ð͔+y͔Ð͕)e024",
    ///     "+(+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-Y͕z͔+Z͕y͔)e032",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(-X͕y͔+Y͕x͔)e021",
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
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔)e0",
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔+v͕ð͔)e4",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔)e234",
    ///     "+(+b͕ð͔-d͕z͔+f͕x͔)e314",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔)e123",
    ///     "+(-W͔d͕+X͕ð͔+Y͕z͔-Z͕y͔+x͔Ð͕)e014",
    ///     "+(-W͔e͕-X͕z͔+Y͕ð͔+Z͕x͔+y͔Ð͕)e024",
    ///     "+(-W͔f͕+X͕y͔-Y͕x͔+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-W͔a͕-X͕ð͔-Y͕z͔+Z͕y͔+x͔Ð͕)e032",
    ///     "+(-W͔b͕+X͕z͔-Y͕ð͔-Z͕x͔+y͔Ð͕)e013",
    ///     "+(-W͔c͕-X͕y͔+Y͕x͔-Z͕ð͔+z͔Ð͕)e021",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔)I",
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
    ///     "+xx+yy+zz+ðð+øø",
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
    /// The multivector of normal $`v^4_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4 + ø\e_5`$.
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
    /// v_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{40} + Ø\e_{05}
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
    /// \ell_0 \equiv x\e_{2345} + y\e_{3145} + z\e_{1245} + ð\e_{1235} + ø\e_{1234}
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
    /// P_\infty \equiv X\e_{03245} + Y\e_{01345} + Z\e_{02145} + Ð\e_{01235} + Ø\e_{01243}
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
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔ð͕-x͕ð͔)e41",
    ///     "+(+y͔ð͕-y͕ð͔)e42",
    ///     "+(+z͔ð͕-z͕ð͔)e43",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(-ð͔ø͕+ð͕ø͔)e45",
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
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕+v͔v͕",
    ///     "+(+a͔v͕+a͕v͔-b͔c͕+b͕c͔-e͔f͕+e͕f͔-h͔i͕+h͕i͔)e23",
    ///     "+(+a͔c͕-a͕c͔+b͔v͕+b͕v͔+d͔f͕-d͕f͔+g͔i͕-g͕i͔)e31",
    ///     "+(-a͔b͕+a͕b͔+c͔v͕+c͕v͔-d͔e͕+d͕e͔-g͔h͕+g͕h͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔+d͔v͕+d͕v͔+g͔j͕-g͕j͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔+e͔v͕+e͕v͔+h͔j͕-h͕j͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔+f͔v͕+f͕v͔+i͔j͕-i͕j͔)e43",
    ///     "+(-b͔i͕+b͕i͔+c͔h͕-c͕h͔-d͔j͕+d͕j͔+g͔v͕+g͕v͔)e15",
    ///     "+(+a͔i͕-a͕i͔-c͔g͕+c͕g͔-e͔j͕+e͕j͔+h͔v͕+h͕v͔)e25",
    ///     "+(-a͔h͕+a͕h͔+b͔g͕-b͕g͔-f͔j͕+f͕j͔+i͔v͕+i͕v͔)e35",
    ///     "+(+d͔g͕-d͕g͔+e͔h͕-e͕h͔+f͔i͕-f͕i͔+j͔v͕+j͕v͔)e45",
    ///     "+(+a͔j͕+a͕j͔+e͔i͕+e͕i͔-f͔h͕-f͕h͔)e2345",
    ///     "+(+b͔j͕+b͕j͔-d͔i͕-d͕i͔+f͔g͕+f͕g͔)e3145",
    ///     "+(+c͔j͕+c͕j͔+d͔h͕+d͕h͔-e͔g͕-e͕g͔)e1245",
    ///     "+(+a͔g͕+a͕g͔+b͔h͕+b͕h͔+c͔i͕+c͕i͔)e1235",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-h͔i͕+h͕i͔)e23",
    ///     "+(+a͔c͕-a͕c͔+d͔f͕-d͕f͔+g͔i͕-g͕i͔)e31",
    ///     "+(-a͔b͕+a͕b͔-d͔e͕+d͕e͔-g͔h͕+g͕h͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔+g͔j͕-g͕j͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔+h͔j͕-h͕j͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔+i͔j͕-i͕j͔)e43",
    ///     "+(-b͔i͕+b͕i͔+c͔h͕-c͕h͔-d͔j͕+d͕j͔)e15",
    ///     "+(+a͔i͕-a͕i͔-c͔g͕+c͕g͔-e͔j͕+e͕j͔)e25",
    ///     "+(-a͔h͕+a͕h͔+b͔g͕-b͕g͔-f͔j͕+f͕j͔)e35",
    ///     "+(+d͔g͕-d͕g͔+e͔h͕-e͕h͔+f͔i͕-f͕i͔)e45",
    ///     "+(+a͔j͕+a͕j͔+e͔i͕+e͕i͔-f͔h͕-f͕h͔)e2345",
    ///     "+(+b͔j͕+b͕j͔-d͔i͕-d͕i͔+f͔g͕+f͕g͔)e3145",
    ///     "+(+c͔j͕+c͕j͔+d͔h͕+d͕h͔-e͔g͕-e͕g͔)e1245",
    ///     "+(+a͔g͕+a͕g͔+b͔h͕+b͕h͔+c͔i͕+c͕i͔)e1235",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(+b͔c͕-b͕c͔+f͔g͕-f͕g͔-i͔j͕+i͕j͔)e23",
    ///     "+(+a͔c͕-a͕c͔+e͔g͕-e͕g͔+h͔j͕-h͕j͔)e31",
    ///     "+(+a͔b͕-a͕b͔+e͔f͕-e͕f͔-h͔i͕+h͕i͔)e12",
    ///     "+(-a͔d͕+a͕d͔+f͔j͕-f͕j͔+g͔i͕-g͕i͔)e41",
    ///     "+(+b͔d͕-b͕d͔+e͔j͕-e͕j͔-g͔h͕+g͕h͔)e42",
    ///     "+(-c͔d͕+c͕d͔-e͔i͕+e͕i͔-f͔h͕+f͕h͔)e43",
    ///     "+(-b͔j͕+b͕j͔-c͔i͕+c͕i͔+d͔e͕-d͕e͔)e15",
    ///     "+(-a͔j͕+a͕j͔+c͔h͕-c͕h͔-d͔f͕+d͕f͔)e25",
    ///     "+(+a͔i͕-a͕i͔+b͔h͕-b͕h͔+d͔g͕-d͕g͔)e35",
    ///     "+(+a͔e͕-a͕e͔+b͔f͕-b͕f͔+c͔g͕-c͕g͔)e45",
    ///     "+(-b͔g͕-b͕g͔+c͔f͕+c͕f͔+d͔h͕+d͕h͔)e2345",
    ///     "+(-a͔g͕-a͕g͔+c͔e͕+c͕e͔+d͔i͕+d͕i͔)e3145",
    ///     "+(-a͔f͕-a͕f͔+b͔e͕+b͕e͔+d͔j͕+d͕j͔)e1245",
    ///     "+(-a͔h͕-a͕h͔+b͔i͕+b͕i͔-c͔j͕-c͕j͔)e1235",
    ///     "+(-e͔h͕-e͕h͔+f͔i͕+f͕i͔-g͔j͕-g͕j͔)e1234",
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
    ///     "+w͔w͕",
    ///     "+(-X͔w͕+X͕w͔)e01",
    ///     "+(-Y͔w͕+Y͕w͔)e02",
    ///     "+(-Z͔w͕+Z͕w͔)e03",
    ///     "+(-w͔Ð͕+w͕Ð͔)e40",
    ///     "+(+w͔Ø͕-w͕Ø͔)e05",
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
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(-W͔ð͕+W͕ð͔)e40",
    ///     "+(+W͔ø͕-W͕ø͔)e05",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(-x͔ð͕+x͕ð͔)e41",
    ///     "+(-y͔ð͕+y͕ð͔)e42",
    ///     "+(-z͔ð͕+z͕ð͔)e43",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
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
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕c͔-Z͕b͔+d͔Ð͕+g͔Ø͕)e01",
    ///     "+(-X͕c͔+Y͕v͔+Z͕a͔+e͔Ð͕+h͔Ø͕)e02",
    ///     "+(+X͕b͔-Y͕a͔+Z͕v͔+f͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕d͔-Y͕e͔-Z͕f͔-j͔Ø͕+v͔Ð͕)e40",
    ///     "+(-X͕g͔-Y͕h͔-Z͕i͔+j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e23",
    ///     "+b͔v͕e31",
    ///     "+c͔v͕e12",
    ///     "+d͔v͕e41",
    ///     "+e͔v͕e42",
    ///     "+f͔v͕e43",
    ///     "+g͔v͕e15",
    ///     "+h͔v͕e25",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+(+X͕j͔-d͔Ø͕+g͔Ð͕)e0145",
    ///     "+(+Y͕j͔-e͔Ø͕+h͔Ð͕)e0245",
    ///     "+(+Z͕j͔-f͔Ø͕+i͔Ð͕)e0345",
    ///     "+(-Y͕i͔+Z͕h͔-a͔Ø͕)e0325",
    ///     "+(+X͕i͔-Z͕g͔-b͔Ø͕)e0135",
    ///     "+(-X͕h͔+Y͕g͔-c͔Ø͕)e0215",
    ///     "+(+Y͕f͔-Z͕e͔+a͔Ð͕)e0324",
    ///     "+(-X͕f͔+Z͕d͔+b͔Ð͕)e0134",
    ///     "+(+X͕e͔-Y͕d͔+c͔Ð͕)e0214",
    ///     "+(+X͕a͔+Y͕b͔+Z͕c͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+B͔z͕-B͕z͔-C͔y͕+C͕y͔-D͔ð͕+D͕ð͔-G͔ø͕+G͕ø͔)e01",
    ///     "+(-A͔z͕+A͕z͔+C͔x͕-C͕x͔-E͔ð͕+E͕ð͔-H͔ø͕+H͕ø͔)e02",
    ///     "+(+A͔y͕-A͕y͔-B͔x͕+B͕x͔-F͔ð͕+F͕ð͔-I͔ø͕+I͕ø͔)e03",
    ///     "+(+D͔x͕-D͕x͔+E͔y͕-E͕y͔+F͔z͕-F͕z͔+J͔ø͕-J͕ø͔)e40",
    ///     "+(+G͔x͕-G͕x͔+H͔y͕-H͕y͔+I͔z͕-I͕z͔-J͔ð͕+J͕ð͔)e05",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(-x͔z͕+x͕z͔)e31",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔ð͕-x͕ð͔)e41",
    ///     "+(+y͔ð͕-y͕ð͔)e42",
    ///     "+(+z͔ð͕-z͕ð͔)e43",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(-ð͔ø͕+ð͕ø͔)e45",
    ///     "+(+D͔ø͕+D͕ø͔-G͔ð͕-G͕ð͔-J͔x͕-J͕x͔)e0145",
    ///     "+(+E͔ø͕+E͕ø͔-H͔ð͕-H͕ð͔-J͔y͕-J͕y͔)e0245",
    ///     "+(+F͔ø͕+F͕ø͔-I͔ð͕-I͕ð͔-J͔z͕-J͕z͔)e0345",
    ///     "+(+A͔ø͕+A͕ø͔-H͔z͕-H͕z͔+I͔y͕+I͕y͔)e0325",
    ///     "+(+B͔ø͕+B͕ø͔+G͔z͕+G͕z͔-I͔x͕-I͕x͔)e0135",
    ///     "+(+C͔ø͕+C͕ø͔-G͔y͕-G͕y͔+H͔x͕+H͕x͔)e0215",
    ///     "+(-A͔ð͕-A͕ð͔+E͔z͕+E͕z͔-F͔y͕-F͕y͔)e0324",
    ///     "+(-B͔ð͕-B͕ð͔-D͔z͕-D͕z͔+F͔x͕+F͕x͔)e0134",
    ///     "+(-C͔ð͕-C͕ð͔+D͔y͕+D͕y͔-E͔x͕-E͕x͔)e0214",
    ///     "+(-A͔x͕-A͕x͔-B͔y͕-B͕y͔-C͔z͕-C͕z͔)e0123",
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
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-Y͔c͕+Y͕c͔+Z͔b͕-Z͕b͔+d͔Ð͕-d͕Ð͔+g͔Ø͕-g͕Ø͔)e01",
    ///     "+(+X͔c͕-X͕c͔-Z͔a͕+Z͕a͔+e͔Ð͕-e͕Ð͔+h͔Ø͕-h͕Ø͔)e02",
    ///     "+(-X͔b͕+X͕b͔+Y͔a͕-Y͕a͔+f͔Ð͕-f͕Ð͔+i͔Ø͕-i͕Ø͔)e03",
    ///     "+(+X͔d͕-X͕d͔+Y͔e͕-Y͕e͔+Z͔f͕-Z͕f͔-j͔Ø͕+j͕Ø͔)e40",
    ///     "+(+X͔g͕-X͕g͔+Y͔h͕-Y͕h͔+Z͔i͕-Z͕i͔+j͔Ð͕-j͕Ð͔)e05",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-h͔i͕+h͕i͔)e23",
    ///     "+(+a͔c͕-a͕c͔+d͔f͕-d͕f͔+g͔i͕-g͕i͔)e31",
    ///     "+(-a͔b͕+a͕b͔-d͔e͕+d͕e͔-g͔h͕+g͕h͔)e12",
    ///     "+(-b͔f͕+b͕f͔+c͔e͕-c͕e͔+g͔j͕-g͕j͔)e41",
    ///     "+(+a͔f͕-a͕f͔-c͔d͕+c͕d͔+h͔j͕-h͕j͔)e42",
    ///     "+(-a͔e͕+a͕e͔+b͔d͕-b͕d͔+i͔j͕-i͕j͔)e43",
    ///     "+(-b͔i͕+b͕i͔+c͔h͕-c͕h͔-d͔j͕+d͕j͔)e15",
    ///     "+(+a͔i͕-a͕i͔-c͔g͕+c͕g͔-e͔j͕+e͕j͔)e25",
    ///     "+(-a͔h͕+a͕h͔+b͔g͕-b͕g͔-f͔j͕+f͕j͔)e35",
    ///     "+(+d͔g͕-d͕g͔+e͔h͕-e͕h͔+f͔i͕-f͕i͔)e45",
    ///     "+(+a͔j͕+a͕j͔+e͔i͕+e͕i͔-f͔h͕-f͕h͔)e2345",
    ///     "+(+b͔j͕+b͕j͔-d͔i͕-d͕i͔+f͔g͕+f͕g͔)e3145",
    ///     "+(+c͔j͕+c͕j͔+d͔h͕+d͕h͔-e͔g͕-e͕g͔)e1245",
    ///     "+(+a͔g͕+a͕g͔+b͔h͕+b͕h͔+c͔i͕+c͕i͔)e1235",
    ///     "+(-a͔d͕-a͕d͔-b͔e͕-b͕e͔-c͔f͕-c͕f͔)e1234",
    ///     "+(+X͔j͕+X͕j͔-d͔Ø͕-d͕Ø͔+g͔Ð͕+g͕Ð͔)e0145",
    ///     "+(+Y͔j͕+Y͕j͔-e͔Ø͕-e͕Ø͔+h͔Ð͕+h͕Ð͔)e0245",
    ///     "+(+Z͔j͕+Z͕j͔-f͔Ø͕-f͕Ø͔+i͔Ð͕+i͕Ð͔)e0345",
    ///     "+(-Y͔i͕-Y͕i͔+Z͔h͕+Z͕h͔-a͔Ø͕-a͕Ø͔)e0325",
    ///     "+(+X͔i͕+X͕i͔-Z͔g͕-Z͕g͔-b͔Ø͕-b͕Ø͔)e0135",
    ///     "+(-X͔h͕-X͕h͔+Y͔g͕+Y͕g͔-c͔Ø͕-c͕Ø͔)e0215",
    ///     "+(+Y͔f͕+Y͕f͔-Z͔e͕-Z͕e͔+a͔Ð͕+a͕Ð͔)e0324",
    ///     "+(-X͔f͕-X͕f͔+Z͔d͕+Z͕d͔+b͔Ð͕+b͕Ð͔)e0134",
    ///     "+(+X͔e͕+X͕e͔-Y͔d͕-Y͕d͔+c͔Ð͕+c͕Ð͔)e0214",
    ///     "+(+X͔a͕+X͕a͔+Y͔b͕+Y͕b͔+Z͔c͕+Z͕c͔)e0123",
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
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕c͔-Z͕b͔+d͔Ð͕+g͔Ø͕)e01",
    ///     "+(-X͕c͔+Y͕v͔+Z͕a͔+e͔Ð͕+h͔Ø͕)e02",
    ///     "+(+X͕b͔-Y͕a͔+Z͕v͔+f͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕d͔-Y͕e͔-Z͕f͔-j͔Ø͕+v͔Ð͕)e40",
    ///     "+(-X͕g͔-Y͕h͔-Z͕i͔+j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e23",
    ///     "+b͔v͕e31",
    ///     "+c͔v͕e12",
    ///     "+d͔v͕e41",
    ///     "+e͔v͕e42",
    ///     "+f͔v͕e43",
    ///     "+g͔v͕e15",
    ///     "+h͔v͕e25",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+v͕x͔e2345",
    ///     "+v͕y͔e3145",
    ///     "+v͕z͔e1245",
    ///     "+v͕ð͔e1235",
    ///     "+v͕ø͔e1234",
    ///     "+(+X͕j͔+Y͕z͔-Z͕y͔-d͔Ø͕+g͔Ð͕)e0145",
    ///     "+(-X͕z͔+Y͕j͔+Z͕x͔-e͔Ø͕+h͔Ð͕)e0245",
    ///     "+(+X͕y͔-Y͕x͔+Z͕j͔-f͔Ø͕+i͔Ð͕)e0345",
    ///     "+(+X͕ð͔-Y͕i͔+Z͕h͔-a͔Ø͕-x͔Ð͕)e0325",
    ///     "+(+X͕i͔+Y͕ð͔-Z͕g͔-b͔Ø͕-y͔Ð͕)e0135",
    ///     "+(-X͕h͔+Y͕g͔+Z͕ð͔-c͔Ø͕-z͔Ð͕)e0215",
    ///     "+(+X͕ø͔+Y͕f͔-Z͕e͔+a͔Ð͕-x͔Ø͕)e0324",
    ///     "+(-X͕f͔+Y͕ø͔+Z͕d͔+b͔Ð͕-y͔Ø͕)e0134",
    ///     "+(+X͕e͔-Y͕d͔+Z͕ø͔+c͔Ð͕-z͔Ø͕)e0214",
    ///     "+(+X͕a͔+Y͕b͔+Z͕c͔-Ð͕ø͔+Ø͕ð͔)e0123",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔+Ø͕ø͔)I",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-B͔g͕+B͕g͔+C͔f͕-C͕f͔+D͔h͕-D͕h͔+F͔c͕-F͕c͔-G͔b͕+G͕b͔+H͔d͕-H͕d͔)e01",
    ///     "+(-A͔g͕+A͕g͔+C͔e͕-C͕e͔+D͔i͕-D͕i͔+E͔c͕-E͕c͔-G͔a͕+G͕a͔+I͔d͕-I͕d͔)e02",
    ///     "+(-A͔f͕+A͕f͔+B͔e͕-B͕e͔+D͔j͕-D͕j͔+E͔b͕-E͕b͔-F͔a͕+F͕a͔+J͔d͕-J͕d͔)e03",
    ///     "+(-A͔h͕+A͕h͔+B͔i͕-B͕i͔-C͔j͕+C͕j͔-H͔a͕+H͕a͔+I͔b͕-I͕b͔-J͔c͕+J͕c͔)e40",
    ///     "+(-E͔h͕+E͕h͔+F͔i͕-F͕i͔-G͔j͕+G͕j͔-H͔e͕+H͕e͔+I͔f͕-I͕f͔-J͔g͕+J͕g͔)e05",
    ///     "+(+b͔c͕-b͕c͔+f͔g͕-f͕g͔-i͔j͕+i͕j͔)e23",
    ///     "+(+a͔c͕-a͕c͔+e͔g͕-e͕g͔+h͔j͕-h͕j͔)e31",
    ///     "+(+a͔b͕-a͕b͔+e͔f͕-e͕f͔-h͔i͕+h͕i͔)e12",
    ///     "+(-a͔d͕+a͕d͔+f͔j͕-f͕j͔+g͔i͕-g͕i͔)e41",
    ///     "+(+b͔d͕-b͕d͔+e͔j͕-e͕j͔-g͔h͕+g͕h͔)e42",
    ///     "+(-c͔d͕+c͕d͔-e͔i͕+e͕i͔-f͔h͕+f͕h͔)e43",
    ///     "+(-b͔j͕+b͕j͔-c͔i͕+c͕i͔+d͔e͕-d͕e͔)e15",
    ///     "+(-a͔j͕+a͕j͔+c͔h͕-c͕h͔-d͔f͕+d͕f͔)e25",
    ///     "+(+a͔i͕-a͕i͔+b͔h͕-b͕h͔+d͔g͕-d͕g͔)e35",
    ///     "+(+a͔e͕-a͕e͔+b͔f͕-b͕f͔+c͔g͕-c͕g͔)e45",
    ///     "+(-b͔g͕-b͕g͔+c͔f͕+c͕f͔+d͔h͕+d͕h͔)e2345",
    ///     "+(-a͔g͕-a͕g͔+c͔e͕+c͕e͔+d͔i͕+d͕i͔)e3145",
    ///     "+(-a͔f͕-a͕f͔+b͔e͕+b͕e͔+d͔j͕+d͕j͔)e1245",
    ///     "+(-a͔h͕-a͕h͔+b͔i͕+b͕i͔-c͔j͕-c͕j͔)e1235",
    ///     "+(-e͔h͕-e͕h͔+f͔i͕+f͕i͔-g͔j͕-g͕j͔)e1234",
    ///     "+(-B͔c͕-B͕c͔+C͔b͕+C͕b͔-F͔g͕-F͕g͔+G͔f͕+G͕f͔+I͔j͕+I͕j͔-J͔i͕-J͕i͔)e0145",
    ///     "+(-A͔c͕-A͕c͔+C͔a͕+C͕a͔-E͔g͕-E͕g͔+G͔e͕+G͕e͔-H͔j͕-H͕j͔+J͔h͕+J͕h͔)e0245",
    ///     "+(-A͔b͕-A͕b͔+B͔a͕+B͕a͔-E͔f͕-E͕f͔+F͔e͕+F͕e͔+H͔i͕+H͕i͔-I͔h͕-I͕h͔)e0345",
    ///     "+(+A͔d͕+A͕d͔-D͔a͕-D͕a͔-F͔j͕-F͕j͔-G͔i͕-G͕i͔+I͔g͕+I͕g͔+J͔f͕+J͕f͔)e0325",
    ///     "+(-B͔d͕-B͕d͔+D͔b͕+D͕b͔-E͔j͕-E͕j͔+G͔h͕+G͕h͔-H͔g͕-H͕g͔+J͔e͕+J͕e͔)e0135",
    ///     "+(+C͔d͕+C͕d͔-D͔c͕-D͕c͔+E͔i͕+E͕i͔+F͔h͕+F͕h͔-H͔f͕-H͕f͔-I͔e͕-I͕e͔)e0215",
    ///     "+(+B͔j͕+B͕j͔+C͔i͕+C͕i͔-D͔e͕-D͕e͔+E͔d͕+E͕d͔-I͔c͕-I͕c͔-J͔b͕-J͕b͔)e0324",
    ///     "+(+A͔j͕+A͕j͔-C͔h͕-C͕h͔+D͔f͕+D͕f͔-F͔d͕-F͕d͔+H͔c͕+H͕c͔-J͔a͕-J͕a͔)e0134",
    ///     "+(-A͔i͕-A͕i͔-B͔h͕-B͕h͔-D͔g͕-D͕g͔+G͔d͕+G͕d͔+H͔b͕+H͕b͔+I͔a͕+I͕a͔)e0214",
    ///     "+(-A͔e͕-A͕e͔-B͔f͕-B͕f͔-C͔g͕-C͕g͔+E͔a͕+E͕a͔+F͔b͕+F͕b͔+G͔c͕+G͕c͔)e0123",
    ///     "+(-A͔a͕+A͕a͔-B͔b͕+B͕b͔-C͔c͕+C͕c͔-D͔d͕+D͕d͔-E͔e͕+E͕e͔-F͔f͕+F͕f͔-G͔g͕+G͕g͔-H͔h͕+H͕h͔-I͔i͕+I͕i͔-J͔j͕+J͕j͔)I",
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
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔)e345",
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
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔+x͔ø͕-x͕ø͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔+z͔ø͕-z͕ø͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔-ð͔ø͕-ð͕ø͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔-x͔ð͕-x͕ð͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔+y͔ð͕+y͕ð͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔-z͔ð͕-z͕ð͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔+x͔z͕-x͕z͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕-ð͔ð͕+ø͔ø͕)e12345",
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
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+v͕ø͔e5",
    ///     "+(+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-Y͕ø͔+y͔Ø͕)e052",
    ///     "+(+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ð͕ø͔+Ø͕ð͔)e054",
    ///     "+(+X͕ð͔+x͔Ð͕)e014",
    ///     "+(-Y͕ð͔-y͔Ð͕)e042",
    ///     "+(+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-Y͕z͔+Z͕y͔)e032",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(-X͕y͔+Y͕x͔)e021",
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
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+W͔g͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-W͔h͕-Y͕ø͔+y͔Ø͕)e052",
    ///     "+(+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-W͔j͕+Ð͕ø͔+Ø͕ð͔)e054",
    ///     "+(-W͔d͕+X͕ð͔+x͔Ð͕)e014",
    ///     "+(+W͔e͕-Y͕ð͔-y͔Ð͕)e042",
    ///     "+(-W͔f͕+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-W͔a͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-W͔c͕-X͕y͔+Y͕x͔)e021",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔)e345",
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
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+A͕ð͔+E͕z͔-F͕y͔+W͔g͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-B͕ð͔+D͕z͔-F͕x͔-W͔h͕-Y͕ø͔+y͔Ø͕)e052",
    ///     "+(+C͕ð͔+D͕y͔-E͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+A͕x͔+B͕y͔+C͕z͔-W͔j͕+Ð͕ø͔+Ø͕ð͔)e054",
    ///     "+(-A͕ø͔+H͕z͔-I͕y͔-W͔d͕+X͕ð͔+x͔Ð͕)e014",
    ///     "+(+B͕ø͔+G͕z͔-I͕x͔+W͔e͕-Y͕ð͔-y͔Ð͕)e042",
    ///     "+(-C͕ø͔+G͕y͔-H͕x͔-W͔f͕+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-D͕ø͔-G͕ð͔+J͕x͔-W͔a͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-E͕ø͔-H͕ð͔+J͕y͔-W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-F͕ø͔-I͕ð͔+J͕z͔-W͔c͕-X͕y͔+Y͕x͔)e021",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-B͕z͔+C͕y͔-D͕ð͔+G͕ø͔)e03245",
    ///     "+(+A͕z͔-C͕x͔-E͕ð͔+H͕ø͔)e01345",
    ///     "+(-A͕y͔+B͕x͔-F͕ð͔+I͕ø͔)e02145",
    ///     "+(+D͕x͔+E͕y͔+F͕z͔+J͕ø͔)e01235",
    ///     "+(-G͕x͔-H͕y͔-I͕z͔-J͕ð͔)e01243",
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
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+A͕ð͔+E͕z͔-F͕y͔+W͔g͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-B͕ð͔+D͕z͔-F͕x͔-W͔h͕-Y͕ø͔+y͔Ø͕)e052",
    ///     "+(+C͕ð͔+D͕y͔-E͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+A͕x͔+B͕y͔+C͕z͔-W͔j͕+Ð͕ø͔+Ø͕ð͔)e054",
    ///     "+(-A͕ø͔+H͕z͔-I͕y͔-W͔d͕+X͕ð͔+x͔Ð͕)e014",
    ///     "+(+B͕ø͔+G͕z͔-I͕x͔+W͔e͕-Y͕ð͔-y͔Ð͕)e042",
    ///     "+(-C͕ø͔+G͕y͔-H͕x͔-W͔f͕+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-D͕ø͔-G͕ð͔+J͕x͔-W͔a͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-E͕ø͔-H͕ð͔+J͕y͔-W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-F͕ø͔-I͕ð͔+J͕z͔-W͔c͕-X͕y͔+Y͕x͔)e021",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔+x͔ø͕-x͕ø͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔+z͔ø͕-z͕ø͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔-ð͔ø͕-ð͕ø͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔-x͔ð͕-x͕ð͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔+y͔ð͕+y͕ð͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔-z͔ð͕-z͕ð͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔+x͔z͕-x͕z͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕-ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(-B͕z͔+C͕y͔-D͕ð͔+G͕ø͔-W͔x͕)e03245",
    ///     "+(+A͕z͔-C͕x͔-E͕ð͔+H͕ø͔-W͔y͕)e01345",
    ///     "+(-A͕y͔+B͕x͔-F͕ð͔+I͕ø͔-W͔z͕)e02145",
    ///     "+(+D͕x͔+E͕y͔+F͕z͔+J͕ø͔+W͔ð͕)e01235",
    ///     "+(-G͕x͔-H͕y͔-I͕z͔-J͕ð͔-W͔ø͕)e01243",
    /// ]);
    ///
    /// let double_flector = Vee::volume4().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// format_eq!(double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔+Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(+b͕z͔-c͕y͔+d͕ð͔-g͕ø͔+v͕x͔)e1",
    ///     "+(-a͕z͔+c͕x͔+e͕ð͔-h͕ø͔+v͕y͔)e2",
    ///     "+(+a͕y͔-b͕x͔+f͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(-d͕x͔-e͕y͔-f͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+g͕x͔+h͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+A͕ð͔+E͕z͔-F͕y͔+W͔g͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-B͕ð͔+D͕z͔-F͕x͔-W͔h͕-Y͕ø͔+y͔Ø͕)e052",
    ///     "+(+C͕ð͔+D͕y͔-E͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+A͕x͔+B͕y͔+C͕z͔-W͔j͕+Ð͕ø͔+Ø͕ð͔)e054",
    ///     "+(-A͕ø͔+H͕z͔-I͕y͔-W͔d͕+X͕ð͔+x͔Ð͕)e014",
    ///     "+(+B͕ø͔+G͕z͔-I͕x͔+W͔e͕-Y͕ð͔-y͔Ð͕)e042",
    ///     "+(-C͕ø͔+G͕y͔-H͕x͔-W͔f͕+Z͕ð͔+z͔Ð͕)e034",
    ///     "+(-D͕ø͔-G͕ð͔+J͕x͔-W͔a͕-Y͕z͔+Z͕y͔)e032",
    ///     "+(-E͕ø͔-H͕ð͔+J͕y͔-W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-F͕ø͔-I͕ð͔+J͕z͔-W͔c͕-X͕y͔+Y͕x͔)e021",
    ///     "+(+a͕ð͔+e͕z͔-f͕y͔+x͔ø͕-x͕ø͔)e234",
    ///     "+(-b͕ð͔+d͕z͔-f͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+c͕ð͔+d͕y͔-e͕x͔+z͔ø͕-z͕ø͔)e124",
    ///     "+(+a͕x͔+b͕y͔+c͕z͔-ð͔ø͕-ð͕ø͔)e123",
    ///     "+(-a͕ø͔+h͕z͔-i͕y͔-x͔ð͕-x͕ð͔)e253",
    ///     "+(+b͕ø͔+g͕z͔-i͕x͔+y͔ð͕+y͕ð͔)e315",
    ///     "+(-c͕ø͔+g͕y͔-h͕x͔-z͔ð͕-z͕ð͔)e152",
    ///     "+(-d͕ø͔-g͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-e͕ø͔-h͕ð͔+j͕y͔+x͔z͕-x͕z͔)e245",
    ///     "+(-f͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕-ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(-B͕z͔+C͕y͔-D͕ð͔+G͕ø͔+V͕x͔-W͔x͕)e03245",
    ///     "+(+A͕z͔-C͕x͔-E͕ð͔+H͕ø͔+V͕y͔-W͔y͕)e01345",
    ///     "+(-A͕y͔+B͕x͔-F͕ð͔+I͕ø͔+V͕z͔-W͔z͕)e02145",
    ///     "+(+D͕x͔+E͕y͔+F͕z͔+J͕ø͔+V͕ð͔+W͔ð͕)e01235",
    ///     "+(-G͕x͔-H͕y͔-I͕z͔-J͕ð͔+V͕ø͔-W͔ø͕)e01243",
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
    /// The multivector of bias $`v^5_\infty \equiv W\e_0`$.
    #[must_use]
    #[inline]
    pub fn bias() -> Self {
        Self::e0()
    }
    /// The multivector of normal $`v^5_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4 + ø\e_5 + þ\e_6`$.
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
    /// The multivector of $`4`$-volume moment $`v^4_\infty`$.
    #[must_use]
    #[inline]
    pub fn volume4_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e40() + Self::e05() + Self::e60()
    }
    /// The multivector of $`4`$-volume displacement $`v^4_0`$.
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
    /// The multivector of volume moment $`v_\infty`$.
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
    /// The multivector of volume displacement $`v_0`$.
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
    /// The multivector of plane moment $`p_\infty`$.
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
    /// The multivector of plane displacement $`p_0`$.
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
    /// The multivector of line moment $`\ell_\infty`$.
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
    /// The multivector of line displacement $`\ell_0`$.
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
    /// The multivector of direction $`P_\infty`$.
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
    /// The multivector of weight $`P_0`$.
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
    /// The multivector of normal
    /// $`v^6_0 \equiv x\e_1 + y\e_2 + z\e_3 + ð\e_4 + ø\e_5 + þ\e_6 + œ\e_7`$.
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
    /// The multivector of $`5`$-volume moment $`v^5_\infty`$.
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
    /// The multivector of $`5`$-volume displacement $`v^5_0`$.
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
    /// The multivector of $`4`$-volume moment $`v^4_\infty`$.
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
    /// The multivector of $`4`$-volume displacement $`v^4_0`$.
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
    /// The multivector of volume moment $`v_\infty`$.
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
    /// The multivector of volume displacement $`v_0`$.
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
    /// The multivector of plane moment $`p_\infty`$.
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
    /// The multivector of plane displacement $`p_0`$.
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
    /// The multivector of line moment $`\ell_\infty`$.
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
    /// The multivector of line displacement $`\ell_0`$.
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
    /// The multivector of direction $`P_\infty`$.
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
    /// The multivector of weight $`P_0`$.
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
#[allow(clippy::cognitive_complexity)]
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
            assert_eq!(table, Ok(text));
        } else {
            write(&path, table.unwrap()).unwrap();
        }
    }
}
