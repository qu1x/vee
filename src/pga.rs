// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Plane-Based Pistachio Flavor -- Projective Geometric Algebra (PGA)

use crate::Symbol;

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
/// Basis blade of Elliptic 4D PGA.
pub type PgaE4 = Pga<1, 4>;
/// Basis blade of Elliptic 5D PGA (exploratory).
pub type PgaE5 = Pga<1, 5>;
/// Basis blade of Elliptic 6D PGA (exploratory).
pub type PgaE6 = Pga<1, 6>;
/// Basis blade of Elliptic 7D PGA (exploratory).
pub type PgaE7 = Pga<1, 7>;

/// Basis blade of Hyperbolic 0D PGA.
pub type PgaH0 = Pga<-1, 0>;
/// Basis blade of Hyperbolic 1D PGA.
pub type PgaH1 = Pga<-1, 1>;
/// Basis blade of Hyperbolic 2D PGA.
pub type PgaH2 = Pga<-1, 2>;
/// Basis blade of Hyperbolic 3D PGA.
pub type PgaH3 = Pga<-1, 3>;
/// Basis blade of Hyperbolic 4D PGA.
pub type PgaH4 = Pga<-1, 4>;
/// Basis blade of Hyperbolic 5D PGA (exploratory).
pub type PgaH5 = Pga<-1, 5>;
/// Basis blade of Hyperbolic 6D PGA (exploratory).
pub type PgaH6 = Pga<-1, 6>;
/// Basis blade of Hyperbolic 7D PGA (exploratory).
pub type PgaH7 = Pga<-1, 7>;

/// Basis blade of Parabolic (Euclidean) 0D PGA.
pub type PgaP0 = Pga<0, 0>;
/// Basis blade of Parabolic (Euclidean) 1D PGA.
pub type PgaP1 = Pga<0, 1>;
/// Basis blade of Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Pga<0, 2>;
/// Basis blade of Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Pga<0, 3>;
/// Basis blade of Parabolic (Euclidean) 4D PGA.
pub type PgaP4 = Pga<0, 4>;
/// Basis blade of Parabolic (Euclidean) 5D PGA (exploratory).
pub type PgaP5 = Pga<0, 5>;
/// Basis blade of Parabolic (Euclidean) 6D PGA (exploratory).
pub type PgaP6 = Pga<0, 6>;
/// Basis blade of Parabolic (Euclidean) 7D PGA (exploratory).
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
    idx: u8,
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
    #[inline]
    pub const fn new(sym: &str) -> Self {
        Self {
            idx: BasisBlade::new(sym).idx,
        }
    }
    /// Creates basis blade for pseudoscalar.
    #[must_use]
    #[inline]
    pub const fn pss() -> Self {
        Self {
            idx: u8::MAX >> (u8::BITS - (N + 1)),
        }
    }
    /// Counts swaps until `self * other` is ordered.
    #[must_use]
    pub fn cnt(self, other: Self) -> u32 {
        ((1..=N).fold(0, |p, n| p ^ (self.idx >> n)) & other.idx).count_ones()
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

impl<const M: i8, const N: u32> From<Pga<M, N>> for Symbol {
    #[inline]
    fn from(b: Pga<M, N>) -> Self {
        let s = Self::new(Self::VEC, LUT[N as usize][b.idx as usize].sym);
        if b.idx == Pga::<M, N>::pss().idx {
            s.alt()
        } else {
            s
        }
    }
}

impl<const M: i8, const N: u32> TryFrom<Symbol> for Pga<M, N> {
    type Error = Symbol;

    #[inline]
    fn try_from(s: Symbol) -> Result<Self, Symbol> {
        let b = Self::new(s.lab);
        if s.lab == LUT[N as usize][b.idx as usize].sym {
            Ok(b)
        } else {
            Err(s)
        }
    }
}

impl<const M: i8, const N: u32> Algebra for Pga<M, N> {
    const N: u32 = N;

    #[inline]
    fn scalar() -> Self {
        Self::default()
    }
    #[inline]
    fn pseudoscalar() -> Self {
        Self::pss()
    }
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
        let cnt = self.cnt(other)
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
        let not = Self {
            idx: !self.idx & Self::pss().idx,
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
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        match self.idx {
            0 if !fmt.alternate() => write!(fmt, "1"),
            idx if !fmt.alternate() && idx == Self::pss().idx => write!(fmt, "I"),
            _ => Display::fmt(LUT[N as usize][self.idx as usize].sym, fmt),
        }
    }
}

#[derive(Debug, Copy, Clone)]
struct BasisBlade<'a> {
    sym: &'a str,
    cnt: u8,
    idx: u8,
}

impl<'a> BasisBlade<'a> {
    const fn new(sym: &'a str) -> Self {
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
        Self { sym, cnt, idx }
    }
    const fn tab<const LEN: usize>(sym: [&'a str; LEN]) -> [Self; LEN] {
        let mut tab = [Self {
            sym: "",
            cnt: 0,
            idx: 0,
        }; LEN];
        let mut i = 0;
        while i < LEN {
            tab[i] = Self::new(sym[i]);
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
    ($t:ident, $u:ident, $n:tt, [$(($s:tt, $b:tt),)*]) => {
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
                    Self::new([
                        ((stringify!($s), stringify!($b)), const { Pga::new(stringify!($b)) }),
                    ])
                }
            )*
        }
        const $t: [BasisBlade; 1 << ($n + 1)] = BasisBlade::tab([$(stringify!($b),)*]);
        const $u: [BasisBlade; 1 << ($n + 1)] = BasisBlade::lut($t);
    };
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
    (ð, e4),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e04),
    (a, e12),
    (b, e13),
    (c, e14),
    (d, e23),
    (e, e24),
    (f, e34),
    (F, e012),
    (E, e031),
    (D, e014),
    (C, e023),
    (B, e042),
    (A, e034),
    (ð, e132),
    (z, e124),
    (y, e143),
    (x, e234),
    (Ð, e0123),
    (Z, e0142),
    (Y, e0134),
    (X, e0243),
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
    (ð, e4),
    (ø, e5),
    (X, e01),
    (Y, e02),
    (Z, e03),
    (Ð, e04),
    (Ø, e05),
    (a, e12),
    (b, e13),
    (c, e14),
    (d, e15),
    (e, e23),
    (f, e24),
    (g, e25),
    (h, e34),
    (i, e35),
    (j, e45),
    (A, e012),
    (B, e013),
    (C, e014),
    (D, e015),
    (E, e023),
    (F, e024),
    (G, e025),
    (H, e034),
    (I, e035),
    (J, e045),
    (j, e123),
    (i, e142),
    (h, e125),
    (g, e134),
    (f, e153),
    (e, e145),
    (d, e243),
    (c, e235),
    (b, e254),
    (a, e345),
    (J, e0123),
    (I, e0142),
    (H, e0125),
    (G, e0134),
    (F, e0153),
    (E, e0145),
    (D, e0243),
    (C, e0235),
    (B, e0254),
    (A, e0345),
    (ø, e1234),
    (ð, e1253),
    (z, e1245),
    (y, e1354),
    (x, e2345),
    (Ø, e01243),
    (Ð, e01235),
    (Z, e01254),
    (Y, e01345),
    (X, e02354),
    (w, e12345),
    (V, e012345),
]);
#[rustfmt::skip]
basis!(TAB6, LUT6, 6, [
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
    (Ð, e04),
    (Ø, e05),
    (Þ, e06),
    // 15
    (α, e12),
    (β, e13),
    (γ, e14),
    (δ, e15),
    (ε, e16),
    (ζ, e23),
    (η, e24),
    (θ, e25),
    (ι, e26),
    (κ, e34),
    (λ, e35),
    (μ, e36),
    (ν, e45),
    (ξ, e46),
    (ο, e56),
    // 15
    (Α, e012),
    (Β, e013),
    (Γ, e014),
    (Δ, e015),
    (Ε, e016),
    (Ζ, e023),
    (Η, e024),
    (Θ, e025),
    (Ι, e026),
    (Κ, e034),
    (Λ, e035),
    (Μ, e036),
    (Ν, e045),
    (Ξ, e046),
    (Ο, e056),
    // 20
    (a, e123),
    (b, e124),
    (c, e125),
    (d, e126),
    (e, e134),
    (f, e135),
    (g, e136),
    (h, e145),
    (i, e146),
    (j, e156),
    (k, e234),
    (l, e235),
    (m, e236),
    (n, e245),
    (o, e246),
    (p, e256),
    (q, e345),
    (r, e346),
    (s, e356),
    (t, e456),
    // 20
    (T, e0123),
    (S, e0142),
    (R, e0125),
    (Q, e0162),
    (P, e0134),
    (O, e0153),
    (N, e0136),
    (M, e0145),
    (L, e0164),
    (K, e0156),
    (J, e0243),
    (I, e0235),
    (H, e0263),
    (G, e0254),
    (F, e0246),
    (E, e0265),
    (D, e0345),
    (C, e0364),
    (B, e0356),
    (A, e0465),
    // 15
    (ο, e1234),
    (ξ, e1253),
    (ν, e1236),
    (μ, e1245),
    (λ, e1264),
    (κ, e1256),
    (ι, e1354),
    (θ, e1346),
    (η, e1365),
    (ζ, e1456),
    (ε, e2345),
    (δ, e2364),
    (γ, e2356),
    (β, e2465),
    (α, e3456),
    // 15
    (Ο, e01234),
    (Ξ, e01253),
    (Ν, e01236),
    (Μ, e01245),
    (Λ, e01264),
    (Κ, e01256),
    (Ι, e01354),
    (Θ, e01346),
    (Η, e01365),
    (Ζ, e01456),
    (Ε, e02345),
    (Δ, e02364),
    (Γ, e02356),
    (Β, e02465),
    (Α, e03456),
    // 6
    (þ, e12354),
    (ø, e12346),
    (ð, e12365),
    (z, e12456),
    (y, e13465),
    (x, e23456),
    // 6
    (Þ, e012345),
    (Ø, e012364),
    (Ð, e012356),
    (Z, e012465),
    (Y, e013456),
    (X, e023465),
    // 1
    (w, e123456),
    // 1
    (V, e0123456),
]);
#[rustfmt::skip]
basis!(TAB7, LUT7, 7, [
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
    (α, e12),
    (β, e13),
    (γ, e14),
    (δ, e15),
    (ε, e16),
    (ζ, e17),
    (η, e23),
    (θ, e24),
    (ι, e25),
    (κ, e26),
    (λ, e27),
    (μ, e34),
    (ν, e35),
    (ξ, e36),
    (ο, e37),
    (π, e45),
    (ρ, e46),
    (σ, e47),
    (τ, e56),
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
    (Τ, e012347),
    (Σ, e012356),
    (Ρ, e012375),
    (Π, e012367),
    (Ο, e012465),
    (Ξ, e012457),
    (Ν, e012476),
    (Μ, e012567),
    (Λ, e013456),
    (Κ, e013475),
    (Ι, e013467),
    (Θ, e013576),
    (Η, e014567),
    (Ζ, e023465),
    (Ε, e023457),
    (Δ, e023476),
    (Γ, e023567),
    (Β, e024576),
    (Α, e034567),
    // 7
    (œ, e123456),
    (þ, e123475),
    (ø, e123467),
    (ð, e123576),
    (z, e124567),
    (y, e134576),
    (x, e234567),
    // 7
    (Œ, e0123465),
    (Þ, e0123457),
    (Ø, e0123476),
    (Ð, e0123567),
    (Z, e0124576),
    (Y, e0134567),
    (X, e0234576),
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
    ///     "+(+vvW͓+2Vvw͓)e0",
    ///     "+vvw͓e1",
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
    ///     "+(+vv+ww)w͓e12",
    ///     "+(+(+vv-ww)X͓+2vwY͓+2(+Xw-Yv)w͓)e20",
    ///     "+(-2vwX͓+(+vv-ww)Y͓+2(+Xv+Yw)w͓)e01",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of rotoreflector $`f_r \equiv \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP2 as Vee};
    ///
    /// let rotoreflector = Vee::displacement().lhs() * Vee::rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// format_eq!(rotoreflector, [
    ///     "+(+v͕x͔-w͕y͔)e1",
    ///     "+(+v͕y͔+w͕x͔)e2",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(-xx-yy)w͓e12",
    ///     "+(+(+xx-yy)X͓+2xyY͓)e20",
    ///     "+(+2xyX͓+(-xx+yy)Y͓)e01",
    /// ]);
    /// ```
    #[must_use]
    pub fn rotoreflector() -> Self {
        Self::displacement()
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
    ///     "+(-xx-yy)w͓e12",
    ///     "+(+(+xx-yy)X͓+2xyY͓+2(+Vy+Wx)w͓)e20",
    ///     "+(+2xyX͓+(-xx+yy)Y͓+2(-Vx+Wy)w͓)e01",
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
    /// let norm_squared = Vee::line().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+xx+yy+zz",
    ///     "+2(-Xx-Yy-Zz)I",
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
    /// let norm_squared = motor.norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+x͔x͔x͕x͕+x͔x͔y͕y͕+x͔x͔z͕z͕+x͕x͕y͔y͔+x͕x͕z͔z͔+y͔y͔y͕y͕+y͔y͔z͕z͕+y͕y͕z͔z͔+z͔z͔z͕z͕",
    ///     "+2(-X͔x͔x͕x͕-X͔x͔y͕y͕-X͔x͔z͕z͕\
    ///         -X͕x͔x͔x͕-X͕x͕y͔y͔-X͕x͕z͔z͔\
    ///         -Y͔x͕x͕y͔-Y͔y͔y͕y͕-Y͔y͔z͕z͕\
    ///         -Y͕x͔x͔y͕-Y͕y͔y͔y͕-Y͕y͕z͔z͔\
    ///         -Z͔x͕x͕z͔-Z͔y͕y͕z͔-Z͔z͔z͕z͕\
    ///         -Z͕x͔x͔z͕-Z͕y͔y͔z͕-Z͕z͔z͔z͕)I",
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
    ///
    /// let norm_squared = motor.norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().vector(0).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+v͔v͔v͕v͕+v͕v͕x͔x͔+v͕v͕y͔y͔+v͕v͕z͔z͔",
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

/// The named entities of the PGA with embedded dimension $`N = 4`$.
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
    /// let quadvector_norm_squared = Vee::point().norm_squared();
    ///
    /// assert_eq!(quadvector_norm_squared.basis_blades(), Vee::scalar().basis_blades());
    /// format_eq!(quadvector_norm_squared, ["+ww"]);
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
    /// $`p_0 \equiv a\e_{12} + b\e_{13} + c\e_{14} + d\e_{23} + e\e_{24} + f\e_{34}`$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        (Self::e12() + Self::e13() + Self::e14() + Self::e23() + Self::e24() + Self::e34()).alt()
    }
    /// The multivector of plane moment
    /// $`p_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{04}`$.
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        (Self::e01() + Self::e02() + Self::e03() + Self::e04()).alt()
    }
    /// The multivector of plane $`p \equiv p_0 + p_\infty`$.
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::plane_displacement() + Self::plane_moment()
    }
    /// The multivector of line displacement
    /// $`\ell_0 \equiv x\e_{132} + y\e_{124} + z\e_{143} + ð\e_{234}`$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        (Self::e132() + Self::e124() + Self::e143() + Self::e234()).alt()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{012} + B\e_{031} + C\e_{014} + D\e_{023} + E\e_{042} + F\e_{034}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        (Self::e012() + Self::e031() + Self::e014() + Self::e023() + Self::e042() + Self::e034())
            .alt()
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
    /// $`P_\infty \equiv X\e_{0123} + Y\e_{0142} + Z\e_{0134} + Ð\e_{0243}`$.
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e0123() + Self::e0142() + Self::e0134() + Self::e0243()
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
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "-ẋ͔ẋ͕-ẏ͔ẏ͕-ż͔ż͕-ð͔̇ð͕̇",
    ///     "+(-ẋ͔ẏ͕+ẋ͕ẏ͔)e12",
    ///     "+(-ẋ͔ż͕+ẋ͕ż͔)e13",
    ///     "+(-ẋ͔ð͕̇+ẋ͕ð͔̇)e14",
    ///     "+(-ẏ͔ż͕+ẏ͕ż͔)e23",
    ///     "+(-ẏ͔ð͕̇+ẏ͕ð͔̇)e24",
    ///     "+(-ż͔ð͕̇+ż͕ð͔̇)e34",
    /// ]);
    ///
    /// let norm_squared = Vee::single_rotator().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv",
    ///     "+2(-ȧḟ+ḃė-ċḋ)e1234",
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
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕+v͔v͕",
    ///     "+(+ȧ͔v͕+ȧ͕v͔-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔+ḃ͔v͕+ḃ͕v͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔+ċ͔v͕+ċ͕v͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔+ḋ͔v͕+ḋ͕v͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔+ė͔v͕+ė͕v͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔+ḟ͔v͕+ḟ͕v͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕",
    ///     "+(-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    /// ]);
    ///
    /// let norm_squared = Vee::double_rotator().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv+ww",
    ///     "+2(-ȧḟ+ḃė-ċḋ+vw)e1234",
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
    ///     "+(+w͔Ð͕-w͕Ð͔)e04",
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
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    /// ]);
    ///
    /// let norm_squared = Vee::simple_motor().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     // Scalar condition.
    ///     "+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv",
    ///     // Point condition.
    ///     "+2(-ȧḟ+ḃė-ċḋ)e1234", // Weight condition.
    ///     "+2(+Ẏḟ-Żė+ḋÐ̇)e0243", // Direction condition.
    ///     "+2(-Ẋḟ+Żċ-ḃÐ̇)e0134", // Direction condition.
    ///     "+2(+Ẋė-Ẏċ+ȧÐ̇)e0142", // Direction condition.
    ///     "+2(-Ẋḋ+Ẏḃ-Żȧ)e0123", // Direction condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(-ȧḟ+ḃė-ċḋ)w͓", // Vanishes with weight condition.
    ///     "+(+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv)w͓e1234",
    ///     "+(+(-ȧȧ-ḃḃ-ċċ+ḋḋ+ėė+ḟḟ+vv)X͓+2(+ȧv-ḃḋ-ċė)Y͓+2(+ȧḋ+ḃv-ċḟ)Z͓+2(-Ẋv-Ẏȧ-Żḃ-ċÐ̇)w͓\
    ///        +2(+ȧė+ḃḟ+ċv)Ð͓)e0243",
    ///     "+(+2(-ȧv-ḃḋ-ċė)X͓+(-ȧȧ+ḃḃ+ċċ-ḋḋ-ėė+ḟḟ+vv)Y͓+2(-ȧḃ+ḋv-ėḟ)Z͓+2(+Ẋȧ-Ẏv-Żḋ-ėÐ̇)w͓\
    ///        +2(-ȧċ+ḋḟ+ėv)Ð͓)e0134",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)X͓+2(-ȧḃ-ḋv-ėḟ)Y͓+(+ȧȧ-ḃḃ+ċċ-ḋḋ+ėė-ḟḟ+vv)Z͓+2(+Ẋḃ+Ẏḋ-Żv-ḟÐ̇)w͓\
    ///        +2(-ḃċ-ḋė+ḟv)Ð͓)e0142",
    ///     "+(+2(+ȧė+ḃḟ-ċv)X͓+2(-ȧċ+ḋḟ-ėv)Y͓+2(-ḃċ-ḋė-ḟv)Z͓+2(+Ẋċ+Ẏė+Żḟ-vÐ̇)w͓\
    ///      +(+ȧȧ+ḃḃ-ċċ+ḋḋ-ėė-ḟḟ+vv)Ð͓)e0123",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor().unit();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+w͓e1234",
    ///     "+(+(+1-2ȧȧ-2ḃḃ-2ċċ)X͓+2(+ȧv-ḃḋ-ċė)Y͓+2(+ȧḋ+ḃv-ċḟ)Z͓+2(-Ẋv-Ẏȧ-Żḃ-ċÐ̇)w͓+2(+ȧė+ḃḟ+ċv)Ð͓)e0243",
    ///     "+(+2(-ȧv-ḃḋ-ċė)X͓+(+1-2ȧȧ-2ḋḋ-2ėė)Y͓+2(-ȧḃ+ḋv-ėḟ)Z͓+2(+Ẋȧ-Ẏv-Żḋ-ėÐ̇)w͓+2(-ȧċ+ḋḟ+ėv)Ð͓)e0134",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)X͓+2(-ȧḃ-ḋv-ėḟ)Y͓+(+1-2ḃḃ-2ḋḋ-2ḟḟ)Z͓+2(+Ẋḃ+Ẏḋ-Żv-ḟÐ̇)w͓+2(-ḃċ-ḋė+ḟv)Ð͓)e0142",
    ///     "+(+2(+ȧė+ḃḟ-ċv)X͓+2(-ȧċ+ḋḟ-ėv)Y͓+2(-ḃċ-ḋė-ḟv)Z͓+2(+Ẋċ+Ẏė+Żḟ-vÐ̇)w͓+(+1-2ċċ-2ėė-2ḟḟ)Ð͓)e0123",
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
    ///     "+(+(+ȧȧ+ḃḃ+ċċ+ḋḋ+ėė+ḟḟ+vv)W͓+2(+Ẋv-Ẏȧ-Żḃ-ċÐ̇)x͓+2(+Ẋȧ+Ẏv-Żḋ-ėÐ̇)y͓+2(+Ẋḃ+Ẏḋ+Żv-ḟÐ̇)z͓\
    ///      +2(+Ẋċ+Ẏė+Żḟ+vÐ̇)ð͓)e0",
    ///     "+(+(-ȧȧ-ḃḃ-ċċ+ḋḋ+ėė+ḟḟ+vv)x͓+2(+ȧv-ḃḋ-ċė)y͓+2(+ȧḋ+ḃv-ċḟ)z͓+2(+ȧė+ḃḟ+ċv)ð͓)e1",
    ///     "+(+2(-ȧv-ḃḋ-ċė)x͓+(-ȧȧ+ḃḃ+ċċ-ḋḋ-ėė+ḟḟ+vv)y͓+2(-ȧḃ+ḋv-ėḟ)z͓+2(-ȧċ+ḋḟ+ėv)ð͓)e2",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)x͓+2(-ȧḃ-ḋv-ėḟ)y͓+(+ȧȧ-ḃḃ+ċċ-ḋḋ+ėė-ḟḟ+vv)z͓+2(-ḃċ-ḋė+ḟv)ð͓)e3",
    ///     "+(+2(+ȧė+ḃḟ-ċv)x͓+2(-ȧċ+ḋḟ-ėv)y͓+2(-ḃċ-ḋė-ḟv)z͓+(+ȧȧ+ḃḃ-ċċ+ḋḋ-ėė-ḟḟ+vv)ð͓)e4",
    ///     // Vanishes with point condition.
    ///     "+2(+(-ȧḟ+ḃė-ċḋ)W͓+(+Ẏḟ-Żė+ḋÐ̇)x͓+(-Ẋḟ+Żċ-ḃÐ̇)y͓+(+Ẋė-Ẏċ+ȧÐ̇)z͓+(-Ẋḋ+Ẏḃ-Żȧ)ð͓)I",
    /// ]);
    ///
    /// let volume = Vee::volume().pin() << Vee::simple_motor().unit();
    ///
    /// assert_eq!(volume.basis_blades(), Vee::volume().basis_blades());
    /// format_eq!(volume, [
    ///     "+(+W͓+2(+Ẋv-Ẏȧ-Żḃ-ċÐ̇)x͓+2(+Ẋȧ+Ẏv-Żḋ-ėÐ̇)y͓+2(+Ẋḃ+Ẏḋ+Żv-ḟÐ̇)z͓+2(+Ẋċ+Ẏė+Żḟ+vÐ̇)ð͓)e0",
    ///     "+(+(+1-2ȧȧ-2ḃḃ-2ċċ)x͓+2(+ȧv-ḃḋ-ċė)y͓+2(+ȧḋ+ḃv-ċḟ)z͓+2(+ȧė+ḃḟ+ċv)ð͓)e1",
    ///     "+(+2(-ȧv-ḃḋ-ċė)x͓+(+1-2ȧȧ-2ḋḋ-2ėė)y͓+2(-ȧḃ+ḋv-ėḟ)z͓+2(-ȧċ+ḋḟ+ėv)ð͓)e2",
    ///     "+(+2(+ȧḋ-ḃv-ċḟ)x͓+2(-ȧḃ-ḋv-ėḟ)y͓+(+1-2ḃḃ-2ḋḋ-2ḟḟ)z͓+2(-ḃċ-ḋė+ḟv)ð͓)e3",
    ///     "+(+2(+ȧė+ḃḟ-ċv)x͓+2(-ȧċ+ḋḟ-ėv)y͓+2(-ḃċ-ḋė-ḟv)z͓+(+1-2ċċ-2ėė-2ḟḟ)ð͓)e4",
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
    ///     "+(+Ẋ͕v͔+Ẏ͕ȧ͔+Ż͕ḃ͔+ċ͔Ð͕̇)e01",
    ///     "+(-Ẋ͕ȧ͔+Ẏ͕v͔+Ż͕ḋ͔+ė͔Ð͕̇)e02",
    ///     "+(-Ẋ͕ḃ͔-Ẏ͕ḋ͔+Ż͕v͔+ḟ͔Ð͕̇)e03",
    ///     "+(-Ẋ͕ċ͔-Ẏ͕ė͔-Ż͕ḟ͔+v͔Ð͕̇)e04",
    ///     "+ȧ͔v͕e12",
    ///     "+ḃ͔v͕e13",
    ///     "+ċ͔v͕e14",
    ///     "+ḋ͔v͕e23",
    ///     "+ė͔v͕e24",
    ///     "+ḟ͔v͕e34",
    ///     "+(-Ẏ͕ḟ͔+Ż͕ė͔-ḋ͔Ð͕̇)e0243",
    ///     "+(+Ẋ͕ḟ͔-Ż͕ċ͔+ḃ͔Ð͕̇)e0134",
    ///     "+(-Ẋ͕ė͔+Ẏ͕ċ͔-ȧ͔Ð͕̇)e0142",
    ///     "+(+Ẋ͕ḋ͔-Ẏ͕ḃ͔+Ż͕ȧ͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "-ẋ͔ẋ͕-ẏ͔ẏ͕-ż͔ż͕-ð͔̇ð͕̇",
    ///     "+(+Ȧ͔ẏ͕-Ȧ͕ẏ͔+Ḃ͔ż͕-Ḃ͕ż͔+Ċ͔ð͕̇-Ċ͕ð͔̇)e01",
    ///     "+(-Ȧ͔ẋ͕+Ȧ͕ẋ͔+Ḋ͔ż͕-Ḋ͕ż͔+Ė͔ð͕̇-Ė͕ð͔̇)e02",
    ///     "+(-Ḃ͔ẋ͕+Ḃ͕ẋ͔-Ḋ͔ẏ͕+Ḋ͕ẏ͔+Ḟ͔ð͕̇-Ḟ͕ð͔̇)e03",
    ///     "+(-Ċ͔ẋ͕+Ċ͕ẋ͔-Ė͔ẏ͕+Ė͕ẏ͔-Ḟ͔ż͕+Ḟ͕ż͔)e04",
    ///     "+(-ẋ͔ẏ͕+ẋ͕ẏ͔)e12",
    ///     "+(-ẋ͔ż͕+ẋ͕ż͔)e13",
    ///     "+(-ẋ͔ð͕̇+ẋ͕ð͔̇)e14",
    ///     "+(-ẏ͔ż͕+ẏ͕ż͔)e23",
    ///     "+(-ẏ͔ð͕̇+ẏ͕ð͔̇)e24",
    ///     "+(-ż͔ð͕̇+ż͕ð͔̇)e34",
    ///     "+(-Ḋ͔ð͕̇-Ḋ͕ð͔̇+Ė͔ż͕+Ė͕ż͔-Ḟ͔ẏ͕-Ḟ͕ẏ͔)e0243",
    ///     "+(+Ḃ͔ð͕̇+Ḃ͕ð͔̇-Ċ͔ż͕-Ċ͕ż͔+Ḟ͔ẋ͕+Ḟ͕ẋ͔)e0134",
    ///     "+(-Ȧ͔ð͕̇-Ȧ͕ð͔̇+Ċ͔ẏ͕+Ċ͕ẏ͔-Ė͔ẋ͕-Ė͕ẋ͔)e0142",
    ///     "+(+Ȧ͔ż͕+Ȧ͕ż͔-Ḃ͔ẏ͕-Ḃ͕ẏ͔+Ḋ͔ẋ͕+Ḋ͕ẋ͔)e0123",
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
    ///     "+(+Ẋ͕v͔+Ẏ͕ȧ͔+Ż͕ḃ͔+ċ͔Ð͕̇)e01",
    ///     "+(-Ẋ͕ȧ͔+Ẏ͕v͔+Ż͕ḋ͔+ė͔Ð͕̇)e02",
    ///     "+(-Ẋ͕ḃ͔-Ẏ͕ḋ͔+Ż͕v͔+ḟ͔Ð͕̇)e03",
    ///     "+(-Ẋ͕ċ͔-Ẏ͕ė͔-Ż͕ḟ͔+v͔Ð͕̇)e04",
    ///     "+ȧ͔v͕e12",
    ///     "+ḃ͔v͕e13",
    ///     "+ċ͔v͕e14",
    ///     "+ḋ͔v͕e23",
    ///     "+ė͔v͕e24",
    ///     "+ḟ͔v͕e34",
    ///     "+v͕w͔e1234",
    ///     "+(+Ẋ͕w͔-Ẏ͕ḟ͔+Ż͕ė͔-ḋ͔Ð͕̇)e0243",
    ///     "+(+Ẋ͕ḟ͔+Ẏ͕w͔-Ż͕ċ͔+ḃ͔Ð͕̇)e0134",
    ///     "+(-Ẋ͕ė͔+Ẏ͕ċ͔+Ż͕w͔-ȧ͔Ð͕̇)e0142",
    ///     "+(+Ẋ͕ḋ͔-Ẏ͕ḃ͔+Ż͕ȧ͔+w͔Ð͕̇)e0123",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-ȧ͔ȧ͕-ḃ͔ḃ͕-ċ͔ċ͕-ḋ͔ḋ͕-ė͔ė͕-ḟ͔ḟ͕",
    ///     "+(-Ẏ͔ȧ͕+Ẏ͕ȧ͔-Ż͔ḃ͕+Ż͕ḃ͔+ċ͔Ð͕̇-ċ͕Ð͔̇)e01",
    ///     "+(+Ẋ͔ȧ͕-Ẋ͕ȧ͔-Ż͔ḋ͕+Ż͕ḋ͔+ė͔Ð͕̇-ė͕Ð͔̇)e02",
    ///     "+(+Ẋ͔ḃ͕-Ẋ͕ḃ͔+Ẏ͔ḋ͕-Ẏ͕ḋ͔+ḟ͔Ð͕̇-ḟ͕Ð͔̇)e03",
    ///     "+(+Ẋ͔ċ͕-Ẋ͕ċ͔+Ẏ͔ė͕-Ẏ͕ė͔+Ż͔ḟ͕-Ż͕ḟ͔)e04",
    ///     "+(-ḃ͔ḋ͕+ḃ͕ḋ͔-ċ͔ė͕+ċ͕ė͔)e12",
    ///     "+(+ȧ͔ḋ͕-ȧ͕ḋ͔-ċ͔ḟ͕+ċ͕ḟ͔)e13",
    ///     "+(+ȧ͔ė͕-ȧ͕ė͔+ḃ͔ḟ͕-ḃ͕ḟ͔)e14",
    ///     "+(-ȧ͔ḃ͕+ȧ͕ḃ͔-ė͔ḟ͕+ė͕ḟ͔)e23",
    ///     "+(-ȧ͔ċ͕+ȧ͕ċ͔+ḋ͔ḟ͕-ḋ͕ḟ͔)e24",
    ///     "+(-ḃ͔ċ͕+ḃ͕ċ͔-ḋ͔ė͕+ḋ͕ė͔)e34",
    ///     "+(+ȧ͔ḟ͕+ȧ͕ḟ͔-ḃ͔ė͕-ḃ͕ė͔+ċ͔ḋ͕+ċ͕ḋ͔)e1234",
    ///     "+(-Ẏ͔ḟ͕-Ẏ͕ḟ͔+Ż͔ė͕+Ż͕ė͔-ḋ͔Ð͕̇-ḋ͕Ð͔̇)e0243",
    ///     "+(+Ẋ͔ḟ͕+Ẋ͕ḟ͔-Ż͔ċ͕-Ż͕ċ͔+ḃ͔Ð͕̇+ḃ͕Ð͔̇)e0134",
    ///     "+(-Ẋ͔ė͕-Ẋ͕ė͔+Ẏ͔ċ͕+Ẏ͕ċ͔-ȧ͔Ð͕̇-ȧ͕Ð͔̇)e0142",
    ///     "+(+Ẋ͔ḋ͕+Ẋ͕ḋ͔-Ẏ͔ḃ͕-Ẏ͕ḃ͔+Ż͔ȧ͕+Ż͕ȧ͔)e0123",
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
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    /// ]);
    ///
    /// let norm_squared = Vee::rotoreflector().norm_squared();
    /// assert_eq!(norm_squared.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+xx+ẋẋ+yy+ẏẏ+zz+żż+ðð+ð̇ð̇",
    ///     "+2(-xẋ-yẏ-zż-ðð̇)e1234", // Weight condition.
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(+xẋ+yẏ+zż+ðð̇)w͓", // Vanishes with weight condition.
    ///     "+(-xx-ẋẋ-yy-ẏẏ-zz-żż-ðð-ð̇ð̇)w͓e1234",
    ///     "+(+(+xx-ẋẋ-yy+ẏẏ-zz+żż-ðð+ð̇ð̇)X͓+2(+xy-ẋẏ+zð̇-żð)Y͓+2(+xz-ẋż-yð̇+ẏð)Z͓\
    ///      +2(+xð-ẋð̇+yż-ẏz)Ð͓)e0243",
    ///     "+(+2(+xy-ẋẏ-zð̇+żð)X͓+(-xx+ẋẋ+yy-ẏẏ-zz+żż-ðð+ð̇ð̇)Y͓+2(+xð̇-ẋð+yz-ẏż)Z͓\
    ///      +2(-xż+ẋz+yð-ẏð̇)Ð͓)e0134",
    ///     "+(+2(+xz-ẋż+yð̇-ẏð)X͓+2(-xð̇+ẋð+yz-ẏż)Y͓+(-xx+ẋẋ-yy+ẏẏ+zz-żż-ðð+ð̇ð̇)Z͓\
    ///      +2(+xẏ-ẋy+zð-żð̇)Ð͓)e0142",
    ///     "+(+2(+xð-ẋð̇-yż+ẏz)X͓+2(+xż-ẋz+yð-ẏð̇)Y͓+2(-xẏ+ẋy+zð-żð̇)Z͓\
    ///      +(-xx+ẋẋ-yy+ẏẏ-zz+żż+ðð-ð̇ð̇)Ð͓)e0123",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector().unit();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// format_eq!(point, [
    ///     "+(+1-2xx-2ẋẋ-2yy-2ẏẏ-2zz-2żż-2ðð-2ð̇ð̇)w͓e1234",
    ///     "+(+(+1-2ẋẋ-2yy-2zz-2ðð)X͓+2(+xy-ẋẏ+zð̇-żð)Y͓+2(+xz-ẋż-yð̇+ẏð)Z͓+2(+xð-ẋð̇+yż-ẏz)Ð͓)e0243",
    ///     "+(+2(+xy-ẋẏ-zð̇+żð)X͓+(+1-2xx-2ẏẏ-2zz-2ðð)Y͓+2(+xð̇-ẋð+yz-ẏż)Z͓+2(-xż+ẋz+yð-ẏð̇)Ð͓)e0134",
    ///     "+(+2(+xz-ẋż+yð̇-ẏð)X͓+2(-xð̇+ẋð+yz-ẏż)Y͓+(+1-2xx-2yy-2żż-2ðð)Z͓+2(+xẏ-ẋy+zð-żð̇)Ð͓)e0142",
    ///     "+(+2(+xð-ẋð̇-yż+ẏz)X͓+2(+xż-ẋz+yð-ẏð̇)Y͓+2(-xẏ+ẋy+zð-żð̇)Z͓+(+1-2xx-2yy-2zz-2ð̇ð̇)Ð͓)e0123",
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
    ///     "+(-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+(+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-Ẏ͕ð͔+y͔Ð͕̇)e042",
    ///     "+(+Ẏ͕z͔-Ż͕y͔)e023",
    ///     "+(+Ẋ͕ð͔-x͔Ð͕̇)e014",
    ///     "+(-Ẋ͕z͔+Ż͕x͔)e031",
    ///     "+(+Ẋ͕y͔-Ẏ͕x͔)e012",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume() + Self::line_moment()
    }
    /// The multivector of flector $`f_s \equiv v + \ell`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP4 as Vee};
    ///
    /// let flector = Vee::volume().lhs() * Vee::simple_motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::simple_flector().basis_blades());
    /// format_eq!(flector, [
    ///     "+(+W͔v͕-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    ///     "+(+W͔ḟ͕+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-W͔ė͕-Ẏ͕ð͔+y͔Ð͕̇)e042",
    ///     "+(+W͔ḋ͕+Ẏ͕z͔-Ż͕y͔)e023",
    ///     "+(+W͔ċ͕+Ẋ͕ð͔-x͔Ð͕̇)e014",
    ///     "+(-W͔ḃ͕-Ẋ͕z͔+Ż͕x͔)e031",
    ///     "+(+W͔ȧ͕+Ẋ͕y͔-Ẏ͕x͔)e012",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_flector() -> Self {
        Self::volume() + Self::line()
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
    ///     "+(+W͔v͕-Ẋ͕x͔-Ẏ͕y͔-Ż͕z͔-Ð͕̇ð͔)e0",
    ///     "+(-ȧ͕y͔-ḃ͕z͔-ċ͕ð͔+v͕x͔)e1",
    ///     "+(+ȧ͕x͔-ḋ͕z͔-ė͕ð͔+v͕y͔)e2",
    ///     "+(+ḃ͕x͔+ḋ͕y͔-ḟ͕ð͔+v͕z͔)e3",
    ///     "+(+ċ͕x͔+ė͕y͔+ḟ͕z͔+v͕ð͔)e4",
    ///     "+(+ḋ͕ð͔-ė͕z͔+ḟ͕y͔)e234",
    ///     "+(-ḃ͕ð͔+ċ͕z͔-ḟ͕x͔)e143",
    ///     "+(+ȧ͕ð͔-ċ͕y͔+ė͕x͔)e124",
    ///     "+(-ȧ͕z͔+ḃ͕y͔-ḋ͕x͔)e132",
    ///     "+(+W͔ḟ͕+X͕y͔-Y͕x͔+Ż͕ð͔-z͔Ð͕̇)e034",
    ///     "+(-W͔ė͕+X͕z͔-Ẏ͕ð͔-Z͕x͔+y͔Ð͕̇)e042",
    ///     "+(+W͔ḋ͕+X͕ð͔+Ẏ͕z͔-Ż͕y͔-x͔Ð͕)e023",
    ///     "+(+W͔ċ͕+Ẋ͕ð͔+Y͕z͔-Z͕y͔-x͔Ð͕̇)e014",
    ///     "+(-W͔ḃ͕-Ẋ͕z͔+Y͕ð͔+Ż͕x͔-y͔Ð͕)e031",
    ///     "+(+W͔ȧ͕+Ẋ͕y͔-Ẏ͕x͔+Z͕ð͔-z͔Ð͕)e012",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::volume() + Self::line() + Self::pseudoscalar()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 5`$ (exploratory).
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
    /// let quadvector_norm_squared = Vee::line().norm_squared();
    ///
    /// assert_eq!(quadvector_norm_squared.basis_blades(),
    ///     (Vee::scalar() + Vee::line_moment()).basis_blades());
    /// format_eq!(quadvector_norm_squared, [
    ///     "+xx+yy+zz+ðð+øø",
    ///     "+2(-Hø+Ið-Jz)e0345",
    ///     "+2(+Fø-Gð+Jy)e0254",
    ///     "+2(-Eø+Gz-Iy)e0235",
    ///     "+2(+Eð-Fz+Hy)e0243",
    ///     "+2(-Cø+Dð-Jx)e0145",
    ///     "+2(+Bø-Dz+Ix)e0153",
    ///     "+2(-Bð+Cz-Hx)e0134",
    ///     "+2(-Aø+Dy-Gx)e0125",
    ///     "+2(+Að-Cy+Fx)e0142",
    ///     "+2(-Az+By-Ex)e0123",
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
    /// v_0 \equiv a\e_{12} + b\e_{13} + c\e_{14} + d\e_{15} + e\e_{23}
    ///          + f\e_{24} + g\e_{25} + h\e_{34} + i\e_{35} + j\e_{45}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e12()
            + Self::e13()
            + Self::e14()
            + Self::e15()
            + Self::e23()
            + Self::e24()
            + Self::e25()
            + Self::e34()
            + Self::e35()
            + Self::e45()
    }
    /// The multivector of volume moment.
    ///
    /// ```math
    /// v_\infty \equiv X\e_{01} + Y\e_{02} + Z\e_{03} + Ð\e_{04} + Ø\e_{05}
    /// ```
    #[must_use]
    #[inline]
    pub fn volume_moment() -> Self {
        Self::e01() + Self::e02() + Self::e03() + Self::e04() + Self::e05()
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
    /// p_0 \equiv a\e_{123} + b\e_{142} + c\e_{125} + d\e_{134} + e\e_{153}
    ///          + f\e_{145} + g\e_{243} + h\e_{235} + i\e_{254} + j\e_{345}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        Self::e123()
            + Self::e142()
            + Self::e125()
            + Self::e134()
            + Self::e153()
            + Self::e145()
            + Self::e243()
            + Self::e235()
            + Self::e254()
            + Self::e345()
    }
    /// The multivector of plane moment.
    ///
    /// ```math
    /// p_\infty \equiv A\e_{012} + B\e_{013} + C\e_{014} + D\e_{015} + E\e_{023}
    ///               + F\e_{024} + G\e_{025} + H\e_{034} + I\e_{035} + J\e_{045}
    /// ```
    #[must_use]
    #[inline]
    pub fn plane_moment() -> Self {
        Self::e012()
            + Self::e013()
            + Self::e014()
            + Self::e015()
            + Self::e023()
            + Self::e024()
            + Self::e025()
            + Self::e034()
            + Self::e035()
            + Self::e045()
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
    /// \ell_0 \equiv x\e_{1234} + y\e_{1253} + z\e_{1245} + ð\e_{1354} + ø\e_{2345}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        Self::e1234() + Self::e1253() + Self::e1245() + Self::e1354() + Self::e2345()
    }
    /// The multivector of line moment.
    ///
    /// ```math
    /// \ell_\infty \equiv A\e_{0123} + B\e_{0142} + C\e_{0125} + D\e_{0134} + E\e_{0153}
    ///                  + F\e_{0145} + G\e_{0243} + H\e_{0235} + I\e_{0254} + J\e_{0345}
    /// ```
    #[must_use]
    #[inline]
    pub fn line_moment() -> Self {
        Self::e0123()
            + Self::e0142()
            + Self::e0125()
            + Self::e0134()
            + Self::e0153()
            + Self::e0145()
            + Self::e0243()
            + Self::e0235()
            + Self::e0254()
            + Self::e0345()
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
    /// P_\infty \equiv X\e_{01243} + Y\e_{01235} + Z\e_{01254} + Ð\e_{01345} + Ø\e_{02354}
    /// ```
    #[must_use]
    #[inline]
    pub fn direction() -> Self {
        Self::e01243() + Self::e01235() + Self::e01254() + Self::e01345() + Self::e02354()
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
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
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
    ///     "+(+a͔v͕+a͕v͔-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔+b͔v͕+b͕v͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔+c͔v͕+c͕v͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔+d͔v͕+d͕v͔)e15",
    ///     "+(-a͔b͕+a͕b͔+e͔v͕+e͕v͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔+f͔v͕+f͕v͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔+g͔v͕+g͕v͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔+h͔v͕+h͕v͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔+i͔v͕+i͕v͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔+j͔v͕+j͕v͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
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
    ///     "+(+w͔Ð͕-w͕Ð͔)e04",
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
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+W͔ø͕-W͕ø͔)e05",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
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
    ///     "+(+X͕v͔+Y͕a͔+Z͕b͔+c͔Ð͕+d͔Ø͕)e01",
    ///     "+(-X͕a͔+Y͕v͔+Z͕e͔+f͔Ð͕+g͔Ø͕)e02",
    ///     "+(-X͕b͔-Y͕e͔+Z͕v͔+h͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕c͔-Y͕f͔-Z͕h͔+j͔Ø͕+v͔Ð͕)e04",
    ///     "+(-X͕d͔-Y͕g͔-Z͕i͔-j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e12",
    ///     "+b͔v͕e13",
    ///     "+c͔v͕e14",
    ///     "+d͔v͕e15",
    ///     "+e͔v͕e23",
    ///     "+f͔v͕e24",
    ///     "+g͔v͕e25",
    ///     "+h͔v͕e34",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+(+Z͕j͔+h͔Ø͕-i͔Ð͕)e0345",
    ///     "+(-Y͕j͔-f͔Ø͕+g͔Ð͕)e0254",
    ///     "+(+Y͕i͔-Z͕g͔+e͔Ø͕)e0235",
    ///     "+(-Y͕h͔+Z͕f͔-e͔Ð͕)e0243",
    ///     "+(+X͕j͔+c͔Ø͕-d͔Ð͕)e0145",
    ///     "+(-X͕i͔+Z͕d͔-b͔Ø͕)e0153",
    ///     "+(+X͕h͔-Z͕c͔+b͔Ð͕)e0134",
    ///     "+(+X͕g͔-Y͕d͔+a͔Ø͕)e0125",
    ///     "+(-X͕f͔+Y͕c͔-a͔Ð͕)e0142",
    ///     "+(+X͕e͔-Y͕b͔+Z͕a͔)e0123",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕",
    ///     "+(-A͔y͕+A͕y͔-B͔z͕+B͕z͔-C͔ð͕+C͕ð͔-D͔ø͕+D͕ø͔)e01",
    ///     "+(+A͔x͕-A͕x͔-E͔z͕+E͕z͔-F͔ð͕+F͕ð͔-G͔ø͕+G͕ø͔)e02",
    ///     "+(+B͔x͕-B͕x͔+E͔y͕-E͕y͔-H͔ð͕+H͕ð͔-I͔ø͕+I͕ø͔)e03",
    ///     "+(+C͔x͕-C͕x͔+F͔y͕-F͕y͔+H͔z͕-H͕z͔-J͔ø͕+J͕ø͔)e04",
    ///     "+(+D͔x͕-D͕x͔+G͔y͕-G͕y͔+I͔z͕-I͕z͔+J͔ð͕-J͕ð͔)e05",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(-H͔ø͕-H͕ø͔+I͔ð͕+I͕ð͔-J͔z͕-J͕z͔)e0345",
    ///     "+(+F͔ø͕+F͕ø͔-G͔ð͕-G͕ð͔+J͔y͕+J͕y͔)e0254",
    ///     "+(-E͔ø͕-E͕ø͔+G͔z͕+G͕z͔-I͔y͕-I͕y͔)e0235",
    ///     "+(+E͔ð͕+E͕ð͔-F͔z͕-F͕z͔+H͔y͕+H͕y͔)e0243",
    ///     "+(-C͔ø͕-C͕ø͔+D͔ð͕+D͕ð͔-J͔x͕-J͕x͔)e0145",
    ///     "+(+B͔ø͕+B͕ø͔-D͔z͕-D͕z͔+I͔x͕+I͕x͔)e0153",
    ///     "+(-B͔ð͕-B͕ð͔+C͔z͕+C͕z͔-H͔x͕-H͕x͔)e0134",
    ///     "+(-A͔ø͕-A͕ø͔+D͔y͕+D͕y͔-G͔x͕-G͕x͔)e0125",
    ///     "+(+A͔ð͕+A͕ð͔-C͔y͕-C͕y͔+F͔x͕+F͕x͔)e0142",
    ///     "+(-A͔z͕-A͕z͔+B͔y͕+B͕y͔-E͔x͕-E͕x͔)e0123",
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
    ///     "+(-Y͔a͕+Y͕a͔-Z͔b͕+Z͕b͔+c͔Ð͕-c͕Ð͔+d͔Ø͕-d͕Ø͔)e01",
    ///     "+(+X͔a͕-X͕a͔-Z͔e͕+Z͕e͔+f͔Ð͕-f͕Ð͔+g͔Ø͕-g͕Ø͔)e02",
    ///     "+(+X͔b͕-X͕b͔+Y͔e͕-Y͕e͔+h͔Ð͕-h͕Ð͔+i͔Ø͕-i͕Ø͔)e03",
    ///     "+(+X͔c͕-X͕c͔+Y͔f͕-Y͕f͔+Z͔h͕-Z͕h͔+j͔Ø͕-j͕Ø͔)e04",
    ///     "+(+X͔d͕-X͕d͔+Y͔g͕-Y͕g͔+Z͔i͕-Z͕i͔-j͔Ð͕+j͕Ð͔)e05",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    ///     "+(+Z͔j͕+Z͕j͔+h͔Ø͕+h͕Ø͔-i͔Ð͕-i͕Ð͔)e0345",
    ///     "+(-Y͔j͕-Y͕j͔-f͔Ø͕-f͕Ø͔+g͔Ð͕+g͕Ð͔)e0254",
    ///     "+(+Y͔i͕+Y͕i͔-Z͔g͕-Z͕g͔+e͔Ø͕+e͕Ø͔)e0235",
    ///     "+(-Y͔h͕-Y͕h͔+Z͔f͕+Z͕f͔-e͔Ð͕-e͕Ð͔)e0243",
    ///     "+(+X͔j͕+X͕j͔+c͔Ø͕+c͕Ø͔-d͔Ð͕-d͕Ð͔)e0145",
    ///     "+(-X͔i͕-X͕i͔+Z͔d͕+Z͕d͔-b͔Ø͕-b͕Ø͔)e0153",
    ///     "+(+X͔h͕+X͕h͔-Z͔c͕-Z͕c͔+b͔Ð͕+b͕Ð͔)e0134",
    ///     "+(+X͔g͕+X͕g͔-Y͔d͕-Y͕d͔+a͔Ø͕+a͕Ø͔)e0125",
    ///     "+(-X͔f͕-X͕f͔+Y͔c͕+Y͕c͔-a͔Ð͕-a͕Ð͔)e0142",
    ///     "+(+X͔e͕+X͕e͔-Y͔b͕-Y͕b͔+Z͔a͕+Z͕a͔)e0123",
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
    ///     "+(+X͕v͔+Y͕a͔+Z͕b͔+c͔Ð͕+d͔Ø͕)e01",
    ///     "+(-X͕a͔+Y͕v͔+Z͕e͔+f͔Ð͕+g͔Ø͕)e02",
    ///     "+(-X͕b͔-Y͕e͔+Z͕v͔+h͔Ð͕+i͔Ø͕)e03",
    ///     "+(-X͕c͔-Y͕f͔-Z͕h͔+j͔Ø͕+v͔Ð͕)e04",
    ///     "+(-X͕d͔-Y͕g͔-Z͕i͔-j͔Ð͕+v͔Ø͕)e05",
    ///     "+a͔v͕e12",
    ///     "+b͔v͕e13",
    ///     "+c͔v͕e14",
    ///     "+d͔v͕e15",
    ///     "+e͔v͕e23",
    ///     "+f͔v͕e24",
    ///     "+g͔v͕e25",
    ///     "+h͔v͕e34",
    ///     "+i͔v͕e35",
    ///     "+j͔v͕e45",
    ///     "+v͕x͔e2345",
    ///     "+v͕y͔e1354",
    ///     "+v͕z͔e1245",
    ///     "+v͕ð͔e1253",
    ///     "+v͕ø͔e1234",
    ///     "+(+X͕y͔-Y͕x͔+Z͕j͔+h͔Ø͕-i͔Ð͕)e0345",
    ///     "+(+X͕z͔-Y͕j͔-Z͕x͔-f͔Ø͕+g͔Ð͕)e0254",
    ///     "+(+X͕ð͔+Y͕i͔-Z͕g͔+e͔Ø͕-x͔Ð͕)e0235",
    ///     "+(+X͕ø͔-Y͕h͔+Z͕f͔-e͔Ð͕-x͔Ø͕)e0243",
    ///     "+(+X͕j͔+Y͕z͔-Z͕y͔+c͔Ø͕-d͔Ð͕)e0145",
    ///     "+(-X͕i͔+Y͕ð͔+Z͕d͔-b͔Ø͕-y͔Ð͕)e0153",
    ///     "+(+X͕h͔+Y͕ø͔-Z͕c͔+b͔Ð͕-y͔Ø͕)e0134",
    ///     "+(+X͕g͔-Y͕d͔+Z͕ð͔+a͔Ø͕-z͔Ð͕)e0125",
    ///     "+(-X͕f͔+Y͕c͔+Z͕ø͔-a͔Ð͕-z͔Ø͕)e0142",
    ///     "+(+X͕e͔-Y͕b͔+Z͕a͔+Ð͕ø͔-Ø͕ð͔)e0123",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔+Ø͕ø͔)I",
    /// ]);
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕",
    ///     "+(-E͔j͕+E͕j͔+F͔i͕-F͕i͔-G͔h͕+G͕h͔-H͔g͕+H͕g͔+I͔f͕-I͕f͔-J͔e͕+J͕e͔)e01",
    ///     "+(+B͔j͕-B͕j͔-C͔i͕+C͕i͔+D͔h͕-D͕h͔+H͔d͕-H͕d͔-I͔c͕+I͕c͔+J͔b͕-J͕b͔)e02",
    ///     "+(-A͔j͕+A͕j͔+C͔g͕-C͕g͔-D͔f͕+D͕f͔-F͔d͕+F͕d͔+G͔c͕-G͕c͔-J͔a͕+J͕a͔)e03",
    ///     "+(+A͔i͕-A͕i͔-B͔g͕+B͕g͔+D͔e͕-D͕e͔+E͔d͕-E͕d͔-G͔b͕+G͕b͔+I͔a͕-I͕a͔)e04",
    ///     "+(-A͔h͕+A͕h͔+B͔f͕-B͕f͔-C͔e͕+C͕e͔-E͔c͕+E͕c͔+F͔b͕-F͕b͔-H͔a͕+H͕a͔)e05",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔)e12",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔)e13",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔)e14",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔)e15",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔)e23",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔)e24",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔)e25",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔)e34",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔)e35",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔)e45",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔)e2345",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔)e1354",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔)e1245",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔)e1253",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔)e1234",
    ///     "+(-B͔e͕-B͕e͔-C͔f͕-C͕f͔-D͔g͕-D͕g͔+E͔b͕+E͕b͔+F͔c͕+F͕c͔+G͔d͕+G͕d͔)e0345",
    ///     "+(+A͔e͕+A͕e͔-C͔h͕-C͕h͔-D͔i͕-D͕i͔-E͔a͕-E͕a͔+H͔c͕+H͕c͔+I͔d͕+I͕d͔)e0254",
    ///     "+(+A͔f͕+A͕f͔+B͔h͕+B͕h͔-D͔j͕-D͕j͔-F͔a͕-F͕a͔-H͔b͕-H͕b͔+J͔d͕+J͕d͔)e0235",
    ///     "+(+A͔g͕+A͕g͔+B͔i͕+B͕i͔+C͔j͕+C͕j͔-G͔a͕-G͕a͔-I͔b͕-I͕b͔-J͔c͕-J͕c͔)e0243",
    ///     "+(-A͔b͕-A͕b͔+B͔a͕+B͕a͔-F͔h͕-F͕h͔-G͔i͕-G͕i͔+H͔f͕+H͕f͔+I͔g͕+I͕g͔)e0145",
    ///     "+(-A͔c͕-A͕c͔+C͔a͕+C͕a͔+E͔h͕+E͕h͔-G͔j͕-G͕j͔-H͔e͕-H͕e͔+J͔g͕+J͕g͔)e0153",
    ///     "+(-A͔d͕-A͕d͔+D͔a͕+D͕a͔+E͔i͕+E͕i͔+F͔j͕+F͕j͔-I͔e͕-I͕e͔-J͔f͕-J͕f͔)e0134",
    ///     "+(-B͔c͕-B͕c͔+C͔b͕+C͕b͔-E͔f͕-E͕f͔+F͔e͕+F͕e͔-I͔j͕-I͕j͔+J͔i͕+J͕i͔)e0125",
    ///     "+(-B͔d͕-B͕d͔+D͔b͕+D͕b͔-E͔g͕-E͕g͔+G͔e͕+G͕e͔+H͔j͕+H͕j͔-J͔h͕-J͕h͔)e0142",
    ///     "+(-C͔d͕-C͕d͔+D͔c͕+D͕c͔-F͔g͕-F͕g͔+G͔f͕+G͕f͔-H͔i͕-H͕i͔+I͔h͕+I͕h͔)e0123",
    ///     "+(+A͔a͕-A͕a͔+B͔b͕-B͕b͔+C͔c͕-C͕c͔+D͔d͕-D͕d͔+E͔e͕-E͕e͔+F͔f͕-F͕f͔+G͔g͕-G͕g͔+H͔h͕-H͕h͔+I͔i͕-I͕i͔+J͔j͕-J͕j͔)I",
    /// ]);
    ///
    /// let point = Vee::point().pin() << Vee::double_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::point() + Vee::volume4()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(+(+ay+bz+cð+dø+ej-fi+gh-vx)X͓+(-ax-bj+ci-dh+ez+fð+gø-vy)Y͓\
    ///         +(+aj-bx-cg+df-ey+hð+iø-vz)Z͓+(-Aa-Bb-Cc-Dd-Ee-Ff-Gg-Hh-Ii-Jj+Vv+Xx+Yy+Zz+Ðð+Øø)w͓\
    ///         +(-ai+bg-cx-de-fy-hz+jø-vð)Ð͓+(+ah-bf+ce-dx-gy-iz-jð-vø)Ø͓)e0",
    ///     "+2(+ay+bz+cð+dø-ej+fi-gh+vx)w͓e1",
    ///     "+2(-ax+bj-ci+dh+ez+fð+gø+vy)w͓e2",
    ///     "+2(-aj-bx+cg-df-ey+hð+iø+vz)w͓e3",
    ///     "+2(+ai-bg-cx+de-fy-hz+jø+vð)w͓e4",
    ///     "+2(-ah+bf-ce-dx-gy-iz-jð+vø)w͓e5",
    ///     "+(+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+vv+xx+yy+zz+ðð+øø)w͓e12345",
    ///     "+(+(-aa-bb-cc-dd+ee+ff+gg+hh+ii+jj+vv+xx-yy-zz-ðð-øø)X͓+2(+av-be-cf-dg+hø-ið+jz+xy)Y͓\
    ///        +2(+ae+bv-ch-di-fø+gð-jy+xz)Z͓+2(+Ay+Bz+Cð+Dø-Ej+Fi-Gh-Hg+If-Je-Vx-Xv-Ya-Zb-cÐ-dØ)w͓\
    ///        +2(+af+bh+cv-dj+eø-gz+iy+xð)Ð͓+2(+ag+bi+cj+dv-eð+fz-hy+xø)Ø͓)e02354",
    ///     "+(+2(-av-be-cf-dg-hø+ið-jz+xy)X͓+(-aa+bb+cc+dd-ee-ff-gg+hh+ii+jj+vv-xx+yy-zz-ðð-øø)Y͓\
    ///        +2(-ab+cø-dð+ev-fh-gi+jx+yz)Z͓+2(-Ax+Bj-Ci+Dh+Ez+Fð+Gø+Hd-Ic+Jb-Vy+Xa-Yv-Ze-fÐ-gØ)w͓\
    ///        +2(-ac-bø+dz+eh+fv-gj-ix+yð)Ð͓+2(-ad+bð-cz+ei+fj+gv+hx+yø)Ø͓)e01345",
    ///     "+(+2(+ae-bv-ch-di+fø-gð+jy+xz)X͓+2(-ab-cø+dð-ev-fh-gi-jx+yz)Y͓\
    ///        +(+aa-bb+cc+dd-ee+ff+gg-hh-ii+jj+vv-xx-yy+zz-ðð-øø)Z͓\
    ///        +2(-Aj-Bx+Cg-Df-Ey-Fd+Gc+Hð+Iø-Ja-Vz+Xb+Ye-Zv-hÐ-iØ)w͓+2(+aø-bc-dy-ef+gx+hv-ij+zð)Ð͓\
    ///        +2(-að-bd+cy-eg-fx+hj+iv+zø)Ø͓)e01254",
    ///     "+(+2(+af+bh-cv-dj-eø+gz-iy+xð)X͓+2(-ac+bø-dz+eh-fv-gj+ix+yð)Y͓\
    ///        +2(-aø-bc+dy-ef-gx-hv-ij+zð)Z͓+2(+Ai-Bg-Cx+De+Ed-Fy-Gb-Hz+Ia+Jø-Vð+Xc+Yf+Zh-jØ-vÐ)w͓\
    ///        +(+aa+bb-cc+dd+ee-ff+gg-hh+ii-jj+vv-xx-yy-zz+ðð-øø)Ð͓\
    ///        +2(+az-by-cd+ex-fg-hi+jv+ðø)Ø͓)e01235",
    ///     "+(+2(+ag+bi+cj-dv+eð-fz+hy+xø)X͓+2(-ad-bð+cz+ei+fj-gv-hx+yø)Y͓\
    ///        +2(+að-bd-cy-eg+fx+hj-iv+zø)Z͓+2(-Ah+Bf-Ce-Dx-Ec+Fb-Gy-Ha-Iz-Jð-Vø+Xd+Yg+Zi+jÐ-vØ)w͓\
    ///        +2(-az+by-cd-ex-fg-hi-jv+ðø)Ð͓\
    ///        +(+aa+bb+cc-dd+ee+ff-gg+hh-ii-jj+vv-xx-yy-zz-ðð+øø)Ø͓)e01243",
    /// ]);
    ///
    /// let norm_squared = Vee::double_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+vv+xx+yy+zz+ðð+øø",
    ///     "+2(+ay+bz+cð+dø-ej+fi-gh+vx)e2345",
    ///     "+2(-ax+bj-ci+dh+ez+fð+gø+vy)e1354",
    ///     "+2(-aj-bx+cg-df-ey+hð+iø+vz)e1245",
    ///     "+2(+ai-bg-cx+de-fy-hz+jø+vð)e1253",
    ///     "+2(-ah+bf-ce-dx-gy-iz-jð+vø)e1234",
    ///     "+2(+Av+Be+Cf+Dg-Eb-Fc-Gd-Hø+Ið-Jz+Va-Xy+Yx-Zj-hØ+iÐ)e0345",
    ///     "+2(-Ae+Bv+Ch+Di+Ea+Fø-Gð-Hc-Id+Jy+Vb-Xz+Yj+Zx+fØ-gÐ)e0254",
    ///     "+2(-Af-Bh+Cv+Dj-Eø+Fa+Gz+Hb-Iy-Jd+Vc-Xð-Yi+Zg-eØ+xÐ)e0235",
    ///     "+2(-Ag-Bi-Cj+Dv+Eð-Fz+Ga+Hy+Ib+Jc+Vd-Xø+Yh-Zf+eÐ+xØ)e0243",
    ///     "+2(+Ab-Ba-Cø+Dð+Ev+Fh+Gi-Hf-Ig-Jx+Ve-Xj-Yz+Zy-cØ+dÐ)e0145",
    ///     "+2(+Ac+Bø-Ca-Dz-Eh+Fv+Gj+He+Ix-Jg+Vf+Xi-Yð-Zd+bØ+yÐ)e0153",
    ///     "+2(+Ad-Bð+Cz-Da-Ei-Fj+Gv-Hx+Ie+Jf+Vg-Xh-Yø+Zc-bÐ+yØ)e0134",
    ///     "+2(-Aø+Bc-Cb+Dy+Ef-Fe-Gx+Hv+Ij-Ji+Vh-Xg+Yd-Zð-aØ+zÐ)e0125",
    ///     "+2(+Að+Bd-Cy-Db+Eg+Fx-Ge-Hj+Iv+Jh+Vi+Xf-Yc-Zø+aÐ+zØ)e0142",
    ///     "+2(-Az+By+Cd-Dc-Ex+Fg-Gf+Hi-Ih+Jv+Vj-Xe+Yb-Za-Ðø+Øð)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
    /// The multivector of single rotoreflector $`f_{r1} \equiv v^4_0 + p_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(single_rotoreflector.basis_blades(), Vee::single_rotoreflector().basis_blades());
    /// format_eq!(single_rotoreflector, [
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement()
    }
    /// The multivector of double rotoreflector $`f_{r2} \equiv v^4_0 + p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let double_rotoreflector = Vee::normal().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(double_rotoreflector.basis_blades(), Vee::double_rotoreflector().basis_blades());
    /// format_eq!(double_rotoreflector, [
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv v^4 + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+v͕ø͔e5",
    ///     "+(+X͕y͔-Y͕x͔)e012",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+Y͕z͔-Z͕y͔)e023",
    ///     "+(+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ð͕ø͔-Ø͕ð͔)e045",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume4() + Self::plane_moment()
    }
    /// The multivector of simple single flector $`f_{s1} \equiv v^4 + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let simple_single_flector = Vee::volume4().lhs() * Vee::simple_single_motor().rhs();
    ///
    /// assert_eq!(simple_single_flector.basis_blades(),
    ///     Vee::simple_single_flector().basis_blades());
    /// format_eq!(simple_single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_flector() -> Self {
        Self::volume4() + Self::plane()
    }
    /// The multivector of single flector $`f_1 \equiv v^4 + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let single_flector = Vee::volume4().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(single_flector.basis_blades(), Vee::single_flector().basis_blades());
    /// format_eq!(single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔)e123",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::volume4() + Self::plane() + Self::direction()
    }
    /// The multivector of simple double flector $`f_{s2} \equiv v^4 + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP5 as Vee};
    ///
    /// let simple_double_flector = Vee::volume4().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades());
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔-W͔x͕)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔-W͔y͕)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔-W͔z͕)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔-W͔ð͕)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔-W͔ø͕)e01243",
    /// ]);
    ///
    /// let simple_double_flector = Vee::volume4().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades());
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔)e0",
    ///     "+(-a͕y͔-b͕z͔-c͕ð͔-d͕ø͔+v͕x͔)e1",
    ///     "+(+a͕x͔-e͕z͔-f͕ð͔-g͕ø͔+v͕y͔)e2",
    ///     "+(+b͕x͔+e͕y͔-h͕ð͔-i͕ø͔+v͕z͔)e3",
    ///     "+(+c͕x͔+f͕y͔+h͕z͔-j͕ø͔+v͕ð͔)e4",
    ///     "+(+d͕x͔+g͕y͔+i͕z͔+j͕ð͔+v͕ø͔)e5",
    ///     "+(-H͕ø͔+I͕ð͔-J͕z͔+W͔a͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+F͕ø͔-G͕ð͔+J͕y͔+W͔b͕+X͕z͔-Z͕x͔)e013",
    ///     "+(-E͕ø͔+G͕z͔-I͕y͔+W͔c͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+E͕ð͔-F͕z͔+H͕y͔+W͔d͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(-C͕ø͔+D͕ð͔-J͕x͔+W͔e͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+B͕ø͔-D͕z͔+I͕x͔+W͔f͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(-B͕ð͔+C͕z͔-H͕x͔+W͔g͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-A͕ø͔+D͕y͔-G͕x͔+W͔h͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+A͕ð͔-C͕y͔+F͕x͔+W͔i͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(-A͕z͔+B͕y͔-E͕x͔+W͔j͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+h͕ø͔-i͕ð͔+j͕z͔-x͔y͕+x͕y͔)e345",
    ///     "+(-f͕ø͔+g͕ð͔-j͕y͔-x͔z͕+x͕z͔)e254",
    ///     "+(+e͕ø͔-g͕z͔+i͕y͔-x͔ð͕+x͕ð͔)e235",
    ///     "+(-e͕ð͔+f͕z͔-h͕y͔-x͔ø͕+x͕ø͔)e243",
    ///     "+(+c͕ø͔-d͕ð͔+j͕x͔-y͔z͕+y͕z͔)e145",
    ///     "+(-b͕ø͔+d͕z͔-i͕x͔-y͔ð͕+y͕ð͔)e153",
    ///     "+(+b͕ð͔-c͕z͔+h͕x͔-y͔ø͕+y͕ø͔)e134",
    ///     "+(+a͕ø͔-d͕y͔+g͕x͔-z͔ð͕+z͕ð͔)e125",
    ///     "+(-a͕ð͔+c͕y͔-f͕x͔-z͔ø͕+z͕ø͔)e142",
    ///     "+(+a͕z͔-b͕y͔+e͕x͔-ð͔ø͕+ð͕ø͔)e123",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕)e12345",
    ///     "+(+A͕y͔+B͕z͔+C͕ð͔+D͕ø͔+V͕x͔-W͔x͕)e02354",
    ///     "+(-A͕x͔+E͕z͔+F͕ð͔+G͕ø͔+V͕y͔-W͔y͕)e01345",
    ///     "+(-B͕x͔-E͕y͔+H͕ð͔+I͕ø͔+V͕z͔-W͔z͕)e01254",
    ///     "+(-C͕x͔-F͕y͔-H͕z͔+J͕ø͔+V͕ð͔-W͔ð͕)e01235",
    ///     "+(-D͕x͔-G͕y͔-I͕z͔-J͕ð͔+V͕ø͔-W͔ø͕)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_flector() -> Self {
        Self::volume4() + Self::plane() + Self::point()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 6`$ (exploratory).
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
    /// The multivector of norm $`n \equiv s + p`$.
    ///
    /// Quadvector $`p`$ does not square to a scalar, therefore $`n`$ is **not** a Study number.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let quadvector_norm_squared = Vee::plane().norm_squared();
    ///
    /// assert_eq!(quadvector_norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(quadvector_norm_squared, [
    ///     "+α̇α̇+β̇β̇+γ̇γ̇+δ̇δ̇+ε̇ε̇+ζ̇ζ̇+η̇η̇+θ̇θ̇+ι̇ι̇+κ̇κ̇+λ̇λ̇+μ̇μ̇+ν̇ν̇+ξ̇ξ̇+ο̇ο̇",
    ///     "+2(-κ̇ο̇+λ̇ξ̇-μ̇ν̇)e3456",
    ///     "+2(+η̇ο̇-θ̇ξ̇+ι̇ν̇)e2465",
    ///     "+2(-ζ̇ο̇+θ̇μ̇-ι̇λ̇)e2356",
    ///     "+2(+ζ̇ξ̇-η̇μ̇+ι̇κ̇)e2364",
    ///     "+2(-ζ̇ν̇+η̇λ̇-θ̇κ̇)e2345",
    ///     "+2(-γ̇ο̇+δ̇ξ̇-ε̇ν̇)e1456",
    ///     "+2(+β̇ο̇-δ̇μ̇+ε̇λ̇)e1365",
    ///     "+2(-β̇ξ̇+γ̇μ̇-ε̇κ̇)e1346",
    ///     "+2(+β̇ν̇-γ̇λ̇+δ̇κ̇)e1354",
    ///     "+2(-α̇ο̇+δ̇ι̇-ε̇θ̇)e1256",
    ///     "+2(+α̇ξ̇-γ̇ι̇+ε̇η̇)e1264",
    ///     "+2(-α̇ν̇+γ̇θ̇-δ̇η̇)e1245",
    ///     "+2(-α̇μ̇+β̇ι̇-ε̇ζ̇)e1236",
    ///     "+2(+α̇λ̇-β̇θ̇+δ̇ζ̇)e1253",
    ///     "+2(-α̇κ̇+β̇η̇-γ̇ζ̇)e1234",
    ///     "+2(-Ḣε̇+İδ̇-J̇γ̇-Ṅι̇+Ȯθ̇-Ṗη̇-Q̇μ̇+Ṙλ̇-Ṡκ̇)e0465",
    ///     "+2(+Ḟε̇-Ġδ̇+J̇β̇+L̇ι̇-Ṁθ̇+Ṗζ̇-Q̇ξ̇+Ṙν̇-Ṫκ̇)e0356",
    ///     "+2(-Ėε̇+Ġγ̇-İβ̇-K̇ι̇+Ṁη̇-Ȯζ̇-Q̇ο̇+Ṡν̇-Ṫλ̇)e0364",
    ///     "+2(+Ėδ̇-Ḟγ̇+Ḣβ̇+K̇θ̇-L̇η̇+Ṅζ̇-Ṙο̇+Ṡξ̇-Ṫμ̇)e0345",
    ///     "+2(-Ċε̇+Ḋδ̇-J̇α̇+L̇μ̇-Ṁλ̇+Ṅξ̇-Ȯν̇+Ṡζ̇+Ṫη̇)e0265",
    ///     "+2(+Ḃε̇-Ḋγ̇+İα̇-K̇μ̇+Ṁκ̇+Ṅο̇-Ṗν̇-Ṙζ̇+Ṫθ̇)e0246",
    ///     "+2(-Ḃδ̇+Ċγ̇-Ḣα̇+K̇λ̇-L̇κ̇+Ȯο̇-Ṗξ̇+Q̇ζ̇+Ṫι̇)e0254",
    ///     "+2(-Ȧε̇+Ḋβ̇-Ġα̇-K̇ξ̇-L̇ο̇+Ȯκ̇+Ṗλ̇-Ṙη̇-Ṡθ̇)e0263",
    ///     "+2(+Ȧδ̇-Ċβ̇+Ḟα̇+K̇ν̇-Ṁο̇-Ṅκ̇+Ṗμ̇+Q̇η̇-Ṡι̇)e0235",
    ///     "+2(-Ȧγ̇+Ḃβ̇-Ėα̇+L̇ν̇+Ṁξ̇-Ṅλ̇-Ȯμ̇+Q̇θ̇+Ṙι̇)e0243",
    ///     "+2(-Ċι̇+Ḋθ̇-Ḟμ̇+Ġλ̇-Ḣξ̇+İν̇-Ṗα̇-Ṡβ̇-Ṫγ̇)e0156",
    ///     "+2(+Ḃι̇-Ḋη̇+Ėμ̇-Ġκ̇-Ḣο̇+J̇ν̇+Ȯα̇+Ṙβ̇-Ṫδ̇)e0164",
    ///     "+2(-Ḃθ̇+Ċη̇-Ėλ̇+Ḟκ̇-İο̇+J̇ξ̇-Ṅα̇-Q̇β̇-Ṫε̇)e0145",
    ///     "+2(-Ȧι̇+Ḋζ̇+Ėξ̇+Ḟο̇-İκ̇-J̇λ̇-Ṁα̇+Ṙγ̇+Ṡδ̇)e0136",
    ///     "+2(+Ȧθ̇-Ċζ̇-Ėν̇+Ġο̇+Ḣκ̇-J̇μ̇+L̇α̇-Q̇γ̇+Ṡε̇)e0153",
    ///     "+2(-Ȧη̇+Ḃζ̇-Ḟν̇-Ġξ̇+Ḣλ̇+İμ̇-K̇α̇-Q̇δ̇-Ṙε̇)e0134",
    ///     "+2(-Ȧμ̇-Ḃξ̇-Ċο̇+Ġζ̇+İη̇+J̇θ̇-Ṁβ̇-Ȯγ̇-Ṗδ̇)e0162",
    ///     "+2(+Ȧλ̇+Ḃν̇-Ḋο̇-Ḟζ̇-Ḣη̇+J̇ι̇+L̇β̇+Ṅγ̇-Ṗε̇)e0125",
    ///     "+2(-Ȧκ̇+Ċν̇+Ḋξ̇+Ėζ̇-Ḣθ̇-İι̇-K̇β̇+Ṅδ̇+Ȯε̇)e0142",
    ///     "+2(-Ḃκ̇-Ċλ̇-Ḋμ̇+Ėη̇+Ḟθ̇+Ġι̇-K̇γ̇-L̇δ̇-Ṁε̇)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::plane()
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
        Self::e01() + Self::e02() + Self::e03() + Self::e04() + Self::e05() + Self::e06()
    }
    /// The multivector of $`4`$-volume displacement $`v^4_0`$.
    #[must_use]
    #[inline]
    pub fn volume4_displacement() -> Self {
        Self::e12()
            + Self::e13()
            + Self::e14()
            + Self::e15()
            + Self::e16()
            + Self::e23()
            + Self::e24()
            + Self::e25()
            + Self::e26()
            + Self::e34()
            + Self::e35()
            + Self::e36()
            + Self::e45()
            + Self::e46()
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
        Self::e012()
            + Self::e013()
            + Self::e014()
            + Self::e015()
            + Self::e016()
            + Self::e023()
            + Self::e024()
            + Self::e025()
            + Self::e026()
            + Self::e034()
            + Self::e035()
            + Self::e036()
            + Self::e045()
            + Self::e046()
            + Self::e056()
    }
    /// The multivector of volume displacement $`v_0`$.
    #[must_use]
    #[inline]
    pub fn volume_displacement() -> Self {
        Self::e123()
            + Self::e124()
            + Self::e125()
            + Self::e126()
            + Self::e134()
            + Self::e135()
            + Self::e136()
            + Self::e145()
            + Self::e146()
            + Self::e156()
            + Self::e234()
            + Self::e235()
            + Self::e236()
            + Self::e245()
            + Self::e246()
            + Self::e256()
            + Self::e345()
            + Self::e346()
            + Self::e356()
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
        (Self::e0123()
            + Self::e0142()
            + Self::e0125()
            + Self::e0162()
            + Self::e0134()
            + Self::e0153()
            + Self::e0136()
            + Self::e0145()
            + Self::e0164()
            + Self::e0156()
            + Self::e0243()
            + Self::e0235()
            + Self::e0263()
            + Self::e0254()
            + Self::e0246()
            + Self::e0265()
            + Self::e0345()
            + Self::e0364()
            + Self::e0356()
            + Self::e0465())
        .alt()
    }
    /// The multivector of plane displacement $`p_0`$.
    #[must_use]
    #[inline]
    pub fn plane_displacement() -> Self {
        (Self::e1234()
            + Self::e1253()
            + Self::e1236()
            + Self::e1245()
            + Self::e1264()
            + Self::e1256()
            + Self::e1354()
            + Self::e1346()
            + Self::e1365()
            + Self::e1456()
            + Self::e2345()
            + Self::e2364()
            + Self::e2356()
            + Self::e2465()
            + Self::e3456())
        .alt()
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
        (Self::e01234()
            + Self::e01253()
            + Self::e01236()
            + Self::e01245()
            + Self::e01264()
            + Self::e01256()
            + Self::e01354()
            + Self::e01346()
            + Self::e01365()
            + Self::e01456()
            + Self::e02345()
            + Self::e02364()
            + Self::e02356()
            + Self::e02465()
            + Self::e03456())
        .alt()
    }
    /// The multivector of line displacement $`\ell_0`$.
    #[must_use]
    #[inline]
    pub fn line_displacement() -> Self {
        (Self::e12354()
            + Self::e12346()
            + Self::e12365()
            + Self::e12456()
            + Self::e13465()
            + Self::e23456())
        .alt()
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
        (Self::e012345()
            + Self::e012364()
            + Self::e012356()
            + Self::e012465()
            + Self::e013456()
            + Self::e023465())
        .alt()
    }
    /// The multivector of weight $`P_0`$.
    #[must_use]
    #[inline]
    pub fn weight() -> Self {
        Self::e123456().alt()
    }
    /// The multivector of point $`P \equiv P_0 + P_\infty`$.
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::direction() + Self::weight()
    }
    /// The multivector of single rotator $`r_1 \equiv s + v^4_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+x͔þ͕-x͕þ͔)e16",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+y͔þ͕-y͕þ͔)e26",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+z͔þ͕-z͕þ͔)e36",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(+ð͔þ͕-ð͕þ͔)e46",
    ///     "+(+ø͔þ͕-ø͕þ͔)e56",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+ẋ͔ẋ͕+ẏ͔ẏ͕+ż͔ż͕+ð͔̇ð͕̇+ø͔̇ø͕̇+þ͔̇þ͕̇",
    ///     "+(+ẋ͔ẏ͕-ẋ͕ẏ͔)e12",
    ///     "+(+ẋ͔ż͕-ẋ͕ż͔)e13",
    ///     "+(+ẋ͔ð͕̇-ẋ͕ð͔̇)e14",
    ///     "+(+ẋ͔ø͕̇-ẋ͕ø͔̇)e15",
    ///     "+(+ẋ͔þ͕̇-ẋ͕þ͔̇)e16",
    ///     "+(+ẏ͔ż͕-ẏ͕ż͔)e23",
    ///     "+(+ẏ͔ð͕̇-ẏ͕ð͔̇)e24",
    ///     "+(+ẏ͔ø͕̇-ẏ͕ø͔̇)e25",
    ///     "+(+ẏ͔þ͕̇-ẏ͕þ͔̇)e26",
    ///     "+(+ż͔ð͕̇-ż͕ð͔̇)e34",
    ///     "+(+ż͔ø͕̇-ż͕ø͔̇)e35",
    ///     "+(+ż͔þ͕̇-ż͕þ͔̇)e36",
    ///     "+(+ð͔̇ø͕̇-ð͕̇ø͔̇)e45",
    ///     "+(+ð͔̇þ͕̇-ð͕̇þ͔̇)e46",
    ///     "+(+ø͔̇þ͕̇-ø͕̇þ͔̇)e56",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::volume4_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + v^4_0 + p_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "+v͔v͕-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕",
    ///     "+(+v͔α͕+v͕α͔-β͔ζ͕+β͕ζ͔-γ͔η͕+γ͕η͔-δ͔θ͕+δ͕θ͔-ε͔ι͕+ε͕ι͔)e12",
    ///     "+(+v͔β͕+v͕β͔+α͔ζ͕-α͕ζ͔-γ͔κ͕+γ͕κ͔-δ͔λ͕+δ͕λ͔-ε͔μ͕+ε͕μ͔)e13",
    ///     "+(+v͔γ͕+v͕γ͔+α͔η͕-α͕η͔+β͔κ͕-β͕κ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔)e14",
    ///     "+(+v͔δ͕+v͕δ͔+α͔θ͕-α͕θ͔+β͔λ͕-β͕λ͔+γ͔ν͕-γ͕ν͔-ε͔ο͕+ε͕ο͔)e15",
    ///     "+(+v͔ε͕+v͕ε͔+α͔ι͕-α͕ι͔+β͔μ͕-β͕μ͔+γ͔ξ͕-γ͕ξ͔+δ͔ο͕-δ͕ο͔)e16",
    ///     "+(+v͔ζ͕+v͕ζ͔-α͔β͕+α͕β͔-η͔κ͕+η͕κ͔-θ͔λ͕+θ͕λ͔-ι͔μ͕+ι͕μ͔)e23",
    ///     "+(+v͔η͕+v͕η͔-α͔γ͕+α͕γ͔+ζ͔κ͕-ζ͕κ͔-θ͔ν͕+θ͕ν͔-ι͔ξ͕+ι͕ξ͔)e24",
    ///     "+(+v͔θ͕+v͕θ͔-α͔δ͕+α͕δ͔+ζ͔λ͕-ζ͕λ͔+η͔ν͕-η͕ν͔-ι͔ο͕+ι͕ο͔)e25",
    ///     "+(+v͔ι͕+v͕ι͔-α͔ε͕+α͕ε͔+ζ͔μ͕-ζ͕μ͔+η͔ξ͕-η͕ξ͔+θ͔ο͕-θ͕ο͔)e26",
    ///     "+(+v͔κ͕+v͕κ͔-β͔γ͕+β͕γ͔-ζ͔η͕+ζ͕η͔-λ͔ν͕+λ͕ν͔-μ͔ξ͕+μ͕ξ͔)e34",
    ///     "+(+v͔λ͕+v͕λ͔-β͔δ͕+β͕δ͔-ζ͔θ͕+ζ͕θ͔+κ͔ν͕-κ͕ν͔-μ͔ο͕+μ͕ο͔)e35",
    ///     "+(+v͔μ͕+v͕μ͔-β͔ε͕+β͕ε͔-ζ͔ι͕+ζ͕ι͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e36",
    ///     "+(+v͔ν͕+v͕ν͔-γ͔δ͕+γ͕δ͔-η͔θ͕+η͕θ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔)e45",
    ///     "+(+v͔ξ͕+v͕ξ͔-γ͔ε͕+γ͕ε͔-η͔ι͕+η͕ι͔-κ͔μ͕+κ͕μ͔+ν͔ο͕-ν͕ο͔)e46",
    ///     "+(+v͔ο͕+v͕ο͔-δ͔ε͕+δ͕ε͔-θ͔ι͕+θ͕ι͔-λ͔μ͕+λ͕μ͔-ν͔ξ͕+ν͕ξ͔)e56",
    ///     "+(+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔+μ͔ν͕+μ͕ν͔)e3456",
    ///     "+(-η͔ο͕-η͕ο͔+θ͔ξ͕+θ͕ξ͔-ι͔ν͕-ι͕ν͔)e2465",
    ///     "+(+ζ͔ο͕+ζ͕ο͔-θ͔μ͕-θ͕μ͔+ι͔λ͕+ι͕λ͔)e2356",
    ///     "+(-ζ͔ξ͕-ζ͕ξ͔+η͔μ͕+η͕μ͔-ι͔κ͕-ι͕κ͔)e2364",
    ///     "+(+ζ͔ν͕+ζ͕ν͔-η͔λ͕-η͕λ͔+θ͔κ͕+θ͕κ͔)e2345",
    ///     "+(+γ͔ο͕+γ͕ο͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔)e1456",
    ///     "+(-β͔ο͕-β͕ο͔+δ͔μ͕+δ͕μ͔-ε͔λ͕-ε͕λ͔)e1365",
    ///     "+(+β͔ξ͕+β͕ξ͔-γ͔μ͕-γ͕μ͔+ε͔κ͕+ε͕κ͔)e1346",
    ///     "+(-β͔ν͕-β͕ν͔+γ͔λ͕+γ͕λ͔-δ͔κ͕-δ͕κ͔)e1354",
    ///     "+(+α͔ο͕+α͕ο͔-δ͔ι͕-δ͕ι͔+ε͔θ͕+ε͕θ͔)e1256",
    ///     "+(-α͔ξ͕-α͕ξ͔+γ͔ι͕+γ͕ι͔-ε͔η͕-ε͕η͔)e1264",
    ///     "+(+α͔ν͕+α͕ν͔-γ͔θ͕-γ͕θ͔+δ͔η͕+δ͕η͔)e1245",
    ///     "+(+α͔μ͕+α͕μ͔-β͔ι͕-β͕ι͔+ε͔ζ͕+ε͕ζ͔)e1236",
    ///     "+(-α͔λ͕-α͕λ͔+β͔θ͕+β͕θ͔-δ͔ζ͕-δ͕ζ͔)e1253",
    ///     "+(+α͔κ͕+α͕κ͔-β͔η͕-β͕η͔+γ͔ζ͕+γ͕ζ͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::volume4_displacement().lhs() * Vee::volume4_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕",
    ///     "+(-β͔ζ͕+β͕ζ͔-γ͔η͕+γ͕η͔-δ͔θ͕+δ͕θ͔-ε͔ι͕+ε͕ι͔)e12",
    ///     "+(+α͔ζ͕-α͕ζ͔-γ͔κ͕+γ͕κ͔-δ͔λ͕+δ͕λ͔-ε͔μ͕+ε͕μ͔)e13",
    ///     "+(+α͔η͕-α͕η͔+β͔κ͕-β͕κ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔)e14",
    ///     "+(+α͔θ͕-α͕θ͔+β͔λ͕-β͕λ͔+γ͔ν͕-γ͕ν͔-ε͔ο͕+ε͕ο͔)e15",
    ///     "+(+α͔ι͕-α͕ι͔+β͔μ͕-β͕μ͔+γ͔ξ͕-γ͕ξ͔+δ͔ο͕-δ͕ο͔)e16",
    ///     "+(-α͔β͕+α͕β͔-η͔κ͕+η͕κ͔-θ͔λ͕+θ͕λ͔-ι͔μ͕+ι͕μ͔)e23",
    ///     "+(-α͔γ͕+α͕γ͔+ζ͔κ͕-ζ͕κ͔-θ͔ν͕+θ͕ν͔-ι͔ξ͕+ι͕ξ͔)e24",
    ///     "+(-α͔δ͕+α͕δ͔+ζ͔λ͕-ζ͕λ͔+η͔ν͕-η͕ν͔-ι͔ο͕+ι͕ο͔)e25",
    ///     "+(-α͔ε͕+α͕ε͔+ζ͔μ͕-ζ͕μ͔+η͔ξ͕-η͕ξ͔+θ͔ο͕-θ͕ο͔)e26",
    ///     "+(-β͔γ͕+β͕γ͔-ζ͔η͕+ζ͕η͔-λ͔ν͕+λ͕ν͔-μ͔ξ͕+μ͕ξ͔)e34",
    ///     "+(-β͔δ͕+β͕δ͔-ζ͔θ͕+ζ͕θ͔+κ͔ν͕-κ͕ν͔-μ͔ο͕+μ͕ο͔)e35",
    ///     "+(-β͔ε͕+β͕ε͔-ζ͔ι͕+ζ͕ι͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e36",
    ///     "+(-γ͔δ͕+γ͕δ͔-η͔θ͕+η͕θ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔)e45",
    ///     "+(-γ͔ε͕+γ͕ε͔-η͔ι͕+η͕ι͔-κ͔μ͕+κ͕μ͔+ν͔ο͕-ν͕ο͔)e46",
    ///     "+(-δ͔ε͕+δ͕ε͔-θ͔ι͕+θ͕ι͔-λ͔μ͕+λ͕μ͔-ν͔ξ͕+ν͕ξ͔)e56",
    ///     "+(+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔+μ͔ν͕+μ͕ν͔)e3456",
    ///     "+(-η͔ο͕-η͕ο͔+θ͔ξ͕+θ͕ξ͔-ι͔ν͕-ι͕ν͔)e2465",
    ///     "+(+ζ͔ο͕+ζ͕ο͔-θ͔μ͕-θ͕μ͔+ι͔λ͕+ι͕λ͔)e2356",
    ///     "+(-ζ͔ξ͕-ζ͕ξ͔+η͔μ͕+η͕μ͔-ι͔κ͕-ι͕κ͔)e2364",
    ///     "+(+ζ͔ν͕+ζ͕ν͔-η͔λ͕-η͕λ͔+θ͔κ͕+θ͕κ͔)e2345",
    ///     "+(+γ͔ο͕+γ͕ο͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔)e1456",
    ///     "+(-β͔ο͕-β͕ο͔+δ͔μ͕+δ͕μ͔-ε͔λ͕-ε͕λ͔)e1365",
    ///     "+(+β͔ξ͕+β͕ξ͔-γ͔μ͕-γ͕μ͔+ε͔κ͕+ε͕κ͔)e1346",
    ///     "+(-β͔ν͕-β͕ν͔+γ͔λ͕+γ͕λ͔-δ͔κ͕-δ͕κ͔)e1354",
    ///     "+(+α͔ο͕+α͕ο͔-δ͔ι͕-δ͕ι͔+ε͔θ͕+ε͕θ͔)e1256",
    ///     "+(-α͔ξ͕-α͕ξ͔+γ͔ι͕+γ͕ι͔-ε͔η͕-ε͕η͔)e1264",
    ///     "+(+α͔ν͕+α͕ν͔-γ͔θ͕-γ͕θ͔+δ͔η͕+δ͕η͔)e1245",
    ///     "+(+α͔μ͕+α͕μ͔-β͔ι͕-β͕ι͔+ε͔ζ͕+ε͕ζ͔)e1236",
    ///     "+(-α͔λ͕-α͕λ͔+β͔θ͕+β͕θ͔-δ͔ζ͕-δ͕ζ͔)e1253",
    ///     "+(+α͔κ͕+α͕κ͔-β͔η͕-β͕η͔+γ͔ζ͕+γ͕ζ͔)e1234",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "+α͔̇α͕̇+β͔̇β͕̇+γ͔̇γ͕̇+δ͔̇δ͕̇+ε͔̇ε͕̇+ζ͔̇ζ͕̇+η͔̇η͕̇+θ͔̇θ͕̇+ι͔̇ι͕̇+κ͔̇κ͕̇+λ͔̇λ͕̇+μ͔̇μ͕̇+ν͔̇ν͕̇+ξ͔̇ξ͕̇+ο͔̇ο͕̇",
    ///     "+(+β͔̇ζ͕̇-β͕̇ζ͔̇+γ͔̇η͕̇-γ͕̇η͔̇+δ͔̇θ͕̇-δ͕̇θ͔̇+ε͔̇ι͕̇-ε͕̇ι͔̇)e12",
    ///     "+(-α͔̇ζ͕̇+α͕̇ζ͔̇+γ͔̇κ͕̇-γ͕̇κ͔̇+δ͔̇λ͕̇-δ͕̇λ͔̇+ε͔̇μ͕̇-ε͕̇μ͔̇)e13",
    ///     "+(-α͔̇η͕̇+α͕̇η͔̇-β͔̇κ͕̇+β͕̇κ͔̇+δ͔̇ν͕̇-δ͕̇ν͔̇+ε͔̇ξ͕̇-ε͕̇ξ͔̇)e14",
    ///     "+(-α͔̇θ͕̇+α͕̇θ͔̇-β͔̇λ͕̇+β͕̇λ͔̇-γ͔̇ν͕̇+γ͕̇ν͔̇+ε͔̇ο͕̇-ε͕̇ο͔̇)e15",
    ///     "+(-α͔̇ι͕̇+α͕̇ι͔̇-β͔̇μ͕̇+β͕̇μ͔̇-γ͔̇ξ͕̇+γ͕̇ξ͔̇-δ͔̇ο͕̇+δ͕̇ο͔̇)e16",
    ///     "+(+α͔̇β͕̇-α͕̇β͔̇+η͔̇κ͕̇-η͕̇κ͔̇+θ͔̇λ͕̇-θ͕̇λ͔̇+ι͔̇μ͕̇-ι͕̇μ͔̇)e23",
    ///     "+(+α͔̇γ͕̇-α͕̇γ͔̇-ζ͔̇κ͕̇+ζ͕̇κ͔̇+θ͔̇ν͕̇-θ͕̇ν͔̇+ι͔̇ξ͕̇-ι͕̇ξ͔̇)e24",
    ///     "+(+α͔̇δ͕̇-α͕̇δ͔̇-ζ͔̇λ͕̇+ζ͕̇λ͔̇-η͔̇ν͕̇+η͕̇ν͔̇+ι͔̇ο͕̇-ι͕̇ο͔̇)e25",
    ///     "+(+α͔̇ε͕̇-α͕̇ε͔̇-ζ͔̇μ͕̇+ζ͕̇μ͔̇-η͔̇ξ͕̇+η͕̇ξ͔̇-θ͔̇ο͕̇+θ͕̇ο͔̇)e26",
    ///     "+(+β͔̇γ͕̇-β͕̇γ͔̇+ζ͔̇η͕̇-ζ͕̇η͔̇+λ͔̇ν͕̇-λ͕̇ν͔̇+μ͔̇ξ͕̇-μ͕̇ξ͔̇)e34",
    ///     "+(+β͔̇δ͕̇-β͕̇δ͔̇+ζ͔̇θ͕̇-ζ͕̇θ͔̇-κ͔̇ν͕̇+κ͕̇ν͔̇+μ͔̇ο͕̇-μ͕̇ο͔̇)e35",
    ///     "+(+β͔̇ε͕̇-β͕̇ε͔̇+ζ͔̇ι͕̇-ζ͕̇ι͔̇-κ͔̇ξ͕̇+κ͕̇ξ͔̇-λ͔̇ο͕̇+λ͕̇ο͔̇)e36",
    ///     "+(+γ͔̇δ͕̇-γ͕̇δ͔̇+η͔̇θ͕̇-η͕̇θ͔̇+κ͔̇λ͕̇-κ͕̇λ͔̇+ξ͔̇ο͕̇-ξ͕̇ο͔̇)e45",
    ///     "+(+γ͔̇ε͕̇-γ͕̇ε͔̇+η͔̇ι͕̇-η͕̇ι͔̇+κ͔̇μ͕̇-κ͕̇μ͔̇-ν͔̇ο͕̇+ν͕̇ο͔̇)e46",
    ///     "+(+δ͔̇ε͕̇-δ͕̇ε͔̇+θ͔̇ι͕̇-θ͕̇ι͔̇+λ͔̇μ͕̇-λ͕̇μ͔̇+ν͔̇ξ͕̇-ν͕̇ξ͔̇)e56",
    ///     "+(-κ͔̇ο͕̇-κ͕̇ο͔̇+λ͔̇ξ͕̇+λ͕̇ξ͔̇-μ͔̇ν͕̇-μ͕̇ν͔̇)e3456",
    ///     "+(+η͔̇ο͕̇+η͕̇ο͔̇-θ͔̇ξ͕̇-θ͕̇ξ͔̇+ι͔̇ν͕̇+ι͕̇ν͔̇)e2465",
    ///     "+(-ζ͔̇ο͕̇-ζ͕̇ο͔̇+θ͔̇μ͕̇+θ͕̇μ͔̇-ι͔̇λ͕̇-ι͕̇λ͔̇)e2356",
    ///     "+(+ζ͔̇ξ͕̇+ζ͕̇ξ͔̇-η͔̇μ͕̇-η͕̇μ͔̇+ι͔̇κ͕̇+ι͕̇κ͔̇)e2364",
    ///     "+(-ζ͔̇ν͕̇-ζ͕̇ν͔̇+η͔̇λ͕̇+η͕̇λ͔̇-θ͔̇κ͕̇-θ͕̇κ͔̇)e2345",
    ///     "+(-γ͔̇ο͕̇-γ͕̇ο͔̇+δ͔̇ξ͕̇+δ͕̇ξ͔̇-ε͔̇ν͕̇-ε͕̇ν͔̇)e1456",
    ///     "+(+β͔̇ο͕̇+β͕̇ο͔̇-δ͔̇μ͕̇-δ͕̇μ͔̇+ε͔̇λ͕̇+ε͕̇λ͔̇)e1365",
    ///     "+(-β͔̇ξ͕̇-β͕̇ξ͔̇+γ͔̇μ͕̇+γ͕̇μ͔̇-ε͔̇κ͕̇-ε͕̇κ͔̇)e1346",
    ///     "+(+β͔̇ν͕̇+β͕̇ν͔̇-γ͔̇λ͕̇-γ͕̇λ͔̇+δ͔̇κ͕̇+δ͕̇κ͔̇)e1354",
    ///     "+(-α͔̇ο͕̇-α͕̇ο͔̇+δ͔̇ι͕̇+δ͕̇ι͔̇-ε͔̇θ͕̇-ε͕̇θ͔̇)e1256",
    ///     "+(+α͔̇ξ͕̇+α͕̇ξ͔̇-γ͔̇ι͕̇-γ͕̇ι͔̇+ε͔̇η͕̇+ε͕̇η͔̇)e1264",
    ///     "+(-α͔̇ν͕̇-α͕̇ν͔̇+γ͔̇θ͕̇+γ͕̇θ͔̇-δ͔̇η͕̇-δ͕̇η͔̇)e1245",
    ///     "+(-α͔̇μ͕̇-α͕̇μ͔̇+β͔̇ι͕̇+β͕̇ι͔̇-ε͔̇ζ͕̇-ε͕̇ζ͔̇)e1236",
    ///     "+(+α͔̇λ͕̇+α͕̇λ͔̇-β͔̇θ͕̇-β͕̇θ͔̇+δ͔̇ζ͕̇+δ͕̇ζ͔̇)e1253",
    ///     "+(-α͔̇κ͕̇-α͕̇κ͔̇+β͔̇η͕̇+β͕̇η͔̇-γ͔̇ζ͕̇-γ͕̇ζ͔̇)e1234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::volume4_displacement() + Self::plane_displacement()
    }
    /// The multivector of triple rotator $`r_3 \equiv s + v^4_0 + p_0 + P_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let triple_rotator = Vee::single_rotator().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(triple_rotator.basis_blades(), Vee::triple_rotator().basis_blades());
    /// format_eq!(triple_rotator, [
    ///     "+v͔v͕-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕",
    ///     "+(+v͔α͕+v͕α͔-β͔ζ͕+β͕ζ͔-γ͔η͕+γ͕η͔-δ͔θ͕+δ͕θ͔-ε͔ι͕+ε͕ι͔-κ͔ο͕̇-κ͕̇ο͔+λ͔ξ͕̇+λ͕̇ξ͔-μ͔ν͕̇-μ͕̇ν͔)e12",
    ///     "+(+v͔β͕+v͕β͔+α͔ζ͕-α͕ζ͔-γ͔κ͕+γ͕κ͔-δ͔λ͕+δ͕λ͔-ε͔μ͕+ε͕μ͔+η͔ο͕̇+η͕̇ο͔-θ͔ξ͕̇-θ͕̇ξ͔+ι͔ν͕̇+ι͕̇ν͔)e13",
    ///     "+(+v͔γ͕+v͕γ͔+α͔η͕-α͕η͔+β͔κ͕-β͕κ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔-ζ͔ο͕̇-ζ͕̇ο͔+θ͔μ͕̇+θ͕̇μ͔-ι͔λ͕̇-ι͕̇λ͔)e14",
    ///     "+(+v͔δ͕+v͕δ͔+α͔θ͕-α͕θ͔+β͔λ͕-β͕λ͔+γ͔ν͕-γ͕ν͔-ε͔ο͕+ε͕ο͔+ζ͔ξ͕̇+ζ͕̇ξ͔-η͔μ͕̇-η͕̇μ͔+ι͔κ͕̇+ι͕̇κ͔)e15",
    ///     "+(+v͔ε͕+v͕ε͔+α͔ι͕-α͕ι͔+β͔μ͕-β͕μ͔+γ͔ξ͕-γ͕ξ͔+δ͔ο͕-δ͕ο͔-ζ͔ν͕̇-ζ͕̇ν͔+η͔λ͕̇+η͕̇λ͔-θ͔κ͕̇-θ͕̇κ͔)e16",
    ///     "+(+v͔ζ͕+v͕ζ͔-α͔β͕+α͕β͔-γ͔ο͕̇-γ͕̇ο͔+δ͔ξ͕̇+δ͕̇ξ͔-ε͔ν͕̇-ε͕̇ν͔-η͔κ͕+η͕κ͔-θ͔λ͕+θ͕λ͔-ι͔μ͕+ι͕μ͔)e23",
    ///     "+(+v͔η͕+v͕η͔-α͔γ͕+α͕γ͔+β͔ο͕̇+β͕̇ο͔-δ͔μ͕̇-δ͕̇μ͔+ε͔λ͕̇+ε͕̇λ͔+ζ͔κ͕-ζ͕κ͔-θ͔ν͕+θ͕ν͔-ι͔ξ͕+ι͕ξ͔)e24",
    ///     "+(+v͔θ͕+v͕θ͔-α͔δ͕+α͕δ͔-β͔ξ͕̇-β͕̇ξ͔+γ͔μ͕̇+γ͕̇μ͔-ε͔κ͕̇-ε͕̇κ͔+ζ͔λ͕-ζ͕λ͔+η͔ν͕-η͕ν͔-ι͔ο͕+ι͕ο͔)e25",
    ///     "+(+v͔ι͕+v͕ι͔-α͔ε͕+α͕ε͔+β͔ν͕̇+β͕̇ν͔-γ͔λ͕̇-γ͕̇λ͔+δ͔κ͕̇+δ͕̇κ͔+ζ͔μ͕-ζ͕μ͔+η͔ξ͕-η͕ξ͔+θ͔ο͕-θ͕ο͔)e26",
    ///     "+(+v͔κ͕+v͕κ͔-α͔ο͕̇-α͕̇ο͔-β͔γ͕+β͕γ͔+δ͔ι͕̇+δ͕̇ι͔-ε͔θ͕̇-ε͕̇θ͔-ζ͔η͕+ζ͕η͔-λ͔ν͕+λ͕ν͔-μ͔ξ͕+μ͕ξ͔)e34",
    ///     "+(+v͔λ͕+v͕λ͔+α͔ξ͕̇+α͕̇ξ͔-β͔δ͕+β͕δ͔-γ͔ι͕̇-γ͕̇ι͔+ε͔η͕̇+ε͕̇η͔-ζ͔θ͕+ζ͕θ͔+κ͔ν͕-κ͕ν͔-μ͔ο͕+μ͕ο͔)e35",
    ///     "+(+v͔μ͕+v͕μ͔-α͔ν͕̇-α͕̇ν͔-β͔ε͕+β͕ε͔+γ͔θ͕̇+γ͕̇θ͔-δ͔η͕̇-δ͕̇η͔-ζ͔ι͕+ζ͕ι͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e36",
    ///     "+(+v͔ν͕+v͕ν͔-α͔μ͕̇-α͕̇μ͔+β͔ι͕̇+β͕̇ι͔-γ͔δ͕+γ͕δ͔-ε͔ζ͕̇-ε͕̇ζ͔-η͔θ͕+η͕θ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔)e45",
    ///     "+(+v͔ξ͕+v͕ξ͔+α͔λ͕̇+α͕̇λ͔-β͔θ͕̇-β͕̇θ͔-γ͔ε͕+γ͕ε͔+δ͔ζ͕̇+δ͕̇ζ͔-η͔ι͕+η͕ι͔-κ͔μ͕+κ͕μ͔+ν͔ο͕-ν͕ο͔)e46",
    ///     "+(+v͔ο͕+v͕ο͔-α͔κ͕̇-α͕̇κ͔+β͔η͕̇+β͕̇η͔-γ͔ζ͕̇-γ͕̇ζ͔-δ͔ε͕+δ͕ε͔-θ͔ι͕+θ͕ι͔-λ͔μ͕+λ͕μ͔-ν͔ξ͕+ν͕ξ͔)e56",
    ///     "+(+v͔α͕̇-β͔ζ͕̇+β͕̇ζ͔-γ͔η͕̇+γ͕̇η͔-δ͔θ͕̇+δ͕̇θ͔-ε͔ι͕̇+ε͕̇ι͔+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔+μ͔ν͕+μ͕ν͔)e3456",
    ///     "+(+v͔β͕̇+α͔ζ͕̇-α͕̇ζ͔-γ͔κ͕̇+γ͕̇κ͔-δ͔λ͕̇+δ͕̇λ͔-ε͔μ͕̇+ε͕̇μ͔-η͔ο͕-η͕ο͔+θ͔ξ͕+θ͕ξ͔-ι͔ν͕-ι͕ν͔)e2465",
    ///     "+(+v͔γ͕̇+α͔η͕̇-α͕̇η͔+β͔κ͕̇-β͕̇κ͔-δ͔ν͕̇+δ͕̇ν͔-ε͔ξ͕̇+ε͕̇ξ͔+ζ͔ο͕+ζ͕ο͔-θ͔μ͕-θ͕μ͔+ι͔λ͕+ι͕λ͔)e2356",
    ///     "+(+v͔δ͕̇+α͔θ͕̇-α͕̇θ͔+β͔λ͕̇-β͕̇λ͔+γ͔ν͕̇-γ͕̇ν͔-ε͔ο͕̇+ε͕̇ο͔-ζ͔ξ͕-ζ͕ξ͔+η͔μ͕+η͕μ͔-ι͔κ͕-ι͕κ͔)e2364",
    ///     "+(+v͔ε͕̇+α͔ι͕̇-α͕̇ι͔+β͔μ͕̇-β͕̇μ͔+γ͔ξ͕̇-γ͕̇ξ͔+δ͔ο͕̇-δ͕̇ο͔+ζ͔ν͕+ζ͕ν͔-η͔λ͕-η͕λ͔+θ͔κ͕+θ͕κ͔)e2345",
    ///     "+(+v͔ζ͕̇-α͔β͕̇+α͕̇β͔+γ͔ο͕+γ͕ο͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔-η͔κ͕̇+η͕̇κ͔-θ͔λ͕̇+θ͕̇λ͔-ι͔μ͕̇+ι͕̇μ͔)e1456",
    ///     "+(+v͔η͕̇-α͔γ͕̇+α͕̇γ͔-β͔ο͕-β͕ο͔+δ͔μ͕+δ͕μ͔-ε͔λ͕-ε͕λ͔+ζ͔κ͕̇-ζ͕̇κ͔-θ͔ν͕̇+θ͕̇ν͔-ι͔ξ͕̇+ι͕̇ξ͔)e1365",
    ///     "+(+v͔θ͕̇-α͔δ͕̇+α͕̇δ͔+β͔ξ͕+β͕ξ͔-γ͔μ͕-γ͕μ͔+ε͔κ͕+ε͕κ͔+ζ͔λ͕̇-ζ͕̇λ͔+η͔ν͕̇-η͕̇ν͔-ι͔ο͕̇+ι͕̇ο͔)e1346",
    ///     "+(+v͔ι͕̇-α͔ε͕̇+α͕̇ε͔-β͔ν͕-β͕ν͔+γ͔λ͕+γ͕λ͔-δ͔κ͕-δ͕κ͔+ζ͔μ͕̇-ζ͕̇μ͔+η͔ξ͕̇-η͕̇ξ͔+θ͔ο͕̇-θ͕̇ο͔)e1354",
    ///     "+(+v͔κ͕̇+α͔ο͕+α͕ο͔-β͔γ͕̇+β͕̇γ͔-δ͔ι͕-δ͕ι͔+ε͔θ͕+ε͕θ͔-ζ͔η͕̇+ζ͕̇η͔-λ͔ν͕̇+λ͕̇ν͔-μ͔ξ͕̇+μ͕̇ξ͔)e1256",
    ///     "+(+v͔λ͕̇-α͔ξ͕-α͕ξ͔-β͔δ͕̇+β͕̇δ͔+γ͔ι͕+γ͕ι͔-ε͔η͕-ε͕η͔-ζ͔θ͕̇+ζ͕̇θ͔+κ͔ν͕̇-κ͕̇ν͔-μ͔ο͕̇+μ͕̇ο͔)e1264",
    ///     "+(+v͔μ͕̇+α͔ν͕+α͕ν͔-β͔ε͕̇+β͕̇ε͔-γ͔θ͕-γ͕θ͔+δ͔η͕+δ͕η͔-ζ͔ι͕̇+ζ͕̇ι͔+κ͔ξ͕̇-κ͕̇ξ͔+λ͔ο͕̇-λ͕̇ο͔)e1245",
    ///     "+(+v͔ν͕̇+α͔μ͕+α͕μ͔-β͔ι͕-β͕ι͔-γ͔δ͕̇+γ͕̇δ͔+ε͔ζ͕+ε͕ζ͔-η͔θ͕̇+η͕̇θ͔-κ͔λ͕̇+κ͕̇λ͔-ξ͔ο͕̇+ξ͕̇ο͔)e1236",
    ///     "+(+v͔ξ͕̇-α͔λ͕-α͕λ͔+β͔θ͕+β͕θ͔-γ͔ε͕̇+γ͕̇ε͔-δ͔ζ͕-δ͕ζ͔-η͔ι͕̇+η͕̇ι͔-κ͔μ͕̇+κ͕̇μ͔+ν͔ο͕̇-ν͕̇ο͔)e1253",
    ///     "+(+v͔ο͕̇+α͔κ͕+α͕κ͔-β͔η͕-β͕η͔+γ͔ζ͕+γ͕ζ͔-δ͔ε͕̇+δ͕̇ε͔-θ͔ι͕̇+θ͕̇ι͔-λ͔μ͕̇+λ͕̇μ͔-ν͔ξ͕̇+ν͕̇ξ͔)e1234",
    ///     "+(+α͔α͕̇+β͔β͕̇+γ͔γ͕̇+δ͔δ͕̇+ε͔ε͕̇+ζ͔ζ͕̇+η͔η͕̇+θ͔θ͕̇+ι͔ι͕̇+κ͔κ͕̇+λ͔λ͕̇+μ͔μ͕̇+ν͔ν͕̇+ξ͔ξ͕̇+ο͔ο͕̇)e123456",
    /// ]);
    ///
    /// let triple_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(triple_rotator.basis_blades(), Vee::triple_rotator().basis_blades());
    /// format_eq!(triple_rotator, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕-k͔k͕-l͔l͕-m͔m͕-n͔n͕-o͔o͕-p͔p͕-q͔q͕-r͔r͕-s͔s͕-t͔t͕",
    ///     "+(-e͔k͕+e͕k͔-f͔l͕+f͕l͔-g͔m͕+g͕m͔-h͔n͕+h͕n͔-i͔o͕+i͕o͔-j͔p͕+j͕p͔)e12",
    ///     "+(+b͔k͕-b͕k͔+c͔l͕-c͕l͔+d͔m͕-d͕m͔-h͔q͕+h͕q͔-i͔r͕+i͕r͔-j͔s͕+j͕s͔)e13",
    ///     "+(-a͔k͕+a͕k͔+c͔n͕-c͕n͔+d͔o͕-d͕o͔+f͔q͕-f͕q͔+g͔r͕-g͕r͔-j͔t͕+j͕t͔)e14",
    ///     "+(-a͔l͕+a͕l͔-b͔n͕+b͕n͔+d͔p͕-d͕p͔-e͔q͕+e͕q͔+g͔s͕-g͕s͔+i͔t͕-i͕t͔)e15",
    ///     "+(-a͔m͕+a͕m͔-b͔o͕+b͕o͔-c͔p͕+c͕p͔-e͔r͕+e͕r͔-f͔s͕+f͕s͔-h͔t͕+h͕t͔)e16",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔-n͔q͕+n͕q͔-o͔r͕+o͕r͔-p͔s͕+p͕s͔)e23",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔+l͔q͕-l͕q͔+m͔r͕-m͕r͔-p͔t͕+p͕t͔)e24",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔-k͔q͕+k͕q͔+m͔s͕-m͕s͔+o͔t͕-o͕t͔)e25",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔-k͔r͕+k͕r͔-l͔s͕+l͕s͔-n͔t͕+n͕t͔)e26",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔-l͔n͕+l͕n͔-m͔o͕+m͕o͔-s͔t͕+s͕t͔)e34",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔+k͔n͕-k͕n͔-m͔p͕+m͕p͔+r͔t͕-r͕t͔)e35",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔+k͔o͕-k͕o͔+l͔p͕-l͕p͔-q͔t͕+q͕t͔)e36",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔-k͔l͕+k͕l͔-o͔p͕+o͕p͔-r͔s͕+r͕s͔)e45",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔-k͔m͕+k͕m͔+n͔p͕-n͕p͔+q͔s͕-q͕s͔)e46",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔-l͔m͕+l͕m͔-n͔o͕+n͕o͔-q͔r͕+q͕r͔)e56",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔+k͔p͕+k͕p͔-l͔o͕-l͕o͔+m͔n͕+m͕n͔)e3456",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔+k͔s͕+k͕s͔-l͔r͕-l͕r͔+m͔q͕+m͕q͔)e2465",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔+k͔t͕+k͕t͔-n͔r͕-n͕r͔+o͔q͕+o͕q͔)e2356",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔+l͔t͕+l͕t͔-n͔s͕-n͕s͔+p͔q͕+p͕q͔)e2364",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔+m͔t͕+m͕t͔-o͔s͕-o͕s͔+p͔r͕+p͕r͔)e2345",
    ///     "+(-b͔p͕-b͕p͔+c͔o͕+c͕o͔-d͔n͕-d͕n͔-e͔s͕-e͕s͔+f͔r͕+f͕r͔-g͔q͕-g͕q͔)e1456",
    ///     "+(+a͔p͕+a͕p͔-c͔m͕-c͕m͔+d͔l͕+d͕l͔-e͔t͕-e͕t͔+h͔r͕+h͕r͔-i͔q͕-i͕q͔)e1365",
    ///     "+(-a͔o͕-a͕o͔+b͔m͕+b͕m͔-d͔k͕-d͕k͔-f͔t͕-f͕t͔+h͔s͕+h͕s͔-j͔q͕-j͕q͔)e1346",
    ///     "+(+a͔n͕+a͕n͔-b͔l͕-b͕l͔+c͔k͕+c͕k͔-g͔t͕-g͕t͔+i͔s͕+i͕s͔-j͔r͕-j͕r͔)e1354",
    ///     "+(+a͔s͕+a͕s͔+b͔t͕+b͕t͔-f͔m͕-f͕m͔+g͔l͕+g͕l͔-h͔o͕-h͕o͔+i͔n͕+i͕n͔)e1256",
    ///     "+(-a͔r͕-a͕r͔+c͔t͕+c͕t͔+e͔m͕+e͕m͔-g͔k͕-g͕k͔-h͔p͕-h͕p͔+j͔n͕+j͕n͔)e1264",
    ///     "+(+a͔q͕+a͕q͔+d͔t͕+d͕t͔-e͔l͕-e͕l͔+f͔k͕+f͕k͔-i͔p͕-i͕p͔+j͔o͕+j͕o͔)e1245",
    ///     "+(-b͔r͕-b͕r͔-c͔s͕-c͕s͔+e͔o͕+e͕o͔+f͔p͕+f͕p͔-i͔k͕-i͕k͔-j͔l͕-j͕l͔)e1236",
    ///     "+(+b͔q͕+b͕q͔-d͔s͕-d͕s͔-e͔n͕-e͕n͔+g͔p͕+g͕p͔+h͔k͕+h͕k͔-j͔m͕-j͕m͔)e1253",
    ///     "+(+c͔q͕+c͕q͔+d͔r͕+d͕r͔-f͔n͕-f͕n͔-g͔o͕-g͕o͔+h͔l͕+h͕l͔+i͔m͕+i͕m͔)e1234",
    ///     "+(+a͔t͕-a͕t͔-b͔s͕+b͕s͔+c͔r͕-c͕r͔-d͔q͕+d͕q͔+e͔p͕-e͕p͔-f͔o͕+f͕o͔+g͔n͕-g͕n͔+h͔m͕-h͕m͔-i͔l͕+i͕l͔+j͔k͕-j͕k͔)e123456",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_rotator() -> Self {
        Self::scalar() + Self::volume4_displacement() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + v^4_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-ẇ͔ẇ͕",
    ///     "+(+Ẋ͔ẇ͕-Ẋ͕ẇ͔)e01",
    ///     "+(+Ẏ͔ẇ͕-Ẏ͕ẇ͔)e02",
    ///     "+(+Ż͔ẇ͕-Ż͕ẇ͔)e03",
    ///     "+(-ẇ͔Ð͕̇+ẇ͕Ð͔̇)e04",
    ///     "+(-ẇ͔Ø͕̇+ẇ͕Ø͔̇)e05",
    ///     "+(-ẇ͔Þ͕̇+ẇ͕Þ͔̇)e06",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::volume4_moment()
    }
    /// The multivector of single motor $`m_1 \equiv s + v^4`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(),
    ///     (Vee::single_motor() + Vee::plane_moment()).basis_blades());
    /// format_eq!(single_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕ζ͔+Ð͕η͔+Ø͕θ͔+Þ͕ι͔)e02",
    ///     "+(-X͕β͔-Y͕ζ͔+Z͕v͔+Ð͕κ͔+Ø͕λ͔+Þ͕μ͔)e03",
    ///     "+(-X͕γ͔-Y͕η͔-Z͕κ͔+v͔Ð͕+Ø͕ν͔+Þ͕ξ͔)e04",
    ///     "+(-X͕δ͔-Y͕θ͔-Z͕λ͔+v͔Ø͕-Ð͕ν͔+Þ͕ο͔)e05",
    ///     "+(-X͕ε͔-Y͕ι͔-Z͕μ͔+v͔Þ͕-Ð͕ξ͔-Ø͕ο͔)e06",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e23",
    ///     "+v͕η͔e24",
    ///     "+v͕θ͔e25",
    ///     "+v͕ι͔e26",
    ///     "+v͕κ͔e34",
    ///     "+v͕λ͔e35",
    ///     "+v͕μ͔e36",
    ///     "+v͕ν͔e45",
    ///     "+v͕ξ͔e46",
    ///     "+v͕ο͔e56",
    ///     // Additional orthogonality condition besides `norm_squared` condition:
    ///     "+(-Ð͕ο͔+Ø͕ξ͔-Þ͕ν͔)e0465",
    ///     "+(+Z͕ο͔-Ø͕μ͔+Þ͕λ͔)e0356",
    ///     "+(-Z͕ξ͔+Ð͕μ͔-Þ͕κ͔)e0364",
    ///     "+(+Z͕ν͔-Ð͕λ͔+Ø͕κ͔)e0345",
    ///     "+(-Y͕ο͔+Ø͕ι͔-Þ͕θ͔)e0265",
    ///     "+(+Y͕ξ͔-Ð͕ι͔+Þ͕η͔)e0246",
    ///     "+(-Y͕ν͔+Ð͕θ͔-Ø͕η͔)e0254",
    ///     "+(-Y͕μ͔+Z͕ι͔-Þ͕ζ͔)e0263",
    ///     "+(+Y͕λ͔-Z͕θ͔+Ø͕ζ͔)e0235",
    ///     "+(-Y͕κ͔+Z͕η͔-Ð͕ζ͔)e0243",
    ///     "+(+X͕ο͔-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(-X͕ξ͔+Ð͕ε͔-Þ͕γ͔)e0164",
    ///     "+(+X͕ν͔-Ð͕δ͔+Ø͕γ͔)e0145",
    ///     "+(+X͕μ͔-Z͕ε͔+Þ͕β͔)e0136",
    ///     "+(-X͕λ͔+Z͕δ͔-Ø͕β͔)e0153",
    ///     "+(+X͕κ͔-Z͕γ͔+Ð͕β͔)e0134",
    ///     "+(-X͕ι͔+Y͕ε͔-Þ͕α͔)e0162",
    ///     "+(+X͕θ͔-Y͕δ͔+Ø͕α͔)e0125",
    ///     "+(-X͕η͔+Y͕γ͔-Ð͕α͔)e0142",
    ///     "+(+X͕ζ͔-Y͕β͔+Z͕α͔)e0123",
    /// ]);
    ///
    /// // Orthogonality condition between single rotator and translator bivectors. Without this
    /// // condition, it is a double motor for arbitrary choices of single rotator and translator.
    /// format_eq!(Vee::single_rotator().vector(2).lhs() ^ Vee::translator().vector(2).rhs(), [
    ///     "+(-Ð͕ο͔+Ø͕ξ͔-Þ͕ν͔)e0465",
    ///     "+(+Z͕ο͔-Ø͕μ͔+Þ͕λ͔)e0356",
    ///     "+(-Z͕ξ͔+Ð͕μ͔-Þ͕κ͔)e0364",
    ///     "+(+Z͕ν͔-Ð͕λ͔+Ø͕κ͔)e0345",
    ///     "+(-Y͕ο͔+Ø͕ι͔-Þ͕θ͔)e0265",
    ///     "+(+Y͕ξ͔-Ð͕ι͔+Þ͕η͔)e0246",
    ///     "+(-Y͕ν͔+Ð͕θ͔-Ø͕η͔)e0254",
    ///     "+(-Y͕μ͔+Z͕ι͔-Þ͕ζ͔)e0263",
    ///     "+(+Y͕λ͔-Z͕θ͔+Ø͕ζ͔)e0235",
    ///     "+(-Y͕κ͔+Z͕η͔-Ð͕ζ͔)e0243",
    ///     "+(+X͕ο͔-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(-X͕ξ͔+Ð͕ε͔-Þ͕γ͔)e0164",
    ///     "+(+X͕ν͔-Ð͕δ͔+Ø͕γ͔)e0145",
    ///     "+(+X͕μ͔-Z͕ε͔+Þ͕β͔)e0136",
    ///     "+(-X͕λ͔+Z͕δ͔-Ø͕β͔)e0153",
    ///     "+(+X͕κ͔-Z͕γ͔+Ð͕β͔)e0134",
    ///     "+(-X͕ι͔+Y͕ε͔-Þ͕α͔)e0162",
    ///     "+(+X͕θ͔-Y͕δ͔+Ø͕α͔)e0125",
    ///     "+(-X͕η͔+Y͕γ͔-Ð͕α͔)e0142",
    ///     "+(+X͕ζ͔-Y͕β͔+Z͕α͔)e0123",
    /// ]);
    ///
    /// let norm_squared = single_motor.norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(),
    ///     (Vee::scalar() + Vee::plane_displacement()).basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+v͔v͔v͕v͕+v͕v͕α͔α͔+v͕v͕β͔β͔+v͕v͕γ͔γ͔+v͕v͕δ͔δ͔+v͕v͕ε͔ε͔+v͕v͕ζ͔ζ͔+v͕v͕η͔η͔+v͕v͕θ͔θ͔+v͕v͕ι͔ι͔+v͕v͕κ͔κ͔+v͕v͕λ͔λ͔+v͕v͕μ͔μ͔+v͕v͕ν͔ν͔+v͕v͕ξ͔ξ͔+v͕v͕ο͔ο͔",
    ///     // Orthogonality condition:
    ///     "+2(-v͕v͕κ͔ο͔+v͕v͕λ͔ξ͔-v͕v͕μ͔ν͔)e3456",
    ///     "+2(+v͕v͕η͔ο͔-v͕v͕θ͔ξ͔+v͕v͕ι͔ν͔)e2465",
    ///     "+2(-v͕v͕ζ͔ο͔+v͕v͕θ͔μ͔-v͕v͕ι͔λ͔)e2356",
    ///     "+2(+v͕v͕ζ͔ξ͔-v͕v͕η͔μ͔+v͕v͕ι͔κ͔)e2364",
    ///     "+2(-v͕v͕ζ͔ν͔+v͕v͕η͔λ͔-v͕v͕θ͔κ͔)e2345",
    ///     "+2(-v͕v͕γ͔ο͔+v͕v͕δ͔ξ͔-v͕v͕ε͔ν͔)e1456",
    ///     "+2(+v͕v͕β͔ο͔-v͕v͕δ͔μ͔+v͕v͕ε͔λ͔)e1365",
    ///     "+2(-v͕v͕β͔ξ͔+v͕v͕γ͔μ͔-v͕v͕ε͔κ͔)e1346",
    ///     "+2(+v͕v͕β͔ν͔-v͕v͕γ͔λ͔+v͕v͕δ͔κ͔)e1354",
    ///     "+2(-v͕v͕α͔ο͔+v͕v͕δ͔ι͔-v͕v͕ε͔θ͔)e1256",
    ///     "+2(+v͕v͕α͔ξ͔-v͕v͕γ͔ι͔+v͕v͕ε͔η͔)e1264",
    ///     "+2(-v͕v͕α͔ν͔+v͕v͕γ͔θ͔-v͕v͕δ͔η͔)e1245",
    ///     "+2(-v͕v͕α͔μ͔+v͕v͕β͔ι͔-v͕v͕ε͔ζ͔)e1236",
    ///     "+2(+v͕v͕α͔λ͔-v͕v͕β͔θ͔+v͕v͕δ͔ζ͔)e1253",
    ///     "+2(-v͕v͕α͔κ͔+v͕v͕β͔η͔-v͕v͕γ͔ζ͔)e1234",
    /// ]);
    ///
    /// let single_motor = Vee::volume5().lhs() * Vee::volume5().rhs();
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+W͔ø͕-W͕ø͔)e05",
    ///     "+(+W͔þ͕-W͕þ͔)e06",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+x͔þ͕-x͕þ͔)e16",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+y͔þ͕-y͕þ͔)e26",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+z͔þ͕-z͕þ͔)e36",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(+ð͔þ͕-ð͕þ͔)e46",
    ///     "+(+ø͔þ͕-ø͕þ͔)e56",
    /// ]);
    ///
    /// let norm_squared = Vee::single_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+vv+αα+ββ+γγ+δδ+εε+ζζ+ηη+θθ+ιι+κκ+λλ+μμ+νν+ξξ+οο",
    ///     "+2(-κο+λξ-μν)e3456",
    ///     "+2(+ηο-θξ+ιν)e2465",
    ///     "+2(-ζο+θμ-ιλ)e2356",
    ///     "+2(+ζξ-ημ+ικ)e2364",
    ///     "+2(-ζν+ηλ-θκ)e2345",
    ///     "+2(-γο+δξ-εν)e1456",
    ///     "+2(+βο-δμ+ελ)e1365",
    ///     "+2(-βξ+γμ-εκ)e1346",
    ///     "+2(+βν-γλ+δκ)e1354",
    ///     "+2(-αο+δι-εθ)e1256",
    ///     "+2(+αξ-γι+εη)e1264",
    ///     "+2(-αν+γθ-δη)e1245",
    ///     "+2(-αμ+βι-εζ)e1236",
    ///     "+2(+αλ-βθ+δζ)e1253",
    ///     "+2(-ακ+βη-γζ)e1234",
    ///     "+2(+Ðο-Øξ+Þν)e0465",
    ///     "+2(-Zο+Øμ-Þλ)e0356",
    ///     "+2(+Zξ-Ðμ+Þκ)e0364",
    ///     "+2(-Zν+Ðλ-Øκ)e0345",
    ///     "+2(+Yο-Øι+Þθ)e0265",
    ///     "+2(-Yξ+Ðι-Þη)e0246",
    ///     "+2(+Yν-Ðθ+Øη)e0254",
    ///     "+2(+Yμ-Zι+Þζ)e0263",
    ///     "+2(-Yλ+Zθ-Øζ)e0235",
    ///     "+2(+Yκ-Zη+Ðζ)e0243",
    ///     "+2(-Xο+Øε-Þδ)e0156",
    ///     "+2(+Xξ-Ðε+Þγ)e0164",
    ///     "+2(-Xν+Ðδ-Øγ)e0145",
    ///     "+2(-Xμ+Zε-Þβ)e0136",
    ///     "+2(+Xλ-Zδ+Øβ)e0153",
    ///     "+2(-Xκ+Zγ-Ðβ)e0134",
    ///     "+2(+Xι-Yε+Þα)e0162",
    ///     "+2(-Xθ+Yδ-Øα)e0125",
    ///     "+2(+Xη-Yγ+Ðα)e0142",
    ///     "+2(-Xζ+Yβ-Zα)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::volume4()
    }
    /// The multivector of simple double motor $`m_{s2} \equiv s + v^4 + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let simple_double_motor = Vee::line().lhs() * Vee::line().rhs();
    /// assert_eq!(simple_double_motor.basis_blades(), Vee::simple_double_motor().basis_blades());
    /// format_eq!(simple_double_motor, [
    ///     "+ẋ͔ẋ͕+ẏ͔ẏ͕+ż͔ż͕+ð͔̇ð͕̇+ø͔̇ø͕̇+þ͔̇þ͕̇",
    ///     "+(+ẏ͔Α͕̇-ẏ͕Α͔̇+ż͔Β͕̇-ż͕Β͔̇+ð͔̇Γ͕̇-ð͕̇Γ͔̇+ø͔̇Δ͕̇-ø͕̇Δ͔̇+þ͔̇Ε͕̇-þ͕̇Ε͔̇)e01",
    ///     "+(-ẋ͔Α͕̇+ẋ͕Α͔̇+ż͔Ζ͕̇-ż͕Ζ͔̇+ð͔̇Η͕̇-ð͕̇Η͔̇+ø͔̇Θ͕̇-ø͕̇Θ͔̇+þ͔̇Ι͕̇-þ͕̇Ι͔̇)e02",
    ///     "+(-ẋ͔Β͕̇+ẋ͕Β͔̇-ẏ͔Ζ͕̇+ẏ͕Ζ͔̇+ð͔̇Κ͕̇-ð͕̇Κ͔̇+ø͔̇Λ͕̇-ø͕̇Λ͔̇+þ͔̇Μ͕̇-þ͕̇Μ͔̇)e03",
    ///     "+(-ẋ͔Γ͕̇+ẋ͕Γ͔̇-ẏ͔Η͕̇+ẏ͕Η͔̇-ż͔Κ͕̇+ż͕Κ͔̇+ø͔̇Ν͕̇-ø͕̇Ν͔̇+þ͔̇Ξ͕̇-þ͕̇Ξ͔̇)e04",
    ///     "+(-ẋ͔Δ͕̇+ẋ͕Δ͔̇-ẏ͔Θ͕̇+ẏ͕Θ͔̇-ż͔Λ͕̇+ż͕Λ͔̇-ð͔̇Ν͕̇+ð͕̇Ν͔̇+þ͔̇Ο͕̇-þ͕̇Ο͔̇)e05",
    ///     "+(-ẋ͔Ε͕̇+ẋ͕Ε͔̇-ẏ͔Ι͕̇+ẏ͕Ι͔̇-ż͔Μ͕̇+ż͕Μ͔̇-ð͔̇Ξ͕̇+ð͕̇Ξ͔̇-ø͔̇Ο͕̇+ø͕̇Ο͔̇)e06",
    ///     "+(+ẋ͔ẏ͕-ẋ͕ẏ͔)e12",
    ///     "+(+ẋ͔ż͕-ẋ͕ż͔)e13",
    ///     "+(+ẋ͔ð͕̇-ẋ͕ð͔̇)e14",
    ///     "+(+ẋ͔ø͕̇-ẋ͕ø͔̇)e15",
    ///     "+(+ẋ͔þ͕̇-ẋ͕þ͔̇)e16",
    ///     "+(+ẏ͔ż͕-ẏ͕ż͔)e23",
    ///     "+(+ẏ͔ð͕̇-ẏ͕ð͔̇)e24",
    ///     "+(+ẏ͔ø͕̇-ẏ͕ø͔̇)e25",
    ///     "+(+ẏ͔þ͕̇-ẏ͕þ͔̇)e26",
    ///     "+(+ż͔ð͕̇-ż͕ð͔̇)e34",
    ///     "+(+ż͔ø͕̇-ż͕ø͔̇)e35",
    ///     "+(+ż͔þ͕̇-ż͕þ͔̇)e36",
    ///     "+(+ð͔̇ø͕̇-ð͕̇ø͔̇)e45",
    ///     "+(+ð͔̇þ͕̇-ð͕̇þ͔̇)e46",
    ///     "+(+ø͔̇þ͕̇-ø͕̇þ͔̇)e56",
    ///     "+(+ð͔̇Ο͕̇+ð͕̇Ο͔̇-ø͔̇Ξ͕̇-ø͕̇Ξ͔̇+þ͔̇Ν͕̇+þ͕̇Ν͔̇)e0465",
    ///     "+(-ż͔Ο͕̇-ż͕Ο͔̇+ø͔̇Μ͕̇+ø͕̇Μ͔̇-þ͔̇Λ͕̇-þ͕̇Λ͔̇)e0356",
    ///     "+(+ż͔Ξ͕̇+ż͕Ξ͔̇-ð͔̇Μ͕̇-ð͕̇Μ͔̇+þ͔̇Κ͕̇+þ͕̇Κ͔̇)e0364",
    ///     "+(-ż͔Ν͕̇-ż͕Ν͔̇+ð͔̇Λ͕̇+ð͕̇Λ͔̇-ø͔̇Κ͕̇-ø͕̇Κ͔̇)e0345",
    ///     "+(+ẏ͔Ο͕̇+ẏ͕Ο͔̇-ø͔̇Ι͕̇-ø͕̇Ι͔̇+þ͔̇Θ͕̇+þ͕̇Θ͔̇)e0265",
    ///     "+(-ẏ͔Ξ͕̇-ẏ͕Ξ͔̇+ð͔̇Ι͕̇+ð͕̇Ι͔̇-þ͔̇Η͕̇-þ͕̇Η͔̇)e0246",
    ///     "+(+ẏ͔Ν͕̇+ẏ͕Ν͔̇-ð͔̇Θ͕̇-ð͕̇Θ͔̇+ø͔̇Η͕̇+ø͕̇Η͔̇)e0254",
    ///     "+(+ẏ͔Μ͕̇+ẏ͕Μ͔̇-ż͔Ι͕̇-ż͕Ι͔̇+þ͔̇Ζ͕̇+þ͕̇Ζ͔̇)e0263",
    ///     "+(-ẏ͔Λ͕̇-ẏ͕Λ͔̇+ż͔Θ͕̇+ż͕Θ͔̇-ø͔̇Ζ͕̇-ø͕̇Ζ͔̇)e0235",
    ///     "+(+ẏ͔Κ͕̇+ẏ͕Κ͔̇-ż͔Η͕̇-ż͕Η͔̇+ð͔̇Ζ͕̇+ð͕̇Ζ͔̇)e0243",
    ///     "+(-ẋ͔Ο͕̇-ẋ͕Ο͔̇+ø͔̇Ε͕̇+ø͕̇Ε͔̇-þ͔̇Δ͕̇-þ͕̇Δ͔̇)e0156",
    ///     "+(+ẋ͔Ξ͕̇+ẋ͕Ξ͔̇-ð͔̇Ε͕̇-ð͕̇Ε͔̇+þ͔̇Γ͕̇+þ͕̇Γ͔̇)e0164",
    ///     "+(-ẋ͔Ν͕̇-ẋ͕Ν͔̇+ð͔̇Δ͕̇+ð͕̇Δ͔̇-ø͔̇Γ͕̇-ø͕̇Γ͔̇)e0145",
    ///     "+(-ẋ͔Μ͕̇-ẋ͕Μ͔̇+ż͔Ε͕̇+ż͕Ε͔̇-þ͔̇Β͕̇-þ͕̇Β͔̇)e0136",
    ///     "+(+ẋ͔Λ͕̇+ẋ͕Λ͔̇-ż͔Δ͕̇-ż͕Δ͔̇+ø͔̇Β͕̇+ø͕̇Β͔̇)e0153",
    ///     "+(-ẋ͔Κ͕̇-ẋ͕Κ͔̇+ż͔Γ͕̇+ż͕Γ͔̇-ð͔̇Β͕̇-ð͕̇Β͔̇)e0134",
    ///     "+(+ẋ͔Ι͕̇+ẋ͕Ι͔̇-ẏ͔Ε͕̇-ẏ͕Ε͔̇+þ͔̇Α͕̇+þ͕̇Α͔̇)e0162",
    ///     "+(-ẋ͔Θ͕̇-ẋ͕Θ͔̇+ẏ͔Δ͕̇+ẏ͕Δ͔̇-ø͔̇Α͕̇-ø͕̇Α͔̇)e0125",
    ///     "+(+ẋ͔Η͕̇+ẋ͕Η͔̇-ẏ͔Γ͕̇-ẏ͕Γ͔̇+ð͔̇Α͕̇+ð͕̇Α͔̇)e0142",
    ///     "+(-ẋ͔Ζ͕̇-ẋ͕Ζ͔̇+ẏ͔Β͕̇+ẏ͕Β͔̇-ż͔Α͕̇-ż͕Α͔̇)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_motor() -> Self {
        Self::scalar() + Self::volume4() + Self::plane_moment()
    }
    /// The multivector of double motor $`m_2 \equiv s + v^4 + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(),
    ///     (Vee::double_motor() + Vee::direction()).basis_blades());
    /// format_eq!(double_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕ζ͔+Ð͕η͔+Ø͕θ͔+Þ͕ι͔)e02",
    ///     "+(-X͕β͔-Y͕ζ͔+Z͕v͔+Ð͕κ͔+Ø͕λ͔+Þ͕μ͔)e03",
    ///     "+(-X͕γ͔-Y͕η͔-Z͕κ͔+v͔Ð͕+Ø͕ν͔+Þ͕ξ͔)e04",
    ///     "+(-X͕δ͔-Y͕θ͔-Z͕λ͔+v͔Ø͕-Ð͕ν͔+Þ͕ο͔)e05",
    ///     "+(-X͕ε͔-Y͕ι͔-Z͕μ͔+v͔Þ͕-Ð͕ξ͔-Ø͕ο͔)e06",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e23",
    ///     "+v͕η͔e24",
    ///     "+v͕θ͔e25",
    ///     "+v͕ι͔e26",
    ///     "+v͕κ͔e34",
    ///     "+v͕λ͔e35",
    ///     "+v͕μ͔e36",
    ///     "+v͕ν͔e45",
    ///     "+v͕ξ͔e46",
    ///     "+v͕ο͔e56",
    ///     "+v͕α͔̇e3456",
    ///     "+v͕β͔̇e2465",
    ///     "+v͕γ͔̇e2356",
    ///     "+v͕δ͔̇e2364",
    ///     "+v͕ε͔̇e2345",
    ///     "+v͕ζ͔̇e1456",
    ///     "+v͕η͔̇e1365",
    ///     "+v͕θ͔̇e1346",
    ///     "+v͕ι͔̇e1354",
    ///     "+v͕κ͔̇e1256",
    ///     "+v͕λ͔̇e1264",
    ///     "+v͕μ͔̇e1245",
    ///     "+v͕ν͔̇e1236",
    ///     "+v͕ξ͔̇e1253",
    ///     "+v͕ο͔̇e1234",
    ///     "+(+X͕ζ͔̇-Y͕β͔̇+Z͕α͔̇-Ð͕ο͔+Ø͕ξ͔-Þ͕ν͔)e0465",
    ///     "+(+X͕η͔̇-Y͕γ͔̇+Z͕ο͔+Ð͕α͔̇-Ø͕μ͔+Þ͕λ͔)e0356",
    ///     "+(+X͕θ͔̇-Y͕δ͔̇-Z͕ξ͔+Ð͕μ͔+Ø͕α͔̇-Þ͕κ͔)e0364",
    ///     "+(+X͕ι͔̇-Y͕ε͔̇+Z͕ν͔-Ð͕λ͔+Ø͕κ͔+Þ͕α͔̇)e0345",
    ///     "+(+X͕κ͔̇-Y͕ο͔-Z͕γ͔̇+Ð͕β͔̇+Ø͕ι͔-Þ͕θ͔)e0265",
    ///     "+(+X͕λ͔̇+Y͕ξ͔-Z͕δ͔̇-Ð͕ι͔+Ø͕β͔̇+Þ͕η͔)e0246",
    ///     "+(+X͕μ͔̇-Y͕ν͔-Z͕ε͔̇+Ð͕θ͔-Ø͕η͔+Þ͕β͔̇)e0254",
    ///     "+(+X͕ν͔̇-Y͕μ͔+Z͕ι͔-Ð͕δ͔̇+Ø͕γ͔̇-Þ͕ζ͔)e0263",
    ///     "+(+X͕ξ͔̇+Y͕λ͔-Z͕θ͔-Ð͕ε͔̇+Ø͕ζ͔+Þ͕γ͔̇)e0235",
    ///     "+(+X͕ο͔̇-Y͕κ͔+Z͕η͔-Ð͕ζ͔-Ø͕ε͔̇+Þ͕δ͔̇)e0243",
    ///     "+(+X͕ο͔+Y͕κ͔̇-Z͕η͔̇+Ð͕ζ͔̇-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(-X͕ξ͔+Y͕λ͔̇-Z͕θ͔̇+Ð͕ε͔+Ø͕ζ͔̇-Þ͕γ͔)e0164",
    ///     "+(+X͕ν͔+Y͕μ͔̇-Z͕ι͔̇-Ð͕δ͔+Ø͕γ͔+Þ͕ζ͔̇)e0145",
    ///     "+(+X͕μ͔+Y͕ν͔̇-Z͕ε͔-Ð͕θ͔̇+Ø͕η͔̇+Þ͕β͔)e0136",
    ///     "+(-X͕λ͔+Y͕ξ͔̇+Z͕δ͔-Ð͕ι͔̇-Ø͕β͔+Þ͕η͔̇)e0153",
    ///     "+(+X͕κ͔+Y͕ο͔̇-Z͕γ͔+Ð͕β͔-Ø͕ι͔̇+Þ͕θ͔̇)e0134",
    ///     "+(-X͕ι͔+Y͕ε͔+Z͕ν͔̇-Ð͕λ͔̇+Ø͕κ͔̇-Þ͕α͔)e0162",
    ///     "+(+X͕θ͔-Y͕δ͔+Z͕ξ͔̇-Ð͕μ͔̇+Ø͕α͔+Þ͕κ͔̇)e0125",
    ///     "+(-X͕η͔+Y͕γ͔+Z͕ο͔̇-Ð͕α͔-Ø͕μ͔̇+Þ͕λ͔̇)e0142",
    ///     "+(+X͕ζ͔-Y͕β͔+Z͕α͔+Ð͕ο͔̇-Ø͕ξ͔̇+Þ͕ν͔̇)e0123",
    ///     // Additional orthogonality condition besides `norm_squared` condition:
    ///     "+(-Y͕α͔̇-Z͕β͔̇-Ð͕γ͔̇-Ø͕δ͔̇-Þ͕ε͔̇)e023465",
    ///     "+(+X͕α͔̇-Z͕ζ͔̇-Ð͕η͔̇-Ø͕θ͔̇-Þ͕ι͔̇)e013456",
    ///     "+(+X͕β͔̇+Y͕ζ͔̇-Ð͕κ͔̇-Ø͕λ͔̇-Þ͕μ͔̇)e012465",
    ///     "+(+X͕γ͔̇+Y͕η͔̇+Z͕κ͔̇-Ø͕ν͔̇-Þ͕ξ͔̇)e012356",
    ///     "+(+X͕δ͔̇+Y͕θ͔̇+Z͕λ͔̇+Ð͕ν͔̇-Þ͕ο͔̇)e012364",
    ///     "+(+X͕ε͔̇+Y͕ι͔̇+Z͕μ͔̇+Ð͕ξ͔̇+Ø͕ο͔̇)e012345",
    /// ]);
    ///
    /// // Orthogonality condition between single rotator quadvector and translator bivectors.
    /// // Without this condition, it is a triple motor for arbitrary choices of double rotator and
    /// // translator.
    /// format_eq!(Vee::double_rotator().vector(4).lhs() ^ Vee::translator().vector(2).rhs(), [
    ///     "+(-Y͕α͔̇-Z͕β͔̇-Ð͕γ͔̇-Ø͕δ͔̇-Þ͕ε͔̇)e023465",
    ///     "+(+X͕α͔̇-Z͕ζ͔̇-Ð͕η͔̇-Ø͕θ͔̇-Þ͕ι͔̇)e013456",
    ///     "+(+X͕β͔̇+Y͕ζ͔̇-Ð͕κ͔̇-Ø͕λ͔̇-Þ͕μ͔̇)e012465",
    ///     "+(+X͕γ͔̇+Y͕η͔̇+Z͕κ͔̇-Ø͕ν͔̇-Þ͕ξ͔̇)e012356",
    ///     "+(+X͕δ͔̇+Y͕θ͔̇+Z͕λ͔̇+Ð͕ν͔̇-Þ͕ο͔̇)e012364",
    ///     "+(+X͕ε͔̇+Y͕ι͔̇+Z͕μ͔̇+Ð͕ξ͔̇+Ø͕ο͔̇)e012345",
    /// ]);
    ///
    /// let norm_squared = Vee::double_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+vv+αα+α̇α̇+ββ+β̇β̇+γγ+γ̇γ̇+δδ+δ̇δ̇+εε+ε̇ε̇+ζζ+ζ̇ζ̇+ηη+η̇η̇\
    ///      +θθ+θ̇θ̇+ιι+ι̇ι̇+κκ+κ̇κ̇+λλ+λ̇λ̇+μμ+μ̇μ̇+νν+ν̇ν̇+ξξ+ξ̇ξ̇+οο+ο̇ο̇",
    ///     "+2(+vα̇-βζ̇+β̇ζ-γη̇+γ̇η-δθ̇+δ̇θ-ει̇+ε̇ι-κο-κ̇ο̇+λξ+λ̇ξ̇-μν-μ̇ν̇)e3456",
    ///     "+2(+vβ̇+αζ̇-α̇ζ-γκ̇+γ̇κ-δλ̇+δ̇λ-εμ̇+ε̇μ+ηο+η̇ο̇-θξ-θ̇ξ̇+ιν+ι̇ν̇)e2465",
    ///     "+2(+vγ̇+αη̇-α̇η+βκ̇-β̇κ-δν̇+δ̇ν-εξ̇+ε̇ξ-ζο-ζ̇ο̇+θμ+θ̇μ̇-ιλ-ι̇λ̇)e2356",
    ///     "+2(+vδ̇+αθ̇-α̇θ+βλ̇-β̇λ+γν̇-γ̇ν-εο̇+ε̇ο+ζξ+ζ̇ξ̇-ημ-η̇μ̇+ικ+ι̇κ̇)e2364",
    ///     "+2(+vε̇+αι̇-α̇ι+βμ̇-β̇μ+γξ̇-γ̇ξ+δο̇-δ̇ο-ζν-ζ̇ν̇+ηλ+η̇λ̇-θκ-θ̇κ̇)e2345",
    ///     "+2(+vζ̇-αβ̇+α̇β-γο-γ̇ο̇+δξ+δ̇ξ̇-εν-ε̇ν̇-ηκ̇+η̇κ-θλ̇+θ̇λ-ιμ̇+ι̇μ)e1456",
    ///     "+2(+vη̇-αγ̇+α̇γ+βο+β̇ο̇-δμ-δ̇μ̇+ελ+ε̇λ̇+ζκ̇-ζ̇κ-θν̇+θ̇ν-ιξ̇+ι̇ξ)e1365",
    ///     "+2(+vθ̇-αδ̇+α̇δ-βξ-β̇ξ̇+γμ+γ̇μ̇-εκ-ε̇κ̇+ζλ̇-ζ̇λ+ην̇-η̇ν-ιο̇+ι̇ο)e1346",
    ///     "+2(+vι̇-αε̇+α̇ε+βν+β̇ν̇-γλ-γ̇λ̇+δκ+δ̇κ̇+ζμ̇-ζ̇μ+ηξ̇-η̇ξ+θο̇-θ̇ο)e1354",
    ///     "+2(+vκ̇-αο-α̇ο̇-βγ̇+β̇γ+δι+δ̇ι̇-εθ-ε̇θ̇-ζη̇+ζ̇η-λν̇+λ̇ν-μξ̇+μ̇ξ)e1256",
    ///     "+2(+vλ̇+αξ+α̇ξ̇-βδ̇+β̇δ-γι-γ̇ι̇+εη+ε̇η̇-ζθ̇+ζ̇θ+κν̇-κ̇ν-μο̇+μ̇ο)e1264",
    ///     "+2(+vμ̇-αν-α̇ν̇-βε̇+β̇ε+γθ+γ̇θ̇-δη-δ̇η̇-ζι̇+ζ̇ι+κξ̇-κ̇ξ+λο̇-λ̇ο)e1245",
    ///     "+2(+vν̇-αμ-α̇μ̇+βι+β̇ι̇-γδ̇+γ̇δ-εζ-ε̇ζ̇-ηθ̇+η̇θ-κλ̇+κ̇λ-ξο̇+ξ̇ο)e1236",
    ///     "+2(+vξ̇+αλ+α̇λ̇-βθ-β̇θ̇-γε̇+γ̇ε+δζ+δ̇ζ̇-ηι̇+η̇ι-κμ̇+κ̇μ+νο̇-ν̇ο)e1253",
    ///     "+2(+vο̇-ακ-α̇κ̇+βη+β̇η̇-γζ-γ̇ζ̇-δε̇+δ̇ε-θι̇+θ̇ι-λμ̇+λ̇μ-νξ̇+ν̇ξ)e1234",
    ///     "+2(+Ȧv+Ḃκ+Ċλ+Ḋμ-Ėη-Ḟθ-Ġι-Ḣε̇+İδ̇-J̇γ̇+K̇γ+L̇δ\
    ///         +Ṁε-Ṅι̇+Ȯθ̇-Ṗη̇-Q̇μ̇+Ṙλ̇-Ṡκ̇-Xζ̇+Yβ̇-Zα̇+Ðο-Øξ+Þν)e0465",
    ///     "+2(-Ȧκ+Ḃv+Ċν+Ḋξ+Ėζ+Ḟε̇-Ġδ̇-Ḣθ-İι+J̇β̇-K̇β+L̇ι̇\
    ///         -Ṁθ̇+Ṅδ+Ȯε+Ṗζ̇-Q̇ξ̇+Ṙν̇-Ṫκ̇-Xη̇+Yγ̇-Zο-Ðα̇+Øμ-Þλ)e0356",
    ///     "+2(-Ȧλ-Ḃν+Ċv+Ḋο-Ėε̇+Ḟζ+Ġγ̇+Ḣη-İβ̇-J̇ι-K̇ι̇-L̇β\
    ///         +Ṁη̇-Ṅγ-Ȯζ̇+Ṗε-Q̇ο̇+Ṡν̇-Ṫλ̇-Xθ̇+Yδ̇+Zξ-Ðμ-Øα̇+Þκ)e0364",
    ///     "+2(-Ȧμ-Ḃξ-Ċο+Ḋv+Ėδ̇-Ḟγ̇+Ġζ+Ḣβ̇+İη+J̇θ+K̇θ̇-L̇η̇\
    ///         -Ṁβ+Ṅζ̇-Ȯγ-Ṗδ-Ṙο̇+Ṡξ̇-Ṫμ̇-Xι̇+Yε̇-Zν+Ðλ-Øκ-Þα̇)e0345",
    ///     "+2(+Ȧη-Ḃζ-Ċε̇+Ḋδ̇+Ėv+Ḟν+Ġξ-Ḣλ-İμ-J̇α̇+K̇α+L̇μ̇\
    ///         -Ṁλ̇+Ṅξ̇-Ȯν̇+Q̇δ+Ṙε+Ṡζ̇+Ṫη̇-Xκ̇+Yο+Zγ̇-Ðβ̇-Øι+Þθ)e0265",
    ///     "+2(+Ȧθ+Ḃε̇-Ċζ-Ḋγ̇-Ėν+Ḟv+Ġο+Ḣκ+İα̇-J̇μ-K̇μ̇+L̇α\
    ///         +Ṁκ̇+Ṅο̇-Ṗν̇-Q̇γ-Ṙζ̇+Ṡε+Ṫθ̇-Xλ̇-Yξ+Zδ̇+Ðι-Øβ̇-Þη)e0246",
    ///     "+2(+Ȧι-Ḃδ̇+Ċγ̇-Ḋζ-Ėξ-Ḟο+Ġv-Ḣα̇+İκ+J̇λ+K̇λ̇-L̇κ̇\
    ///         +Ṁα+Ȯο̇-Ṗξ̇+Q̇ζ̇-Ṙγ-Ṡδ+Ṫι̇-Xμ̇+Yν+Zε̇-Ðθ+Øη-Þβ̇)e0254",
    ///     "+2(-Ȧε̇+Ḃθ-Ċη+Ḋβ̇+Ėλ-Ḟκ-Ġα̇+Ḣv+İο-J̇ξ-K̇ξ̇-L̇ο̇\
    ///         +Ṅα+Ȯκ̇+Ṗλ̇+Q̇β-Ṙη̇-Ṡθ̇+Ṫε-Xν̇+Yμ-Zι+Ðδ̇-Øγ̇+Þζ)e0263",
    ///     "+2(+Ȧδ̇+Ḃι-Ċβ̇-Ḋη+Ėμ+Ḟα̇-Ġκ-Ḣο+İv+J̇ν+K̇ν̇-Ṁο̇\
    ///         -Ṅκ̇+Ȯα+Ṗμ̇+Q̇η̇+Ṙβ-Ṡι̇-Ṫδ-Xξ̇-Yλ+Zθ+Ðε̇-Øζ-Þγ̇)e0235",
    ///     "+2(-Ȧγ̇+Ḃβ̇+Ċι-Ḋθ-Ėα̇+Ḟμ-Ġλ+Ḣξ-İν+J̇v+L̇ν̇+Ṁξ̇\
    ///         -Ṅλ̇-Ȯμ̇+Ṗα+Q̇θ̇+Ṙι̇+Ṡβ+Ṫγ-Xο̇+Yκ-Zη+Ðζ+Øε̇-Þδ̇)e0243",
    ///     "+2(-Ȧγ+Ḃβ-Ċι̇+Ḋθ̇-Ėα-Ḟμ̇+Ġλ̇-Ḣξ̇+İν̇+K̇v+L̇ν+Ṁξ\
    ///         -Ṅλ-Ȯμ-Ṗα̇+Q̇θ+Ṙι-Ṡβ̇-Ṫγ̇-Xο-Yκ̇+Zη̇-Ðζ̇+Øε-Þδ)e0156",
    ///     "+2(-Ȧδ+Ḃι̇+Ċβ-Ḋη̇+Ėμ̇-Ḟα-Ġκ̇-Ḣο̇+J̇ν̇-K̇ν+L̇v+Ṁο\
    ///         +Ṅκ+Ȯα̇-Ṗμ-Q̇η+Ṙβ̇+Ṡι-Ṫδ̇+Xξ-Yλ̇+Zθ̇-Ðε-Øζ̇+Þγ)e0164",
    ///     "+2(-Ȧε-Ḃθ̇+Ċη̇+Ḋβ-Ėλ̇+Ḟκ̇-Ġα-İο̇+J̇ξ̇-K̇ξ-L̇ο+Ṁv\
    ///         -Ṅα̇+Ȯκ+Ṗλ-Q̇β̇-Ṙη-Ṡθ-Ṫε̇-Xν-Yμ̇+Zι̇+Ðδ-Øγ-Þζ̇)e0145",
    ///     "+2(-Ȧι̇-Ḃδ+Ċγ+Ḋζ̇+Ėξ̇+Ḟο̇-Ḣα-İκ̇-J̇λ̇+K̇λ-L̇κ-Ṁα̇\
    ///         +Ṅv+Ȯο-Ṗξ+Q̇ζ+Ṙγ̇+Ṡδ̇+Ṫι-Xμ-Yν̇+Zε+Ðθ̇-Øη̇-Þβ)e0136",
    ///     "+2(+Ȧθ̇-Ḃε-Ċζ̇+Ḋγ-Ėν̇+Ġο̇+Ḣκ̇-İα-J̇μ̇+K̇μ+L̇α̇-Ṁκ\
    ///         -Ṅο+Ȯv+Ṗν-Q̇γ̇+Ṙζ+Ṡε̇-Ṫθ+Xλ-Yξ̇-Zδ+Ðι̇+Øβ-Þη̇)e0153",
    ///     "+2(-Ȧη̇+Ḃζ̇-Ċε+Ḋδ-Ḟν̇-Ġξ̇+Ḣλ̇+İμ̇-J̇α-K̇α̇+L̇μ-Ṁλ\
    ///         +Ṅξ-Ȯν+Ṗv-Q̇δ̇-Ṙε̇+Ṡζ+Ṫη-Xκ-Yο̇+Zγ-Ðβ+Øι̇-Þθ̇)e0134",
    ///     "+2(-Ȧμ̇-Ḃξ̇-Ċο̇-Ėδ+Ḟγ+Ġζ̇-Ḣβ+İη̇+J̇θ̇-K̇θ+L̇η-Ṁβ̇\
    ///         -Ṅζ-Ȯγ̇-Ṗδ̇+Q̇v+Ṙο-Ṡξ+Ṫμ+Xι-Yε-Zν̇+Ðλ̇-Øκ̇+Þα)e0162",
    ///     "+2(+Ȧλ̇+Ḃν̇-Ḋο̇-Ėε-Ḟζ̇+Ġγ-Ḣη̇-İβ+J̇ι̇-K̇ι+L̇β̇+Ṁη\
    ///         +Ṅγ̇-Ȯζ-Ṗε̇-Q̇ο+Ṙv+Ṡν-Ṫλ-Xθ+Yδ-Zξ̇+Ðμ̇-Øα-Þκ̇)e0125",
    ///     "+2(-Ȧκ̇+Ċν̇+Ḋξ̇+Ėζ̇-Ḟε+Ġδ-Ḣθ̇-İι̇-J̇β-K̇β̇-L̇ι+Ṁθ\
    ///         +Ṅδ̇+Ȯε̇-Ṗζ+Q̇ξ-Ṙν+Ṡv+Ṫκ+Xη-Yγ-Zο̇+Ðα+Øμ̇-Þλ̇)e0142",
    ///     "+2(-Ḃκ̇-Ċλ̇-Ḋμ̇+Ėη̇+Ḟθ̇+Ġι̇-Ḣε+İδ-J̇γ-K̇γ̇-L̇δ̇-Ṁε̇\
    ///         -Ṅι+Ȯθ-Ṗη-Q̇μ+Ṙλ-Ṡκ+Ṫv-Xζ+Yβ-Zα-Ðο̇+Øξ̇-Þν̇)e0123",
    /// ]);
    ///
    /// let double_motor = Vee::volume4().lhs() * Vee::volume4().rhs();
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕",
    ///     "+(-Y͔α͕+Y͕α͔-Z͔β͕+Z͕β͔-Ð͔γ͕+Ð͕γ͔-Ø͔δ͕+Ø͕δ͔-Þ͔ε͕+Þ͕ε͔)e01",
    ///     "+(+X͔α͕-X͕α͔-Z͔ζ͕+Z͕ζ͔-Ð͔η͕+Ð͕η͔-Ø͔θ͕+Ø͕θ͔-Þ͔ι͕+Þ͕ι͔)e02",
    ///     "+(+X͔β͕-X͕β͔+Y͔ζ͕-Y͕ζ͔-Ð͔κ͕+Ð͕κ͔-Ø͔λ͕+Ø͕λ͔-Þ͔μ͕+Þ͕μ͔)e03",
    ///     "+(+X͔γ͕-X͕γ͔+Y͔η͕-Y͕η͔+Z͔κ͕-Z͕κ͔-Ø͔ν͕+Ø͕ν͔-Þ͔ξ͕+Þ͕ξ͔)e04",
    ///     "+(+X͔δ͕-X͕δ͔+Y͔θ͕-Y͕θ͔+Z͔λ͕-Z͕λ͔+Ð͔ν͕-Ð͕ν͔-Þ͔ο͕+Þ͕ο͔)e05",
    ///     "+(+X͔ε͕-X͕ε͔+Y͔ι͕-Y͕ι͔+Z͔μ͕-Z͕μ͔+Ð͔ξ͕-Ð͕ξ͔+Ø͔ο͕-Ø͕ο͔)e06",
    ///     "+(-β͔ζ͕+β͕ζ͔-γ͔η͕+γ͕η͔-δ͔θ͕+δ͕θ͔-ε͔ι͕+ε͕ι͔)e12",
    ///     "+(+α͔ζ͕-α͕ζ͔-γ͔κ͕+γ͕κ͔-δ͔λ͕+δ͕λ͔-ε͔μ͕+ε͕μ͔)e13",
    ///     "+(+α͔η͕-α͕η͔+β͔κ͕-β͕κ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔)e14",
    ///     "+(+α͔θ͕-α͕θ͔+β͔λ͕-β͕λ͔+γ͔ν͕-γ͕ν͔-ε͔ο͕+ε͕ο͔)e15",
    ///     "+(+α͔ι͕-α͕ι͔+β͔μ͕-β͕μ͔+γ͔ξ͕-γ͕ξ͔+δ͔ο͕-δ͕ο͔)e16",
    ///     "+(-α͔β͕+α͕β͔-η͔κ͕+η͕κ͔-θ͔λ͕+θ͕λ͔-ι͔μ͕+ι͕μ͔)e23",
    ///     "+(-α͔γ͕+α͕γ͔+ζ͔κ͕-ζ͕κ͔-θ͔ν͕+θ͕ν͔-ι͔ξ͕+ι͕ξ͔)e24",
    ///     "+(-α͔δ͕+α͕δ͔+ζ͔λ͕-ζ͕λ͔+η͔ν͕-η͕ν͔-ι͔ο͕+ι͕ο͔)e25",
    ///     "+(-α͔ε͕+α͕ε͔+ζ͔μ͕-ζ͕μ͔+η͔ξ͕-η͕ξ͔+θ͔ο͕-θ͕ο͔)e26",
    ///     "+(-β͔γ͕+β͕γ͔-ζ͔η͕+ζ͕η͔-λ͔ν͕+λ͕ν͔-μ͔ξ͕+μ͕ξ͔)e34",
    ///     "+(-β͔δ͕+β͕δ͔-ζ͔θ͕+ζ͕θ͔+κ͔ν͕-κ͕ν͔-μ͔ο͕+μ͕ο͔)e35",
    ///     "+(-β͔ε͕+β͕ε͔-ζ͔ι͕+ζ͕ι͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e36",
    ///     "+(-γ͔δ͕+γ͕δ͔-η͔θ͕+η͕θ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔)e45",
    ///     "+(-γ͔ε͕+γ͕ε͔-η͔ι͕+η͕ι͔-κ͔μ͕+κ͕μ͔+ν͔ο͕-ν͕ο͔)e46",
    ///     "+(-δ͔ε͕+δ͕ε͔-θ͔ι͕+θ͕ι͔-λ͔μ͕+λ͕μ͔-ν͔ξ͕+ν͕ξ͔)e56",
    ///     "+(+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔+μ͔ν͕+μ͕ν͔)e3456",
    ///     "+(-η͔ο͕-η͕ο͔+θ͔ξ͕+θ͕ξ͔-ι͔ν͕-ι͕ν͔)e2465",
    ///     "+(+ζ͔ο͕+ζ͕ο͔-θ͔μ͕-θ͕μ͔+ι͔λ͕+ι͕λ͔)e2356",
    ///     "+(-ζ͔ξ͕-ζ͕ξ͔+η͔μ͕+η͕μ͔-ι͔κ͕-ι͕κ͔)e2364",
    ///     "+(+ζ͔ν͕+ζ͕ν͔-η͔λ͕-η͕λ͔+θ͔κ͕+θ͕κ͔)e2345",
    ///     "+(+γ͔ο͕+γ͕ο͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔)e1456",
    ///     "+(-β͔ο͕-β͕ο͔+δ͔μ͕+δ͕μ͔-ε͔λ͕-ε͕λ͔)e1365",
    ///     "+(+β͔ξ͕+β͕ξ͔-γ͔μ͕-γ͕μ͔+ε͔κ͕+ε͕κ͔)e1346",
    ///     "+(-β͔ν͕-β͕ν͔+γ͔λ͕+γ͕λ͔-δ͔κ͕-δ͕κ͔)e1354",
    ///     "+(+α͔ο͕+α͕ο͔-δ͔ι͕-δ͕ι͔+ε͔θ͕+ε͕θ͔)e1256",
    ///     "+(-α͔ξ͕-α͕ξ͔+γ͔ι͕+γ͕ι͔-ε͔η͕-ε͕η͔)e1264",
    ///     "+(+α͔ν͕+α͕ν͔-γ͔θ͕-γ͕θ͔+δ͔η͕+δ͕η͔)e1245",
    ///     "+(+α͔μ͕+α͕μ͔-β͔ι͕-β͕ι͔+ε͔ζ͕+ε͕ζ͔)e1236",
    ///     "+(-α͔λ͕-α͕λ͔+β͔θ͕+β͕θ͔-δ͔ζ͕-δ͕ζ͔)e1253",
    ///     "+(+α͔κ͕+α͕κ͔-β͔η͕-β͕η͔+γ͔ζ͕+γ͕ζ͔)e1234",
    ///     "+(-Ð͔ο͕-Ð͕ο͔+Ø͔ξ͕+Ø͕ξ͔-Þ͔ν͕-Þ͕ν͔)e0465",
    ///     "+(+Z͔ο͕+Z͕ο͔-Ø͔μ͕-Ø͕μ͔+Þ͔λ͕+Þ͕λ͔)e0356",
    ///     "+(-Z͔ξ͕-Z͕ξ͔+Ð͔μ͕+Ð͕μ͔-Þ͔κ͕-Þ͕κ͔)e0364",
    ///     "+(+Z͔ν͕+Z͕ν͔-Ð͔λ͕-Ð͕λ͔+Ø͔κ͕+Ø͕κ͔)e0345",
    ///     "+(-Y͔ο͕-Y͕ο͔+Ø͔ι͕+Ø͕ι͔-Þ͔θ͕-Þ͕θ͔)e0265",
    ///     "+(+Y͔ξ͕+Y͕ξ͔-Ð͔ι͕-Ð͕ι͔+Þ͔η͕+Þ͕η͔)e0246",
    ///     "+(-Y͔ν͕-Y͕ν͔+Ð͔θ͕+Ð͕θ͔-Ø͔η͕-Ø͕η͔)e0254",
    ///     "+(-Y͔μ͕-Y͕μ͔+Z͔ι͕+Z͕ι͔-Þ͔ζ͕-Þ͕ζ͔)e0263",
    ///     "+(+Y͔λ͕+Y͕λ͔-Z͔θ͕-Z͕θ͔+Ø͔ζ͕+Ø͕ζ͔)e0235",
    ///     "+(-Y͔κ͕-Y͕κ͔+Z͔η͕+Z͕η͔-Ð͔ζ͕-Ð͕ζ͔)e0243",
    ///     "+(+X͔ο͕+X͕ο͔-Ø͔ε͕-Ø͕ε͔+Þ͔δ͕+Þ͕δ͔)e0156",
    ///     "+(-X͔ξ͕-X͕ξ͔+Ð͔ε͕+Ð͕ε͔-Þ͔γ͕-Þ͕γ͔)e0164",
    ///     "+(+X͔ν͕+X͕ν͔-Ð͔δ͕-Ð͕δ͔+Ø͔γ͕+Ø͕γ͔)e0145",
    ///     "+(+X͔μ͕+X͕μ͔-Z͔ε͕-Z͕ε͔+Þ͔β͕+Þ͕β͔)e0136",
    ///     "+(-X͔λ͕-X͕λ͔+Z͔δ͕+Z͕δ͔-Ø͔β͕-Ø͕β͔)e0153",
    ///     "+(+X͔κ͕+X͕κ͔-Z͔γ͕-Z͕γ͔+Ð͔β͕+Ð͕β͔)e0134",
    ///     "+(-X͔ι͕-X͕ι͔+Y͔ε͕+Y͕ε͔-Þ͔α͕-Þ͕α͔)e0162",
    ///     "+(+X͔θ͕+X͕θ͔-Y͔δ͕-Y͕δ͔+Ø͔α͕+Ø͕α͔)e0125",
    ///     "+(-X͔η͕-X͕η͔+Y͔γ͕+Y͕γ͔-Ð͔α͕-Ð͕α͔)e0142",
    ///     "+(+X͔ζ͕+X͕ζ͔-Y͔β͕-Y͕β͔+Z͔α͕+Z͕α͔)e0123",
    /// ]);
    ///
    /// let norm_squared = Vee::double_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+vv+αα+α̇α̇+ββ+β̇β̇+γγ+γ̇γ̇+δδ+δ̇δ̇+εε+ε̇ε̇+ζζ+ζ̇ζ̇+ηη+η̇η̇\
    ///      +θθ+θ̇θ̇+ιι+ι̇ι̇+κκ+κ̇κ̇+λλ+λ̇λ̇+μμ+μ̇μ̇+νν+ν̇ν̇+ξξ+ξ̇ξ̇+οο+ο̇ο̇",
    ///     "+2(+vα̇-βζ̇+β̇ζ-γη̇+γ̇η-δθ̇+δ̇θ-ει̇+ε̇ι-κο-κ̇ο̇+λξ+λ̇ξ̇-μν-μ̇ν̇)e3456",
    ///     "+2(+vβ̇+αζ̇-α̇ζ-γκ̇+γ̇κ-δλ̇+δ̇λ-εμ̇+ε̇μ+ηο+η̇ο̇-θξ-θ̇ξ̇+ιν+ι̇ν̇)e2465",
    ///     "+2(+vγ̇+αη̇-α̇η+βκ̇-β̇κ-δν̇+δ̇ν-εξ̇+ε̇ξ-ζο-ζ̇ο̇+θμ+θ̇μ̇-ιλ-ι̇λ̇)e2356",
    ///     "+2(+vδ̇+αθ̇-α̇θ+βλ̇-β̇λ+γν̇-γ̇ν-εο̇+ε̇ο+ζξ+ζ̇ξ̇-ημ-η̇μ̇+ικ+ι̇κ̇)e2364",
    ///     "+2(+vε̇+αι̇-α̇ι+βμ̇-β̇μ+γξ̇-γ̇ξ+δο̇-δ̇ο-ζν-ζ̇ν̇+ηλ+η̇λ̇-θκ-θ̇κ̇)e2345",
    ///     "+2(+vζ̇-αβ̇+α̇β-γο-γ̇ο̇+δξ+δ̇ξ̇-εν-ε̇ν̇-ηκ̇+η̇κ-θλ̇+θ̇λ-ιμ̇+ι̇μ)e1456",
    ///     "+2(+vη̇-αγ̇+α̇γ+βο+β̇ο̇-δμ-δ̇μ̇+ελ+ε̇λ̇+ζκ̇-ζ̇κ-θν̇+θ̇ν-ιξ̇+ι̇ξ)e1365",
    ///     "+2(+vθ̇-αδ̇+α̇δ-βξ-β̇ξ̇+γμ+γ̇μ̇-εκ-ε̇κ̇+ζλ̇-ζ̇λ+ην̇-η̇ν-ιο̇+ι̇ο)e1346",
    ///     "+2(+vι̇-αε̇+α̇ε+βν+β̇ν̇-γλ-γ̇λ̇+δκ+δ̇κ̇+ζμ̇-ζ̇μ+ηξ̇-η̇ξ+θο̇-θ̇ο)e1354",
    ///     "+2(+vκ̇-αο-α̇ο̇-βγ̇+β̇γ+δι+δ̇ι̇-εθ-ε̇θ̇-ζη̇+ζ̇η-λν̇+λ̇ν-μξ̇+μ̇ξ)e1256",
    ///     "+2(+vλ̇+αξ+α̇ξ̇-βδ̇+β̇δ-γι-γ̇ι̇+εη+ε̇η̇-ζθ̇+ζ̇θ+κν̇-κ̇ν-μο̇+μ̇ο)e1264",
    ///     "+2(+vμ̇-αν-α̇ν̇-βε̇+β̇ε+γθ+γ̇θ̇-δη-δ̇η̇-ζι̇+ζ̇ι+κξ̇-κ̇ξ+λο̇-λ̇ο)e1245",
    ///     "+2(+vν̇-αμ-α̇μ̇+βι+β̇ι̇-γδ̇+γ̇δ-εζ-ε̇ζ̇-ηθ̇+η̇θ-κλ̇+κ̇λ-ξο̇+ξ̇ο)e1236",
    ///     "+2(+vξ̇+αλ+α̇λ̇-βθ-β̇θ̇-γε̇+γ̇ε+δζ+δ̇ζ̇-ηι̇+η̇ι-κμ̇+κ̇μ+νο̇-ν̇ο)e1253",
    ///     "+2(+vο̇-ακ-α̇κ̇+βη+β̇η̇-γζ-γ̇ζ̇-δε̇+δ̇ε-θι̇+θ̇ι-λμ̇+λ̇μ-νξ̇+ν̇ξ)e1234",
    ///     "+2(+Ȧv+Ḃκ+Ċλ+Ḋμ-Ėη-Ḟθ-Ġι-Ḣε̇+İδ̇-J̇γ̇+K̇γ+L̇δ\
    ///         +Ṁε-Ṅι̇+Ȯθ̇-Ṗη̇-Q̇μ̇+Ṙλ̇-Ṡκ̇-Xζ̇+Yβ̇-Zα̇+Ðο-Øξ+Þν)e0465",
    ///     "+2(-Ȧκ+Ḃv+Ċν+Ḋξ+Ėζ+Ḟε̇-Ġδ̇-Ḣθ-İι+J̇β̇-K̇β+L̇ι̇\
    ///         -Ṁθ̇+Ṅδ+Ȯε+Ṗζ̇-Q̇ξ̇+Ṙν̇-Ṫκ̇-Xη̇+Yγ̇-Zο-Ðα̇+Øμ-Þλ)e0356",
    ///     "+2(-Ȧλ-Ḃν+Ċv+Ḋο-Ėε̇+Ḟζ+Ġγ̇+Ḣη-İβ̇-J̇ι-K̇ι̇-L̇β\
    ///         +Ṁη̇-Ṅγ-Ȯζ̇+Ṗε-Q̇ο̇+Ṡν̇-Ṫλ̇-Xθ̇+Yδ̇+Zξ-Ðμ-Øα̇+Þκ)e0364",
    ///     "+2(-Ȧμ-Ḃξ-Ċο+Ḋv+Ėδ̇-Ḟγ̇+Ġζ+Ḣβ̇+İη+J̇θ+K̇θ̇-L̇η̇\
    ///         -Ṁβ+Ṅζ̇-Ȯγ-Ṗδ-Ṙο̇+Ṡξ̇-Ṫμ̇-Xι̇+Yε̇-Zν+Ðλ-Øκ-Þα̇)e0345",
    ///     "+2(+Ȧη-Ḃζ-Ċε̇+Ḋδ̇+Ėv+Ḟν+Ġξ-Ḣλ-İμ-J̇α̇+K̇α+L̇μ̇\
    ///         -Ṁλ̇+Ṅξ̇-Ȯν̇+Q̇δ+Ṙε+Ṡζ̇+Ṫη̇-Xκ̇+Yο+Zγ̇-Ðβ̇-Øι+Þθ)e0265",
    ///     "+2(+Ȧθ+Ḃε̇-Ċζ-Ḋγ̇-Ėν+Ḟv+Ġο+Ḣκ+İα̇-J̇μ-K̇μ̇+L̇α\
    ///         +Ṁκ̇+Ṅο̇-Ṗν̇-Q̇γ-Ṙζ̇+Ṡε+Ṫθ̇-Xλ̇-Yξ+Zδ̇+Ðι-Øβ̇-Þη)e0246",
    ///     "+2(+Ȧι-Ḃδ̇+Ċγ̇-Ḋζ-Ėξ-Ḟο+Ġv-Ḣα̇+İκ+J̇λ+K̇λ̇-L̇κ̇\
    ///         +Ṁα+Ȯο̇-Ṗξ̇+Q̇ζ̇-Ṙγ-Ṡδ+Ṫι̇-Xμ̇+Yν+Zε̇-Ðθ+Øη-Þβ̇)e0254",
    ///     "+2(-Ȧε̇+Ḃθ-Ċη+Ḋβ̇+Ėλ-Ḟκ-Ġα̇+Ḣv+İο-J̇ξ-K̇ξ̇-L̇ο̇\
    ///         +Ṅα+Ȯκ̇+Ṗλ̇+Q̇β-Ṙη̇-Ṡθ̇+Ṫε-Xν̇+Yμ-Zι+Ðδ̇-Øγ̇+Þζ)e0263",
    ///     "+2(+Ȧδ̇+Ḃι-Ċβ̇-Ḋη+Ėμ+Ḟα̇-Ġκ-Ḣο+İv+J̇ν+K̇ν̇-Ṁο̇\
    ///         -Ṅκ̇+Ȯα+Ṗμ̇+Q̇η̇+Ṙβ-Ṡι̇-Ṫδ-Xξ̇-Yλ+Zθ+Ðε̇-Øζ-Þγ̇)e0235",
    ///     "+2(-Ȧγ̇+Ḃβ̇+Ċι-Ḋθ-Ėα̇+Ḟμ-Ġλ+Ḣξ-İν+J̇v+L̇ν̇+Ṁξ̇\
    ///         -Ṅλ̇-Ȯμ̇+Ṗα+Q̇θ̇+Ṙι̇+Ṡβ+Ṫγ-Xο̇+Yκ-Zη+Ðζ+Øε̇-Þδ̇)e0243",
    ///     "+2(-Ȧγ+Ḃβ-Ċι̇+Ḋθ̇-Ėα-Ḟμ̇+Ġλ̇-Ḣξ̇+İν̇+K̇v+L̇ν+Ṁξ\
    ///         -Ṅλ-Ȯμ-Ṗα̇+Q̇θ+Ṙι-Ṡβ̇-Ṫγ̇-Xο-Yκ̇+Zη̇-Ðζ̇+Øε-Þδ)e0156",
    ///     "+2(-Ȧδ+Ḃι̇+Ċβ-Ḋη̇+Ėμ̇-Ḟα-Ġκ̇-Ḣο̇+J̇ν̇-K̇ν+L̇v+Ṁο\
    ///         +Ṅκ+Ȯα̇-Ṗμ-Q̇η+Ṙβ̇+Ṡι-Ṫδ̇+Xξ-Yλ̇+Zθ̇-Ðε-Øζ̇+Þγ)e0164",
    ///     "+2(-Ȧε-Ḃθ̇+Ċη̇+Ḋβ-Ėλ̇+Ḟκ̇-Ġα-İο̇+J̇ξ̇-K̇ξ-L̇ο+Ṁv\
    ///         -Ṅα̇+Ȯκ+Ṗλ-Q̇β̇-Ṙη-Ṡθ-Ṫε̇-Xν-Yμ̇+Zι̇+Ðδ-Øγ-Þζ̇)e0145",
    ///     "+2(-Ȧι̇-Ḃδ+Ċγ+Ḋζ̇+Ėξ̇+Ḟο̇-Ḣα-İκ̇-J̇λ̇+K̇λ-L̇κ-Ṁα̇\
    ///         +Ṅv+Ȯο-Ṗξ+Q̇ζ+Ṙγ̇+Ṡδ̇+Ṫι-Xμ-Yν̇+Zε+Ðθ̇-Øη̇-Þβ)e0136",
    ///     "+2(+Ȧθ̇-Ḃε-Ċζ̇+Ḋγ-Ėν̇+Ġο̇+Ḣκ̇-İα-J̇μ̇+K̇μ+L̇α̇-Ṁκ\
    ///         -Ṅο+Ȯv+Ṗν-Q̇γ̇+Ṙζ+Ṡε̇-Ṫθ+Xλ-Yξ̇-Zδ+Ðι̇+Øβ-Þη̇)e0153",
    ///     "+2(-Ȧη̇+Ḃζ̇-Ċε+Ḋδ-Ḟν̇-Ġξ̇+Ḣλ̇+İμ̇-J̇α-K̇α̇+L̇μ-Ṁλ\
    ///         +Ṅξ-Ȯν+Ṗv-Q̇δ̇-Ṙε̇+Ṡζ+Ṫη-Xκ-Yο̇+Zγ-Ðβ+Øι̇-Þθ̇)e0134",
    ///     "+2(-Ȧμ̇-Ḃξ̇-Ċο̇-Ėδ+Ḟγ+Ġζ̇-Ḣβ+İη̇+J̇θ̇-K̇θ+L̇η-Ṁβ̇\
    ///         -Ṅζ-Ȯγ̇-Ṗδ̇+Q̇v+Ṙο-Ṡξ+Ṫμ+Xι-Yε-Zν̇+Ðλ̇-Øκ̇+Þα)e0162",
    ///     "+2(+Ȧλ̇+Ḃν̇-Ḋο̇-Ėε-Ḟζ̇+Ġγ-Ḣη̇-İβ+J̇ι̇-K̇ι+L̇β̇+Ṁη\
    ///         +Ṅγ̇-Ȯζ-Ṗε̇-Q̇ο+Ṙv+Ṡν-Ṫλ-Xθ+Yδ-Zξ̇+Ðμ̇-Øα-Þκ̇)e0125",
    ///     "+2(-Ȧκ̇+Ċν̇+Ḋξ̇+Ėζ̇-Ḟε+Ġδ-Ḣθ̇-İι̇-J̇β-K̇β̇-L̇ι+Ṁθ\
    ///         +Ṅδ̇+Ȯε̇-Ṗζ+Q̇ξ-Ṙν+Ṡv+Ṫκ+Xη-Yγ-Zο̇+Ðα+Øμ̇-Þλ̇)e0142",
    ///     "+2(-Ḃκ̇-Ċλ̇-Ḋμ̇+Ėη̇+Ḟθ̇+Ġι̇-Ḣε+İδ-J̇γ-K̇γ̇-L̇δ̇-Ṁε̇\
    ///         -Ṅι+Ȯθ-Ṗη-Q̇μ+Ṙλ-Ṡκ+Ṫv-Xζ+Yβ-Zα-Ðο̇+Øξ̇-Þν̇)e0123",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume4() + Self::plane()
    }
    /// The multivector of simple triple motor $`m_{s3} \equiv s + v^4 + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let simple_triple_motor = Vee::plane().lhs() * Vee::plane().rhs();
    /// assert_eq!(simple_triple_motor.basis_blades(), Vee::simple_triple_motor().basis_blades());
    /// format_eq!(simple_triple_motor, [
    ///     "+α͔̇α͕̇+β͔̇β͕̇+γ͔̇γ͕̇+δ͔̇δ͕̇+ε͔̇ε͕̇+ζ͔̇ζ͕̇+η͔̇η͕̇+θ͔̇θ͕̇+ι͔̇ι͕̇+κ͔̇κ͕̇+λ͔̇λ͕̇+μ͔̇μ͕̇+ν͔̇ν͕̇+ξ͔̇ξ͕̇+ο͔̇ο͕̇",
    ///     "+(-Ȧ͔ζ͕̇+Ȧ͕ζ͔̇-Ḃ͔η͕̇+Ḃ͕η͔̇-Ċ͔θ͕̇+Ċ͕θ͔̇-Ḋ͔ι͕̇+Ḋ͕ι͔̇-Ė͔κ͕̇+Ė͕κ͔̇-Ḟ͔λ͕̇+Ḟ͕λ͔̇-Ġ͔μ͕̇+Ġ͕μ͔̇-Ḣ͔ν͕̇+Ḣ͕ν͔̇-İ͔ξ͕̇+İ͕ξ͔̇-J͔̇ο͕̇+J͕̇ο͔̇)e01",
    ///     "+(+Ȧ͔β͕̇-Ȧ͕β͔̇+Ḃ͔γ͕̇-Ḃ͕γ͔̇+Ċ͔δ͕̇-Ċ͕δ͔̇+Ḋ͔ε͕̇-Ḋ͕ε͔̇-K͔̇κ͕̇+K͕̇κ͔̇-L͔̇λ͕̇+L͕̇λ͔̇-Ṁ͔μ͕̇+Ṁ͕μ͔̇-Ṅ͔ν͕̇+Ṅ͕ν͔̇-Ȯ͔ξ͕̇+Ȯ͕ξ͔̇-Ṗ͔ο͕̇+Ṗ͕ο͔̇)e02",
    ///     "+(-Ȧ͔α͕̇+Ȧ͕α͔̇+Ė͔γ͕̇-Ė͕γ͔̇+Ḟ͔δ͕̇-Ḟ͕δ͔̇+Ġ͔ε͕̇-Ġ͕ε͔̇+K͔̇η͕̇-K͕̇η͔̇+L͔̇θ͕̇-L͕̇θ͔̇+Ṁ͔ι͕̇-Ṁ͕ι͔̇-Q͔̇ν͕̇+Q͕̇ν͔̇-Ṙ͔ξ͕̇+Ṙ͕ξ͔̇-Ṡ͔ο͕̇+Ṡ͕ο͔̇)e03",
    ///     "+(-Ḃ͔α͕̇+Ḃ͕α͔̇-Ė͔β͕̇+Ė͕β͔̇+Ḣ͔δ͕̇-Ḣ͕δ͔̇+İ͔ε͕̇-İ͕ε͔̇-K͔̇ζ͕̇+K͕̇ζ͔̇+Ṅ͔θ͕̇-Ṅ͕θ͔̇+Ȯ͔ι͕̇-Ȯ͕ι͔̇+Q͔̇λ͕̇-Q͕̇λ͔̇+Ṙ͔μ͕̇-Ṙ͕μ͔̇-Ṫ͔ο͕̇+Ṫ͕ο͔̇)e04",
    ///     "+(-Ċ͔α͕̇+Ċ͕α͔̇-Ḟ͔β͕̇+Ḟ͕β͔̇-Ḣ͔γ͕̇+Ḣ͕γ͔̇+J͔̇ε͕̇-J͕̇ε͔̇-L͔̇ζ͕̇+L͕̇ζ͔̇-Ṅ͔η͕̇+Ṅ͕η͔̇+Ṗ͔ι͕̇-Ṗ͕ι͔̇-Q͔̇κ͕̇+Q͕̇κ͔̇+Ṡ͔μ͕̇-Ṡ͕μ͔̇+Ṫ͔ξ͕̇-Ṫ͕ξ͔̇)e05",
    ///     "+(-Ḋ͔α͕̇+Ḋ͕α͔̇-Ġ͔β͕̇+Ġ͕β͔̇-İ͔γ͕̇+İ͕γ͔̇-J͔̇δ͕̇+J͕̇δ͔̇-Ṁ͔ζ͕̇+Ṁ͕ζ͔̇-Ȯ͔η͕̇+Ȯ͕η͔̇-Ṗ͔θ͕̇+Ṗ͕θ͔̇-Ṙ͔κ͕̇+Ṙ͕κ͔̇-Ṡ͔λ͕̇+Ṡ͕λ͔̇-Ṫ͔ν͕̇+Ṫ͕ν͔̇)e06",
    ///     "+(+β͔̇ζ͕̇-β͕̇ζ͔̇+γ͔̇η͕̇-γ͕̇η͔̇+δ͔̇θ͕̇-δ͕̇θ͔̇+ε͔̇ι͕̇-ε͕̇ι͔̇)e12",
    ///     "+(-α͔̇ζ͕̇+α͕̇ζ͔̇+γ͔̇κ͕̇-γ͕̇κ͔̇+δ͔̇λ͕̇-δ͕̇λ͔̇+ε͔̇μ͕̇-ε͕̇μ͔̇)e13",
    ///     "+(-α͔̇η͕̇+α͕̇η͔̇-β͔̇κ͕̇+β͕̇κ͔̇+δ͔̇ν͕̇-δ͕̇ν͔̇+ε͔̇ξ͕̇-ε͕̇ξ͔̇)e14",
    ///     "+(-α͔̇θ͕̇+α͕̇θ͔̇-β͔̇λ͕̇+β͕̇λ͔̇-γ͔̇ν͕̇+γ͕̇ν͔̇+ε͔̇ο͕̇-ε͕̇ο͔̇)e15",
    ///     "+(-α͔̇ι͕̇+α͕̇ι͔̇-β͔̇μ͕̇+β͕̇μ͔̇-γ͔̇ξ͕̇+γ͕̇ξ͔̇-δ͔̇ο͕̇+δ͕̇ο͔̇)e16",
    ///     "+(+α͔̇β͕̇-α͕̇β͔̇+η͔̇κ͕̇-η͕̇κ͔̇+θ͔̇λ͕̇-θ͕̇λ͔̇+ι͔̇μ͕̇-ι͕̇μ͔̇)e23",
    ///     "+(+α͔̇γ͕̇-α͕̇γ͔̇-ζ͔̇κ͕̇+ζ͕̇κ͔̇+θ͔̇ν͕̇-θ͕̇ν͔̇+ι͔̇ξ͕̇-ι͕̇ξ͔̇)e24",
    ///     "+(+α͔̇δ͕̇-α͕̇δ͔̇-ζ͔̇λ͕̇+ζ͕̇λ͔̇-η͔̇ν͕̇+η͕̇ν͔̇+ι͔̇ο͕̇-ι͕̇ο͔̇)e25",
    ///     "+(+α͔̇ε͕̇-α͕̇ε͔̇-ζ͔̇μ͕̇+ζ͕̇μ͔̇-η͔̇ξ͕̇+η͕̇ξ͔̇-θ͔̇ο͕̇+θ͕̇ο͔̇)e26",
    ///     "+(+β͔̇γ͕̇-β͕̇γ͔̇+ζ͔̇η͕̇-ζ͕̇η͔̇+λ͔̇ν͕̇-λ͕̇ν͔̇+μ͔̇ξ͕̇-μ͕̇ξ͔̇)e34",
    ///     "+(+β͔̇δ͕̇-β͕̇δ͔̇+ζ͔̇θ͕̇-ζ͕̇θ͔̇-κ͔̇ν͕̇+κ͕̇ν͔̇+μ͔̇ο͕̇-μ͕̇ο͔̇)e35",
    ///     "+(+β͔̇ε͕̇-β͕̇ε͔̇+ζ͔̇ι͕̇-ζ͕̇ι͔̇-κ͔̇ξ͕̇+κ͕̇ξ͔̇-λ͔̇ο͕̇+λ͕̇ο͔̇)e36",
    ///     "+(+γ͔̇δ͕̇-γ͕̇δ͔̇+η͔̇θ͕̇-η͕̇θ͔̇+κ͔̇λ͕̇-κ͕̇λ͔̇+ξ͔̇ο͕̇-ξ͕̇ο͔̇)e45",
    ///     "+(+γ͔̇ε͕̇-γ͕̇ε͔̇+η͔̇ι͕̇-η͕̇ι͔̇+κ͔̇μ͕̇-κ͕̇μ͔̇-ν͔̇ο͕̇+ν͕̇ο͔̇)e46",
    ///     "+(+δ͔̇ε͕̇-δ͕̇ε͔̇+θ͔̇ι͕̇-θ͕̇ι͔̇+λ͔̇μ͕̇-λ͕̇μ͔̇+ν͔̇ξ͕̇-ν͕̇ξ͔̇)e56",
    ///     "+(-κ͔̇ο͕̇-κ͕̇ο͔̇+λ͔̇ξ͕̇+λ͕̇ξ͔̇-μ͔̇ν͕̇-μ͕̇ν͔̇)e3456",
    ///     "+(+η͔̇ο͕̇+η͕̇ο͔̇-θ͔̇ξ͕̇-θ͕̇ξ͔̇+ι͔̇ν͕̇+ι͕̇ν͔̇)e2465",
    ///     "+(-ζ͔̇ο͕̇-ζ͕̇ο͔̇+θ͔̇μ͕̇+θ͕̇μ͔̇-ι͔̇λ͕̇-ι͕̇λ͔̇)e2356",
    ///     "+(+ζ͔̇ξ͕̇+ζ͕̇ξ͔̇-η͔̇μ͕̇-η͕̇μ͔̇+ι͔̇κ͕̇+ι͕̇κ͔̇)e2364",
    ///     "+(-ζ͔̇ν͕̇-ζ͕̇ν͔̇+η͔̇λ͕̇+η͕̇λ͔̇-θ͔̇κ͕̇-θ͕̇κ͔̇)e2345",
    ///     "+(-γ͔̇ο͕̇-γ͕̇ο͔̇+δ͔̇ξ͕̇+δ͕̇ξ͔̇-ε͔̇ν͕̇-ε͕̇ν͔̇)e1456",
    ///     "+(+β͔̇ο͕̇+β͕̇ο͔̇-δ͔̇μ͕̇-δ͕̇μ͔̇+ε͔̇λ͕̇+ε͕̇λ͔̇)e1365",
    ///     "+(-β͔̇ξ͕̇-β͕̇ξ͔̇+γ͔̇μ͕̇+γ͕̇μ͔̇-ε͔̇κ͕̇-ε͕̇κ͔̇)e1346",
    ///     "+(+β͔̇ν͕̇+β͕̇ν͔̇-γ͔̇λ͕̇-γ͕̇λ͔̇+δ͔̇κ͕̇+δ͕̇κ͔̇)e1354",
    ///     "+(-α͔̇ο͕̇-α͕̇ο͔̇+δ͔̇ι͕̇+δ͕̇ι͔̇-ε͔̇θ͕̇-ε͕̇θ͔̇)e1256",
    ///     "+(+α͔̇ξ͕̇+α͕̇ξ͔̇-γ͔̇ι͕̇-γ͕̇ι͔̇+ε͔̇η͕̇+ε͕̇η͔̇)e1264",
    ///     "+(-α͔̇ν͕̇-α͕̇ν͔̇+γ͔̇θ͕̇+γ͕̇θ͔̇-δ͔̇η͕̇-δ͕̇η͔̇)e1245",
    ///     "+(-α͔̇μ͕̇-α͕̇μ͔̇+β͔̇ι͕̇+β͕̇ι͔̇-ε͔̇ζ͕̇-ε͕̇ζ͔̇)e1236",
    ///     "+(+α͔̇λ͕̇+α͕̇λ͔̇-β͔̇θ͕̇-β͕̇θ͔̇+δ͔̇ζ͕̇+δ͕̇ζ͔̇)e1253",
    ///     "+(-α͔̇κ͕̇-α͕̇κ͔̇+β͔̇η͕̇+β͕̇η͔̇-γ͔̇ζ͕̇-γ͕̇ζ͔̇)e1234",
    ///     "+(-Ḣ͔ε͕̇-Ḣ͕ε͔̇+İ͔δ͕̇+İ͕δ͔̇-J͔̇γ͕̇-J͕̇γ͔̇-Ṅ͔ι͕̇-Ṅ͕ι͔̇+Ȯ͔θ͕̇+Ȯ͕θ͔̇-Ṗ͔η͕̇-Ṗ͕η͔̇-Q͔̇μ͕̇-Q͕̇μ͔̇+Ṙ͔λ͕̇+Ṙ͕λ͔̇-Ṡ͔κ͕̇-Ṡ͕κ͔̇)e0465",
    ///     "+(+Ḟ͔ε͕̇+Ḟ͕ε͔̇-Ġ͔δ͕̇-Ġ͕δ͔̇+J͔̇β͕̇+J͕̇β͔̇+L͔̇ι͕̇+L͕̇ι͔̇-Ṁ͔θ͕̇-Ṁ͕θ͔̇+Ṗ͔ζ͕̇+Ṗ͕ζ͔̇-Q͔̇ξ͕̇-Q͕̇ξ͔̇+Ṙ͔ν͕̇+Ṙ͕ν͔̇-Ṫ͔κ͕̇-Ṫ͕κ͔̇)e0356",
    ///     "+(-Ė͔ε͕̇-Ė͕ε͔̇+Ġ͔γ͕̇+Ġ͕γ͔̇-İ͔β͕̇-İ͕β͔̇-K͔̇ι͕̇-K͕̇ι͔̇+Ṁ͔η͕̇+Ṁ͕η͔̇-Ȯ͔ζ͕̇-Ȯ͕ζ͔̇-Q͔̇ο͕̇-Q͕̇ο͔̇+Ṡ͔ν͕̇+Ṡ͕ν͔̇-Ṫ͔λ͕̇-Ṫ͕λ͔̇)e0364",
    ///     "+(+Ė͔δ͕̇+Ė͕δ͔̇-Ḟ͔γ͕̇-Ḟ͕γ͔̇+Ḣ͔β͕̇+Ḣ͕β͔̇+K͔̇θ͕̇+K͕̇θ͔̇-L͔̇η͕̇-L͕̇η͔̇+Ṅ͔ζ͕̇+Ṅ͕ζ͔̇-Ṙ͔ο͕̇-Ṙ͕ο͔̇+Ṡ͔ξ͕̇+Ṡ͕ξ͔̇-Ṫ͔μ͕̇-Ṫ͕μ͔̇)e0345",
    ///     "+(-Ċ͔ε͕̇-Ċ͕ε͔̇+Ḋ͔δ͕̇+Ḋ͕δ͔̇-J͔̇α͕̇-J͕̇α͔̇+L͔̇μ͕̇+L͕̇μ͔̇-Ṁ͔λ͕̇-Ṁ͕λ͔̇+Ṅ͔ξ͕̇+Ṅ͕ξ͔̇-Ȯ͔ν͕̇-Ȯ͕ν͔̇+Ṡ͔ζ͕̇+Ṡ͕ζ͔̇+Ṫ͔η͕̇+Ṫ͕η͔̇)e0265",
    ///     "+(+Ḃ͔ε͕̇+Ḃ͕ε͔̇-Ḋ͔γ͕̇-Ḋ͕γ͔̇+İ͔α͕̇+İ͕α͔̇-K͔̇μ͕̇-K͕̇μ͔̇+Ṁ͔κ͕̇+Ṁ͕κ͔̇+Ṅ͔ο͕̇+Ṅ͕ο͔̇-Ṗ͔ν͕̇-Ṗ͕ν͔̇-Ṙ͔ζ͕̇-Ṙ͕ζ͔̇+Ṫ͔θ͕̇+Ṫ͕θ͔̇)e0246",
    ///     "+(-Ḃ͔δ͕̇-Ḃ͕δ͔̇+Ċ͔γ͕̇+Ċ͕γ͔̇-Ḣ͔α͕̇-Ḣ͕α͔̇+K͔̇λ͕̇+K͕̇λ͔̇-L͔̇κ͕̇-L͕̇κ͔̇+Ȯ͔ο͕̇+Ȯ͕ο͔̇-Ṗ͔ξ͕̇-Ṗ͕ξ͔̇+Q͔̇ζ͕̇+Q͕̇ζ͔̇+Ṫ͔ι͕̇+Ṫ͕ι͔̇)e0254",
    ///     "+(-Ȧ͔ε͕̇-Ȧ͕ε͔̇+Ḋ͔β͕̇+Ḋ͕β͔̇-Ġ͔α͕̇-Ġ͕α͔̇-K͔̇ξ͕̇-K͕̇ξ͔̇-L͔̇ο͕̇-L͕̇ο͔̇+Ȯ͔κ͕̇+Ȯ͕κ͔̇+Ṗ͔λ͕̇+Ṗ͕λ͔̇-Ṙ͔η͕̇-Ṙ͕η͔̇-Ṡ͔θ͕̇-Ṡ͕θ͔̇)e0263",
    ///     "+(+Ȧ͔δ͕̇+Ȧ͕δ͔̇-Ċ͔β͕̇-Ċ͕β͔̇+Ḟ͔α͕̇+Ḟ͕α͔̇+K͔̇ν͕̇+K͕̇ν͔̇-Ṁ͔ο͕̇-Ṁ͕ο͔̇-Ṅ͔κ͕̇-Ṅ͕κ͔̇+Ṗ͔μ͕̇+Ṗ͕μ͔̇+Q͔̇η͕̇+Q͕̇η͔̇-Ṡ͔ι͕̇-Ṡ͕ι͔̇)e0235",
    ///     "+(-Ȧ͔γ͕̇-Ȧ͕γ͔̇+Ḃ͔β͕̇+Ḃ͕β͔̇-Ė͔α͕̇-Ė͕α͔̇+L͔̇ν͕̇+L͕̇ν͔̇+Ṁ͔ξ͕̇+Ṁ͕ξ͔̇-Ṅ͔λ͕̇-Ṅ͕λ͔̇-Ȯ͔μ͕̇-Ȯ͕μ͔̇+Q͔̇θ͕̇+Q͕̇θ͔̇+Ṙ͔ι͕̇+Ṙ͕ι͔̇)e0243",
    ///     "+(-Ċ͔ι͕̇-Ċ͕ι͔̇+Ḋ͔θ͕̇+Ḋ͕θ͔̇-Ḟ͔μ͕̇-Ḟ͕μ͔̇+Ġ͔λ͕̇+Ġ͕λ͔̇-Ḣ͔ξ͕̇-Ḣ͕ξ͔̇+İ͔ν͕̇+İ͕ν͔̇-Ṗ͔α͕̇-Ṗ͕α͔̇-Ṡ͔β͕̇-Ṡ͕β͔̇-Ṫ͔γ͕̇-Ṫ͕γ͔̇)e0156",
    ///     "+(+Ḃ͔ι͕̇+Ḃ͕ι͔̇-Ḋ͔η͕̇-Ḋ͕η͔̇+Ė͔μ͕̇+Ė͕μ͔̇-Ġ͔κ͕̇-Ġ͕κ͔̇-Ḣ͔ο͕̇-Ḣ͕ο͔̇+J͔̇ν͕̇+J͕̇ν͔̇+Ȯ͔α͕̇+Ȯ͕α͔̇+Ṙ͔β͕̇+Ṙ͕β͔̇-Ṫ͔δ͕̇-Ṫ͕δ͔̇)e0164",
    ///     "+(-Ḃ͔θ͕̇-Ḃ͕θ͔̇+Ċ͔η͕̇+Ċ͕η͔̇-Ė͔λ͕̇-Ė͕λ͔̇+Ḟ͔κ͕̇+Ḟ͕κ͔̇-İ͔ο͕̇-İ͕ο͔̇+J͔̇ξ͕̇+J͕̇ξ͔̇-Ṅ͔α͕̇-Ṅ͕α͔̇-Q͔̇β͕̇-Q͕̇β͔̇-Ṫ͔ε͕̇-Ṫ͕ε͔̇)e0145",
    ///     "+(-Ȧ͔ι͕̇-Ȧ͕ι͔̇+Ḋ͔ζ͕̇+Ḋ͕ζ͔̇+Ė͔ξ͕̇+Ė͕ξ͔̇+Ḟ͔ο͕̇+Ḟ͕ο͔̇-İ͔κ͕̇-İ͕κ͔̇-J͔̇λ͕̇-J͕̇λ͔̇-Ṁ͔α͕̇-Ṁ͕α͔̇+Ṙ͔γ͕̇+Ṙ͕γ͔̇+Ṡ͔δ͕̇+Ṡ͕δ͔̇)e0136",
    ///     "+(+Ȧ͔θ͕̇+Ȧ͕θ͔̇-Ċ͔ζ͕̇-Ċ͕ζ͔̇-Ė͔ν͕̇-Ė͕ν͔̇+Ġ͔ο͕̇+Ġ͕ο͔̇+Ḣ͔κ͕̇+Ḣ͕κ͔̇-J͔̇μ͕̇-J͕̇μ͔̇+L͔̇α͕̇+L͕̇α͔̇-Q͔̇γ͕̇-Q͕̇γ͔̇+Ṡ͔ε͕̇+Ṡ͕ε͔̇)e0153",
    ///     "+(-Ȧ͔η͕̇-Ȧ͕η͔̇+Ḃ͔ζ͕̇+Ḃ͕ζ͔̇-Ḟ͔ν͕̇-Ḟ͕ν͔̇-Ġ͔ξ͕̇-Ġ͕ξ͔̇+Ḣ͔λ͕̇+Ḣ͕λ͔̇+İ͔μ͕̇+İ͕μ͔̇-K͔̇α͕̇-K͕̇α͔̇-Q͔̇δ͕̇-Q͕̇δ͔̇-Ṙ͔ε͕̇-Ṙ͕ε͔̇)e0134",
    ///     "+(-Ȧ͔μ͕̇-Ȧ͕μ͔̇-Ḃ͔ξ͕̇-Ḃ͕ξ͔̇-Ċ͔ο͕̇-Ċ͕ο͔̇+Ġ͔ζ͕̇+Ġ͕ζ͔̇+İ͔η͕̇+İ͕η͔̇+J͔̇θ͕̇+J͕̇θ͔̇-Ṁ͔β͕̇-Ṁ͕β͔̇-Ȯ͔γ͕̇-Ȯ͕γ͔̇-Ṗ͔δ͕̇-Ṗ͕δ͔̇)e0162",
    ///     "+(+Ȧ͔λ͕̇+Ȧ͕λ͔̇+Ḃ͔ν͕̇+Ḃ͕ν͔̇-Ḋ͔ο͕̇-Ḋ͕ο͔̇-Ḟ͔ζ͕̇-Ḟ͕ζ͔̇-Ḣ͔η͕̇-Ḣ͕η͔̇+J͔̇ι͕̇+J͕̇ι͔̇+L͔̇β͕̇+L͕̇β͔̇+Ṅ͔γ͕̇+Ṅ͕γ͔̇-Ṗ͔ε͕̇-Ṗ͕ε͔̇)e0125",
    ///     "+(-Ȧ͔κ͕̇-Ȧ͕κ͔̇+Ċ͔ν͕̇+Ċ͕ν͔̇+Ḋ͔ξ͕̇+Ḋ͕ξ͔̇+Ė͔ζ͕̇+Ė͕ζ͔̇-Ḣ͔θ͕̇-Ḣ͕θ͔̇-İ͔ι͕̇-İ͕ι͔̇-K͔̇β͕̇-K͕̇β͔̇+Ṅ͔δ͕̇+Ṅ͕δ͔̇+Ȯ͔ε͕̇+Ȯ͕ε͔̇)e0142",
    ///     "+(-Ḃ͔κ͕̇-Ḃ͕κ͔̇-Ċ͔λ͕̇-Ċ͕λ͔̇-Ḋ͔μ͕̇-Ḋ͕μ͔̇+Ė͔η͕̇+Ė͕η͔̇+Ḟ͔θ͕̇+Ḟ͕θ͔̇+Ġ͔ι͕̇+Ġ͕ι͔̇-K͔̇γ͕̇-K͕̇γ͔̇-L͔̇δ͕̇-L͕̇δ͔̇-Ṁ͔ε͕̇-Ṁ͕ε͔̇)e0123",
    ///     "+(-K͔̇ο͕̇+K͕̇ο͔̇+L͔̇ξ͕̇-L͕̇ξ͔̇-Ṁ͔ν͕̇+Ṁ͕ν͔̇-Ṅ͔μ͕̇+Ṅ͕μ͔̇+Ȯ͔λ͕̇-Ȯ͕λ͔̇-Ṗ͔κ͕̇+Ṗ͕κ͔̇+Q͔̇ι͕̇-Q͕̇ι͔̇-Ṙ͔θ͕̇+Ṙ͕θ͔̇+Ṡ͔η͕̇-Ṡ͕η͔̇-Ṫ͔ζ͕̇+Ṫ͕ζ͔̇)e023465",
    ///     "+(+Ė͔ο͕̇-Ė͕ο͔̇-Ḟ͔ξ͕̇+Ḟ͕ξ͔̇+Ġ͔ν͕̇-Ġ͕ν͔̇+Ḣ͔μ͕̇-Ḣ͕μ͔̇-İ͔λ͕̇+İ͕λ͔̇+J͔̇κ͕̇-J͕̇κ͔̇-Q͔̇ε͕̇+Q͕̇ε͔̇+Ṙ͔δ͕̇-Ṙ͕δ͔̇-Ṡ͔γ͕̇+Ṡ͕γ͔̇+Ṫ͔β͕̇-Ṫ͕β͔̇)e013456",
    ///     "+(-Ḃ͔ο͕̇+Ḃ͕ο͔̇+Ċ͔ξ͕̇-Ċ͕ξ͔̇-Ḋ͔ν͕̇+Ḋ͕ν͔̇-Ḣ͔ι͕̇+Ḣ͕ι͔̇+İ͔θ͕̇-İ͕θ͔̇-J͔̇η͕̇+J͕̇η͔̇+Ṅ͔ε͕̇-Ṅ͕ε͔̇-Ȯ͔δ͕̇+Ȯ͕δ͔̇+Ṗ͔γ͕̇-Ṗ͕γ͔̇-Ṫ͔α͕̇+Ṫ͕α͔̇)e012465",
    ///     "+(+Ȧ͔ο͕̇-Ȧ͕ο͔̇-Ċ͔μ͕̇+Ċ͕μ͔̇+Ḋ͔λ͕̇-Ḋ͕λ͔̇+Ḟ͔ι͕̇-Ḟ͕ι͔̇-Ġ͔θ͕̇+Ġ͕θ͔̇+J͔̇ζ͕̇-J͕̇ζ͔̇-L͔̇ε͕̇+L͕̇ε͔̇+Ṁ͔δ͕̇-Ṁ͕δ͔̇-Ṗ͔β͕̇+Ṗ͕β͔̇+Ṡ͔α͕̇-Ṡ͕α͔̇)e012356",
    ///     "+(-Ȧ͔ξ͕̇+Ȧ͕ξ͔̇+Ḃ͔μ͕̇-Ḃ͕μ͔̇-Ḋ͔κ͕̇+Ḋ͕κ͔̇-Ė͔ι͕̇+Ė͕ι͔̇+Ġ͔η͕̇-Ġ͕η͔̇-İ͔ζ͕̇+İ͕ζ͔̇+K͔̇ε͕̇-K͕̇ε͔̇-Ṁ͔γ͕̇+Ṁ͕γ͔̇+Ȯ͔β͕̇-Ȯ͕β͔̇-Ṙ͔α͕̇+Ṙ͕α͔̇)e012364",
    ///     "+(+Ȧ͔ν͕̇-Ȧ͕ν͔̇-Ḃ͔λ͕̇+Ḃ͕λ͔̇+Ċ͔κ͕̇-Ċ͕κ͔̇+Ė͔θ͕̇-Ė͕θ͔̇-Ḟ͔η͕̇+Ḟ͕η͔̇+Ḣ͔ζ͕̇-Ḣ͕ζ͔̇-K͔̇δ͕̇+K͕̇δ͔̇+L͔̇γ͕̇-L͕̇γ͔̇-Ṅ͔β͕̇+Ṅ͕β͔̇+Q͔̇α͕̇-Q͕̇α͔̇)e012345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_triple_motor() -> Self {
        Self::scalar() + Self::volume4() + Self::plane() + Self::direction()
    }
    /// The multivector of triple motor $`m_3 \equiv s + v^4 + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let triple_motor = Vee::triple_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(triple_motor.basis_blades(), Vee::triple_motor().basis_blades());
    /// format_eq!(triple_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕ζ͔+Ð͕η͔+Ø͕θ͔+Þ͕ι͔)e02",
    ///     "+(-X͕β͔-Y͕ζ͔+Z͕v͔+Ð͕κ͔+Ø͕λ͔+Þ͕μ͔)e03",
    ///     "+(-X͕γ͔-Y͕η͔-Z͕κ͔+v͔Ð͕+Ø͕ν͔+Þ͕ξ͔)e04",
    ///     "+(-X͕δ͔-Y͕θ͔-Z͕λ͔+v͔Ø͕-Ð͕ν͔+Þ͕ο͔)e05",
    ///     "+(-X͕ε͔-Y͕ι͔-Z͕μ͔+v͔Þ͕-Ð͕ξ͔-Ø͕ο͔)e06",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e23",
    ///     "+v͕η͔e24",
    ///     "+v͕θ͔e25",
    ///     "+v͕ι͔e26",
    ///     "+v͕κ͔e34",
    ///     "+v͕λ͔e35",
    ///     "+v͕μ͔e36",
    ///     "+v͕ν͔e45",
    ///     "+v͕ξ͔e46",
    ///     "+v͕ο͔e56",
    ///     "+v͕α͔̇e3456",
    ///     "+v͕β͔̇e2465",
    ///     "+v͕γ͔̇e2356",
    ///     "+v͕δ͔̇e2364",
    ///     "+v͕ε͔̇e2345",
    ///     "+v͕ζ͔̇e1456",
    ///     "+v͕η͔̇e1365",
    ///     "+v͕θ͔̇e1346",
    ///     "+v͕ι͔̇e1354",
    ///     "+v͕κ͔̇e1256",
    ///     "+v͕λ͔̇e1264",
    ///     "+v͕μ͔̇e1245",
    ///     "+v͕ν͔̇e1236",
    ///     "+v͕ξ͔̇e1253",
    ///     "+v͕ο͔̇e1234",
    ///     "+(+X͕ζ͔̇-Y͕β͔̇+Z͕α͔̇-Ð͕ο͔+Ø͕ξ͔-Þ͕ν͔)e0465",
    ///     "+(+X͕η͔̇-Y͕γ͔̇+Z͕ο͔+Ð͕α͔̇-Ø͕μ͔+Þ͕λ͔)e0356",
    ///     "+(+X͕θ͔̇-Y͕δ͔̇-Z͕ξ͔+Ð͕μ͔+Ø͕α͔̇-Þ͕κ͔)e0364",
    ///     "+(+X͕ι͔̇-Y͕ε͔̇+Z͕ν͔-Ð͕λ͔+Ø͕κ͔+Þ͕α͔̇)e0345",
    ///     "+(+X͕κ͔̇-Y͕ο͔-Z͕γ͔̇+Ð͕β͔̇+Ø͕ι͔-Þ͕θ͔)e0265",
    ///     "+(+X͕λ͔̇+Y͕ξ͔-Z͕δ͔̇-Ð͕ι͔+Ø͕β͔̇+Þ͕η͔)e0246",
    ///     "+(+X͕μ͔̇-Y͕ν͔-Z͕ε͔̇+Ð͕θ͔-Ø͕η͔+Þ͕β͔̇)e0254",
    ///     "+(+X͕ν͔̇-Y͕μ͔+Z͕ι͔-Ð͕δ͔̇+Ø͕γ͔̇-Þ͕ζ͔)e0263",
    ///     "+(+X͕ξ͔̇+Y͕λ͔-Z͕θ͔-Ð͕ε͔̇+Ø͕ζ͔+Þ͕γ͔̇)e0235",
    ///     "+(+X͕ο͔̇-Y͕κ͔+Z͕η͔-Ð͕ζ͔-Ø͕ε͔̇+Þ͕δ͔̇)e0243",
    ///     "+(+X͕ο͔+Y͕κ͔̇-Z͕η͔̇+Ð͕ζ͔̇-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(-X͕ξ͔+Y͕λ͔̇-Z͕θ͔̇+Ð͕ε͔+Ø͕ζ͔̇-Þ͕γ͔)e0164",
    ///     "+(+X͕ν͔+Y͕μ͔̇-Z͕ι͔̇-Ð͕δ͔+Ø͕γ͔+Þ͕ζ͔̇)e0145",
    ///     "+(+X͕μ͔+Y͕ν͔̇-Z͕ε͔-Ð͕θ͔̇+Ø͕η͔̇+Þ͕β͔)e0136",
    ///     "+(-X͕λ͔+Y͕ξ͔̇+Z͕δ͔-Ð͕ι͔̇-Ø͕β͔+Þ͕η͔̇)e0153",
    ///     "+(+X͕κ͔+Y͕ο͔̇-Z͕γ͔+Ð͕β͔-Ø͕ι͔̇+Þ͕θ͔̇)e0134",
    ///     "+(-X͕ι͔+Y͕ε͔+Z͕ν͔̇-Ð͕λ͔̇+Ø͕κ͔̇-Þ͕α͔)e0162",
    ///     "+(+X͕θ͔-Y͕δ͔+Z͕ξ͔̇-Ð͕μ͔̇+Ø͕α͔+Þ͕κ͔̇)e0125",
    ///     "+(-X͕η͔+Y͕γ͔+Z͕ο͔̇-Ð͕α͔-Ø͕μ͔̇+Þ͕λ͔̇)e0142",
    ///     "+(+X͕ζ͔-Y͕β͔+Z͕α͔+Ð͕ο͔̇-Ø͕ξ͔̇+Þ͕ν͔̇)e0123",
    ///     "+v͕ẇ͔e123456",
    ///     "+(+X͕ẇ͔-Y͕α͔̇-Z͕β͔̇-Ð͕γ͔̇-Ø͕δ͔̇-Þ͕ε͔̇)e023465",
    ///     "+(+X͕α͔̇+Y͕ẇ͔-Z͕ζ͔̇-Ð͕η͔̇-Ø͕θ͔̇-Þ͕ι͔̇)e013456",
    ///     "+(+X͕β͔̇+Y͕ζ͔̇+Z͕ẇ͔-Ð͕κ͔̇-Ø͕λ͔̇-Þ͕μ͔̇)e012465",
    ///     "+(+X͕γ͔̇+Y͕η͔̇+Z͕κ͔̇+ẇ͔Ð͕-Ø͕ν͔̇-Þ͕ξ͔̇)e012356",
    ///     "+(+X͕δ͔̇+Y͕θ͔̇+Z͕λ͔̇+ẇ͔Ø͕+Ð͕ν͔̇-Þ͕ο͔̇)e012364",
    ///     "+(+X͕ε͔̇+Y͕ι͔̇+Z͕μ͔̇+ẇ͔Þ͕+Ð͕ξ͔̇+Ø͕ο͔̇)e012345",
    /// ]);
    ///
    /// let norm_squared = Vee::triple_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+vv+ẇẇ+αα+α̇α̇+ββ+β̇β̇+γγ+γ̇γ̇+δδ+δ̇δ̇+εε+ε̇ε̇+ζζ+ζ̇ζ̇+ηη+η̇η̇\
    ///      +θθ+θ̇θ̇+ιι+ι̇ι̇+κκ+κ̇κ̇+λλ+λ̇λ̇+μμ+μ̇μ̇+νν+ν̇ν̇+ξξ+ξ̇ξ̇+οο+ο̇ο̇",
    ///     "+2(+vα̇+ẇα-βζ̇+β̇ζ-γη̇+γ̇η-δθ̇+δ̇θ-ει̇+ε̇ι-κο-κ̇ο̇+λξ+λ̇ξ̇-μν-μ̇ν̇)e3456",
    ///     "+2(+vβ̇+ẇβ+αζ̇-α̇ζ-γκ̇+γ̇κ-δλ̇+δ̇λ-εμ̇+ε̇μ+ηο+η̇ο̇-θξ-θ̇ξ̇+ιν+ι̇ν̇)e2465",
    ///     "+2(+vγ̇+ẇγ+αη̇-α̇η+βκ̇-β̇κ-δν̇+δ̇ν-εξ̇+ε̇ξ-ζο-ζ̇ο̇+θμ+θ̇μ̇-ιλ-ι̇λ̇)e2356",
    ///     "+2(+vδ̇+ẇδ+αθ̇-α̇θ+βλ̇-β̇λ+γν̇-γ̇ν-εο̇+ε̇ο+ζξ+ζ̇ξ̇-ημ-η̇μ̇+ικ+ι̇κ̇)e2364",
    ///     "+2(+vε̇+ẇε+αι̇-α̇ι+βμ̇-β̇μ+γξ̇-γ̇ξ+δο̇-δ̇ο-ζν-ζ̇ν̇+ηλ+η̇λ̇-θκ-θ̇κ̇)e2345",
    ///     "+2(+vζ̇+ẇζ-αβ̇+α̇β-γο-γ̇ο̇+δξ+δ̇ξ̇-εν-ε̇ν̇-ηκ̇+η̇κ-θλ̇+θ̇λ-ιμ̇+ι̇μ)e1456",
    ///     "+2(+vη̇+ẇη-αγ̇+α̇γ+βο+β̇ο̇-δμ-δ̇μ̇+ελ+ε̇λ̇+ζκ̇-ζ̇κ-θν̇+θ̇ν-ιξ̇+ι̇ξ)e1365",
    ///     "+2(+vθ̇+ẇθ-αδ̇+α̇δ-βξ-β̇ξ̇+γμ+γ̇μ̇-εκ-ε̇κ̇+ζλ̇-ζ̇λ+ην̇-η̇ν-ιο̇+ι̇ο)e1346",
    ///     "+2(+vι̇+ẇι-αε̇+α̇ε+βν+β̇ν̇-γλ-γ̇λ̇+δκ+δ̇κ̇+ζμ̇-ζ̇μ+ηξ̇-η̇ξ+θο̇-θ̇ο)e1354",
    ///     "+2(+vκ̇+ẇκ-αο-α̇ο̇-βγ̇+β̇γ+δι+δ̇ι̇-εθ-ε̇θ̇-ζη̇+ζ̇η-λν̇+λ̇ν-μξ̇+μ̇ξ)e1256",
    ///     "+2(+vλ̇+ẇλ+αξ+α̇ξ̇-βδ̇+β̇δ-γι-γ̇ι̇+εη+ε̇η̇-ζθ̇+ζ̇θ+κν̇-κ̇ν-μο̇+μ̇ο)e1264",
    ///     "+2(+vμ̇+ẇμ-αν-α̇ν̇-βε̇+β̇ε+γθ+γ̇θ̇-δη-δ̇η̇-ζι̇+ζ̇ι+κξ̇-κ̇ξ+λο̇-λ̇ο)e1245",
    ///     "+2(+vν̇+ẇν-αμ-α̇μ̇+βι+β̇ι̇-γδ̇+γ̇δ-εζ-ε̇ζ̇-ηθ̇+η̇θ-κλ̇+κ̇λ-ξο̇+ξ̇ο)e1236",
    ///     "+2(+vξ̇+ẇξ+αλ+α̇λ̇-βθ-β̇θ̇-γε̇+γ̇ε+δζ+δ̇ζ̇-ηι̇+η̇ι-κμ̇+κ̇μ+νο̇-ν̇ο)e1253",
    ///     "+2(+vο̇+ẇο-ακ-α̇κ̇+βη+β̇η̇-γζ-γ̇ζ̇-δε̇+δ̇ε-θι̇+θ̇ι-λμ̇+λ̇μ-νξ̇+ν̇ξ)e1234",
    ///     "+2(+Ȧv+Ḃκ+Ċλ+Ḋμ-Ėη-Ḟθ-Ġι-Ḣε̇+İδ̇-J̇γ̇+K̇γ+L̇δ+Ṁε-Ṅι̇+Ȯθ̇-Ṗη̇\
    ///         -Q̇μ̇+Ṙλ̇-Ṡκ̇-Ṫẇ-Xζ̇+Ẋζ+Yβ̇-Ẏβ-Zα̇+Żα+Ðο+Ð̇ο̇-Øξ-Ø̇ξ̇+Þν+Þ̇ν̇)e0465",
    ///     "+2(-Ȧκ+Ḃv+Ċν+Ḋξ+Ėζ+Ḟε̇-Ġδ̇-Ḣθ-İι+J̇β̇-K̇β+L̇ι̇-Ṁθ̇+Ṅδ+Ȯε+Ṗζ̇\
    ///         -Q̇ξ̇+Ṙν̇+Ṡẇ-Ṫκ̇-Xη̇+Ẋη+Yγ̇-Ẏγ-Zο-Żο̇-Ðα̇+Ð̇α+Øμ+Ø̇μ̇-Þλ-Þ̇λ̇)e0356",
    ///     "+2(-Ȧλ-Ḃν+Ċv+Ḋο-Ėε̇+Ḟζ+Ġγ̇+Ḣη-İβ̇-J̇ι-K̇ι̇-L̇β+Ṁη̇-Ṅγ-Ȯζ̇+Ṗε\
    ///         -Q̇ο̇-Ṙẇ+Ṡν̇-Ṫλ̇-Xθ̇+Ẋθ+Yδ̇-Ẏδ+Zξ+Żξ̇-Ðμ-Ð̇μ̇-Øα̇+Ø̇α+Þκ+Þ̇κ̇)e0364",
    ///     "+2(-Ȧμ-Ḃξ-Ċο+Ḋv+Ėδ̇-Ḟγ̇+Ġζ+Ḣβ̇+İη+J̇θ+K̇θ̇-L̇η̇-Ṁβ+Ṅζ̇-Ȯγ-Ṗδ\
    ///         +Q̇ẇ-Ṙο̇+Ṡξ̇-Ṫμ̇-Xι̇+Ẋι+Yε̇-Ẏε-Zν-Żν̇+Ðλ+Ð̇λ̇-Øκ-Ø̇κ̇-Þα̇+Þ̇α)e0345",
    ///     "+2(+Ȧη-Ḃζ-Ċε̇+Ḋδ̇+Ėv+Ḟν+Ġξ-Ḣλ-İμ-J̇α̇+K̇α+L̇μ̇-Ṁλ̇+Ṅξ̇-Ȯν̇-Ṗẇ\
    ///         +Q̇δ+Ṙε+Ṡζ̇+Ṫη̇-Xκ̇+Ẋκ+Yο+Ẏο̇+Zγ̇-Żγ-Ðβ̇+Ð̇β-Øι-Ø̇ι̇+Þθ+Þ̇θ̇)e0265",
    ///     "+2(+Ȧθ+Ḃε̇-Ċζ-Ḋγ̇-Ėν+Ḟv+Ġο+Ḣκ+İα̇-J̇μ-K̇μ̇+L̇α+Ṁκ̇+Ṅο̇+Ȯẇ-Ṗν̇\
    ///         -Q̇γ-Ṙζ̇+Ṡε+Ṫθ̇-Xλ̇+Ẋλ-Yξ-Ẏξ̇+Zδ̇-Żδ+Ðι+Ð̇ι̇-Øβ̇+Ø̇β-Þη-Þ̇η̇)e0246",
    ///     "+2(+Ȧι-Ḃδ̇+Ċγ̇-Ḋζ-Ėξ-Ḟο+Ġv-Ḣα̇+İκ+J̇λ+K̇λ̇-L̇κ̇+Ṁα-Ṅẇ+Ȯο̇-Ṗξ̇\
    ///         +Q̇ζ̇-Ṙγ-Ṡδ+Ṫι̇-Xμ̇+Ẋμ+Yν+Ẏν̇+Zε̇-Żε-Ðθ-Ð̇θ̇+Øη+Ø̇η̇-Þβ̇+Þ̇β)e0254",
    ///     "+2(-Ȧε̇+Ḃθ-Ċη+Ḋβ̇+Ėλ-Ḟκ-Ġα̇+Ḣv+İο-J̇ξ-K̇ξ̇-L̇ο̇-Ṁẇ+Ṅα+Ȯκ̇+Ṗλ̇\
    ///         +Q̇β-Ṙη̇-Ṡθ̇+Ṫε-Xν̇+Ẋν+Yμ+Ẏμ̇-Zι-Żι̇+Ðδ̇-Ð̇δ-Øγ̇+Ø̇γ+Þζ+Þ̇ζ̇)e0263",
    ///     "+2(+Ȧδ̇+Ḃι-Ċβ̇-Ḋη+Ėμ+Ḟα̇-Ġκ-Ḣο+İv+J̇ν+K̇ν̇+L̇ẇ-Ṁο̇-Ṅκ̇+Ȯα+Ṗμ̇\
    ///         +Q̇η̇+Ṙβ-Ṡι̇-Ṫδ-Xξ̇+Ẋξ-Yλ-Ẏλ̇+Zθ+Żθ̇+Ðε̇-Ð̇ε-Øζ-Ø̇ζ̇-Þγ̇+Þ̇γ)e0235",
    ///     "+2(-Ȧγ̇+Ḃβ̇+Ċι-Ḋθ-Ėα̇+Ḟμ-Ġλ+Ḣξ-İν+J̇v-K̇ẇ+L̇ν̇+Ṁξ̇-Ṅλ̇-Ȯμ̇+Ṗα\
    ///         +Q̇θ̇+Ṙι̇+Ṡβ+Ṫγ-Xο̇+Ẋο+Yκ+Ẏκ̇-Zη-Żη̇+Ðζ+Ð̇ζ̇+Øε̇-Ø̇ε-Þδ̇+Þ̇δ)e0243",
    ///     "+2(-Ȧγ+Ḃβ-Ċι̇+Ḋθ̇-Ėα-Ḟμ̇+Ġλ̇-Ḣξ̇+İν̇+J̇ẇ+K̇v+L̇ν+Ṁξ-Ṅλ-Ȯμ-Ṗα̇\
    ///         +Q̇θ+Ṙι-Ṡβ̇-Ṫγ̇-Xο-Ẋο̇-Yκ̇+Ẏκ+Zη̇-Żη-Ðζ̇+Ð̇ζ+Øε+Ø̇ε̇-Þδ-Þ̇δ̇)e0156",
    ///     "+2(-Ȧδ+Ḃι̇+Ċβ-Ḋη̇+Ėμ̇-Ḟα-Ġκ̇-Ḣο̇-İẇ+J̇ν̇-K̇ν+L̇v+Ṁο+Ṅκ+Ȯα̇-Ṗμ\
    ///         -Q̇η+Ṙβ̇+Ṡι-Ṫδ̇+Xξ+Ẋξ̇-Yλ̇+Ẏλ+Zθ̇-Żθ-Ðε-Ð̇ε̇-Øζ̇+Ø̇ζ+Þγ+Þ̇γ̇)e0164",
    ///     "+2(-Ȧε-Ḃθ̇+Ċη̇+Ḋβ-Ėλ̇+Ḟκ̇-Ġα+Ḣẇ-İο̇+J̇ξ̇-K̇ξ-L̇ο+Ṁv-Ṅα̇+Ȯκ+Ṗλ\
    ///         -Q̇β̇-Ṙη-Ṡθ-Ṫε̇-Xν-Ẋν̇-Yμ̇+Ẏμ+Zι̇-Żι+Ðδ+Ð̇δ̇-Øγ-Ø̇γ̇-Þζ̇+Þ̇ζ)e0145",
    ///     "+2(-Ȧι̇-Ḃδ+Ċγ+Ḋζ̇+Ėξ̇+Ḟο̇+Ġẇ-Ḣα-İκ̇-J̇λ̇+K̇λ-L̇κ-Ṁα̇+Ṅv+Ȯο-Ṗξ\
    ///         +Q̇ζ+Ṙγ̇+Ṡδ̇+Ṫι-Xμ-Ẋμ̇-Yν̇+Ẏν+Zε+Żε̇+Ðθ̇-Ð̇θ-Øη̇+Ø̇η-Þβ-Þ̇β̇)e0136",
    ///     "+2(+Ȧθ̇-Ḃε-Ċζ̇+Ḋγ-Ėν̇-Ḟẇ+Ġο̇+Ḣκ̇-İα-J̇μ̇+K̇μ+L̇α̇-Ṁκ-Ṅο+Ȯv+Ṗν\
    ///         -Q̇γ̇+Ṙζ+Ṡε̇-Ṫθ+Xλ+Ẋλ̇-Yξ̇+Ẏξ-Zδ-Żδ̇+Ðι̇-Ð̇ι+Øβ+Ø̇β̇-Þη̇+Þ̇η)e0153",
    ///     "+2(-Ȧη̇+Ḃζ̇-Ċε+Ḋδ+Ėẇ-Ḟν̇-Ġξ̇+Ḣλ̇+İμ̇-J̇α-K̇α̇+L̇μ-Ṁλ+Ṅξ-Ȯν+Ṗv\
    ///         -Q̇δ̇-Ṙε̇+Ṡζ+Ṫη-Xκ-Ẋκ̇-Yο̇+Ẏο+Zγ+Żγ̇-Ðβ-Ð̇β̇+Øι̇-Ø̇ι-Þθ̇+Þ̇θ)e0134",
    ///     "+2(-Ȧμ̇-Ḃξ̇-Ċο̇-Ḋẇ-Ėδ+Ḟγ+Ġζ̇-Ḣβ+İη̇+J̇θ̇-K̇θ+L̇η-Ṁβ̇-Ṅζ-Ȯγ̇-Ṗδ̇\
    ///         +Q̇v+Ṙο-Ṡξ+Ṫμ+Xι+Ẋι̇-Yε-Ẏε̇-Zν̇+Żν+Ðλ̇-Ð̇λ-Øκ̇+Ø̇κ+Þα+Þ̇α̇)e0162",
    ///     "+2(+Ȧλ̇+Ḃν̇+Ċẇ-Ḋο̇-Ėε-Ḟζ̇+Ġγ-Ḣη̇-İβ+J̇ι̇-K̇ι+L̇β̇+Ṁη+Ṅγ̇-Ȯζ-Ṗε̇\
    ///         -Q̇ο+Ṙv+Ṡν-Ṫλ-Xθ-Ẋθ̇+Yδ+Ẏδ̇-Zξ̇+Żξ+Ðμ̇-Ð̇μ-Øα-Ø̇α̇-Þκ̇+Þ̇κ)e0125",
    ///     "+2(-Ȧκ̇-Ḃẇ+Ċν̇+Ḋξ̇+Ėζ̇-Ḟε+Ġδ-Ḣθ̇-İι̇-J̇β-K̇β̇-L̇ι+Ṁθ+Ṅδ̇+Ȯε̇-Ṗζ\
    ///         +Q̇ξ-Ṙν+Ṡv+Ṫκ+Xη+Ẋη̇-Yγ-Ẏγ̇-Zο̇+Żο+Ðα+Ð̇α̇+Øμ̇-Ø̇μ-Þλ̇+Þ̇λ)e0142",
    ///     "+2(+Ȧẇ-Ḃκ̇-Ċλ̇-Ḋμ̇+Ėη̇+Ḟθ̇+Ġι̇-Ḣε+İδ-J̇γ-K̇γ̇-L̇δ̇-Ṁε̇-Ṅι+Ȯθ-Ṗη\
    ///         -Q̇μ+Ṙλ-Ṡκ+Ṫv-Xζ-Ẋζ̇+Yβ+Ẏβ̇-Zα-Żα̇-Ðο̇+Ð̇ο+Øξ̇-Ø̇ξ-Þν̇+Þ̇ν)e0123",
    /// ]);
    ///
    /// let triple_motor = Vee::volume().lhs() * Vee::volume().rhs();
    /// assert_eq!(triple_motor.basis_blades(), Vee::triple_motor().basis_blades());
    /// format_eq!(triple_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕-k͔k͕-l͔l͕-m͔m͕-n͔n͕-o͔o͕-p͔p͕-q͔q͕-r͔r͕-s͔s͕-t͔t͕",
    ///     "+(+a͔Ζ͕-a͕Ζ͔+b͔Η͕-b͕Η͔+c͔Θ͕-c͕Θ͔+d͔Ι͕-d͕Ι͔+e͔Κ͕-e͕Κ͔+f͔Λ͕-f͕Λ͔+g͔Μ͕-g͕Μ͔+h͔Ν͕-h͕Ν͔+i͔Ξ͕-i͕Ξ͔+j͔Ο͕-j͕Ο͔)e01",
    ///     "+(-a͔Β͕+a͕Β͔-b͔Γ͕+b͕Γ͔-c͔Δ͕+c͕Δ͔-d͔Ε͕+d͕Ε͔+k͔Κ͕-k͕Κ͔+l͔Λ͕-l͕Λ͔+m͔Μ͕-m͕Μ͔+n͔Ν͕-n͕Ν͔+o͔Ξ͕-o͕Ξ͔+p͔Ο͕-p͕Ο͔)e02",
    ///     "+(+a͔Α͕-a͕Α͔-e͔Γ͕+e͕Γ͔-f͔Δ͕+f͕Δ͔-g͔Ε͕+g͕Ε͔-k͔Η͕+k͕Η͔-l͔Θ͕+l͕Θ͔-m͔Ι͕+m͕Ι͔+q͔Ν͕-q͕Ν͔+r͔Ξ͕-r͕Ξ͔+s͔Ο͕-s͕Ο͔)e03",
    ///     "+(+b͔Α͕-b͕Α͔+e͔Β͕-e͕Β͔-h͔Δ͕+h͕Δ͔-i͔Ε͕+i͕Ε͔+k͔Ζ͕-k͕Ζ͔-n͔Θ͕+n͕Θ͔-o͔Ι͕+o͕Ι͔-q͔Λ͕+q͕Λ͔-r͔Μ͕+r͕Μ͔+t͔Ο͕-t͕Ο͔)e04",
    ///     "+(+c͔Α͕-c͕Α͔+f͔Β͕-f͕Β͔+h͔Γ͕-h͕Γ͔-j͔Ε͕+j͕Ε͔+l͔Ζ͕-l͕Ζ͔+n͔Η͕-n͕Η͔-p͔Ι͕+p͕Ι͔+q͔Κ͕-q͕Κ͔-s͔Μ͕+s͕Μ͔-t͔Ξ͕+t͕Ξ͔)e05",
    ///     "+(+d͔Α͕-d͕Α͔+g͔Β͕-g͕Β͔+i͔Γ͕-i͕Γ͔+j͔Δ͕-j͕Δ͔+m͔Ζ͕-m͕Ζ͔+o͔Η͕-o͕Η͔+p͔Θ͕-p͕Θ͔+r͔Κ͕-r͕Κ͔+s͔Λ͕-s͕Λ͔+t͔Ν͕-t͕Ν͔)e06",
    ///     "+(-e͔k͕+e͕k͔-f͔l͕+f͕l͔-g͔m͕+g͕m͔-h͔n͕+h͕n͔-i͔o͕+i͕o͔-j͔p͕+j͕p͔)e12",
    ///     "+(+b͔k͕-b͕k͔+c͔l͕-c͕l͔+d͔m͕-d͕m͔-h͔q͕+h͕q͔-i͔r͕+i͕r͔-j͔s͕+j͕s͔)e13",
    ///     "+(-a͔k͕+a͕k͔+c͔n͕-c͕n͔+d͔o͕-d͕o͔+f͔q͕-f͕q͔+g͔r͕-g͕r͔-j͔t͕+j͕t͔)e14",
    ///     "+(-a͔l͕+a͕l͔-b͔n͕+b͕n͔+d͔p͕-d͕p͔-e͔q͕+e͕q͔+g͔s͕-g͕s͔+i͔t͕-i͕t͔)e15",
    ///     "+(-a͔m͕+a͕m͔-b͔o͕+b͕o͔-c͔p͕+c͕p͔-e͔r͕+e͕r͔-f͔s͕+f͕s͔-h͔t͕+h͕t͔)e16",
    ///     "+(-b͔e͕+b͕e͔-c͔f͕+c͕f͔-d͔g͕+d͕g͔-n͔q͕+n͕q͔-o͔r͕+o͕r͔-p͔s͕+p͕s͔)e23",
    ///     "+(+a͔e͕-a͕e͔-c͔h͕+c͕h͔-d͔i͕+d͕i͔+l͔q͕-l͕q͔+m͔r͕-m͕r͔-p͔t͕+p͕t͔)e24",
    ///     "+(+a͔f͕-a͕f͔+b͔h͕-b͕h͔-d͔j͕+d͕j͔-k͔q͕+k͕q͔+m͔s͕-m͕s͔+o͔t͕-o͕t͔)e25",
    ///     "+(+a͔g͕-a͕g͔+b͔i͕-b͕i͔+c͔j͕-c͕j͔-k͔r͕+k͕r͔-l͔s͕+l͕s͔-n͔t͕+n͕t͔)e26",
    ///     "+(-a͔b͕+a͕b͔-f͔h͕+f͕h͔-g͔i͕+g͕i͔-l͔n͕+l͕n͔-m͔o͕+m͕o͔-s͔t͕+s͕t͔)e34",
    ///     "+(-a͔c͕+a͕c͔+e͔h͕-e͕h͔-g͔j͕+g͕j͔+k͔n͕-k͕n͔-m͔p͕+m͕p͔+r͔t͕-r͕t͔)e35",
    ///     "+(-a͔d͕+a͕d͔+e͔i͕-e͕i͔+f͔j͕-f͕j͔+k͔o͕-k͕o͔+l͔p͕-l͕p͔-q͔t͕+q͕t͔)e36",
    ///     "+(-b͔c͕+b͕c͔-e͔f͕+e͕f͔-i͔j͕+i͕j͔-k͔l͕+k͕l͔-o͔p͕+o͕p͔-r͔s͕+r͕s͔)e45",
    ///     "+(-b͔d͕+b͕d͔-e͔g͕+e͕g͔+h͔j͕-h͕j͔-k͔m͕+k͕m͔+n͔p͕-n͕p͔+q͔s͕-q͕s͔)e46",
    ///     "+(-c͔d͕+c͕d͔-f͔g͕+f͕g͔-h͔i͕+h͕i͔-l͔m͕+l͕m͔-n͔o͕+n͕o͔-q͔r͕+q͕r͔)e56",
    ///     "+(+e͔j͕+e͕j͔-f͔i͕-f͕i͔+g͔h͕+g͕h͔+k͔p͕+k͕p͔-l͔o͕-l͕o͔+m͔n͕+m͕n͔)e3456",
    ///     "+(-b͔j͕-b͕j͔+c͔i͕+c͕i͔-d͔h͕-d͕h͔+k͔s͕+k͕s͔-l͔r͕-l͕r͔+m͔q͕+m͕q͔)e2465",
    ///     "+(+a͔j͕+a͕j͔-c͔g͕-c͕g͔+d͔f͕+d͕f͔+k͔t͕+k͕t͔-n͔r͕-n͕r͔+o͔q͕+o͕q͔)e2356",
    ///     "+(-a͔i͕-a͕i͔+b͔g͕+b͕g͔-d͔e͕-d͕e͔+l͔t͕+l͕t͔-n͔s͕-n͕s͔+p͔q͕+p͕q͔)e2364",
    ///     "+(+a͔h͕+a͕h͔-b͔f͕-b͕f͔+c͔e͕+c͕e͔+m͔t͕+m͕t͔-o͔s͕-o͕s͔+p͔r͕+p͕r͔)e2345",
    ///     "+(-b͔p͕-b͕p͔+c͔o͕+c͕o͔-d͔n͕-d͕n͔-e͔s͕-e͕s͔+f͔r͕+f͕r͔-g͔q͕-g͕q͔)e1456",
    ///     "+(+a͔p͕+a͕p͔-c͔m͕-c͕m͔+d͔l͕+d͕l͔-e͔t͕-e͕t͔+h͔r͕+h͕r͔-i͔q͕-i͕q͔)e1365",
    ///     "+(-a͔o͕-a͕o͔+b͔m͕+b͕m͔-d͔k͕-d͕k͔-f͔t͕-f͕t͔+h͔s͕+h͕s͔-j͔q͕-j͕q͔)e1346",
    ///     "+(+a͔n͕+a͕n͔-b͔l͕-b͕l͔+c͔k͕+c͕k͔-g͔t͕-g͕t͔+i͔s͕+i͕s͔-j͔r͕-j͕r͔)e1354",
    ///     "+(+a͔s͕+a͕s͔+b͔t͕+b͕t͔-f͔m͕-f͕m͔+g͔l͕+g͕l͔-h͔o͕-h͕o͔+i͔n͕+i͕n͔)e1256",
    ///     "+(-a͔r͕-a͕r͔+c͔t͕+c͕t͔+e͔m͕+e͕m͔-g͔k͕-g͕k͔-h͔p͕-h͕p͔+j͔n͕+j͕n͔)e1264",
    ///     "+(+a͔q͕+a͕q͔+d͔t͕+d͕t͔-e͔l͕-e͕l͔+f͔k͕+f͕k͔-i͔p͕-i͕p͔+j͔o͕+j͕o͔)e1245",
    ///     "+(-b͔r͕-b͕r͔-c͔s͕-c͕s͔+e͔o͕+e͕o͔+f͔p͕+f͕p͔-i͔k͕-i͕k͔-j͔l͕-j͕l͔)e1236",
    ///     "+(+b͔q͕+b͕q͔-d͔s͕-d͕s͔-e͔n͕-e͕n͔+g͔p͕+g͕p͔+h͔k͕+h͕k͔-j͔m͕-j͕m͔)e1253",
    ///     "+(+c͔q͕+c͕q͔+d͔r͕+d͕r͔-f͔n͕-f͕n͔-g͔o͕-g͕o͔+h͔l͕+h͕l͔+i͔m͕+i͕m͔)e1234",
    ///     "+(+h͔Ε͕+h͕Ε͔-i͔Δ͕-i͕Δ͔+j͔Γ͕+j͕Γ͔+n͔Ι͕+n͕Ι͔-o͔Θ͕-o͕Θ͔+p͔Η͕+p͕Η͔+q͔Μ͕+q͕Μ͔-r͔Λ͕-r͕Λ͔+s͔Κ͕+s͕Κ͔)e0465",
    ///     "+(-f͔Ε͕-f͕Ε͔+g͔Δ͕+g͕Δ͔-j͔Β͕-j͕Β͔-l͔Ι͕-l͕Ι͔+m͔Θ͕+m͕Θ͔-p͔Ζ͕-p͕Ζ͔+q͔Ξ͕+q͕Ξ͔-r͔Ν͕-r͕Ν͔+t͔Κ͕+t͕Κ͔)e0356",
    ///     "+(+e͔Ε͕+e͕Ε͔-g͔Γ͕-g͕Γ͔+i͔Β͕+i͕Β͔+k͔Ι͕+k͕Ι͔-m͔Η͕-m͕Η͔+o͔Ζ͕+o͕Ζ͔+q͔Ο͕+q͕Ο͔-s͔Ν͕-s͕Ν͔+t͔Λ͕+t͕Λ͔)e0364",
    ///     "+(-e͔Δ͕-e͕Δ͔+f͔Γ͕+f͕Γ͔-h͔Β͕-h͕Β͔-k͔Θ͕-k͕Θ͔+l͔Η͕+l͕Η͔-n͔Ζ͕-n͕Ζ͔+r͔Ο͕+r͕Ο͔-s͔Ξ͕-s͕Ξ͔+t͔Μ͕+t͕Μ͔)e0345",
    ///     "+(+c͔Ε͕+c͕Ε͔-d͔Δ͕-d͕Δ͔+j͔Α͕+j͕Α͔-l͔Μ͕-l͕Μ͔+m͔Λ͕+m͕Λ͔-n͔Ξ͕-n͕Ξ͔+o͔Ν͕+o͕Ν͔-s͔Ζ͕-s͕Ζ͔-t͔Η͕-t͕Η͔)e0265",
    ///     "+(-b͔Ε͕-b͕Ε͔+d͔Γ͕+d͕Γ͔-i͔Α͕-i͕Α͔+k͔Μ͕+k͕Μ͔-m͔Κ͕-m͕Κ͔-n͔Ο͕-n͕Ο͔+p͔Ν͕+p͕Ν͔+r͔Ζ͕+r͕Ζ͔-t͔Θ͕-t͕Θ͔)e0246",
    ///     "+(+b͔Δ͕+b͕Δ͔-c͔Γ͕-c͕Γ͔+h͔Α͕+h͕Α͔-k͔Λ͕-k͕Λ͔+l͔Κ͕+l͕Κ͔-o͔Ο͕-o͕Ο͔+p͔Ξ͕+p͕Ξ͔-q͔Ζ͕-q͕Ζ͔-t͔Ι͕-t͕Ι͔)e0254",
    ///     "+(+a͔Ε͕+a͕Ε͔-d͔Β͕-d͕Β͔+g͔Α͕+g͕Α͔+k͔Ξ͕+k͕Ξ͔+l͔Ο͕+l͕Ο͔-o͔Κ͕-o͕Κ͔-p͔Λ͕-p͕Λ͔+r͔Η͕+r͕Η͔+s͔Θ͕+s͕Θ͔)e0263",
    ///     "+(-a͔Δ͕-a͕Δ͔+c͔Β͕+c͕Β͔-f͔Α͕-f͕Α͔-k͔Ν͕-k͕Ν͔+m͔Ο͕+m͕Ο͔+n͔Κ͕+n͕Κ͔-p͔Μ͕-p͕Μ͔-q͔Η͕-q͕Η͔+s͔Ι͕+s͕Ι͔)e0235",
    ///     "+(+a͔Γ͕+a͕Γ͔-b͔Β͕-b͕Β͔+e͔Α͕+e͕Α͔-l͔Ν͕-l͕Ν͔-m͔Ξ͕-m͕Ξ͔+n͔Λ͕+n͕Λ͔+o͔Μ͕+o͕Μ͔-q͔Θ͕-q͕Θ͔-r͔Ι͕-r͕Ι͔)e0243",
    ///     "+(+c͔Ι͕+c͕Ι͔-d͔Θ͕-d͕Θ͔+f͔Μ͕+f͕Μ͔-g͔Λ͕-g͕Λ͔+h͔Ξ͕+h͕Ξ͔-i͔Ν͕-i͕Ν͔+p͔Α͕+p͕Α͔+s͔Β͕+s͕Β͔+t͔Γ͕+t͕Γ͔)e0156",
    ///     "+(-b͔Ι͕-b͕Ι͔+d͔Η͕+d͕Η͔-e͔Μ͕-e͕Μ͔+g͔Κ͕+g͕Κ͔+h͔Ο͕+h͕Ο͔-j͔Ν͕-j͕Ν͔-o͔Α͕-o͕Α͔-r͔Β͕-r͕Β͔+t͔Δ͕+t͕Δ͔)e0164",
    ///     "+(+b͔Θ͕+b͕Θ͔-c͔Η͕-c͕Η͔+e͔Λ͕+e͕Λ͔-f͔Κ͕-f͕Κ͔+i͔Ο͕+i͕Ο͔-j͔Ξ͕-j͕Ξ͔+n͔Α͕+n͕Α͔+q͔Β͕+q͕Β͔+t͔Ε͕+t͕Ε͔)e0145",
    ///     "+(+a͔Ι͕+a͕Ι͔-d͔Ζ͕-d͕Ζ͔-e͔Ξ͕-e͕Ξ͔-f͔Ο͕-f͕Ο͔+i͔Κ͕+i͕Κ͔+j͔Λ͕+j͕Λ͔+m͔Α͕+m͕Α͔-r͔Γ͕-r͕Γ͔-s͔Δ͕-s͕Δ͔)e0136",
    ///     "+(-a͔Θ͕-a͕Θ͔+c͔Ζ͕+c͕Ζ͔+e͔Ν͕+e͕Ν͔-g͔Ο͕-g͕Ο͔-h͔Κ͕-h͕Κ͔+j͔Μ͕+j͕Μ͔-l͔Α͕-l͕Α͔+q͔Γ͕+q͕Γ͔-s͔Ε͕-s͕Ε͔)e0153",
    ///     "+(+a͔Η͕+a͕Η͔-b͔Ζ͕-b͕Ζ͔+f͔Ν͕+f͕Ν͔+g͔Ξ͕+g͕Ξ͔-h͔Λ͕-h͕Λ͔-i͔Μ͕-i͕Μ͔+k͔Α͕+k͕Α͔+q͔Δ͕+q͕Δ͔+r͔Ε͕+r͕Ε͔)e0134",
    ///     "+(+a͔Μ͕+a͕Μ͔+b͔Ξ͕+b͕Ξ͔+c͔Ο͕+c͕Ο͔-g͔Ζ͕-g͕Ζ͔-i͔Η͕-i͕Η͔-j͔Θ͕-j͕Θ͔+m͔Β͕+m͕Β͔+o͔Γ͕+o͕Γ͔+p͔Δ͕+p͕Δ͔)e0162",
    ///     "+(-a͔Λ͕-a͕Λ͔-b͔Ν͕-b͕Ν͔+d͔Ο͕+d͕Ο͔+f͔Ζ͕+f͕Ζ͔+h͔Η͕+h͕Η͔-j͔Ι͕-j͕Ι͔-l͔Β͕-l͕Β͔-n͔Γ͕-n͕Γ͔+p͔Ε͕+p͕Ε͔)e0125",
    ///     "+(+a͔Κ͕+a͕Κ͔-c͔Ν͕-c͕Ν͔-d͔Ξ͕-d͕Ξ͔-e͔Ζ͕-e͕Ζ͔+h͔Θ͕+h͕Θ͔+i͔Ι͕+i͕Ι͔+k͔Β͕+k͕Β͔-n͔Δ͕-n͕Δ͔-o͔Ε͕-o͕Ε͔)e0142",
    ///     "+(+b͔Κ͕+b͕Κ͔+c͔Λ͕+c͕Λ͔+d͔Μ͕+d͕Μ͔-e͔Η͕-e͕Η͔-f͔Θ͕-f͕Θ͔-g͔Ι͕-g͕Ι͔+k͔Γ͕+k͕Γ͔+l͔Δ͕+l͕Δ͔+m͔Ε͕+m͕Ε͔)e0123",
    ///     "+(+a͔t͕-a͕t͔-b͔s͕+b͕s͔+c͔r͕-c͕r͔-d͔q͕+d͕q͔+e͔p͕-e͕p͔-f͔o͕+f͕o͔+g͔n͕-g͕n͔+h͔m͕-h͕m͔-i͔l͕+i͕l͔+j͔k͕-j͕k͔)e123456",
    ///     "+(+k͔Ο͕-k͕Ο͔-l͔Ξ͕+l͕Ξ͔+m͔Ν͕-m͕Ν͔+n͔Μ͕-n͕Μ͔-o͔Λ͕+o͕Λ͔+p͔Κ͕-p͕Κ͔-q͔Ι͕+q͕Ι͔+r͔Θ͕-r͕Θ͔-s͔Η͕+s͕Η͔+t͔Ζ͕-t͕Ζ͔)e023465",
    ///     "+(-e͔Ο͕+e͕Ο͔+f͔Ξ͕-f͕Ξ͔-g͔Ν͕+g͕Ν͔-h͔Μ͕+h͕Μ͔+i͔Λ͕-i͕Λ͔-j͔Κ͕+j͕Κ͔+q͔Ε͕-q͕Ε͔-r͔Δ͕+r͕Δ͔+s͔Γ͕-s͕Γ͔-t͔Β͕+t͕Β͔)e013456",
    ///     "+(+b͔Ο͕-b͕Ο͔-c͔Ξ͕+c͕Ξ͔+d͔Ν͕-d͕Ν͔+h͔Ι͕-h͕Ι͔-i͔Θ͕+i͕Θ͔+j͔Η͕-j͕Η͔-n͔Ε͕+n͕Ε͔+o͔Δ͕-o͕Δ͔-p͔Γ͕+p͕Γ͔+t͔Α͕-t͕Α͔)e012465",
    ///     "+(-a͔Ο͕+a͕Ο͔+c͔Μ͕-c͕Μ͔-d͔Λ͕+d͕Λ͔-f͔Ι͕+f͕Ι͔+g͔Θ͕-g͕Θ͔-j͔Ζ͕+j͕Ζ͔+l͔Ε͕-l͕Ε͔-m͔Δ͕+m͕Δ͔+p͔Β͕-p͕Β͔-s͔Α͕+s͕Α͔)e012356",
    ///     "+(+a͔Ξ͕-a͕Ξ͔-b͔Μ͕+b͕Μ͔+d͔Κ͕-d͕Κ͔+e͔Ι͕-e͕Ι͔-g͔Η͕+g͕Η͔+i͔Ζ͕-i͕Ζ͔-k͔Ε͕+k͕Ε͔+m͔Γ͕-m͕Γ͔-o͔Β͕+o͕Β͔+r͔Α͕-r͕Α͔)e012364",
    ///     "+(-a͔Ν͕+a͕Ν͔+b͔Λ͕-b͕Λ͔-c͔Κ͕+c͕Κ͔-e͔Θ͕+e͕Θ͔+f͔Η͕-f͕Η͔-h͔Ζ͕+h͕Ζ͔+k͔Δ͕-k͕Δ͔-l͔Γ͕+l͕Γ͔+n͔Β͕-n͕Β͔-q͔Α͕+q͕Α͔)e012345",
    /// ]);
    ///
    /// let point = Vee::point() << Vee::triple_motor();
    /// assert_eq!(point.basis_blades(), (Vee::point() + Vee::volume4()).basis_blades());
    /// format_eq!(point, [
    ///     "+2(-Ȧẇζ-Ḃẇη-Ċẇθ-Ḋẇι-Ėẇκ-Ḟẇλ-Ġẇμ-Ḣẇν-İẇξ-J̇ẇο-K̇ẇο̇+L̇ẇξ̇\
    ///         -Ṁẇν̇-Ṅẇμ̇+Ȯẇλ̇-Ṗẇκ̇+Q̇ẇι̇-Ṙẇθ̇+Ṡẇη̇-Ṫẇζ̇+Xẇẇ-Ẋαα̇-Ẋββ̇-Ẋγγ̇\
    ///         -Ẋδδ̇-Ẋεε̇+Ẋζζ̇+Ẋηη̇+Ẋθθ̇+Ẋιι̇+Ẋκκ̇+Ẋλλ̇+Ẋμμ̇+Ẋνν̇+Ẋξξ̇+Ẋοο̇\
    ///         -Yẇα̇+Ẏvα̇-Ẏβζ̇-Ẏβ̇ζ-Ẏγη̇-Ẏγ̇η-Ẏδθ̇-Ẏδ̇θ-Ẏει̇-Ẏε̇ι-Ẏκο+Ẏκ̇ο̇\
    ///         +Ẏλξ-Ẏλ̇ξ̇-Ẏμν+Ẏμ̇ν̇-Zẇβ̇+Żvβ̇+Żαζ̇+Żα̇ζ-Żγκ̇-Żγ̇κ-Żδλ̇-Żδ̇λ\
    ///         -Żεμ̇-Żε̇μ+Żηο-Żη̇ο̇-Żθξ+Żθ̇ξ̇+Żιν-Żι̇ν̇+vÐ̇γ̇+vØ̇δ̇+vÞ̇ε̇-ẇÐγ̇\
    ///         -ẇØδ̇-ẇÞε̇+Ð̇αη̇+Ð̇α̇η+Ð̇βκ̇+Ð̇β̇κ-Ð̇δν̇-Ð̇δ̇ν-Ð̇εξ̇-Ð̇ε̇ξ-Ð̇ζο+Ð̇ζ̇ο̇\
    ///         +Ð̇θμ-Ð̇θ̇μ̇-Ð̇ιλ+Ð̇ι̇λ̇+Ø̇αθ̇+Ø̇α̇θ+Ø̇βλ̇+Ø̇β̇λ+Ø̇γν̇+Ø̇γ̇ν-Ø̇εο̇-Ø̇ε̇ο\
    ///         +Ø̇ζξ-Ø̇ζ̇ξ̇-Ø̇ημ+Ø̇η̇μ̇+Ø̇ικ-Ø̇ι̇κ̇+Þ̇αι̇+Þ̇α̇ι+Þ̇βμ̇+Þ̇β̇μ+Þ̇γξ̇+Þ̇γ̇ξ\
    ///         +Þ̇δο̇+Þ̇δ̇ο-Þ̇ζν+Þ̇ζ̇ν̇+Þ̇ηλ-Þ̇η̇λ̇-Þ̇θκ+Þ̇θ̇κ̇)e01",
    ///     "+2(+Ȧẇβ+Ḃẇγ+Ċẇδ+Ḋẇε+Ėẇο̇-Ḟẇξ̇+Ġẇν̇+Ḣẇμ̇-İẇλ̇+J̇ẇκ̇-K̇ẇκ-L̇ẇλ\
    ///         -Ṁẇμ-Ṅẇν-Ȯẇξ-Ṗẇο-Q̇ẇε̇+Ṙẇδ̇-Ṡẇγ̇+Ṫẇβ̇+Xẇα̇-Ẋvα̇-Ẋβζ̇-Ẋβ̇ζ\
    ///         -Ẋγη̇-Ẋγ̇η-Ẋδθ̇-Ẋδ̇θ-Ẋει̇-Ẋε̇ι+Ẋκο-Ẋκ̇ο̇-Ẋλξ+Ẋλ̇ξ̇+Ẋμν-Ẋμ̇ν̇\
    ///         +Yẇẇ-Ẏαα̇+Ẏββ̇+Ẏγγ̇+Ẏδδ̇+Ẏεε̇-Ẏζζ̇-Ẏηη̇-Ẏθθ̇-Ẏιι̇+Ẏκκ̇+Ẏλλ̇\
    ///         +Ẏμμ̇+Ẏνν̇+Ẏξξ̇+Ẏοο̇-Zẇζ̇+Żvζ̇-Żαβ̇-Żα̇β-Żγο+Żγ̇ο̇+Żδξ-Żδ̇ξ̇\
    ///         -Żεν+Żε̇ν̇-Żηκ̇-Żη̇κ-Żθλ̇-Żθ̇λ-Żιμ̇-Żι̇μ+vÐ̇η̇+vØ̇θ̇+vÞ̇ι̇-ẇÐη̇\
    ///         -ẇØθ̇-ẇÞι̇-Ð̇αγ̇-Ð̇α̇γ+Ð̇βο-Ð̇β̇ο̇-Ð̇δμ+Ð̇δ̇μ̇+Ð̇ελ-Ð̇ε̇λ̇+Ð̇ζκ̇+Ð̇ζ̇κ\
    ///         -Ð̇θν̇-Ð̇θ̇ν-Ð̇ιξ̇-Ð̇ι̇ξ-Ø̇αδ̇-Ø̇α̇δ-Ø̇βξ+Ø̇β̇ξ̇+Ø̇γμ-Ø̇γ̇μ̇-Ø̇εκ+Ø̇ε̇κ̇\
    ///         +Ø̇ζλ̇+Ø̇ζ̇λ+Ø̇ην̇+Ø̇η̇ν-Ø̇ιο̇-Ø̇ι̇ο-Þ̇αε̇-Þ̇α̇ε+Þ̇βν-Þ̇β̇ν̇-Þ̇γλ+Þ̇γ̇λ̇\
    ///         +Þ̇δκ-Þ̇δ̇κ̇+Þ̇ζμ̇+Þ̇ζ̇μ+Þ̇ηξ̇+Þ̇η̇ξ+Þ̇θο̇+Þ̇θ̇ο)e02",
    ///     "+2(-Ȧẇα-Ḃẇο̇+Ċẇξ̇-Ḋẇν̇+Ėẇγ+Ḟẇδ+Ġẇε-Ḣẇι̇+İẇθ̇-J̇ẇη̇+K̇ẇη+L̇ẇθ\
    ///         +Ṁẇι+Ṅẇε̇-Ȯẇδ̇+Ṗẇγ̇-Q̇ẇν-Ṙẇξ-Ṡẇο-Ṫẇα̇+Xẇβ̇-Ẋvβ̇+Ẋαζ̇+Ẋα̇ζ\
    ///         -Ẋγκ̇-Ẋγ̇κ-Ẋδλ̇-Ẋδ̇λ-Ẋεμ̇-Ẋε̇μ-Ẋηο+Ẋη̇ο̇+Ẋθξ-Ẋθ̇ξ̇-Ẋιν+Ẋι̇ν̇\
    ///         +Yẇζ̇-Ẏvζ̇-Ẏαβ̇-Ẏα̇β+Ẏγο-Ẏγ̇ο̇-Ẏδξ+Ẏδ̇ξ̇+Ẏεν-Ẏε̇ν̇-Ẏηκ̇-Ẏη̇κ\
    ///         -Ẏθλ̇-Ẏθ̇λ-Ẏιμ̇-Ẏι̇μ+Zẇẇ+Żαα̇-Żββ̇+Żγγ̇+Żδδ̇+Żεε̇-Żζζ̇+Żηη̇\
    ///         +Żθθ̇+Żιι̇-Żκκ̇-Żλλ̇-Żμμ̇+Żνν̇+Żξξ̇+Żοο̇+vÐ̇κ̇+vØ̇λ̇+vÞ̇μ̇-ẇÐκ̇\
    ///         -ẇØλ̇-ẇÞμ̇-Ð̇αο+Ð̇α̇ο̇-Ð̇βγ̇-Ð̇β̇γ+Ð̇δι-Ð̇δ̇ι̇-Ð̇εθ+Ð̇ε̇θ̇-Ð̇ζη̇-Ð̇ζ̇η\
    ///         -Ð̇λν̇-Ð̇λ̇ν-Ð̇μξ̇-Ð̇μ̇ξ+Ø̇αξ-Ø̇α̇ξ̇-Ø̇βδ̇-Ø̇β̇δ-Ø̇γι+Ø̇γ̇ι̇+Ø̇εη-Ø̇ε̇η̇\
    ///         -Ø̇ζθ̇-Ø̇ζ̇θ+Ø̇κν̇+Ø̇κ̇ν-Ø̇μο̇-Ø̇μ̇ο-Þ̇αν+Þ̇α̇ν̇-Þ̇βε̇-Þ̇β̇ε+Þ̇γθ-Þ̇γ̇θ̇\
    ///         -Þ̇δη+Þ̇δ̇η̇-Þ̇ζι̇-Þ̇ζ̇ι+Þ̇κξ̇+Þ̇κ̇ξ+Þ̇λο̇+Þ̇λ̇ο)e03",
    ///     "+2(+Ȧẇο̇-Ḃẇα-Ċẇμ̇+Ḋẇλ̇-Ėẇβ+Ḟẇι̇-Ġẇθ̇+Ḣẇδ+İẇε+J̇ẇζ̇-K̇ẇζ-L̇ẇε̇\
    ///         +Ṁẇδ̇+Ṅẇθ+Ȯẇι-Ṗẇβ̇+Q̇ẇλ+Ṙẇμ+Ṡẇα̇-Ṫẇο+Xẇγ̇-Ẋvγ̇+Ẋαη̇+Ẋα̇η\
    ///         +Ẋβκ̇+Ẋβ̇κ-Ẋδν̇-Ẋδ̇ν-Ẋεξ̇-Ẋε̇ξ+Ẋζο-Ẋζ̇ο̇-Ẋθμ+Ẋθ̇μ̇+Ẋιλ-Ẋι̇λ̇\
    ///         +Yẇη̇-Ẏvη̇-Ẏαγ̇-Ẏα̇γ-Ẏβο+Ẏβ̇ο̇+Ẏδμ-Ẏδ̇μ̇-Ẏελ+Ẏε̇λ̇+Ẏζκ̇+Ẏζ̇κ\
    ///         -Ẏθν̇-Ẏθ̇ν-Ẏιξ̇-Ẏι̇ξ+Zẇκ̇-Żvκ̇+Żαο-Żα̇ο̇-Żβγ̇-Żβ̇γ-Żδι+Żδ̇ι̇\
    ///         +Żεθ-Żε̇θ̇-Żζη̇-Żζ̇η-Żλν̇-Żλ̇ν-Żμξ̇-Żμ̇ξ+vØ̇ν̇+vÞ̇ξ̇-ẇØν̇-ẇÞξ̇\
    ///         +ẇẇÐ+Ð̇αα̇+Ð̇ββ̇-Ð̇γγ̇+Ð̇δδ̇+Ð̇εε̇+Ð̇ζζ̇-Ð̇ηη̇+Ð̇θθ̇+Ð̇ιι̇-Ð̇κκ̇+Ð̇λλ̇\
    ///         +Ð̇μμ̇-Ð̇νν̇-Ð̇ξξ̇+Ð̇οο̇-Ø̇αμ+Ø̇α̇μ̇+Ø̇βι-Ø̇β̇ι̇-Ø̇γδ̇-Ø̇γ̇δ-Ø̇εζ+Ø̇ε̇ζ̇\
    ///         -Ø̇ηθ̇-Ø̇η̇θ-Ø̇κλ̇-Ø̇κ̇λ-Ø̇ξο̇-Ø̇ξ̇ο+Þ̇αλ-Þ̇α̇λ̇-Þ̇βθ+Þ̇β̇θ̇-Þ̇γε̇-Þ̇γ̇ε\
    ///         +Þ̇δζ-Þ̇δ̇ζ̇-Þ̇ηι̇-Þ̇η̇ι-Þ̇κμ̇-Þ̇κ̇μ+Þ̇νο̇+Þ̇ν̇ο)e04",
    ///     "+2(-Ȧẇξ̇+Ḃẇμ̇-Ċẇα-Ḋẇκ̇-Ėẇι̇-Ḟẇβ+Ġẇη̇-Ḣẇγ-İẇζ̇+J̇ẇε+K̇ẇε̇-L̇ẇζ\
    ///         -Ṁẇγ̇-Ṅẇη+Ȯẇβ̇+Ṗẇι-Q̇ẇκ-Ṙẇα̇+Ṡẇμ+Ṫẇξ+Xẇδ̇-Ẋvδ̇+Ẋαθ̇+Ẋα̇θ\
    ///         +Ẋβλ̇+Ẋβ̇λ+Ẋγν̇+Ẋγ̇ν-Ẋεο̇-Ẋε̇ο-Ẋζξ+Ẋζ̇ξ̇+Ẋημ-Ẋη̇μ̇-Ẋικ+Ẋι̇κ̇\
    ///         +Yẇθ̇-Ẏvθ̇-Ẏαδ̇-Ẏα̇δ+Ẏβξ-Ẏβ̇ξ̇-Ẏγμ+Ẏγ̇μ̇+Ẏεκ-Ẏε̇κ̇+Ẏζλ̇+Ẏζ̇λ\
    ///         +Ẏην̇+Ẏη̇ν-Ẏιο̇-Ẏι̇ο+Zẇλ̇-Żvλ̇-Żαξ+Żα̇ξ̇-Żβδ̇-Żβ̇δ+Żγι-Żγ̇ι̇\
    ///         -Żεη+Żε̇η̇-Żζθ̇-Żζ̇θ+Żκν̇+Żκ̇ν-Żμο̇-Żμ̇ο-vÐ̇ν̇+vÞ̇ο̇+ẇÐν̇-ẇÞο̇\
    ///         +ẇẇØ+Ð̇αμ-Ð̇α̇μ̇-Ð̇βι+Ð̇β̇ι̇-Ð̇γδ̇-Ð̇γ̇δ+Ð̇εζ-Ð̇ε̇ζ̇-Ð̇ηθ̇-Ð̇η̇θ-Ð̇κλ̇\
    ///         -Ð̇κ̇λ-Ð̇ξο̇-Ð̇ξ̇ο+Ø̇αα̇+Ø̇ββ̇+Ø̇γγ̇-Ø̇δδ̇+Ø̇εε̇+Ø̇ζζ̇+Ø̇ηη̇-Ø̇θθ̇+Ø̇ιι̇\
    ///         +Ø̇κκ̇-Ø̇λλ̇+Ø̇μμ̇-Ø̇νν̇+Ø̇ξξ̇-Ø̇οο̇-Þ̇ακ+Þ̇α̇κ̇+Þ̇βη-Þ̇β̇η̇-Þ̇γζ+Þ̇γ̇ζ̇\
    ///         -Þ̇δε̇-Þ̇δ̇ε-Þ̇θι̇-Þ̇θ̇ι-Þ̇λμ̇-Þ̇λ̇μ-Þ̇νξ̇-Þ̇ν̇ξ)e05",
    ///     "+2(+Ȧẇν̇-Ḃẇλ̇+Ċẇκ̇-Ḋẇα+Ėẇθ̇-Ḟẇη̇-Ġẇβ+Ḣẇζ̇-İẇγ-J̇ẇδ-K̇ẇδ̇+L̇ẇγ̇\
    ///         -Ṁẇζ-Ṅẇβ̇-Ȯẇη-Ṗẇθ+Q̇ẇα̇-Ṙẇκ-Ṡẇλ-Ṫẇν+Xẇε̇-Ẋvε̇+Ẋαι̇+Ẋα̇ι\
    ///         +Ẋβμ̇+Ẋβ̇μ+Ẋγξ̇+Ẋγ̇ξ+Ẋδο̇+Ẋδ̇ο+Ẋζν-Ẋζ̇ν̇-Ẋηλ+Ẋη̇λ̇+Ẋθκ-Ẋθ̇κ̇\
    ///         +Yẇι̇-Ẏvι̇-Ẏαε̇-Ẏα̇ε-Ẏβν+Ẏβ̇ν̇+Ẏγλ-Ẏγ̇λ̇-Ẏδκ+Ẏδ̇κ̇+Ẏζμ̇+Ẏζ̇μ\
    ///         +Ẏηξ̇+Ẏη̇ξ+Ẏθο̇+Ẏθ̇ο+Zẇμ̇-Żvμ̇+Żαν-Żα̇ν̇-Żβε̇-Żβ̇ε-Żγθ+Żγ̇θ̇\
    ///         +Żδη-Żδ̇η̇-Żζι̇-Żζ̇ι+Żκξ̇+Żκ̇ξ+Żλο̇+Żλ̇ο-vÐ̇ξ̇-vØ̇ο̇+ẇÐξ̇+ẇØο̇\
    ///         +ẇẇÞ-Ð̇αλ+Ð̇α̇λ̇+Ð̇βθ-Ð̇β̇θ̇-Ð̇γε̇-Ð̇γ̇ε-Ð̇δζ+Ð̇δ̇ζ̇-Ð̇ηι̇-Ð̇η̇ι-Ð̇κμ̇\
    ///         -Ð̇κ̇μ+Ð̇νο̇+Ð̇ν̇ο+Ø̇ακ-Ø̇α̇κ̇-Ø̇βη+Ø̇β̇η̇+Ø̇γζ-Ø̇γ̇ζ̇-Ø̇δε̇-Ø̇δ̇ε-Ø̇θι̇\
    ///         -Ø̇θ̇ι-Ø̇λμ̇-Ø̇λ̇μ-Ø̇νξ̇-Ø̇ν̇ξ+Þ̇αα̇+Þ̇ββ̇+Þ̇γγ̇+Þ̇δδ̇-Þ̇εε̇+Þ̇ζζ̇+Þ̇ηη̇\
    ///         +Þ̇θθ̇-Þ̇ιι̇+Þ̇κκ̇+Þ̇λλ̇-Þ̇μμ̇+Þ̇νν̇-Þ̇ξξ̇-Þ̇οο̇)e06",
    ///     "+2(+vẇα̇-ẇβζ̇+ẇβ̇ζ-ẇγη̇+ẇγ̇η-ẇδθ̇+ẇδ̇θ-ẇει̇\
    ///         +ẇε̇ι-ẇκο-ẇκ̇ο̇+ẇλξ+ẇλ̇ξ̇-ẇμν-ẇμ̇ν̇+ẇẇα)e12",
    ///     "+2(+vẇβ̇+ẇαζ̇-ẇα̇ζ-ẇγκ̇+ẇγ̇κ-ẇδλ̇+ẇδ̇λ-ẇεμ̇\
    ///         +ẇε̇μ+ẇηο+ẇη̇ο̇-ẇθξ-ẇθ̇ξ̇+ẇιν+ẇι̇ν̇+ẇẇβ)e13",
    ///     "+2(+vẇγ̇+ẇαη̇-ẇα̇η+ẇβκ̇-ẇβ̇κ-ẇδν̇+ẇδ̇ν-ẇεξ̇\
    ///         +ẇε̇ξ-ẇζο-ẇζ̇ο̇+ẇθμ+ẇθ̇μ̇-ẇιλ-ẇι̇λ̇+ẇẇγ)e14",
    ///     "+2(+vẇδ̇+ẇαθ̇-ẇα̇θ+ẇβλ̇-ẇβ̇λ+ẇγν̇-ẇγ̇ν-ẇεο̇\
    ///         +ẇε̇ο+ẇζξ+ẇζ̇ξ̇-ẇημ-ẇη̇μ̇+ẇικ+ẇι̇κ̇+ẇẇδ)e15",
    ///     "+2(+vẇε̇+ẇαι̇-ẇα̇ι+ẇβμ̇-ẇβ̇μ+ẇγξ̇-ẇγ̇ξ+ẇδο̇\
    ///         -ẇδ̇ο-ẇζν-ẇζ̇ν̇+ẇηλ+ẇη̇λ̇-ẇθκ-ẇθ̇κ̇+ẇẇε)e16",
    ///     "+2(+vẇζ̇-ẇαβ̇+ẇα̇β-ẇγο-ẇγ̇ο̇+ẇδξ+ẇδ̇ξ̇-ẇεν\
    ///         -ẇε̇ν̇-ẇηκ̇+ẇη̇κ-ẇθλ̇+ẇθ̇λ-ẇιμ̇+ẇι̇μ+ẇẇζ)e23",
    ///     "+2(+vẇη̇-ẇαγ̇+ẇα̇γ+ẇβο+ẇβ̇ο̇-ẇδμ-ẇδ̇μ̇+ẇελ\
    ///         +ẇε̇λ̇+ẇζκ̇-ẇζ̇κ-ẇθν̇+ẇθ̇ν-ẇιξ̇+ẇι̇ξ+ẇẇη)e24",
    ///     "+2(+vẇθ̇-ẇαδ̇+ẇα̇δ-ẇβξ-ẇβ̇ξ̇+ẇγμ+ẇγ̇μ̇-ẇεκ\
    ///         -ẇε̇κ̇+ẇζλ̇-ẇζ̇λ+ẇην̇-ẇη̇ν-ẇιο̇+ẇι̇ο+ẇẇθ)e25",
    ///     "+2(+vẇι̇-ẇαε̇+ẇα̇ε+ẇβν+ẇβ̇ν̇-ẇγλ-ẇγ̇λ̇+ẇδκ\
    ///         +ẇδ̇κ̇+ẇζμ̇-ẇζ̇μ+ẇηξ̇-ẇη̇ξ+ẇθο̇-ẇθ̇ο+ẇẇι)e26",
    ///     "+2(+vẇκ̇-ẇαο-ẇα̇ο̇-ẇβγ̇+ẇβ̇γ+ẇδι+ẇδ̇ι̇-ẇεθ\
    ///         -ẇε̇θ̇-ẇζη̇+ẇζ̇η-ẇλν̇+ẇλ̇ν-ẇμξ̇+ẇμ̇ξ+ẇẇκ)e34",
    ///     "+2(+vẇλ̇+ẇαξ+ẇα̇ξ̇-ẇβδ̇+ẇβ̇δ-ẇγι-ẇγ̇ι̇+ẇεη\
    ///         +ẇε̇η̇-ẇζθ̇+ẇζ̇θ+ẇκν̇-ẇκ̇ν-ẇμο̇+ẇμ̇ο+ẇẇλ)e35",
    ///     "+2(+vẇμ̇-ẇαν-ẇα̇ν̇-ẇβε̇+ẇβ̇ε+ẇγθ+ẇγ̇θ̇-ẇδη\
    ///         -ẇδ̇η̇-ẇζι̇+ẇζ̇ι+ẇκξ̇-ẇκ̇ξ+ẇλο̇-ẇλ̇ο+ẇẇμ)e36",
    ///     "+2(+vẇν̇-ẇαμ-ẇα̇μ̇+ẇβι+ẇβ̇ι̇-ẇγδ̇+ẇγ̇δ-ẇεζ\
    ///         -ẇε̇ζ̇-ẇηθ̇+ẇη̇θ-ẇκλ̇+ẇκ̇λ-ẇξο̇+ẇξ̇ο+ẇẇν)e45",
    ///     "+2(+vẇξ̇+ẇαλ+ẇα̇λ̇-ẇβθ-ẇβ̇θ̇-ẇγε̇+ẇγ̇ε+ẇδζ\
    ///         +ẇδ̇ζ̇-ẇηι̇+ẇη̇ι-ẇκμ̇+ẇκ̇μ+ẇνο̇-ẇν̇ο+ẇẇξ)e46",
    ///     "+2(+vẇο̇-ẇακ-ẇα̇κ̇+ẇβη+ẇβ̇η̇-ẇγζ-ẇγ̇ζ̇-ẇδε̇\
    ///         +ẇδ̇ε-ẇθι̇+ẇθ̇ι-ẇλμ̇+ẇλ̇μ-ẇνξ̇+ẇν̇ξ+ẇẇο)e56",
    ///     "+(+vvẇ+ẇαα+ẇα̇α̇+ẇββ+ẇβ̇β̇+ẇγγ+ẇγ̇γ̇+ẇδδ+ẇδ̇δ̇+ẇεε+ẇε̇ε̇+ẇζζ+ẇζ̇ζ̇+ẇηη+ẇη̇η̇+ẇθθ\
    ///        +ẇθ̇θ̇+ẇιι+ẇι̇ι̇+ẇκκ+ẇκ̇κ̇+ẇλλ+ẇλ̇λ̇+ẇμμ+ẇμ̇μ̇+ẇνν+ẇν̇ν̇+ẇξξ+ẇξ̇ξ̇+ẇοο+ẇο̇ο̇+ẇẇẇ)e123456",
    ///     "+(+2Ȧẇζ̇+2Ḃẇη̇+2Ċẇθ̇+2Ḋẇι̇+2Ėẇκ̇+2Ḟẇλ̇+2Ġẇμ̇+2Ḣẇν̇+2İẇξ̇+2J̇ẇο̇-2K̇ẇο+2L̇ẇξ-2Ṁẇν\
    ///        -2Ṅẇμ+2Ȯẇλ-2Ṗẇκ+2Q̇ẇι-2Ṙẇθ+2Ṡẇη-2Ṫẇζ-2Xvẇ+Ẋvv+Ẋẇẇ-Ẋαα+Ẋα̇α̇-Ẋββ+Ẋβ̇β̇\
    ///        -Ẋγγ+Ẋγ̇γ̇-Ẋδδ+Ẋδ̇δ̇-Ẋεε+Ẋε̇ε̇+Ẋζζ-Ẋζ̇ζ̇+Ẋηη-Ẋη̇η̇+Ẋθθ-Ẋθ̇θ̇+Ẋιι-Ẋι̇ι̇+Ẋκκ-Ẋκ̇κ̇\
    ///        +Ẋλλ-Ẋλ̇λ̇+Ẋμμ-Ẋμ̇μ̇+Ẋνν-Ẋν̇ν̇+Ẋξξ-Ẋξ̇ξ̇+Ẋοο-Ẋο̇ο̇-2Yẇα+2Ẏvα-2Ẏβζ+2Ẏβ̇ζ̇-2Ẏγη\
    ///        +2Ẏγ̇η̇-2Ẏδθ+2Ẏδ̇θ̇-2Ẏει+2Ẏε̇ι̇+2Ẏκο̇+2Ẏκ̇ο-2Ẏλξ̇-2Ẏλ̇ξ+2Ẏμν̇+2Ẏμ̇ν-2Zẇβ+2Żvβ\
    ///        +2Żαζ-2Żα̇ζ̇-2Żγκ+2Żγ̇κ̇-2Żδλ+2Żδ̇λ̇-2Żεμ+2Żε̇μ̇-2Żηο̇-2Żη̇ο+2Żθξ̇+2Żθ̇ξ-2Żιν̇\
    ///        -2Żι̇ν+2vÐ̇γ+2vØ̇δ+2vÞ̇ε-2ẇÐγ-2ẇØδ-2ẇÞε+2Ð̇αη-2Ð̇α̇η̇+2Ð̇βκ-2Ð̇β̇κ̇-2Ð̇δν+2Ð̇δ̇ν̇\
    ///        -2Ð̇εξ+2Ð̇ε̇ξ̇+2Ð̇ζο̇+2Ð̇ζ̇ο-2Ð̇θμ̇-2Ð̇θ̇μ+2Ð̇ιλ̇+2Ð̇ι̇λ+2Ø̇αθ-2Ø̇α̇θ̇+2Ø̇βλ-2Ø̇β̇λ̇+2Ø̇γν\
    ///        -2Ø̇γ̇ν̇-2Ø̇εο+2Ø̇ε̇ο̇-2Ø̇ζξ̇-2Ø̇ζ̇ξ+2Ø̇ημ̇+2Ø̇η̇μ-2Ø̇ικ̇-2Ø̇ι̇κ+2Þ̇αι-2Þ̇α̇ι̇+2Þ̇βμ-2Þ̇β̇μ̇\
    ///        +2Þ̇γξ-2Þ̇γ̇ξ̇+2Þ̇δο-2Þ̇δ̇ο̇+2Þ̇ζν̇+2Þ̇ζ̇ν-2Þ̇ηλ̇-2Þ̇η̇λ+2Þ̇θκ̇+2Þ̇θ̇κ)e023465",
    ///     "+(-2Ȧẇβ̇-2Ḃẇγ̇-2Ċẇδ̇-2Ḋẇε̇+2Ėẇο-2Ḟẇξ+2Ġẇν+2Ḣẇμ-2İẇλ+2J̇ẇκ+2K̇ẇκ̇+2L̇ẇλ̇+2Ṁẇμ̇\
    ///        +2Ṅẇν̇+2Ȯẇξ̇+2Ṗẇο̇-2Q̇ẇε+2Ṙẇδ-2Ṡẇγ+2Ṫẇβ+2Xẇα-2Ẋvα-2Ẋβζ+2Ẋβ̇ζ̇-2Ẋγη+2Ẋγ̇\
    ///        η̇-2Ẋδθ+2Ẋδ̇θ̇-2Ẋει+2Ẋε̇ι̇-2Ẋκο̇-2Ẋκ̇ο+2Ẋλξ̇+2Ẋλ̇ξ-2Ẋμν̇-2Ẋμ̇ν-2Yvẇ+Ẏvv+Ẏẇẇ\
    ///        -Ẏαα+Ẏα̇α̇+Ẏββ-Ẏβ̇β̇+Ẏγγ-Ẏγ̇γ̇+Ẏδδ-Ẏδ̇δ̇+Ẏεε-Ẏε̇ε̇-Ẏζζ+Ẏζ̇ζ̇-Ẏηη+Ẏη̇η̇-Ẏθθ+Ẏθ̇θ̇-\
    ///        Ẏιι+Ẏι̇ι̇+Ẏκκ-Ẏκ̇κ̇+Ẏλλ-Ẏλ̇λ̇+Ẏμμ-Ẏμ̇μ̇+Ẏνν-Ẏν̇ν̇+Ẏξξ-Ẏξ̇ξ̇+Ẏοο-Ẏο̇ο̇-2Zẇζ+2Żvζ\
    ///        -2Żαβ+2Żα̇β̇+2Żγο̇+2Żγ̇ο-2Żδξ̇-2Żδ̇ξ+2Żεν̇+2Żε̇ν-2Żηκ+2Żη̇κ̇-2Żθλ+2Żθ̇λ̇-2Żιμ\
    ///        +2Żι̇μ̇+2vÐ̇η+2vØ̇θ+2vÞ̇ι-2ẇÐη-2ẇØθ-2ẇÞι-2Ð̇αγ+2Ð̇α̇γ̇-2Ð̇βο̇-2Ð̇β̇ο+2Ð̇δμ̇+2Ð̇δ̇μ\
    ///        -2Ð̇ελ̇-2Ð̇ε̇λ+2Ð̇ζκ-2Ð̇ζ̇κ̇-2Ð̇θν+2Ð̇θ̇ν̇-2Ð̇ιξ+2Ð̇ι̇ξ̇-2Ø̇αδ+2Ø̇α̇δ̇+2Ø̇βξ̇+2Ø̇β̇ξ-2Ø̇γμ̇\
    ///        -2Ø̇γ̇μ+2Ø̇εκ̇+2Ø̇ε̇κ+2Ø̇ζλ-2Ø̇ζ̇λ̇+2Ø̇ην-2Ø̇η̇ν̇-2Ø̇ιο+2Ø̇ι̇ο̇-2Þ̇αε+2Þ̇α̇ε̇-2Þ̇βν̇-2Þ̇β̇ν\
    ///        +2Þ̇γλ̇+2Þ̇γ̇λ-2Þ̇δκ̇-2Þ̇δ̇κ+2Þ̇ζμ-2Þ̇ζ̇μ̇+2Þ̇ηξ-2Þ̇η̇ξ̇+2Þ̇θο-2Þ̇θ̇ο̇)e013456",
    ///     "+(+2Ȧẇα̇-2Ḃẇο+2Ċẇξ-2Ḋẇν-2Ėẇγ̇-2Ḟẇδ̇-2Ġẇε̇-2Ḣẇι+2İẇθ-2J̇ẇη-2K̇ẇη̇-2L̇ẇθ̇-2Ṁẇι̇\
    ///        +2Ṅẇε-2Ȯẇδ+2Ṗẇγ+2Q̇ẇν̇+2Ṙẇξ̇+2Ṡẇο̇-2Ṫẇα+2Xẇβ-2Ẋvβ+2Ẋαζ-2Ẋα̇ζ̇-2Ẋγκ+2Ẋγ̇\
    ///        κ̇-2Ẋδλ+2Ẋδ̇λ̇-2Ẋεμ+2Ẋε̇μ̇+2Ẋηο̇+2Ẋη̇ο-2Ẋθξ̇-2Ẋθ̇ξ+2Ẋιν̇+2Ẋι̇ν+2Yẇζ-2Ẏvζ-2Ẏ\
    ///        αβ+2Ẏα̇β̇-2Ẏγο̇-2Ẏγ̇ο+2Ẏδξ̇+2Ẏδ̇ξ-2Ẏεν̇-2Ẏε̇ν-2Ẏηκ+2Ẏη̇κ̇-2Ẏθλ+2Ẏθ̇λ̇-2Ẏιμ+2Ẏ\
    ///        ι̇μ̇-2Zvẇ+Żvv+Żẇẇ+Żαα-Żα̇α̇-Żββ+Żβ̇β̇+Żγγ-Żγ̇γ̇+Żδδ-Żδ̇δ̇+Żεε-Żε̇ε̇-Żζζ+Żζ̇ζ̇+Ż\
    ///        ηη-Żη̇η̇+Żθθ-Żθ̇θ̇+Żιι-Żι̇ι̇-Żκκ+Żκ̇κ̇-Żλλ+Żλ̇λ̇-Żμμ+Żμ̇μ̇+Żνν-Żν̇ν̇+Żξξ-Żξ̇ξ̇+Żο\
    ///        ο-Żο̇ο̇+2vÐ̇κ+2vØ̇λ+2vÞ̇μ-2ẇÐκ-2ẇØλ-2ẇÞμ+2Ð̇αο̇+2Ð̇α̇ο-2Ð̇βγ+2Ð̇β̇γ̇-2Ð̇δι̇-2Ð̇δ̇ι\
    ///        +2Ð̇εθ̇+2Ð̇ε̇θ-2Ð̇ζη+2Ð̇ζ̇η̇-2Ð̇λν+2Ð̇λ̇ν̇-2Ð̇μξ+2Ð̇μ̇ξ̇-2Ø̇αξ̇-2Ø̇α̇ξ-2Ø̇βδ+2Ø̇β̇δ̇+2Ø̇γι̇\
    ///        +2Ø̇γ̇ι-2Ø̇εη̇-2Ø̇ε̇η-2Ø̇ζθ+2Ø̇ζ̇θ̇+2Ø̇κν-2Ø̇κ̇ν̇-2Ø̇μο+2Ø̇μ̇ο̇+2Þ̇αν̇+2Þ̇α̇ν-2Þ̇βε+2Þ̇β̇ε̇\
    ///        -2Þ̇γθ̇-2Þ̇γ̇θ+2Þ̇δη̇+2Þ̇δ̇η-2Þ̇ζι+2Þ̇ζ̇ι̇+2Þ̇κξ-2Þ̇κ̇ξ̇+2Þ̇λο-2Þ̇λ̇ο̇)e012465",
    ///     "+(+2Ȧẇο+2Ḃẇα̇-2Ċẇμ+2Ḋẇλ+2Ėẇβ̇+2Ḟẇι-2Ġẇθ-2Ḣẇδ̇-2İẇε̇+2J̇ẇζ+2K̇ẇζ̇-2L̇ẇε+2Ṁẇδ\
    ///        -2Ṅẇθ̇-2Ȯẇι̇-2Ṗẇβ-2Q̇ẇλ̇-2Ṙẇμ̇+2Ṡẇα+2Ṫẇο̇+2Xẇγ-2Ẋvγ+2Ẋαη-2Ẋα̇η̇+2Ẋβκ-2Ẋβ̇\
    ///        κ̇-2Ẋδν+2Ẋδ̇ν̇-2Ẋεξ+2Ẋε̇ξ̇-2Ẋζο̇-2Ẋζ̇ο+2Ẋθμ̇+2Ẋθ̇μ-2Ẋιλ̇-2Ẋι̇λ+2Yẇη-2Ẏvη-2Ẏ\
    ///        αγ+2Ẏα̇γ̇+2Ẏβο̇+2Ẏβ̇ο-2Ẏδμ̇-2Ẏδ̇μ+2Ẏελ̇+2Ẏε̇λ+2Ẏζκ-2Ẏζ̇κ̇-2Ẏθν+2Ẏθ̇ν̇-2Ẏιξ+2Ẏ\
    ///        ι̇ξ̇+2Zẇκ-2Żvκ-2Żαο̇-2Żα̇ο-2Żβγ+2Żβ̇γ̇+2Żδι̇+2Żδ̇ι-2Żεθ̇-2Żε̇θ-2Żζη+2Żζ̇η̇-2Ż\
    ///        λν+2Żλ̇ν̇-2Żμξ+2Żμ̇ξ̇-2vẇÐ+2vØ̇ν+2vÞ̇ξ+vvÐ̇-2ẇØν-2ẇÞξ+ẇẇÐ̇+Ð̇αα-Ð̇α̇α̇+Ð̇ββ-Ð̇β̇\
    ///        β̇-Ð̇γγ+Ð̇γ̇γ̇+Ð̇δδ-Ð̇δ̇δ̇+Ð̇εε-Ð̇ε̇ε̇+Ð̇ζζ-Ð̇ζ̇ζ̇-Ð̇ηη+Ð̇η̇η̇+Ð̇θθ-Ð̇θ̇θ̇+Ð̇ιι-Ð̇ι̇ι̇-Ð̇κκ+Ð̇κ̇κ̇\
    ///        +Ð̇λλ-Ð̇λ̇λ̇+Ð̇μμ-Ð̇μ̇μ̇-Ð̇νν+Ð̇ν̇ν̇-Ð̇ξξ+Ð̇ξ̇ξ̇+Ð̇οο-Ð̇ο̇ο̇+2Ø̇αμ̇+2Ø̇α̇μ-2Ø̇βι̇-2Ø̇β̇ι-2Ø̇γδ\
    ///        +2Ø̇γ̇δ̇+2Ø̇εζ̇+2Ø̇ε̇ζ-2Ø̇ηθ+2Ø̇η̇θ̇-2Ø̇κλ+2Ø̇κ̇λ̇-2Ø̇ξο+2Ø̇ξ̇ο̇-2Þ̇αλ̇-2Þ̇α̇λ+2Þ̇βθ̇+2Þ̇β̇θ\
    ///        -2Þ̇γε+2Þ̇γ̇ε̇-2Þ̇δζ̇-2Þ̇δ̇ζ-2Þ̇ηι+2Þ̇η̇ι̇-2Þ̇κμ+2Þ̇κ̇μ̇+2Þ̇νο-2Þ̇ν̇ο̇)e012356",
    ///     "+(-2Ȧẇξ+2Ḃẇμ+2Ċẇα̇-2Ḋẇκ-2Ėẇι+2Ḟẇβ̇+2Ġẇη+2Ḣẇγ̇-2İẇζ-2J̇ẇε̇+2K̇ẇε+2L̇ẇζ̇-2Ṁẇγ\
    ///        +2Ṅẇη̇+2Ȯẇβ-2Ṗẇι̇+2Q̇ẇκ̇-2Ṙẇα-2Ṡẇμ̇-2Ṫẇξ̇+2Xẇδ-2Ẋvδ+2Ẋαθ-2Ẋα̇θ̇+2Ẋβλ-2Ẋβ̇\
    ///        λ̇+2Ẋγν-2Ẋγ̇ν̇-2Ẋεο+2Ẋε̇ο̇+2Ẋζξ̇+2Ẋζ̇ξ-2Ẋημ̇-2Ẋη̇μ+2Ẋικ̇+2Ẋι̇κ+2Yẇθ-2Ẏvθ-2Ẏ\
    ///        αδ+2Ẏα̇δ̇-2Ẏβξ̇-2Ẏβ̇ξ+2Ẏγμ̇+2Ẏγ̇μ-2Ẏεκ̇-2Ẏε̇κ+2Ẏζλ-2Ẏζ̇λ̇+2Ẏην-2Ẏη̇ν̇-2Ẏιο+2Ẏ\
    ///        ι̇ο̇+2Zẇλ-2Żvλ+2Żαξ̇+2Żα̇ξ-2Żβδ+2Żβ̇δ̇-2Żγι̇-2Żγ̇ι+2Żεη̇+2Żε̇η-2Żζθ+2Żζ̇θ̇+2Ż\
    ///        κν-2Żκ̇ν̇-2Żμο+2Żμ̇ο̇-2vẇØ-2vÐ̇ν+2vÞ̇ο+vvØ̇+2ẇÐν-2ẇÞο+ẇẇØ̇-2Ð̇αμ̇-2Ð̇α̇μ+2Ð̇βι̇\
    ///        +2Ð̇β̇ι-2Ð̇γδ+2Ð̇γ̇δ̇-2Ð̇εζ̇-2Ð̇ε̇ζ-2Ð̇ηθ+2Ð̇η̇θ̇-2Ð̇κλ+2Ð̇κ̇λ̇-2Ð̇ξο+2Ð̇ξ̇ο̇+Ø̇αα-Ø̇α̇α̇+Ø̇\
    ///        ββ-Ø̇β̇β̇+Ø̇γγ-Ø̇γ̇γ̇-Ø̇δδ+Ø̇δ̇δ̇+Ø̇εε-Ø̇ε̇ε̇+Ø̇ζζ-Ø̇ζ̇ζ̇+Ø̇ηη-Ø̇η̇η̇-Ø̇θθ+Ø̇θ̇θ̇+Ø̇ιι-Ø̇ι̇ι̇+Ø̇κ\
    ///        κ-Ø̇κ̇κ̇-Ø̇λλ+Ø̇λ̇λ̇+Ø̇μμ-Ø̇μ̇μ̇-Ø̇νν+Ø̇ν̇ν̇+Ø̇ξξ-Ø̇ξ̇ξ̇-Ø̇οο+Ø̇ο̇ο̇+2Þ̇ακ̇+2Þ̇α̇κ-2Þ̇βη̇-2Þ̇β̇η\
    ///        +2Þ̇γζ̇+2Þ̇γ̇ζ-2Þ̇δε+2Þ̇δ̇ε̇-2Þ̇θι+2Þ̇θ̇ι̇-2Þ̇λμ+2Þ̇λ̇μ̇-2Þ̇νξ+2Þ̇ν̇ξ̇)e012364",
    ///     "+(+2Ȧẇν-2Ḃẇλ+2Ċẇκ+2Ḋẇα̇+2Ėẇθ-2Ḟẇη+2Ġẇβ̇+2Ḣẇζ+2İẇγ̇+2J̇ẇδ̇-2K̇ẇδ+2L̇ẇγ+2Ṁẇζ̇\
    ///        -2Ṅẇβ+2Ȯẇη̇+2Ṗẇθ̇+2Q̇ẇα+2Ṙẇκ̇+2Ṡẇλ̇+2Ṫẇν̇+2Xẇε-2Ẋvε+2Ẋαι-2Ẋα̇ι̇+2Ẋβμ-2Ẋβ̇\
    ///        μ̇+2Ẋγξ-2Ẋγ̇ξ̇+2Ẋδο-2Ẋδ̇ο̇-2Ẋζν̇-2Ẋζ̇ν+2Ẋηλ̇+2Ẋη̇λ-2Ẋθκ̇-2Ẋθ̇κ+2Yẇι-2Ẏvι-2Ẏ\
    ///        αε+2Ẏα̇ε̇+2Ẏβν̇+2Ẏβ̇ν-2Ẏγλ̇-2Ẏγ̇λ+2Ẏδκ̇+2Ẏδ̇κ+2Ẏζμ-2Ẏζ̇μ̇+2Ẏηξ-2Ẏη̇ξ̇+2Ẏθο-2Ẏ\
    ///        θ̇ο̇+2Zẇμ-2Żvμ-2Żαν̇-2Żα̇ν-2Żβε+2Żβ̇ε̇+2Żγθ̇+2Żγ̇θ-2Żδη̇-2Żδ̇η-2Żζι+2Żζ̇ι̇+2Ż\
    ///        κξ-2Żκ̇ξ̇+2Żλο-2Żλ̇ο̇-2vẇÞ-2vÐ̇ξ-2vØ̇ο+vvÞ̇+2ẇÐξ+2ẇØο+ẇẇÞ̇+2Ð̇αλ̇+2Ð̇α̇λ-2Ð̇βθ̇\
    ///        -2Ð̇β̇θ-2Ð̇γε+2Ð̇γ̇ε̇+2Ð̇δζ̇+2Ð̇δ̇ζ-2Ð̇ηι+2Ð̇η̇ι̇-2Ð̇κμ+2Ð̇κ̇μ̇+2Ð̇νο-2Ð̇ν̇ο̇-2Ø̇ακ̇-2Ø̇α̇κ\
    ///        +2Ø̇βη̇+2Ø̇β̇η-2Ø̇γζ̇-2Ø̇γ̇ζ-2Ø̇δε+2Ø̇δ̇ε̇-2Ø̇θι+2Ø̇θ̇ι̇-2Ø̇λμ+2Ø̇λ̇μ̇-2Ø̇νξ+2Ø̇ν̇ξ̇+Þ̇αα-\
    ///        Þ̇α̇α̇+Þ̇ββ-Þ̇β̇β̇+Þ̇γγ-Þ̇γ̇γ̇+Þ̇δδ-Þ̇δ̇δ̇-Þ̇εε+Þ̇ε̇ε̇+Þ̇ζζ-Þ̇ζ̇ζ̇+Þ̇ηη-Þ̇η̇η̇+Þ̇θθ-Þ̇θ̇θ̇-Þ̇ιι+Þ̇\
    ///        ι̇ι̇+Þ̇κκ-Þ̇κ̇κ̇+Þ̇λλ-Þ̇λ̇λ̇-Þ̇μμ+Þ̇μ̇μ̇+Þ̇νν-Þ̇ν̇ν̇-Þ̇ξξ+Þ̇ξ̇ξ̇-Þ̇οο+Þ̇ο̇ο̇)e012345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_motor() -> Self {
        Self::scalar() + Self::volume4() + Self::plane() + Self::point()
    }
    /// The multivector of single rotoreflector $`f_{r1} \equiv v^5_0 + v_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let single_rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(single_rotoreflector.basis_blades(), Vee::single_rotoreflector().basis_blades());
    /// format_eq!(single_rotoreflector, [
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔η͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔θ͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔ι͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔κ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔λ͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔μ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ν͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ξ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+y͔κ͕-z͔η͕+ð͔ζ͕)e234",
    ///     "+(+y͔λ͕-z͔θ͕+ø͔ζ͕)e235",
    ///     "+(+y͔μ͕-z͔ι͕+þ͔ζ͕)e236",
    ///     "+(+y͔ν͕-ð͔θ͕+ø͔η͕)e245",
    ///     "+(+y͔ξ͕-ð͔ι͕+þ͔η͕)e246",
    ///     "+(+y͔ο͕-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(+z͔ν͕-ð͔λ͕+ø͔κ͕)e345",
    ///     "+(+z͔ξ͕-ð͔μ͕+þ͔κ͕)e346",
    ///     "+(+z͔ο͕-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::normal() + Self::volume_displacement()
    }
    /// The multivector of double rotoreflector $`f_{r2} \equiv v^5_0 + v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let double_rotoreflector = Vee::normal().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(double_rotoreflector.basis_blades(), Vee::double_rotoreflector().basis_blades());
    /// format_eq!(double_rotoreflector, [
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕-ð͔ο͕̇+ø͔ξ͕̇-þ͔ν͕̇)e123",
    ///     "+(+x͔η͕-y͔γ͕+z͔ο͕̇+ð͔α͕-ø͔μ͕̇+þ͔λ͕̇)e124",
    ///     "+(+x͔θ͕-y͔δ͕-z͔ξ͕̇+ð͔μ͕̇+ø͔α͕-þ͔κ͕̇)e125",
    ///     "+(+x͔ι͕-y͔ε͕+z͔ν͕̇-ð͔λ͕̇+ø͔κ͕̇+þ͔α͕)e126",
    ///     "+(+x͔κ͕-y͔ο͕̇-z͔γ͕+ð͔β͕+ø͔ι͕̇-þ͔θ͕̇)e134",
    ///     "+(+x͔λ͕+y͔ξ͕̇-z͔δ͕-ð͔ι͕̇+ø͔β͕+þ͔η͕̇)e135",
    ///     "+(+x͔μ͕-y͔ν͕̇-z͔ε͕+ð͔θ͕̇-ø͔η͕̇+þ͔β͕)e136",
    ///     "+(+x͔ν͕-y͔μ͕̇+z͔ι͕̇-ð͔δ͕+ø͔γ͕-þ͔ζ͕̇)e145",
    ///     "+(+x͔ξ͕+y͔λ͕̇-z͔θ͕̇-ð͔ε͕+ø͔ζ͕̇+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-y͔κ͕̇+z͔η͕̇-ð͔ζ͕̇-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔ο͕̇+y͔κ͕-z͔η͕+ð͔ζ͕-ø͔ε͕̇+þ͔δ͕̇)e234",
    ///     "+(-x͔ξ͕̇+y͔λ͕-z͔θ͕+ð͔ε͕̇+ø͔ζ͕-þ͔γ͕̇)e235",
    ///     "+(+x͔ν͕̇+y͔μ͕-z͔ι͕-ð͔δ͕̇+ø͔γ͕̇+þ͔ζ͕)e236",
    ///     "+(+x͔μ͕̇+y͔ν͕-z͔ε͕̇-ð͔θ͕+ø͔η͕+þ͔β͕̇)e245",
    ///     "+(-x͔λ͕̇+y͔ξ͕+z͔δ͕̇-ð͔ι͕-ø͔β͕̇+þ͔η͕)e246",
    ///     "+(+x͔κ͕̇+y͔ο͕-z͔γ͕̇+ð͔β͕̇-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(-x͔ι͕̇+y͔ε͕̇+z͔ν͕-ð͔λ͕+ø͔κ͕-þ͔α͕̇)e345",
    ///     "+(+x͔θ͕̇-y͔δ͕̇+z͔ξ͕-ð͔μ͕+ø͔α͕̇+þ͔κ͕)e346",
    ///     "+(-x͔η͕̇+y͔γ͕̇+z͔ο͕-ð͔α͕̇-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+x͔ζ͕̇-y͔β͕̇+z͔α͕̇+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    ///     "+(+y͔α͕̇+z͔β͕̇+ð͔γ͕̇+ø͔δ͕̇+þ͔ε͕̇)e23456",
    ///     "+(-x͔α͕̇+z͔ζ͕̇+ð͔η͕̇+ø͔θ͕̇+þ͔ι͕̇)e13465",
    ///     "+(-x͔β͕̇-y͔ζ͕̇+ð͔κ͕̇+ø͔λ͕̇+þ͔μ͕̇)e12456",
    ///     "+(-x͔γ͕̇-y͔η͕̇-z͔κ͕̇+ø͔ν͕̇+þ͔ξ͕̇)e12365",
    ///     "+(-x͔δ͕̇-y͔θ͕̇-z͔λ͕̇-ð͔ν͕̇+þ͔ο͕̇)e12346",
    ///     "+(-x͔ε͕̇-y͔ι͕̇-z͔μ͕̇-ð͔ξ͕̇-ø͔ο͕̇)e12354",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::normal() + Self::volume_displacement() + Self::line_displacement()
    }
    /// The multivector of transflector $`f_t \equiv v^5 + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+v͕ø͔e5",
    ///     "+v͕þ͔e6",
    ///     "+(+X͕y͔-Y͕x͔)e012",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+Y͕z͔-Z͕y͔)e023",
    ///     "+(+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+Ø͕þ͔-Þ͕ø͔)e056",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume5() + Self::volume_moment()
    }
    /// The multivector of single flector $`f_1 \equiv v^5 + v`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let single_flector = Vee::volume5().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(single_flector.basis_blades(), Vee::single_flector().basis_blades());
    /// format_eq!(single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+W͔ζ͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+W͔η͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+W͔θ͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+W͔ι͕+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+W͔κ͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+W͔λ͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+W͔μ͕+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+W͔ν͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+W͔ξ͕+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+W͔ο͕+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔η͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔θ͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔ι͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔κ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔λ͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔μ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ν͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ξ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+y͔κ͕-z͔η͕+ð͔ζ͕)e234",
    ///     "+(+y͔λ͕-z͔θ͕+ø͔ζ͕)e235",
    ///     "+(+y͔μ͕-z͔ι͕+þ͔ζ͕)e236",
    ///     "+(+y͔ν͕-ð͔θ͕+ø͔η͕)e245",
    ///     "+(+y͔ξ͕-ð͔ι͕+þ͔η͕)e246",
    ///     "+(+y͔ο͕-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(+z͔ν͕-ð͔λ͕+ø͔κ͕)e345",
    ///     "+(+z͔ξ͕-ð͔μ͕+þ͔κ͕)e346",
    ///     "+(+z͔ο͕-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::volume5() + Self::volume()
    }
    /// The multivector of simple double flector $`f_{s2} \equiv v^5 + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let simple_double_flector = Vee::volume5().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades());
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+Q͕̇þ͔-Ṙ͕ø͔+Ṡ͕ð͔-Ṫ͕z͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(-Ṅ͕þ͔+Ȯ͕ø͔-Ṗ͕ð͔+Ṫ͕y͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+L͕̇þ͔-Ṁ͕ø͔+Ṗ͕z͔-Ṡ͕y͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(-K͕̇þ͔+Ṁ͕ð͔-Ȯ͕z͔+Ṙ͕y͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+K͕̇ø͔-L͕̇ð͔+Ṅ͕z͔-Q͕̇y͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+Ḣ͕þ͔-İ͕ø͔+J͕̇ð͔-Ṫ͕x͔+W͔ζ͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-Ḟ͕þ͔+Ġ͕ø͔-J͕̇z͔+Ṡ͕x͔+W͔η͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Ė͕þ͔-Ġ͕ð͔+İ͕z͔-Ṙ͕x͔+W͔θ͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-Ė͕ø͔+Ḟ͕ð͔-Ḣ͕z͔+Q͕̇x͔+W͔ι͕+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+Ċ͕þ͔-Ḋ͕ø͔+J͕̇y͔-Ṗ͕x͔+W͔κ͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(-Ḃ͕þ͔+Ḋ͕ð͔-İ͕y͔+Ȯ͕x͔+W͔λ͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ḃ͕ø͔-Ċ͕ð͔+Ḣ͕y͔-Ṅ͕x͔+W͔μ͕+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+Ȧ͕þ͔-Ḋ͕z͔+Ġ͕y͔-Ṁ͕x͔+W͔ν͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(-Ȧ͕ø͔+Ċ͕z͔-Ḟ͕y͔+L͕̇x͔+W͔ξ͕+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+Ȧ͕ð͔-Ḃ͕z͔+Ė͕y͔-K͕̇x͔+W͔ο͕+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔η͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔θ͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔ι͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔κ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔λ͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔μ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ν͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ξ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+y͔κ͕-z͔η͕+ð͔ζ͕)e234",
    ///     "+(+y͔λ͕-z͔θ͕+ø͔ζ͕)e235",
    ///     "+(+y͔μ͕-z͔ι͕+þ͔ζ͕)e236",
    ///     "+(+y͔ν͕-ð͔θ͕+ø͔η͕)e245",
    ///     "+(+y͔ξ͕-ð͔ι͕+þ͔η͕)e246",
    ///     "+(+y͔ο͕-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(+z͔ν͕-ð͔λ͕+ø͔κ͕)e345",
    ///     "+(+z͔ξ͕-ð͔μ͕+þ͔κ͕)e346",
    ///     "+(+z͔ο͕-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    ///     "+(+Ȧ͕z͔+Ḃ͕ð͔+Ċ͕ø͔+Ḋ͕þ͔)e03456",
    ///     "+(-Ȧ͕y͔+Ė͕ð͔+Ḟ͕ø͔+Ġ͕þ͔)e02465",
    ///     "+(-Ḃ͕y͔-Ė͕z͔+Ḣ͕ø͔+İ͕þ͔)e02356",
    ///     "+(-Ċ͕y͔-Ḟ͕z͔-Ḣ͕ð͔+J͕̇þ͔)e02364",
    ///     "+(-Ḋ͕y͔-Ġ͕z͔-İ͕ð͔-J͕̇ø͔)e02345",
    ///     "+(+Ȧ͕x͔+K͕̇ð͔+L͕̇ø͔+Ṁ͕þ͔)e01456",
    ///     "+(+Ḃ͕x͔-K͕̇z͔+Ṅ͕ø͔+Ȯ͕þ͔)e01365",
    ///     "+(+Ċ͕x͔-L͕̇z͔-Ṅ͕ð͔+Ṗ͕þ͔)e01346",
    ///     "+(+Ḋ͕x͔-Ṁ͕z͔-Ȯ͕ð͔-Ṗ͕ø͔)e01354",
    ///     "+(+Ė͕x͔+K͕̇y͔+Q͕̇ø͔+Ṙ͕þ͔)e01256",
    ///     "+(+Ḟ͕x͔+L͕̇y͔-Q͕̇ð͔+Ṡ͕þ͔)e01264",
    ///     "+(+Ġ͕x͔+Ṁ͕y͔-Ṙ͕ð͔-Ṡ͕ø͔)e01245",
    ///     "+(+Ḣ͕x͔+Ṅ͕y͔+Q͕̇z͔+Ṫ͕þ͔)e01236",
    ///     "+(+İ͕x͔+Ȯ͕y͔+Ṙ͕z͔-Ṫ͕ø͔)e01253",
    ///     "+(+J͕̇x͔+Ṗ͕y͔+Ṡ͕z͔+Ṫ͕ð͔)e01234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_flector() -> Self {
        Self::volume5() + Self::volume() + Self::line_moment()
    }
    /// The multivector of double flector $`f_2 \equiv v^5 + v + \ell`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let double_flector = Vee::volume5().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// format_eq!(double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+Q͕̇þ͔-Ṙ͕ø͔+Ṡ͕ð͔-Ṫ͕z͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(-Ṅ͕þ͔+Ȯ͕ø͔-Ṗ͕ð͔+Ṫ͕y͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+L͕̇þ͔-Ṁ͕ø͔+Ṗ͕z͔-Ṡ͕y͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(-K͕̇þ͔+Ṁ͕ð͔-Ȯ͕z͔+Ṙ͕y͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+K͕̇ø͔-L͕̇ð͔+Ṅ͕z͔-Q͕̇y͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+Ḣ͕þ͔-İ͕ø͔+J͕̇ð͔-Ṫ͕x͔+W͔ζ͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-Ḟ͕þ͔+Ġ͕ø͔-J͕̇z͔+Ṡ͕x͔+W͔η͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Ė͕þ͔-Ġ͕ð͔+İ͕z͔-Ṙ͕x͔+W͔θ͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-Ė͕ø͔+Ḟ͕ð͔-Ḣ͕z͔+Q͕̇x͔+W͔ι͕+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+Ċ͕þ͔-Ḋ͕ø͔+J͕̇y͔-Ṗ͕x͔+W͔κ͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(-Ḃ͕þ͔+Ḋ͕ð͔-İ͕y͔+Ȯ͕x͔+W͔λ͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ḃ͕ø͔-Ċ͕ð͔+Ḣ͕y͔-Ṅ͕x͔+W͔μ͕+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+Ȧ͕þ͔-Ḋ͕z͔+Ġ͕y͔-Ṁ͕x͔+W͔ν͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(-Ȧ͕ø͔+Ċ͕z͔-Ḟ͕y͔+L͕̇x͔+W͔ξ͕+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+Ȧ͕ð͔-Ḃ͕z͔+Ė͕y͔-K͕̇x͔+W͔ο͕+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕-ð͔ο͕̇+ø͔ξ͕̇-þ͔ν͕̇)e123",
    ///     "+(+x͔η͕-y͔γ͕+z͔ο͕̇+ð͔α͕-ø͔μ͕̇+þ͔λ͕̇)e124",
    ///     "+(+x͔θ͕-y͔δ͕-z͔ξ͕̇+ð͔μ͕̇+ø͔α͕-þ͔κ͕̇)e125",
    ///     "+(+x͔ι͕-y͔ε͕+z͔ν͕̇-ð͔λ͕̇+ø͔κ͕̇+þ͔α͕)e126",
    ///     "+(+x͔κ͕-y͔ο͕̇-z͔γ͕+ð͔β͕+ø͔ι͕̇-þ͔θ͕̇)e134",
    ///     "+(+x͔λ͕+y͔ξ͕̇-z͔δ͕-ð͔ι͕̇+ø͔β͕+þ͔η͕̇)e135",
    ///     "+(+x͔μ͕-y͔ν͕̇-z͔ε͕+ð͔θ͕̇-ø͔η͕̇+þ͔β͕)e136",
    ///     "+(+x͔ν͕-y͔μ͕̇+z͔ι͕̇-ð͔δ͕+ø͔γ͕-þ͔ζ͕̇)e145",
    ///     "+(+x͔ξ͕+y͔λ͕̇-z͔θ͕̇-ð͔ε͕+ø͔ζ͕̇+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-y͔κ͕̇+z͔η͕̇-ð͔ζ͕̇-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔ο͕̇+y͔κ͕-z͔η͕+ð͔ζ͕-ø͔ε͕̇+þ͔δ͕̇)e234",
    ///     "+(-x͔ξ͕̇+y͔λ͕-z͔θ͕+ð͔ε͕̇+ø͔ζ͕-þ͔γ͕̇)e235",
    ///     "+(+x͔ν͕̇+y͔μ͕-z͔ι͕-ð͔δ͕̇+ø͔γ͕̇+þ͔ζ͕)e236",
    ///     "+(+x͔μ͕̇+y͔ν͕-z͔ε͕̇-ð͔θ͕+ø͔η͕+þ͔β͕̇)e245",
    ///     "+(-x͔λ͕̇+y͔ξ͕+z͔δ͕̇-ð͔ι͕-ø͔β͕̇+þ͔η͕)e246",
    ///     "+(+x͔κ͕̇+y͔ο͕-z͔γ͕̇+ð͔β͕̇-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(-x͔ι͕̇+y͔ε͕̇+z͔ν͕-ð͔λ͕+ø͔κ͕-þ͔α͕̇)e345",
    ///     "+(+x͔θ͕̇-y͔δ͕̇+z͔ξ͕-ð͔μ͕+ø͔α͕̇+þ͔κ͕)e346",
    ///     "+(-x͔η͕̇+y͔γ͕̇+z͔ο͕-ð͔α͕̇-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+x͔ζ͕̇-y͔β͕̇+z͔α͕̇+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    ///     "+(+y͔α͕̇+z͔β͕̇+ð͔γ͕̇+ø͔δ͕̇+þ͔ε͕̇)e23456",
    ///     "+(-x͔α͕̇+z͔ζ͕̇+ð͔η͕̇+ø͔θ͕̇+þ͔ι͕̇)e13465",
    ///     "+(-x͔β͕̇-y͔ζ͕̇+ð͔κ͕̇+ø͔λ͕̇+þ͔μ͕̇)e12456",
    ///     "+(-x͔γ͕̇-y͔η͕̇-z͔κ͕̇+ø͔ν͕̇+þ͔ξ͕̇)e12365",
    ///     "+(-x͔δ͕̇-y͔θ͕̇-z͔λ͕̇-ð͔ν͕̇+þ͔ο͕̇)e12346",
    ///     "+(-x͔ε͕̇-y͔ι͕̇-z͔μ͕̇-ð͔ξ͕̇-ø͔ο͕̇)e12354",
    ///     "+(+Ȧ͕z͔+Ḃ͕ð͔+Ċ͕ø͔+Ḋ͕þ͔+W͔α͕̇)e03456",
    ///     "+(-Ȧ͕y͔+Ė͕ð͔+Ḟ͕ø͔+Ġ͕þ͔+W͔β͕̇)e02465",
    ///     "+(-Ḃ͕y͔-Ė͕z͔+Ḣ͕ø͔+İ͕þ͔+W͔γ͕̇)e02356",
    ///     "+(-Ċ͕y͔-Ḟ͕z͔-Ḣ͕ð͔+J͕̇þ͔+W͔δ͕̇)e02364",
    ///     "+(-Ḋ͕y͔-Ġ͕z͔-İ͕ð͔-J͕̇ø͔+W͔ε͕̇)e02345",
    ///     "+(+Ȧ͕x͔+K͕̇ð͔+L͕̇ø͔+Ṁ͕þ͔+W͔ζ͕̇)e01456",
    ///     "+(+Ḃ͕x͔-K͕̇z͔+Ṅ͕ø͔+Ȯ͕þ͔+W͔η͕̇)e01365",
    ///     "+(+Ċ͕x͔-L͕̇z͔-Ṅ͕ð͔+Ṗ͕þ͔+W͔θ͕̇)e01346",
    ///     "+(+Ḋ͕x͔-Ṁ͕z͔-Ȯ͕ð͔-Ṗ͕ø͔+W͔ι͕̇)e01354",
    ///     "+(+Ė͕x͔+K͕̇y͔+Q͕̇ø͔+Ṙ͕þ͔+W͔κ͕̇)e01256",
    ///     "+(+Ḟ͕x͔+L͕̇y͔-Q͕̇ð͔+Ṡ͕þ͔+W͔λ͕̇)e01264",
    ///     "+(+Ġ͕x͔+Ṁ͕y͔-Ṙ͕ð͔-Ṡ͕ø͔+W͔μ͕̇)e01245",
    ///     "+(+Ḣ͕x͔+Ṅ͕y͔+Q͕̇z͔+Ṫ͕þ͔+W͔ν͕̇)e01236",
    ///     "+(+İ͕x͔+Ȯ͕y͔+Ṙ͕z͔-Ṫ͕ø͔+W͔ξ͕̇)e01253",
    ///     "+(+J͕̇x͔+Ṗ͕y͔+Ṡ͕z͔+Ṫ͕ð͔+W͔ο͕̇)e01234",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_flector() -> Self {
        Self::volume5() + Self::volume() + Self::line()
    }
    /// The multivector of triple flector $`f_3 \equiv v^5 + v + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP6 as Vee};
    ///
    /// let triple_flector = Vee::volume5().lhs() * Vee::triple_motor().rhs();
    ///
    /// assert_eq!(triple_flector.basis_blades(), Vee::triple_flector().basis_blades());
    /// format_eq!(triple_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔ζ͕-ð͔η͕-ø͔θ͕-þ͔ι͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔ζ͕-ð͔κ͕-ø͔λ͕-þ͔μ͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔η͕+z͔κ͕-ø͔ν͕-þ͔ξ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔θ͕+z͔λ͕+ð͔ν͕-þ͔ο͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔ι͕+z͔μ͕+ð͔ξ͕+ø͔ο͕)e6",
    ///     "+(+Q͕̇þ͔-Ṙ͕ø͔+Ṡ͕ð͔-Ṫ͕z͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(-Ṅ͕þ͔+Ȯ͕ø͔-Ṗ͕ð͔+Ṫ͕y͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+L͕̇þ͔-Ṁ͕ø͔+Ṗ͕z͔-Ṡ͕y͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(-K͕̇þ͔+Ṁ͕ð͔-Ȯ͕z͔+Ṙ͕y͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+K͕̇ø͔-L͕̇ð͔+Ṅ͕z͔-Q͕̇y͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+Ḣ͕þ͔-İ͕ø͔+J͕̇ð͔-Ṫ͕x͔+W͔ζ͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-Ḟ͕þ͔+Ġ͕ø͔-J͕̇z͔+Ṡ͕x͔+W͔η͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Ė͕þ͔-Ġ͕ð͔+İ͕z͔-Ṙ͕x͔+W͔θ͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(-Ė͕ø͔+Ḟ͕ð͔-Ḣ͕z͔+Q͕̇x͔+W͔ι͕+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+Ċ͕þ͔-Ḋ͕ø͔+J͕̇y͔-Ṗ͕x͔+W͔κ͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(-Ḃ͕þ͔+Ḋ͕ð͔-İ͕y͔+Ȯ͕x͔+W͔λ͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Ḃ͕ø͔-Ċ͕ð͔+Ḣ͕y͔-Ṅ͕x͔+W͔μ͕+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+Ȧ͕þ͔-Ḋ͕z͔+Ġ͕y͔-Ṁ͕x͔+W͔ν͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(-Ȧ͕ø͔+Ċ͕z͔-Ḟ͕y͔+L͕̇x͔+W͔ξ͕+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+Ȧ͕ð͔-Ḃ͕z͔+Ė͕y͔-K͕̇x͔+W͔ο͕+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+x͔ζ͕-y͔β͕+z͔α͕-ð͔ο͕̇+ø͔ξ͕̇-þ͔ν͕̇)e123",
    ///     "+(+x͔η͕-y͔γ͕+z͔ο͕̇+ð͔α͕-ø͔μ͕̇+þ͔λ͕̇)e124",
    ///     "+(+x͔θ͕-y͔δ͕-z͔ξ͕̇+ð͔μ͕̇+ø͔α͕-þ͔κ͕̇)e125",
    ///     "+(+x͔ι͕-y͔ε͕+z͔ν͕̇-ð͔λ͕̇+ø͔κ͕̇+þ͔α͕)e126",
    ///     "+(+x͔κ͕-y͔ο͕̇-z͔γ͕+ð͔β͕+ø͔ι͕̇-þ͔θ͕̇)e134",
    ///     "+(+x͔λ͕+y͔ξ͕̇-z͔δ͕-ð͔ι͕̇+ø͔β͕+þ͔η͕̇)e135",
    ///     "+(+x͔μ͕-y͔ν͕̇-z͔ε͕+ð͔θ͕̇-ø͔η͕̇+þ͔β͕)e136",
    ///     "+(+x͔ν͕-y͔μ͕̇+z͔ι͕̇-ð͔δ͕+ø͔γ͕-þ͔ζ͕̇)e145",
    ///     "+(+x͔ξ͕+y͔λ͕̇-z͔θ͕̇-ð͔ε͕+ø͔ζ͕̇+þ͔γ͕)e146",
    ///     "+(+x͔ο͕-y͔κ͕̇+z͔η͕̇-ð͔ζ͕̇-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔ο͕̇+y͔κ͕-z͔η͕+ð͔ζ͕-ø͔ε͕̇+þ͔δ͕̇)e234",
    ///     "+(-x͔ξ͕̇+y͔λ͕-z͔θ͕+ð͔ε͕̇+ø͔ζ͕-þ͔γ͕̇)e235",
    ///     "+(+x͔ν͕̇+y͔μ͕-z͔ι͕-ð͔δ͕̇+ø͔γ͕̇+þ͔ζ͕)e236",
    ///     "+(+x͔μ͕̇+y͔ν͕-z͔ε͕̇-ð͔θ͕+ø͔η͕+þ͔β͕̇)e245",
    ///     "+(-x͔λ͕̇+y͔ξ͕+z͔δ͕̇-ð͔ι͕-ø͔β͕̇+þ͔η͕)e246",
    ///     "+(+x͔κ͕̇+y͔ο͕-z͔γ͕̇+ð͔β͕̇-ø͔ι͕+þ͔θ͕)e256",
    ///     "+(-x͔ι͕̇+y͔ε͕̇+z͔ν͕-ð͔λ͕+ø͔κ͕-þ͔α͕̇)e345",
    ///     "+(+x͔θ͕̇-y͔δ͕̇+z͔ξ͕-ð͔μ͕+ø͔α͕̇+þ͔κ͕)e346",
    ///     "+(-x͔η͕̇+y͔γ͕̇+z͔ο͕-ð͔α͕̇-ø͔μ͕+þ͔λ͕)e356",
    ///     "+(+x͔ζ͕̇-y͔β͕̇+z͔α͕̇+ð͔ο͕-ø͔ξ͕+þ͔ν͕)e456",
    ///     "+(+ẇ͕x͔+y͔α͕̇+z͔β͕̇+ð͔γ͕̇+ø͔δ͕̇+þ͔ε͕̇)e23456",
    ///     "+(+ẇ͕y͔-x͔α͕̇+z͔ζ͕̇+ð͔η͕̇+ø͔θ͕̇+þ͔ι͕̇)e13465",
    ///     "+(+ẇ͕z͔-x͔β͕̇-y͔ζ͕̇+ð͔κ͕̇+ø͔λ͕̇+þ͔μ͕̇)e12456",
    ///     "+(+ẇ͕ð͔-x͔γ͕̇-y͔η͕̇-z͔κ͕̇+ø͔ν͕̇+þ͔ξ͕̇)e12365",
    ///     "+(+ẇ͕ø͔-x͔δ͕̇-y͔θ͕̇-z͔λ͕̇-ð͔ν͕̇+þ͔ο͕̇)e12346",
    ///     "+(+ẇ͕þ͔-x͔ε͕̇-y͔ι͕̇-z͔μ͕̇-ð͔ξ͕̇-ø͔ο͕̇)e12354",
    ///     "+(+Ȧ͕z͔+Ḃ͕ð͔+Ċ͕ø͔+Ḋ͕þ͔+W͔α͕̇+Ẋ͕y͔-Ẏ͕x͔)e03456",
    ///     "+(-Ȧ͕y͔+Ė͕ð͔+Ḟ͕ø͔+Ġ͕þ͔+W͔β͕̇+Ẋ͕z͔-Ż͕x͔)e02465",
    ///     "+(-Ḃ͕y͔-Ė͕z͔+Ḣ͕ø͔+İ͕þ͔+W͔γ͕̇+Ẋ͕ð͔-x͔Ð͕̇)e02356",
    ///     "+(-Ċ͕y͔-Ḟ͕z͔-Ḣ͕ð͔+J͕̇þ͔+W͔δ͕̇+Ẋ͕ø͔-x͔Ø͕̇)e02364",
    ///     "+(-Ḋ͕y͔-Ġ͕z͔-İ͕ð͔-J͕̇ø͔+W͔ε͕̇+Ẋ͕þ͔-x͔Þ͕̇)e02345",
    ///     "+(+Ȧ͕x͔+K͕̇ð͔+L͕̇ø͔+Ṁ͕þ͔+W͔ζ͕̇+Ẏ͕z͔-Ż͕y͔)e01456",
    ///     "+(+Ḃ͕x͔-K͕̇z͔+Ṅ͕ø͔+Ȯ͕þ͔+W͔η͕̇+Ẏ͕ð͔-y͔Ð͕̇)e01365",
    ///     "+(+Ċ͕x͔-L͕̇z͔-Ṅ͕ð͔+Ṗ͕þ͔+W͔θ͕̇+Ẏ͕ø͔-y͔Ø͕̇)e01346",
    ///     "+(+Ḋ͕x͔-Ṁ͕z͔-Ȯ͕ð͔-Ṗ͕ø͔+W͔ι͕̇+Ẏ͕þ͔-y͔Þ͕̇)e01354",
    ///     "+(+Ė͕x͔+K͕̇y͔+Q͕̇ø͔+Ṙ͕þ͔+W͔κ͕̇+Ż͕ð͔-z͔Ð͕̇)e01256",
    ///     "+(+Ḟ͕x͔+L͕̇y͔-Q͕̇ð͔+Ṡ͕þ͔+W͔λ͕̇+Ż͕ø͔-z͔Ø͕̇)e01264",
    ///     "+(+Ġ͕x͔+Ṁ͕y͔-Ṙ͕ð͔-Ṡ͕ø͔+W͔μ͕̇+Ż͕þ͔-z͔Þ͕̇)e01245",
    ///     "+(+Ḣ͕x͔+Ṅ͕y͔+Q͕̇z͔+Ṫ͕þ͔+W͔ν͕̇+Ð͕̇ø͔-Ø͕̇ð͔)e01236",
    ///     "+(+İ͕x͔+Ȯ͕y͔+Ṙ͕z͔-Ṫ͕ø͔+W͔ξ͕̇+Ð͕̇þ͔-Þ͕̇ð͔)e01253",
    ///     "+(+J͕̇x͔+Ṗ͕y͔+Ṡ͕z͔+Ṫ͕ð͔+W͔ο͕̇+Ø͕̇þ͔-Þ͕̇ø͔)e01234",
    ///     "+(+W͔ẇ͕+Ẋ͕x͔+Ẏ͕y͔+Ż͕z͔+Ð͕̇ð͔+Ø͕̇ø͔+Þ͕̇þ͔)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_flector() -> Self {
        Self::volume5() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
}

/// The named entities of the PGA with embedded dimension $`N = 7`$ (exploratory).
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
    /// The multivector of norm $`n \equiv s + v + S`$.
    ///
    /// Quadvector $`v`$ does not square to a scalar, therefore $`n`$ is **not** a Study number.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let quadvector_norm_squared = Vee::volume().norm_squared();
    ///
    /// assert_eq!(quadvector_norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(quadvector_norm_squared, [
    ///     "+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk+ll+mm+nn+oo+pp+qq\
    ///      +rr+ss+tt+uu+áá+ää+åå+ææ+çç+éé+ëë+íí+ïï+ññ+óó+öö+úú+üü",
    ///     "+2(-cç-dé-eë+gt+hu+iá-jq-kr-ls)e1234",
    ///     "+2(-Jo+Kn-Lm-Ml+Nk-Oj-Tæ+Uå-tÆ+uÅ-Áä-Äá-Çñ+Éï-Ëí-Íë+Ïé-Ñç)e0123",
    ///     "+2(+Go-Hn+Im+Mi-Nh+Og+Qæ-Rå+Sä+qÆ-rÅ+sÄ-Çú+Éö-Ëó-Óë+Öé-Úç)e0124",
    ///     "+2(-Fo+Hl-Ik-Ki+Lh-Of-Pæ+Rá-Su-Us-pÆ+rÁ-Çü+Íö-Ïó-Óï+Öí-Üç)e0125",
    ///     "+2(+Fn-Gl+Ij+Ji-Lg+Nf+På-Qá+St+Ts+pÅ-qÁ-Éü+Íú-Ñó-Óñ+Úí-Üé)e0126",
    ///     "+2(-Fm+Gk-Hj-Jh+Kg-Mf-Pä+Qu-Rt-Tr+Uq-pÄ-Ëü+Ïú-Ñö-Öñ+Úï-Üë)e0127",
    ///     "+2(-Co+Dn-Em-Me+Nd-Oc+Qñ-Rï+Sí+Tú-Uö+qÑ-rÏ+sÍ+tÚ-uÖ+Áó+Óá)e0134",
    ///     "+2(+Bo-Dl+Ek+Ke-Ld+Ob-Pñ+Rë-Sé+Tü-pÑ+rË-sÉ+tÜ-Äö+Åó+Óå-Öä)e0135",
    ///     "+2(-Bn+Cl-Ej-Je+Lc-Nb+Pï-Që+Sç+Uü+pÏ-qË+sÇ+uÜ-Äú+Æó+Óæ-Úä)e0136",
    ///     "+2(+Bm-Ck+Dj+Jd-Kc+Mb-Pí+Qé-Rç-pÍ+qÉ-rÇ+Áü-Åú+Æö+Öæ-Úå+Üá)e0137",
    ///     "+2(-Ao+Di-Eh-He+Id-Oa-Pú-Qü+Uë-pÚ-qÜ+uË-Áé+Äï-Åí-Éá-Íå+Ïä)e0145",
    ///     "+2(+An-Ci+Eg+Ge-Ic+Na+Pö-Rü-Të+pÖ-rÜ-tË+Áç+Äñ-Æí+Çá-Íæ+Ñä)e0146",
    ///     "+2(-Am+Ch-Dg-Gd+Hc-Ma-Pó-Sü+Té-Uç-pÓ-sÜ+tÉ-uÇ+Åñ-Æï-Ïæ+Ñå)e0147",
    ///     "+2(-Al+Bi-Ef-Fe+Ib-La+Qö+Rú-Tï-Uñ+qÖ+rÚ-tÏ-uÑ+Åç+Æé+Çå+Éæ)e0156",
    ///     "+2(+Ak-Bh+Df+Fd-Hb+Ka-Qó+Sú+Tí-qÓ+sÚ+tÍ-Áñ-Äç+Æë-Çä+Ëæ-Ñá)e0157",
    ///     "+2(-Aj+Bg-Cf-Fc+Gb-Ja-Ró-Sö+Uí-rÓ-sÖ+uÍ+Áï-Äé-Åë-Éä-Ëå+Ïá)e0167",
    ///     "+2(-Cæ+Då-Eä-Gñ+Hï-Ií-Jú+Kö-Ló-cÆ+dÅ-eÄ-gÑ+hÏ-iÍ-jÚ+kÖ-lÓ)e0234",
    ///     "+2(+Bæ-Dá+Eu+Fñ-Hë+Ié-Jü+Mö-Nó+Ue+bÆ-dÁ+fÑ-hË+iÉ-jÜ+mÖ-nÓ)e0235",
    ///     "+2(-Bå+Cá-Et-Fï+Gë-Iç-Kü+Mú-Oó-Te-bÅ+cÁ-fÏ+gË-iÇ-kÜ+mÚ-oÓ)e0236",
    ///     "+2(+Bä-Cu+Dt+Fí-Gé+Hç-Lü+Nú-Oö+Td-Uc+bÄ+fÍ-gÉ+hÇ-lÜ+nÚ-oÖ)e0237",
    ///     "+2(-Aæ+Ds-Er+Fú+Gü-Kë+Lé-Mï+Ní-Re+Sd-aÆ+fÚ+gÜ-kË+lÉ-mÏ+nÍ)e0245",
    ///     "+2(+Aå-Cs+Eq-Fö+Hü+Jë-Lç-Mñ+Oí+Qe-Sc+aÅ-fÖ+hÜ+jË-lÇ-mÑ+oÍ)e0246",
    ///     "+2(-Aä+Cr-Dq+Fó+Iü-Jé+Kç-Nñ+Oï-Qd+Rc-aÄ+fÓ+iÜ-jÉ+kÇ-nÑ+oÏ)e0247",
    ///     "+2(-Aá+Bs-Ep-Gö-Hú+Jï+Kñ-Nç-Oé-Pe+Sb-aÁ-gÖ-hÚ+jÏ+kÑ-nÇ-oÉ)e0256",
    ///     "+2(+Au-Br+Dp+Gó-Iú-Jí+Lñ+Mç-Oë+Pd-Rb+Ua+gÓ-iÚ-jÍ+lÑ+mÇ-oË)e0257",
    ///     "+2(-At+Bq-Cp+Hó+Iö-Kí-Lï+Mé+Në-Pc+Qb-Ta+hÓ+iÖ-kÍ-lÏ+mÉ+nË)e0267",
    ///     "+2(-Añ-Bú-Cü+Hs-Ir+Ká-Lu+Må-Nä-Ri+Sh-Ul-aÑ-bÚ-cÜ+kÁ+mÅ-nÄ)e0345",
    ///     "+2(+Aï+Bö-Dü-Gs+Iq-Já+Lt+Mæ-Oä+Qi-Sg+Tl+aÏ+bÖ-dÜ-jÁ+mÆ-oÄ)e0346",
    ///     "+2(-Aí-Bó-Eü+Gr-Hq+Ju-Kt+Næ-Oå-Qh+Rg-Tk+Uj-aÍ-bÓ-eÜ+nÆ-oÅ)e0347",
    ///     "+2(-Aë+Cö+Dú+Fs-Ip-Jå-Kæ+Nt+Ou-Pi+Sf+Tn+Uo-aË+cÖ+dÚ-jÅ-kÆ)e0356",
    ///     "+2(+Aé-Có+Eú-Fr+Hp+Jä-Læ-Mt+Oá+Ph-Rf-Tm+aÉ-cÓ+eÚ+jÄ-lÆ+oÁ)e0357",
    ///     "+2(-Aç-Dó-Eö+Fq-Gp+Kä+Lå-Mu-Ná-Pg+Qf-Um-aÇ-dÓ-eÖ+kÄ+lÅ-nÁ)e0367",
    ///     "+2(-Bë-Cï-Dñ+Fá+Gå+Hæ-Lp-Nq-Or-Pl-Qn-Ro-bË-cÏ-dÑ+fÁ+gÅ+hÆ)e0456",
    ///     "+2(+Bé+Cí-Eñ-Fu-Gä+Iæ+Kp+Mq-Os+Pk+Qm-So-Uf+bÉ+cÍ-eÑ-gÄ+iÆ)e0457",
    ///     "+2(-Bç+Dí+Eï+Ft-Hä-Iå-Jp+Mr+Ns-Pj+Rm+Sn+Tf-bÇ+dÍ+eÏ-hÄ-iÅ)e0467",
    ///     "+2(-Cç-Dé-Eë+Gt+Hu+Iá-Jq-Kr-Ls-Qj-Rk-Sl+Tg+Uh-cÇ-dÉ-eË+iÁ)e0567",
    ///     "+2(-jo+kn-lm-tæ+uå-áä-çñ+éï-ëí)e4567",
    ///     "+2(+go-hn+im+qæ-rå+sä-çú+éö-ëó)e3576",
    ///     "+2(-fo+hl-ik-pæ+rá-su-çü+íö-ïó)e3467",
    ///     "+2(+fn-gl+ij+på-qá+st-éü+íú-ñó)e3475",
    ///     "+2(-fm+gk-hj-pä+qu-rt-ëü+ïú-ñö)e3456",
    ///     "+2(-co+dn-em+qñ-rï+sí+tú-uö+áó)e2567",
    ///     "+2(+bo-dl+ek-pñ+rë-sé+tü-äö+åó)e2476",
    ///     "+2(-bn+cl-ej+pï-që+sç+uü-äú+æó)e2457",
    ///     "+2(+bm-ck+dj-pí+qé-rç+áü-åú+æö)e2465",
    ///     "+2(-ao+di-eh-pú-qü+uë-áé+äï-åí)e2367",
    ///     "+2(+an-ci+eg+pö-rü-të+áç+äñ-æí)e2375",
    ///     "+2(-am+ch-dg-pó-sü+té-uç+åñ-æï)e2356",
    ///     "+2(-al+bi-ef+qö+rú-tï-uñ+åç+æé)e2347",
    ///     "+2(+ak-bh+df-qó+sú+tí-áñ-äç+æë)e2364",
    ///     "+2(-aj+bg-cf-ró-sö+uí+áï-äé-åë)e2345",
    ///     "+2(-cæ+då-eä-gñ+hï-ií-jú+kö-ló)e1576",
    ///     "+2(+bæ-dá+eu+fñ-hë+ié-jü+mö-nó)e1467",
    ///     "+2(-bå+cá-et-fï+gë-iç-kü+mú-oó)e1475",
    ///     "+2(+bä-cu+dt+fí-gé+hç-lü+nú-oö)e1456",
    ///     "+2(-aæ+ds-er+fú+gü-kë+lé-mï+ní)e1376",
    ///     "+2(+aå-cs+eq-fö+hü+jë-lç-mñ+oí)e1357",
    ///     "+2(-aä+cr-dq+fó+iü-jé+kç-nñ+oï)e1365",
    ///     "+2(-aá+bs-ep-gö-hú+jï+kñ-nç-oé)e1374",
    ///     "+2(+au-br+dp+gó-iú-jí+lñ+mç-oë)e1346",
    ///     "+2(-at+bq-cp+hó+iö-kí-lï+mé+në)e1354",
    ///     "+2(-añ-bú-cü+hs-ir+ká-lu+må-nä)e1267",
    ///     "+2(+aï+bö-dü-gs+iq-já+lt+mæ-oä)e1275",
    ///     "+2(-aí-bó-eü+gr-hq+ju-kt+næ-oå)e1256",
    ///     "+2(-aë+cö+dú+fs-ip-jå-kæ+nt+ou)e1247",
    ///     "+2(+aé-có+eú-fr+hp+jä-læ-mt+oá)e1264",
    ///     "+2(-aç-dó-eö+fq-gp+kä+lå-mu-ná)e1245",
    ///     "+2(-bë-cï-dñ+fá+gå+hæ-lp-nq-or)e1273",
    ///     "+2(+bé+cí-eñ-fu-gä+iæ+kp+mq-os)e1236",
    ///     "+2(-bç+dí+eï+ft-hä-iå-jp+mr+ns)e1253",
    ///     "+2(+Aa+Bb+Cc+Dd+Ee+Ff+Gg+Hh+Ii+Jj+Kk+Ll+Mm+Nn+Oo+Pp+Qq\
    ///         +Rr+Ss+Tt+Uu+Áá+Ää+Åå+Ææ+Çç+Éé+Ëë+Íí+Ïï+Ññ+Óó+Öö+Úú+Üü)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::volume() + Self::pseudoscalar()
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
            + Self::e13()
            + Self::e12()
            + Self::e14()
            + Self::e24()
            + Self::e34()
            + Self::e15()
            + Self::e25()
            + Self::e35()
            + Self::e45()
            + Self::e16()
            + Self::e26()
            + Self::e36()
            + Self::e46()
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
            + Self::e012375()
            + Self::e012367()
            + Self::e012465()
            + Self::e012457()
            + Self::e012476()
            + Self::e012567()
            + Self::e013456()
            + Self::e013475()
            + Self::e013467()
            + Self::e013576()
            + Self::e014567()
            + Self::e023465()
            + Self::e023457()
            + Self::e023476()
            + Self::e023567()
            + Self::e024576()
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
            + Self::e134576()
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
        Self::e0123465()
            + Self::e0123457()
            + Self::e0123476()
            + Self::e0123567()
            + Self::e0124576()
            + Self::e0134567()
            + Self::e0234576()
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
    /// The multivector of single rotator $`r_1 \equiv s + v^5_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕+œ͔œ͕",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+x͔þ͕-x͕þ͔)e16",
    ///     "+(+x͔œ͕-x͕œ͔)e17",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+y͔þ͕-y͕þ͔)e26",
    ///     "+(+y͔œ͕-y͕œ͔)e27",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+z͔þ͕-z͕þ͔)e36",
    ///     "+(+z͔œ͕-z͕œ͔)e37",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(+ð͔þ͕-ð͕þ͔)e46",
    ///     "+(+ð͔œ͕-ð͕œ͔)e47",
    ///     "+(+ø͔þ͕-ø͕þ͔)e56",
    ///     "+(+ø͔œ͕-ø͕œ͔)e57",
    ///     "+(+þ͔œ͕-þ͕œ͔)e67",
    /// ]);
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// format_eq!(single_rotator, [
    ///     "-x͔x͕-y͔y͕-z͔z͕-ð͔ð͕-ø͔ø͕-þ͔þ͕-œ͔œ͕",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(-x͔z͕+x͕z͔)e13",
    ///     "+(-x͔ð͕+x͕ð͔)e14",
    ///     "+(-x͔ø͕+x͕ø͔)e15",
    ///     "+(-x͔þ͕+x͕þ͔)e16",
    ///     "+(-x͔œ͕+x͕œ͔)e17",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(-y͔ð͕+y͕ð͔)e24",
    ///     "+(-y͔ø͕+y͕ø͔)e25",
    ///     "+(-y͔þ͕+y͕þ͔)e26",
    ///     "+(-y͔œ͕+y͕œ͔)e27",
    ///     "+(-z͔ð͕+z͕ð͔)e34",
    ///     "+(-z͔ø͕+z͕ø͔)e35",
    ///     "+(-z͔þ͕+z͕þ͔)e36",
    ///     "+(-z͔œ͕+z͕œ͔)e37",
    ///     "+(-ð͔ø͕+ð͕ø͔)e45",
    ///     "+(-ð͔þ͕+ð͕þ͔)e46",
    ///     "+(-ð͔œ͕+ð͕œ͔)e47",
    ///     "+(-ø͔þ͕+ø͕þ͔)e56",
    ///     "+(-ø͔œ͕+ø͕œ͔)e57",
    ///     "+(-þ͔œ͕+þ͕œ͔)e67",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::volume5_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + v^5_0 + v_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "+v͔v͕-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕-π͔π͕-ρ͔ρ͕-σ͔σ͕-τ͔τ͕-υ͔υ͕-φ͔φ͕",
    ///     "+(+v͔α͕+v͕α͔-β͔η͕+β͕η͔-γ͔θ͕+γ͕θ͔-δ͔ι͕+δ͕ι͔-ε͔κ͕+ε͕κ͔-ζ͔λ͕+ζ͕λ͔)e12",
    ///     "+(+v͔β͕+v͕β͔+α͔η͕-α͕η͔-γ͔μ͕+γ͕μ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔-ζ͔ο͕+ζ͕ο͔)e13",
    ///     "+(+v͔γ͕+v͕γ͔+α͔θ͕-α͕θ͔+β͔μ͕-β͕μ͔-δ͔π͕+δ͕π͔-ε͔ρ͕+ε͕ρ͔-ζ͔σ͕+ζ͕σ͔)e14",
    ///     "+(+v͔δ͕+v͕δ͔+α͔ι͕-α͕ι͔+β͔ν͕-β͕ν͔+γ͔π͕-γ͕π͔-ε͔τ͕+ε͕τ͔-ζ͔υ͕+ζ͕υ͔)e15",
    ///     "+(+v͔ε͕+v͕ε͔+α͔κ͕-α͕κ͔+β͔ξ͕-β͕ξ͔+γ͔ρ͕-γ͕ρ͔+δ͔τ͕-δ͕τ͔-ζ͔φ͕+ζ͕φ͔)e16",
    ///     "+(+v͔ζ͕+v͕ζ͔+α͔λ͕-α͕λ͔+β͔ο͕-β͕ο͔+γ͔σ͕-γ͕σ͔+δ͔υ͕-δ͕υ͔+ε͔φ͕-ε͕φ͔)e17",
    ///     "+(+v͔η͕+v͕η͔-α͔β͕+α͕β͔-θ͔μ͕+θ͕μ͔-ι͔ν͕+ι͕ν͔-κ͔ξ͕+κ͕ξ͔-λ͔ο͕+λ͕ο͔)e23",
    ///     "+(+v͔θ͕+v͕θ͔-α͔γ͕+α͕γ͔+η͔μ͕-η͕μ͔-ι͔π͕+ι͕π͔-κ͔ρ͕+κ͕ρ͔-λ͔σ͕+λ͕σ͔)e24",
    ///     "+(+v͔ι͕+v͕ι͔-α͔δ͕+α͕δ͔+η͔ν͕-η͕ν͔+θ͔π͕-θ͕π͔-κ͔τ͕+κ͕τ͔-λ͔υ͕+λ͕υ͔)e25",
    ///     "+(+v͔κ͕+v͕κ͔-α͔ε͕+α͕ε͔+η͔ξ͕-η͕ξ͔+θ͔ρ͕-θ͕ρ͔+ι͔τ͕-ι͕τ͔-λ͔φ͕+λ͕φ͔)e26",
    ///     "+(+v͔λ͕+v͕λ͔-α͔ζ͕+α͕ζ͔+η͔ο͕-η͕ο͔+θ͔σ͕-θ͕σ͔+ι͔υ͕-ι͕υ͔+κ͔φ͕-κ͕φ͔)e27",
    ///     "+(+v͔μ͕+v͕μ͔-β͔γ͕+β͕γ͔-η͔θ͕+η͕θ͔-ν͔π͕+ν͕π͔-ξ͔ρ͕+ξ͕ρ͔-ο͔σ͕+ο͕σ͔)e34",
    ///     "+(+v͔ν͕+v͕ν͔-β͔δ͕+β͕δ͔-η͔ι͕+η͕ι͔+μ͔π͕-μ͕π͔-ξ͔τ͕+ξ͕τ͔-ο͔υ͕+ο͕υ͔)e35",
    ///     "+(+v͔ξ͕+v͕ξ͔-β͔ε͕+β͕ε͔-η͔κ͕+η͕κ͔+μ͔ρ͕-μ͕ρ͔+ν͔τ͕-ν͕τ͔-ο͔φ͕+ο͕φ͔)e36",
    ///     "+(+v͔ο͕+v͕ο͔-β͔ζ͕+β͕ζ͔-η͔λ͕+η͕λ͔+μ͔σ͕-μ͕σ͔+ν͔υ͕-ν͕υ͔+ξ͔φ͕-ξ͕φ͔)e37",
    ///     "+(+v͔π͕+v͕π͔-γ͔δ͕+γ͕δ͔-θ͔ι͕+θ͕ι͔-μ͔ν͕+μ͕ν͔-ρ͔τ͕+ρ͕τ͔-σ͔υ͕+σ͕υ͔)e45",
    ///     "+(+v͔ρ͕+v͕ρ͔-γ͔ε͕+γ͕ε͔-θ͔κ͕+θ͕κ͔-μ͔ξ͕+μ͕ξ͔+π͔τ͕-π͕τ͔-σ͔φ͕+σ͕φ͔)e46",
    ///     "+(+v͔σ͕+v͕σ͔-γ͔ζ͕+γ͕ζ͔-θ͔λ͕+θ͕λ͔-μ͔ο͕+μ͕ο͔+π͔υ͕-π͕υ͔+ρ͔φ͕-ρ͕φ͔)e47",
    ///     "+(+v͔τ͕+v͕τ͔-δ͔ε͕+δ͕ε͔-ι͔κ͕+ι͕κ͔-ν͔ξ͕+ν͕ξ͔-π͔ρ͕+π͕ρ͔-υ͔φ͕+υ͕φ͔)e56",
    ///     "+(+v͔υ͕+v͕υ͔-δ͔ζ͕+δ͕ζ͔-ι͔λ͕+ι͕λ͔-ν͔ο͕+ν͕ο͔-π͔σ͕+π͕σ͔+τ͔φ͕-τ͕φ͔)e57",
    ///     "+(+v͔φ͕+v͕φ͔-ε͔ζ͕+ε͕ζ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔-ρ͔σ͕+ρ͕σ͔-τ͔υ͕+τ͕υ͔)e67",
    ///     "+(+α͔μ͕+α͕μ͔-β͔θ͕-β͕θ͔+γ͔η͕+γ͕η͔)e1234",
    ///     "+(+π͔φ͕+π͕φ͔-ρ͔υ͕-ρ͕υ͔+σ͔τ͕+σ͕τ͔)e4567",
    ///     "+(-ν͔φ͕-ν͕φ͔+ξ͔υ͕+ξ͕υ͔-ο͔τ͕-ο͕τ͔)e3576",
    ///     "+(+μ͔φ͕+μ͕φ͔-ξ͔σ͕-ξ͕σ͔+ο͔ρ͕+ο͕ρ͔)e3467",
    ///     "+(-μ͔υ͕-μ͕υ͔+ν͔σ͕+ν͕σ͔-ο͔π͕-ο͕π͔)e3475",
    ///     "+(+μ͔τ͕+μ͕τ͔-ν͔ρ͕-ν͕ρ͔+ξ͔π͕+ξ͕π͔)e3456",
    ///     "+(+ι͔φ͕+ι͕φ͔-κ͔υ͕-κ͕υ͔+λ͔τ͕+λ͕τ͔)e2567",
    ///     "+(-θ͔φ͕-θ͕φ͔+κ͔σ͕+κ͕σ͔-λ͔ρ͕-λ͕ρ͔)e2476",
    ///     "+(+θ͔υ͕+θ͕υ͔-ι͔σ͕-ι͕σ͔+λ͔π͕+λ͕π͔)e2457",
    ///     "+(-θ͔τ͕-θ͕τ͔+ι͔ρ͕+ι͕ρ͔-κ͔π͕-κ͕π͔)e2465",
    ///     "+(+η͔φ͕+η͕φ͔-κ͔ο͕-κ͕ο͔+λ͔ξ͕+λ͕ξ͔)e2367",
    ///     "+(-η͔υ͕-η͕υ͔+ι͔ο͕+ι͕ο͔-λ͔ν͕-λ͕ν͔)e2375",
    ///     "+(+η͔τ͕+η͕τ͔-ι͔ξ͕-ι͕ξ͔+κ͔ν͕+κ͕ν͔)e2356",
    ///     "+(+η͔σ͕+η͕σ͔-θ͔ο͕-θ͕ο͔+λ͔μ͕+λ͕μ͔)e2347",
    ///     "+(-η͔ρ͕-η͕ρ͔+θ͔ξ͕+θ͕ξ͔-κ͔μ͕-κ͕μ͔)e2364",
    ///     "+(+η͔π͕+η͕π͔-θ͔ν͕-θ͕ν͔+ι͔μ͕+ι͕μ͔)e2345",
    ///     "+(-δ͔φ͕-δ͕φ͔+ε͔υ͕+ε͕υ͔-ζ͔τ͕-ζ͕τ͔)e1576",
    ///     "+(+γ͔φ͕+γ͕φ͔-ε͔σ͕-ε͕σ͔+ζ͔ρ͕+ζ͕ρ͔)e1467",
    ///     "+(-γ͔υ͕-γ͕υ͔+δ͔σ͕+δ͕σ͔-ζ͔π͕-ζ͕π͔)e1475",
    ///     "+(+γ͔τ͕+γ͕τ͔-δ͔ρ͕-δ͕ρ͔+ε͔π͕+ε͕π͔)e1456",
    ///     "+(-β͔φ͕-β͕φ͔+ε͔ο͕+ε͕ο͔-ζ͔ξ͕-ζ͕ξ͔)e1376",
    ///     "+(+β͔υ͕+β͕υ͔-δ͔ο͕-δ͕ο͔+ζ͔ν͕+ζ͕ν͔)e1357",
    ///     "+(-β͔τ͕-β͕τ͔+δ͔ξ͕+δ͕ξ͔-ε͔ν͕-ε͕ν͔)e1365",
    ///     "+(-β͔σ͕-β͕σ͔+γ͔ο͕+γ͕ο͔-ζ͔μ͕-ζ͕μ͔)e1374",
    ///     "+(+β͔ρ͕+β͕ρ͔-γ͔ξ͕-γ͕ξ͔+ε͔μ͕+ε͕μ͔)e1346",
    ///     "+(-β͔π͕-β͕π͔+γ͔ν͕+γ͕ν͔-δ͔μ͕-δ͕μ͔)e1354",
    ///     "+(+α͔φ͕+α͕φ͔-ε͔λ͕-ε͕λ͔+ζ͔κ͕+ζ͕κ͔)e1267",
    ///     "+(-α͔υ͕-α͕υ͔+δ͔λ͕+δ͕λ͔-ζ͔ι͕-ζ͕ι͔)e1275",
    ///     "+(+α͔τ͕+α͕τ͔-δ͔κ͕-δ͕κ͔+ε͔ι͕+ε͕ι͔)e1256",
    ///     "+(+α͔σ͕+α͕σ͔-γ͔λ͕-γ͕λ͔+ζ͔θ͕+ζ͕θ͔)e1247",
    ///     "+(-α͔ρ͕-α͕ρ͔+γ͔κ͕+γ͕κ͔-ε͔θ͕-ε͕θ͔)e1264",
    ///     "+(+α͔π͕+α͕π͔-γ͔ι͕-γ͕ι͔+δ͔θ͕+δ͕θ͔)e1245",
    ///     "+(-α͔ο͕-α͕ο͔+β͔λ͕+β͕λ͔-ζ͔η͕-ζ͕η͔)e1273",
    ///     "+(+α͔ξ͕+α͕ξ͔-β͔κ͕-β͕κ͔+ε͔η͕+ε͕η͔)e1236",
    ///     "+(-α͔ν͕-α͕ν͔+β͔ι͕+β͕ι͔-δ͔η͕-δ͕η͔)e1253",
    /// ]);
    ///
    /// let double_rotator = Vee::volume5_displacement().lhs() * Vee::volume5_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕-π͔π͕-ρ͔ρ͕-σ͔σ͕-τ͔τ͕-υ͔υ͕-φ͔φ͕",
    ///     "+(-β͔η͕+β͕η͔-γ͔θ͕+γ͕θ͔-δ͔ι͕+δ͕ι͔-ε͔κ͕+ε͕κ͔-ζ͔λ͕+ζ͕λ͔)e12",
    ///     "+(+α͔η͕-α͕η͔-γ͔μ͕+γ͕μ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔-ζ͔ο͕+ζ͕ο͔)e13",
    ///     "+(+α͔θ͕-α͕θ͔+β͔μ͕-β͕μ͔-δ͔π͕+δ͕π͔-ε͔ρ͕+ε͕ρ͔-ζ͔σ͕+ζ͕σ͔)e14",
    ///     "+(+α͔ι͕-α͕ι͔+β͔ν͕-β͕ν͔+γ͔π͕-γ͕π͔-ε͔τ͕+ε͕τ͔-ζ͔υ͕+ζ͕υ͔)e15",
    ///     "+(+α͔κ͕-α͕κ͔+β͔ξ͕-β͕ξ͔+γ͔ρ͕-γ͕ρ͔+δ͔τ͕-δ͕τ͔-ζ͔φ͕+ζ͕φ͔)e16",
    ///     "+(+α͔λ͕-α͕λ͔+β͔ο͕-β͕ο͔+γ͔σ͕-γ͕σ͔+δ͔υ͕-δ͕υ͔+ε͔φ͕-ε͕φ͔)e17",
    ///     "+(-α͔β͕+α͕β͔-θ͔μ͕+θ͕μ͔-ι͔ν͕+ι͕ν͔-κ͔ξ͕+κ͕ξ͔-λ͔ο͕+λ͕ο͔)e23",
    ///     "+(-α͔γ͕+α͕γ͔+η͔μ͕-η͕μ͔-ι͔π͕+ι͕π͔-κ͔ρ͕+κ͕ρ͔-λ͔σ͕+λ͕σ͔)e24",
    ///     "+(-α͔δ͕+α͕δ͔+η͔ν͕-η͕ν͔+θ͔π͕-θ͕π͔-κ͔τ͕+κ͕τ͔-λ͔υ͕+λ͕υ͔)e25",
    ///     "+(-α͔ε͕+α͕ε͔+η͔ξ͕-η͕ξ͔+θ͔ρ͕-θ͕ρ͔+ι͔τ͕-ι͕τ͔-λ͔φ͕+λ͕φ͔)e26",
    ///     "+(-α͔ζ͕+α͕ζ͔+η͔ο͕-η͕ο͔+θ͔σ͕-θ͕σ͔+ι͔υ͕-ι͕υ͔+κ͔φ͕-κ͕φ͔)e27",
    ///     "+(-β͔γ͕+β͕γ͔-η͔θ͕+η͕θ͔-ν͔π͕+ν͕π͔-ξ͔ρ͕+ξ͕ρ͔-ο͔σ͕+ο͕σ͔)e34",
    ///     "+(-β͔δ͕+β͕δ͔-η͔ι͕+η͕ι͔+μ͔π͕-μ͕π͔-ξ͔τ͕+ξ͕τ͔-ο͔υ͕+ο͕υ͔)e35",
    ///     "+(-β͔ε͕+β͕ε͔-η͔κ͕+η͕κ͔+μ͔ρ͕-μ͕ρ͔+ν͔τ͕-ν͕τ͔-ο͔φ͕+ο͕φ͔)e36",
    ///     "+(-β͔ζ͕+β͕ζ͔-η͔λ͕+η͕λ͔+μ͔σ͕-μ͕σ͔+ν͔υ͕-ν͕υ͔+ξ͔φ͕-ξ͕φ͔)e37",
    ///     "+(-γ͔δ͕+γ͕δ͔-θ͔ι͕+θ͕ι͔-μ͔ν͕+μ͕ν͔-ρ͔τ͕+ρ͕τ͔-σ͔υ͕+σ͕υ͔)e45",
    ///     "+(-γ͔ε͕+γ͕ε͔-θ͔κ͕+θ͕κ͔-μ͔ξ͕+μ͕ξ͔+π͔τ͕-π͕τ͔-σ͔φ͕+σ͕φ͔)e46",
    ///     "+(-γ͔ζ͕+γ͕ζ͔-θ͔λ͕+θ͕λ͔-μ͔ο͕+μ͕ο͔+π͔υ͕-π͕υ͔+ρ͔φ͕-ρ͕φ͔)e47",
    ///     "+(-δ͔ε͕+δ͕ε͔-ι͔κ͕+ι͕κ͔-ν͔ξ͕+ν͕ξ͔-π͔ρ͕+π͕ρ͔-υ͔φ͕+υ͕φ͔)e56",
    ///     "+(-δ͔ζ͕+δ͕ζ͔-ι͔λ͕+ι͕λ͔-ν͔ο͕+ν͕ο͔-π͔σ͕+π͕σ͔+τ͔φ͕-τ͕φ͔)e57",
    ///     "+(-ε͔ζ͕+ε͕ζ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔-ρ͔σ͕+ρ͕σ͔-τ͔υ͕+τ͕υ͔)e67",
    ///     "+(+α͔μ͕+α͕μ͔-β͔θ͕-β͕θ͔+γ͔η͕+γ͕η͔)e1234",
    ///     "+(+π͔φ͕+π͕φ͔-ρ͔υ͕-ρ͕υ͔+σ͔τ͕+σ͕τ͔)e4567",
    ///     "+(-ν͔φ͕-ν͕φ͔+ξ͔υ͕+ξ͕υ͔-ο͔τ͕-ο͕τ͔)e3576",
    ///     "+(+μ͔φ͕+μ͕φ͔-ξ͔σ͕-ξ͕σ͔+ο͔ρ͕+ο͕ρ͔)e3467",
    ///     "+(-μ͔υ͕-μ͕υ͔+ν͔σ͕+ν͕σ͔-ο͔π͕-ο͕π͔)e3475",
    ///     "+(+μ͔τ͕+μ͕τ͔-ν͔ρ͕-ν͕ρ͔+ξ͔π͕+ξ͕π͔)e3456",
    ///     "+(+ι͔φ͕+ι͕φ͔-κ͔υ͕-κ͕υ͔+λ͔τ͕+λ͕τ͔)e2567",
    ///     "+(-θ͔φ͕-θ͕φ͔+κ͔σ͕+κ͕σ͔-λ͔ρ͕-λ͕ρ͔)e2476",
    ///     "+(+θ͔υ͕+θ͕υ͔-ι͔σ͕-ι͕σ͔+λ͔π͕+λ͕π͔)e2457",
    ///     "+(-θ͔τ͕-θ͕τ͔+ι͔ρ͕+ι͕ρ͔-κ͔π͕-κ͕π͔)e2465",
    ///     "+(+η͔φ͕+η͕φ͔-κ͔ο͕-κ͕ο͔+λ͔ξ͕+λ͕ξ͔)e2367",
    ///     "+(-η͔υ͕-η͕υ͔+ι͔ο͕+ι͕ο͔-λ͔ν͕-λ͕ν͔)e2375",
    ///     "+(+η͔τ͕+η͕τ͔-ι͔ξ͕-ι͕ξ͔+κ͔ν͕+κ͕ν͔)e2356",
    ///     "+(+η͔σ͕+η͕σ͔-θ͔ο͕-θ͕ο͔+λ͔μ͕+λ͕μ͔)e2347",
    ///     "+(-η͔ρ͕-η͕ρ͔+θ͔ξ͕+θ͕ξ͔-κ͔μ͕-κ͕μ͔)e2364",
    ///     "+(+η͔π͕+η͕π͔-θ͔ν͕-θ͕ν͔+ι͔μ͕+ι͕μ͔)e2345",
    ///     "+(-δ͔φ͕-δ͕φ͔+ε͔υ͕+ε͕υ͔-ζ͔τ͕-ζ͕τ͔)e1576",
    ///     "+(+γ͔φ͕+γ͕φ͔-ε͔σ͕-ε͕σ͔+ζ͔ρ͕+ζ͕ρ͔)e1467",
    ///     "+(-γ͔υ͕-γ͕υ͔+δ͔σ͕+δ͕σ͔-ζ͔π͕-ζ͕π͔)e1475",
    ///     "+(+γ͔τ͕+γ͕τ͔-δ͔ρ͕-δ͕ρ͔+ε͔π͕+ε͕π͔)e1456",
    ///     "+(-β͔φ͕-β͕φ͔+ε͔ο͕+ε͕ο͔-ζ͔ξ͕-ζ͕ξ͔)e1376",
    ///     "+(+β͔υ͕+β͕υ͔-δ͔ο͕-δ͕ο͔+ζ͔ν͕+ζ͕ν͔)e1357",
    ///     "+(-β͔τ͕-β͕τ͔+δ͔ξ͕+δ͕ξ͔-ε͔ν͕-ε͕ν͔)e1365",
    ///     "+(-β͔σ͕-β͕σ͔+γ͔ο͕+γ͕ο͔-ζ͔μ͕-ζ͕μ͔)e1374",
    ///     "+(+β͔ρ͕+β͕ρ͔-γ͔ξ͕-γ͕ξ͔+ε͔μ͕+ε͕μ͔)e1346",
    ///     "+(-β͔π͕-β͕π͔+γ͔ν͕+γ͕ν͔-δ͔μ͕-δ͕μ͔)e1354",
    ///     "+(+α͔φ͕+α͕φ͔-ε͔λ͕-ε͕λ͔+ζ͔κ͕+ζ͕κ͔)e1267",
    ///     "+(-α͔υ͕-α͕υ͔+δ͔λ͕+δ͕λ͔-ζ͔ι͕-ζ͕ι͔)e1275",
    ///     "+(+α͔τ͕+α͕τ͔-δ͔κ͕-δ͕κ͔+ε͔ι͕+ε͕ι͔)e1256",
    ///     "+(+α͔σ͕+α͕σ͔-γ͔λ͕-γ͕λ͔+ζ͔θ͕+ζ͕θ͔)e1247",
    ///     "+(-α͔ρ͕-α͕ρ͔+γ͔κ͕+γ͕κ͔-ε͔θ͕-ε͕θ͔)e1264",
    ///     "+(+α͔π͕+α͕π͔-γ͔ι͕-γ͕ι͔+δ͔θ͕+δ͕θ͔)e1245",
    ///     "+(-α͔ο͕-α͕ο͔+β͔λ͕+β͕λ͔-ζ͔η͕-ζ͕η͔)e1273",
    ///     "+(+α͔ξ͕+α͕ξ͔-β͔κ͕-β͕κ͔+ε͔η͕+ε͕η͔)e1236",
    ///     "+(-α͔ν͕-α͕ν͔+β͔ι͕+β͕ι͔-δ͔η͕-δ͕η͔)e1253",
    /// ]);
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// format_eq!(double_rotator, [
    ///     "+α͔α͕+β͔β͕+γ͔γ͕+δ͔δ͕+ε͔ε͕+ζ͔ζ͕+η͔η͕+θ͔θ͕+ι͔ι͕+κ͔κ͕+λ͔λ͕+μ͔μ͕+ν͔ν͕+ξ͔ξ͕+ο͔ο͕+π͔π͕+ρ͔ρ͕+σ͔σ͕+τ͔τ͕+υ͔υ͕+φ͔φ͕",
    ///     "+(+β͔η͕-β͕η͔+γ͔θ͕-γ͕θ͔+δ͔ι͕-δ͕ι͔+ε͔κ͕-ε͕κ͔+ζ͔λ͕-ζ͕λ͔)e12",
    ///     "+(-α͔η͕+α͕η͔+γ͔μ͕-γ͕μ͔+δ͔ν͕-δ͕ν͔+ε͔ξ͕-ε͕ξ͔+ζ͔ο͕-ζ͕ο͔)e13",
    ///     "+(-α͔θ͕+α͕θ͔-β͔μ͕+β͕μ͔+δ͔π͕-δ͕π͔+ε͔ρ͕-ε͕ρ͔+ζ͔σ͕-ζ͕σ͔)e14",
    ///     "+(-α͔ι͕+α͕ι͔-β͔ν͕+β͕ν͔-γ͔π͕+γ͕π͔+ε͔τ͕-ε͕τ͔+ζ͔υ͕-ζ͕υ͔)e15",
    ///     "+(-α͔κ͕+α͕κ͔-β͔ξ͕+β͕ξ͔-γ͔ρ͕+γ͕ρ͔-δ͔τ͕+δ͕τ͔+ζ͔φ͕-ζ͕φ͔)e16",
    ///     "+(-α͔λ͕+α͕λ͔-β͔ο͕+β͕ο͔-γ͔σ͕+γ͕σ͔-δ͔υ͕+δ͕υ͔-ε͔φ͕+ε͕φ͔)e17",
    ///     "+(+α͔β͕-α͕β͔+θ͔μ͕-θ͕μ͔+ι͔ν͕-ι͕ν͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e23",
    ///     "+(+α͔γ͕-α͕γ͔-η͔μ͕+η͕μ͔+ι͔π͕-ι͕π͔+κ͔ρ͕-κ͕ρ͔+λ͔σ͕-λ͕σ͔)e24",
    ///     "+(+α͔δ͕-α͕δ͔-η͔ν͕+η͕ν͔-θ͔π͕+θ͕π͔+κ͔τ͕-κ͕τ͔+λ͔υ͕-λ͕υ͔)e25",
    ///     "+(+α͔ε͕-α͕ε͔-η͔ξ͕+η͕ξ͔-θ͔ρ͕+θ͕ρ͔-ι͔τ͕+ι͕τ͔+λ͔φ͕-λ͕φ͔)e26",
    ///     "+(+α͔ζ͕-α͕ζ͔-η͔ο͕+η͕ο͔-θ͔σ͕+θ͕σ͔-ι͔υ͕+ι͕υ͔-κ͔φ͕+κ͕φ͔)e27",
    ///     "+(+β͔γ͕-β͕γ͔+η͔θ͕-η͕θ͔+ν͔π͕-ν͕π͔+ξ͔ρ͕-ξ͕ρ͔+ο͔σ͕-ο͕σ͔)e34",
    ///     "+(+β͔δ͕-β͕δ͔+η͔ι͕-η͕ι͔-μ͔π͕+μ͕π͔+ξ͔τ͕-ξ͕τ͔+ο͔υ͕-ο͕υ͔)e35",
    ///     "+(+β͔ε͕-β͕ε͔+η͔κ͕-η͕κ͔-μ͔ρ͕+μ͕ρ͔-ν͔τ͕+ν͕τ͔+ο͔φ͕-ο͕φ͔)e36",
    ///     "+(+β͔ζ͕-β͕ζ͔+η͔λ͕-η͕λ͔-μ͔σ͕+μ͕σ͔-ν͔υ͕+ν͕υ͔-ξ͔φ͕+ξ͕φ͔)e37",
    ///     "+(+γ͔δ͕-γ͕δ͔+θ͔ι͕-θ͕ι͔+μ͔ν͕-μ͕ν͔+ρ͔τ͕-ρ͕τ͔+σ͔υ͕-σ͕υ͔)e45",
    ///     "+(+γ͔ε͕-γ͕ε͔+θ͔κ͕-θ͕κ͔+μ͔ξ͕-μ͕ξ͔-π͔τ͕+π͕τ͔+σ͔φ͕-σ͕φ͔)e46",
    ///     "+(+γ͔ζ͕-γ͕ζ͔+θ͔λ͕-θ͕λ͔+μ͔ο͕-μ͕ο͔-π͔υ͕+π͕υ͔-ρ͔φ͕+ρ͕φ͔)e47",
    ///     "+(+δ͔ε͕-δ͕ε͔+ι͔κ͕-ι͕κ͔+ν͔ξ͕-ν͕ξ͔+π͔ρ͕-π͕ρ͔+υ͔φ͕-υ͕φ͔)e56",
    ///     "+(+δ͔ζ͕-δ͕ζ͔+ι͔λ͕-ι͕λ͔+ν͔ο͕-ν͕ο͔+π͔σ͕-π͕σ͔-τ͔φ͕+τ͕φ͔)e57",
    ///     "+(+ε͔ζ͕-ε͕ζ͔+κ͔λ͕-κ͕λ͔+ξ͔ο͕-ξ͕ο͔+ρ͔σ͕-ρ͕σ͔+τ͔υ͕-τ͕υ͔)e67",
    ///     "+(-α͔μ͕-α͕μ͔+β͔θ͕+β͕θ͔-γ͔η͕-γ͕η͔)e1234",
    ///     "+(-π͔φ͕-π͕φ͔+ρ͔υ͕+ρ͕υ͔-σ͔τ͕-σ͕τ͔)e4567",
    ///     "+(+ν͔φ͕+ν͕φ͔-ξ͔υ͕-ξ͕υ͔+ο͔τ͕+ο͕τ͔)e3576",
    ///     "+(-μ͔φ͕-μ͕φ͔+ξ͔σ͕+ξ͕σ͔-ο͔ρ͕-ο͕ρ͔)e3467",
    ///     "+(+μ͔υ͕+μ͕υ͔-ν͔σ͕-ν͕σ͔+ο͔π͕+ο͕π͔)e3475",
    ///     "+(-μ͔τ͕-μ͕τ͔+ν͔ρ͕+ν͕ρ͔-ξ͔π͕-ξ͕π͔)e3456",
    ///     "+(-ι͔φ͕-ι͕φ͔+κ͔υ͕+κ͕υ͔-λ͔τ͕-λ͕τ͔)e2567",
    ///     "+(+θ͔φ͕+θ͕φ͔-κ͔σ͕-κ͕σ͔+λ͔ρ͕+λ͕ρ͔)e2476",
    ///     "+(-θ͔υ͕-θ͕υ͔+ι͔σ͕+ι͕σ͔-λ͔π͕-λ͕π͔)e2457",
    ///     "+(+θ͔τ͕+θ͕τ͔-ι͔ρ͕-ι͕ρ͔+κ͔π͕+κ͕π͔)e2465",
    ///     "+(-η͔φ͕-η͕φ͔+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔)e2367",
    ///     "+(+η͔υ͕+η͕υ͔-ι͔ο͕-ι͕ο͔+λ͔ν͕+λ͕ν͔)e2375",
    ///     "+(-η͔τ͕-η͕τ͔+ι͔ξ͕+ι͕ξ͔-κ͔ν͕-κ͕ν͔)e2356",
    ///     "+(-η͔σ͕-η͕σ͔+θ͔ο͕+θ͕ο͔-λ͔μ͕-λ͕μ͔)e2347",
    ///     "+(+η͔ρ͕+η͕ρ͔-θ͔ξ͕-θ͕ξ͔+κ͔μ͕+κ͕μ͔)e2364",
    ///     "+(-η͔π͕-η͕π͔+θ͔ν͕+θ͕ν͔-ι͔μ͕-ι͕μ͔)e2345",
    ///     "+(+δ͔φ͕+δ͕φ͔-ε͔υ͕-ε͕υ͔+ζ͔τ͕+ζ͕τ͔)e1576",
    ///     "+(-γ͔φ͕-γ͕φ͔+ε͔σ͕+ε͕σ͔-ζ͔ρ͕-ζ͕ρ͔)e1467",
    ///     "+(+γ͔υ͕+γ͕υ͔-δ͔σ͕-δ͕σ͔+ζ͔π͕+ζ͕π͔)e1475",
    ///     "+(-γ͔τ͕-γ͕τ͔+δ͔ρ͕+δ͕ρ͔-ε͔π͕-ε͕π͔)e1456",
    ///     "+(+β͔φ͕+β͕φ͔-ε͔ο͕-ε͕ο͔+ζ͔ξ͕+ζ͕ξ͔)e1376",
    ///     "+(-β͔υ͕-β͕υ͔+δ͔ο͕+δ͕ο͔-ζ͔ν͕-ζ͕ν͔)e1357",
    ///     "+(+β͔τ͕+β͕τ͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔)e1365",
    ///     "+(+β͔σ͕+β͕σ͔-γ͔ο͕-γ͕ο͔+ζ͔μ͕+ζ͕μ͔)e1374",
    ///     "+(-β͔ρ͕-β͕ρ͔+γ͔ξ͕+γ͕ξ͔-ε͔μ͕-ε͕μ͔)e1346",
    ///     "+(+β͔π͕+β͕π͔-γ͔ν͕-γ͕ν͔+δ͔μ͕+δ͕μ͔)e1354",
    ///     "+(-α͔φ͕-α͕φ͔+ε͔λ͕+ε͕λ͔-ζ͔κ͕-ζ͕κ͔)e1267",
    ///     "+(+α͔υ͕+α͕υ͔-δ͔λ͕-δ͕λ͔+ζ͔ι͕+ζ͕ι͔)e1275",
    ///     "+(-α͔τ͕-α͕τ͔+δ͔κ͕+δ͕κ͔-ε͔ι͕-ε͕ι͔)e1256",
    ///     "+(-α͔σ͕-α͕σ͔+γ͔λ͕+γ͕λ͔-ζ͔θ͕-ζ͕θ͔)e1247",
    ///     "+(+α͔ρ͕+α͕ρ͔-γ͔κ͕-γ͕κ͔+ε͔θ͕+ε͕θ͔)e1264",
    ///     "+(-α͔π͕-α͕π͔+γ͔ι͕+γ͕ι͔-δ͔θ͕-δ͕θ͔)e1245",
    ///     "+(+α͔ο͕+α͕ο͔-β͔λ͕-β͕λ͔+ζ͔η͕+ζ͕η͔)e1273",
    ///     "+(-α͔ξ͕-α͕ξ͔+β͔κ͕+β͕κ͔-ε͔η͕-ε͕η͔)e1236",
    ///     "+(+α͔ν͕+α͕ν͔-β͔ι͕-β͕ι͔+δ͔η͕+δ͕η͔)e1253",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::volume5_displacement() + Self::volume_displacement()
    }
    /// The multivector of triple rotator $`r_3 \equiv s + v^5_0 + v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let triple_rotator = Vee::single_rotator().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(triple_rotator.basis_blades(), Vee::triple_rotator().basis_blades());
    /// format_eq!(triple_rotator, [
    ///     "+v͔v͕-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕-π͔π͕-ρ͔ρ͕-σ͔σ͕-τ͔τ͕-υ͔υ͕-φ͔φ͕",
    ///     "+(+v͔α͕+v͕α͔-ç͕φ͔+é͕υ͔-ë͕τ͔-í͕σ͔+ï͕ρ͔-ñ͕π͔+ó͕ο͔-ö͕ξ͔+ú͕ν͔-ü͕μ͔-β͔η͕+β͕η͔-γ͔θ͕+γ͕θ͔-δ͔ι͕+δ͕ι͔-ε͔κ͕+ε͕κ͔-ζ͔λ͕+ζ͕λ͔)e12",
    ///     "+(+t͕φ͔-u͕υ͔+v͔β͕+v͕β͔+á͕τ͔+ä͕σ͔-å͕ρ͔+æ͕π͔-ó͕λ͔+ö͕κ͔-ú͕ι͔+ü͕θ͔+α͔η͕-α͕η͔-γ͔μ͕+γ͕μ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔-ζ͔ο͕+ζ͕ο͔)e13",
    ///     "+(-q͕φ͔+r͕υ͔-s͕τ͔+v͔γ͕+v͕γ͔-ä͕ο͔+å͕ξ͔-æ͕ν͔+í͕λ͔-ï͕κ͔+ñ͕ι͔-ü͕η͔+α͔θ͕-α͕θ͔+β͔μ͕-β͕μ͔-δ͔π͕+δ͕π͔-ε͔ρ͕+ε͕ρ͔-ζ͔σ͕+ζ͕σ͔)e14",
    ///     "+(+p͕φ͔-r͕σ͔+s͕ρ͔+u͕ο͔+v͔δ͕+v͕δ͔-á͕ξ͔+æ͕μ͔-é͕λ͔+ë͕κ͔-ñ͕θ͔+ú͕η͔+α͔ι͕-α͕ι͔+β͔ν͕-β͕ν͔+γ͔π͕-γ͕π͔-ε͔τ͕+ε͕τ͔-ζ͔υ͕+ζ͕υ͔)e15",
    ///     "+(-p͕υ͔+q͕σ͔-s͕π͔-t͕ο͔+v͔ε͕+v͕ε͔+á͕ν͔-å͕μ͔+ç͕λ͔-ë͕ι͔+ï͕θ͔-ö͕η͔+α͔κ͕-α͕κ͔+β͔ξ͕-β͕ξ͔+γ͔ρ͕-γ͕ρ͔+δ͔τ͕-δ͕τ͔-ζ͔φ͕+ζ͕φ͔)e16",
    ///     "+(+p͕τ͔-q͕ρ͔+r͕π͔+t͕ξ͔-u͕ν͔+v͔ζ͕+v͕ζ͔+ä͕μ͔-ç͕κ͔+é͕ι͔-í͕θ͔+ó͕η͔+α͔λ͕-α͕λ͔+β͔ο͕-β͕ο͔+γ͔σ͕-γ͕σ͔+δ͔υ͕-δ͕υ͔+ε͔φ͕-ε͕φ͔)e17",
    ///     "+(-j͕φ͔+k͕υ͔-l͕τ͔-m͕σ͔+n͕ρ͔-o͕π͔+v͔η͕+v͕η͔+ó͕ζ͔-ö͕ε͔+ú͕δ͔-ü͕γ͔-α͔β͕+α͕β͔-θ͔μ͕+θ͕μ͔-ι͔ν͕+ι͕ν͔-κ͔ξ͕+κ͕ξ͔-λ͔ο͕+λ͕ο͔)e23",
    ///     "+(+g͕φ͔-h͕υ͔+i͕τ͔+m͕ο͔-n͕ξ͔+o͕ν͔+v͔θ͕+v͕θ͔-í͕ζ͔+ï͕ε͔-ñ͕δ͔+ü͕β͔-α͔γ͕+α͕γ͔+η͔μ͕-η͕μ͔-ι͔π͕+ι͕π͔-κ͔ρ͕+κ͕ρ͔-λ͔σ͕+λ͕σ͔)e24",
    ///     "+(-f͕φ͔+h͕σ͔-i͕ρ͔-k͕ο͔+l͕ξ͔-o͕μ͔+v͔ι͕+v͕ι͔+é͕ζ͔-ë͕ε͔+ñ͕γ͔-ú͕β͔-α͔δ͕+α͕δ͔+η͔ν͕-η͕ν͔+θ͔π͕-θ͕π͔-κ͔τ͕+κ͕τ͔-λ͔υ͕+λ͕υ͔)e25",
    ///     "+(+f͕υ͔-g͕σ͔+i͕π͔+j͕ο͔-l͕ν͔+n͕μ͔+v͔κ͕+v͕κ͔-ç͕ζ͔+ë͕δ͔-ï͕γ͔+ö͕β͔-α͔ε͕+α͕ε͔+η͔ξ͕-η͕ξ͔+θ͔ρ͕-θ͕ρ͔+ι͔τ͕-ι͕τ͔-λ͔φ͕+λ͕φ͔)e26",
    ///     "+(-f͕τ͔+g͕ρ͔-h͕π͔-j͕ξ͔+k͕ν͔-m͕μ͔+v͔λ͕+v͕λ͔+ç͕ε͔-é͕δ͔+í͕γ͔-ó͕β͔-α͔ζ͕+α͕ζ͔+η͔ο͕-η͕ο͔+θ͔σ͕-θ͕σ͔+ι͔υ͕-ι͕υ͔+κ͔φ͕-κ͕φ͔)e27",
    ///     "+(-c͕φ͔+d͕υ͔-e͕τ͔-m͕λ͔+n͕κ͔-o͕ι͔+v͔μ͕+v͕μ͔+ä͕ζ͔-å͕ε͔+æ͕δ͔-ü͕α͔-β͔γ͕+β͕γ͔-η͔θ͕+η͕θ͔-ν͔π͕+ν͕π͔-ξ͔ρ͕+ξ͕ρ͔-ο͔σ͕+ο͕σ͔)e34",
    ///     "+(+b͕φ͔-d͕σ͔+e͕ρ͔+k͕λ͔-l͕κ͔+o͕θ͔-u͕ζ͔+v͔ν͕+v͕ν͔+á͕ε͔-æ͕γ͔+ú͕α͔-β͔δ͕+β͕δ͔-η͔ι͕+η͕ι͔+μ͔π͕-μ͕π͔-ξ͔τ͕+ξ͕τ͔-ο͔υ͕+ο͕υ͔)e35",
    ///     "+(-b͕υ͔+c͕σ͔-e͕π͔-j͕λ͔+l͕ι͔-n͕θ͔+t͕ζ͔+v͔ξ͕+v͕ξ͔-á͕δ͔+å͕γ͔-ö͕α͔-β͔ε͕+β͕ε͔-η͔κ͕+η͕κ͔+μ͔ρ͕-μ͕ρ͔+ν͔τ͕-ν͕τ͔-ο͔φ͕+ο͕φ͔)e36",
    ///     "+(+b͕τ͔-c͕ρ͔+d͕π͔+j͕κ͔-k͕ι͔+m͕θ͔-t͕ε͔+u͕δ͔+v͔ο͕+v͕ο͔-ä͕γ͔+ó͕α͔-β͔ζ͕+β͕ζ͔-η͔λ͕+η͕λ͔+μ͔σ͕-μ͕σ͔+ν͔υ͕-ν͕υ͔+ξ͔φ͕-ξ͕φ͔)e37",
    ///     "+(-a͕φ͔+d͕ο͔-e͕ξ͔-h͕λ͔+i͕κ͔-o͕η͔+r͕ζ͔-s͕ε͔+v͔π͕+v͕π͔+æ͕β͔-ñ͕α͔-γ͔δ͕+γ͕δ͔-θ͔ι͕+θ͕ι͔-μ͔ν͕+μ͕ν͔-ρ͔τ͕+ρ͕τ͔-σ͔υ͕+σ͕υ͔)e45",
    ///     "+(+a͕υ͔-c͕ο͔+e͕ν͔+g͕λ͔-i͕ι͔+n͕η͔-q͕ζ͔+s͕δ͔+v͔ρ͕+v͕ρ͔-å͕β͔+ï͕α͔-γ͔ε͕+γ͕ε͔-θ͔κ͕+θ͕κ͔-μ͔ξ͕+μ͕ξ͔+π͔τ͕-π͕τ͔-σ͔φ͕+σ͕φ͔)e46",
    ///     "+(-a͕τ͔+c͕ξ͔-d͕ν͔-g͕κ͔+h͕ι͔-m͕η͔+q͕ε͔-r͕δ͔+v͔σ͕+v͕σ͔+ä͕β͔-í͕α͔-γ͔ζ͕+γ͕ζ͔-θ͔λ͕+θ͕λ͔-μ͔ο͕+μ͕ο͔+π͔υ͕-π͕υ͔+ρ͔φ͕-ρ͕φ͔)e47",
    ///     "+(-a͕σ͔+b͕ο͔-e͕μ͔-f͕λ͔+i͕θ͔-l͕η͔+p͕ζ͔-s͕γ͔+v͔τ͕+v͕τ͔+á͕β͔-ë͕α͔-δ͔ε͕+δ͕ε͔-ι͔κ͕+ι͕κ͔-ν͔ξ͕+ν͕ξ͔-π͔ρ͕+π͕ρ͔-υ͔φ͕+υ͕φ͔)e56",
    ///     "+(+a͕ρ͔-b͕ξ͔+d͕μ͔+f͕κ͔-h͕θ͔+k͕η͔-p͕ε͔+r͕γ͔-u͕β͔+v͔υ͕+v͕υ͔+é͕α͔-δ͔ζ͕+δ͕ζ͔-ι͔λ͕+ι͕λ͔-ν͔ο͕+ν͕ο͔-π͔σ͕+π͕σ͔+τ͔φ͕-τ͕φ͔)e57",
    ///     "+(-a͕π͔+b͕ν͔-c͕μ͔-f͕ι͔+g͕θ͔-j͕η͔+p͕δ͔-q͕γ͔+t͕β͔+v͔φ͕+v͕φ͔-ç͕α͔-ε͔ζ͕+ε͕ζ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔-ρ͔σ͕+ρ͕σ͔-τ͔υ͕+τ͕υ͔)e67",
    ///     "+(-m͕ζ͔+n͕ε͔-o͕δ͔+v͔ü͕-ä͕λ͔+å͕κ͔-æ͕ι͔-í͕ο͔+ï͕ξ͔-ñ͕ν͔-ó͕σ͔+ö͕ρ͔-ú͕π͔+α͔μ͕+α͕μ͔-β͔θ͕-β͕θ͔+γ͔η͕+γ͕η͔)e1234",
    ///     "+(+a͕v͔+b͕μ͔+c͕ν͔+d͕ξ͔+e͕ο͔-f͕θ͔-g͕ι͔-h͕κ͔-i͕λ͔+p͕γ͔+q͕δ͔+r͕ε͔+s͕ζ͔+π͔φ͕+π͕φ͔-ρ͔υ͕-ρ͕υ͔+σ͔τ͕+σ͕τ͔)e4567",
    ///     "+(-a͕μ͔+b͕v͔+c͕π͔+d͕ρ͔+e͕σ͔+f͕η͔-j͕ι͔-k͕κ͔-l͕λ͔-p͕β͔+t͕δ͔+u͕ε͔+á͕ζ͔-ν͔φ͕-ν͕φ͔+ξ͔υ͕+ξ͕υ͔-ο͔τ͕-ο͕τ͔)e3576",
    ///     "+(-a͕ν͔-b͕π͔+c͕v͔+d͕τ͔+e͕υ͔+g͕η͔+j͕θ͔-m͕κ͔-n͕λ͔-q͕β͔-t͕γ͔+ä͕ε͔+å͕ζ͔+μ͔φ͕+μ͕φ͔-ξ͔σ͕-ξ͕σ͔+ο͔ρ͕+ο͕ρ͔)e3467",
    ///     "+(-a͕ξ͔-b͕ρ͔-c͕τ͔+d͕v͔+e͕φ͔+h͕η͔+k͕θ͔+m͕ι͔-o͕λ͔-r͕β͔-u͕γ͔-ä͕δ͔+æ͕ζ͔-μ͔υ͕-μ͕υ͔+ν͔σ͕+ν͕σ͔-ο͔π͕-ο͕π͔)e3475",
    ///     "+(-a͕ο͔-b͕σ͔-c͕υ͔-d͕φ͔+e͕v͔+i͕η͔+l͕θ͔+n͕ι͔+o͕κ͔-s͕β͔-á͕γ͔-å͕δ͔-æ͕ε͔+μ͔τ͕+μ͕τ͔-ν͔ρ͕-ν͕ρ͔+ξ͔π͕+ξ͕π͔)e3456",
    ///     "+(+a͕θ͔-b͕η͔+f͕v͔+g͕π͔+h͕ρ͔+i͕σ͔-j͕ν͔-k͕ξ͔-l͕ο͔+p͕α͔+ç͕δ͔+é͕ε͔+ë͕ζ͔+ι͔φ͕+ι͕φ͔-κ͔υ͕-κ͕υ͔+λ͔τ͕+λ͕τ͔)e2567",
    ///     "+(+a͕ι͔-c͕η͔-f͕π͔+g͕v͔+h͕τ͔+i͕υ͔+j͕μ͔-m͕ξ͔-n͕ο͔+q͕α͔-ç͕γ͔+í͕ε͔+ï͕ζ͔-θ͔φ͕-θ͕φ͔+κ͔σ͕+κ͕σ͔-λ͔ρ͕-λ͕ρ͔)e2476",
    ///     "+(+a͕κ͔-d͕η͔-f͕ρ͔-g͕τ͔+h͕v͔+i͕φ͔+k͕μ͔+m͕ν͔-o͕ο͔+r͕α͔-é͕γ͔-í͕δ͔+ñ͕ζ͔+θ͔υ͕+θ͕υ͔-ι͔σ͕-ι͕σ͔+λ͔π͕+λ͕π͔)e2457",
    ///     "+(+a͕λ͔-e͕η͔-f͕σ͔-g͕υ͔-h͕φ͔+i͕v͔+l͕μ͔+n͕ν͔+o͕ξ͔+s͕α͔-ë͕γ͔-ï͕δ͔-ñ͕ε͔-θ͔τ͕-θ͕τ͔+ι͔ρ͕+ι͕ρ͔-κ͔π͕-κ͕π͔)e2465",
    ///     "+(+b͕ι͔-c͕θ͔+f͕ν͔-g͕μ͔+j͕v͔+k͕τ͔+l͕υ͔-m͕ρ͔-n͕σ͔+t͕α͔+ç͕β͔+ó͕ε͔+ö͕ζ͔+η͔φ͕+η͕φ͔-κ͔ο͕-κ͕ο͔+λ͔ξ͕+λ͕ξ͔)e2367",
    ///     "+(+b͕κ͔-d͕θ͔+f͕ξ͔-h͕μ͔-j͕τ͔+k͕v͔+l͕φ͔+m͕π͔-o͕σ͔+u͕α͔+é͕β͔-ó͕δ͔+ú͕ζ͔-η͔υ͕-η͕υ͔+ι͔ο͕+ι͕ο͔-λ͔ν͕-λ͕ν͔)e2375",
    ///     "+(+b͕λ͔-e͕θ͔+f͕ο͔-i͕μ͔-j͕υ͔-k͕φ͔+l͕v͔+n͕π͔+o͕ρ͔+á͕α͔+ë͕β͔-ö͕δ͔-ú͕ε͔+η͔τ͕+η͕τ͔-ι͔ξ͕-ι͕ξ͔+κ͔ν͕+κ͕ν͔)e2356",
    ///     "+(+c͕κ͔-d͕ι͔+g͕ξ͔-h͕ν͔+j͕ρ͔-k͕π͔+m͕v͔+n͕φ͔-o͕υ͔+ä͕α͔+í͕β͔+ó͕γ͔+ü͕ζ͔+η͔σ͕+η͕σ͔-θ͔ο͕-θ͕ο͔+λ͔μ͕+λ͕μ͔)e2347",
    ///     "+(+c͕λ͔-e͕ι͔+g͕ο͔-i͕ν͔+j͕σ͔-l͕π͔-m͕φ͔+n͕v͔+o͕τ͔+å͕α͔+ï͕β͔+ö͕γ͔-ü͕ε͔-η͔ρ͕-η͕ρ͔+θ͔ξ͕+θ͕ξ͔-κ͔μ͕-κ͕μ͔)e2364",
    ///     "+(+d͕λ͔-e͕κ͔+h͕ο͔-i͕ξ͔+k͕σ͔-l͕ρ͔+m͕υ͔-n͕τ͔+o͕v͔+æ͕α͔+ñ͕β͔+ú͕γ͔+ü͕δ͔+η͔π͕+η͕π͔-θ͔ν͕-θ͕ν͔+ι͔μ͕+ι͕μ͔)e2345",
    ///     "+(-a͕γ͔+b͕β͔-f͕α͔+p͕v͔+q͕π͔+r͕ρ͔+s͕σ͔-t͕ν͔-u͕ξ͔-á͕ο͔+ç͕ι͔+é͕κ͔+ë͕λ͔-δ͔φ͕-δ͕φ͔+ε͔υ͕+ε͕υ͔-ζ͔τ͕-ζ͕τ͔)e1576",
    ///     "+(-a͕δ͔+c͕β͔-g͕α͔-p͕π͔+q͕v͔+r͕τ͔+s͕υ͔+t͕μ͔-ä͕ξ͔-å͕ο͔-ç͕θ͔+í͕κ͔+ï͕λ͔+γ͔φ͕+γ͕φ͔-ε͔σ͕-ε͕σ͔+ζ͔ρ͕+ζ͕ρ͔)e1467",
    ///     "+(-a͕ε͔+d͕β͔-h͕α͔-p͕ρ͔-q͕τ͔+r͕v͔+s͕φ͔+u͕μ͔+ä͕ν͔-æ͕ο͔-é͕θ͔-í͕ι͔+ñ͕λ͔-γ͔υ͕-γ͕υ͔+δ͔σ͕+δ͕σ͔-ζ͔π͕-ζ͕π͔)e1475",
    ///     "+(-a͕ζ͔+e͕β͔-i͕α͔-p͕σ͔-q͕υ͔-r͕φ͔+s͕v͔+á͕μ͔+å͕ν͔+æ͕ξ͔-ë͕θ͔-ï͕ι͔-ñ͕κ͔+γ͔τ͕+γ͕τ͔-δ͔ρ͕-δ͕ρ͔+ε͔π͕+ε͕π͔)e1456",
    ///     "+(-b͕δ͔+c͕γ͔-j͕α͔+p͕ν͔-q͕μ͔+t͕v͔+u͕τ͔+á͕υ͔-ä͕ρ͔-å͕σ͔+ç͕η͔+ó͕κ͔+ö͕λ͔-β͔φ͕-β͕φ͔+ε͔ο͕+ε͕ο͔-ζ͔ξ͕-ζ͕ξ͔)e1376",
    ///     "+(-b͕ε͔+d͕γ͔-k͕α͔+p͕ξ͔-r͕μ͔-t͕τ͔+u͕v͔+á͕φ͔+ä͕π͔-æ͕σ͔+é͕η͔-ó͕ι͔+ú͕λ͔+β͔υ͕+β͕υ͔-δ͔ο͕-δ͕ο͔+ζ͔ν͕+ζ͕ν͔)e1357",
    ///     "+(-b͕ζ͔+e͕γ͔-l͕α͔+p͕ο͔-s͕μ͔-t͕υ͔-u͕φ͔+v͔á͕+å͕π͔+æ͕ρ͔+ë͕η͔-ö͕ι͔-ú͕κ͔-β͔τ͕-β͕τ͔+δ͔ξ͕+δ͕ξ͔-ε͔ν͕-ε͕ν͔)e1365",
    ///     "+(-c͕ε͔+d͕δ͔-m͕α͔+q͕ξ͔-r͕ν͔+t͕ρ͔-u͕π͔+v͔ä͕+å͕φ͔-æ͕υ͔+í͕η͔+ó͕θ͔+ü͕λ͔-β͔σ͕-β͕σ͔+γ͔ο͕+γ͕ο͔-ζ͔μ͕-ζ͕μ͔)e1374",
    ///     "+(-c͕ζ͔+e͕δ͔-n͕α͔+q͕ο͔-s͕ν͔+t͕σ͔+v͔å͕-á͕π͔-ä͕φ͔+æ͕τ͔+ï͕η͔+ö͕θ͔-ü͕κ͔+β͔ρ͕+β͕ρ͔-γ͔ξ͕-γ͕ξ͔+ε͔μ͕+ε͕μ͔)e1346",
    ///     "+(-d͕ζ͔+e͕ε͔-o͕α͔+r͕ο͔-s͕ξ͔+u͕σ͔+v͔æ͕-á͕ρ͔+ä͕υ͔-å͕τ͔+ñ͕η͔+ú͕θ͔+ü͕ι͔-β͔π͕-β͕π͔+γ͔ν͕+γ͕ν͔-δ͔μ͕-δ͕μ͔)e1354",
    ///     "+(-f͕δ͔+g͕γ͔-j͕β͔-p͕ι͔+q͕θ͔-t͕η͔+v͔ç͕+é͕τ͔+ë͕υ͔-í͕ρ͔-ï͕σ͔+ó͕ξ͔+ö͕ο͔+α͔φ͕+α͕φ͔-ε͔λ͕-ε͕λ͔+ζ͔κ͕+ζ͕κ͔)e1267",
    ///     "+(-f͕ε͔+h͕γ͔-k͕β͔-p͕κ͔+r͕θ͔-u͕η͔+v͔é͕-ç͕τ͔+ë͕φ͔+í͕π͔-ñ͕σ͔-ó͕ν͔+ú͕ο͔-α͔υ͕-α͕υ͔+δ͔λ͕+δ͕λ͔-ζ͔ι͕-ζ͕ι͔)e1275",
    ///     "+(-f͕ζ͔+i͕γ͔-l͕β͔-p͕λ͔+s͕θ͔+v͔ë͕-á͕η͔-ç͕υ͔-é͕φ͔+ï͕π͔+ñ͕ρ͔-ö͕ν͔-ú͕ξ͔+α͔τ͕+α͕τ͔-δ͔κ͕-δ͕κ͔+ε͔ι͕+ε͕ι͔)e1256",
    ///     "+(-g͕ε͔+h͕δ͔-m͕β͔-q͕κ͔+r͕ι͔+v͔í͕-ä͕η͔+ç͕ρ͔-é͕π͔+ï͕φ͔-ñ͕υ͔+ó͕μ͔+ü͕ο͔+α͔σ͕+α͕σ͔-γ͔λ͕-γ͕λ͔+ζ͔θ͕+ζ͕θ͔)e1247",
    ///     "+(-g͕ζ͔+i͕δ͔-n͕β͔-q͕λ͔+s͕ι͔+v͔ï͕-å͕η͔+ç͕σ͔-ë͕π͔-í͕φ͔+ñ͕τ͔+ö͕μ͔-ü͕ξ͔-α͔ρ͕-α͕ρ͔+γ͔κ͕+γ͕κ͔-ε͔θ͕-ε͕θ͔)e1264",
    ///     "+(-h͕ζ͔+i͕ε͔-o͕β͔-r͕λ͔+s͕κ͔+v͔ñ͕-æ͕η͔+é͕σ͔-ë͕ρ͔+í͕υ͔-ï͕τ͔+ú͕μ͔+ü͕ν͔+α͔π͕+α͕π͔-γ͔ι͕-γ͕ι͔+δ͔θ͕+δ͕θ͔)e1245",
    ///     "+(-j͕ε͔+k͕δ͔-m͕γ͔-t͕κ͔+u͕ι͔+v͔ó͕-ä͕θ͔-ç͕ξ͔+é͕ν͔-í͕μ͔+ö͕φ͔-ú͕υ͔+ü͕σ͔-α͔ο͕-α͕ο͔+β͔λ͕+β͕λ͔-ζ͔η͕-ζ͕η͔)e1273",
    ///     "+(-j͕ζ͔+l͕δ͔-n͕γ͔-t͕λ͔+v͔ö͕+á͕ι͔-å͕θ͔-ç͕ο͔+ë͕ν͔-ï͕μ͔-ó͕φ͔+ú͕τ͔-ü͕ρ͔+α͔ξ͕+α͕ξ͔-β͔κ͕-β͕κ͔+ε͔η͕+ε͕η͔)e1236",
    ///     "+(-k͕ζ͔+l͕ε͔-o͕γ͔-u͕λ͔+v͔ú͕+á͕κ͔-æ͕θ͔-é͕ο͔+ë͕ξ͔-ñ͕μ͔+ó͕υ͔-ö͕τ͔+ü͕π͔-α͔ν͕-α͕ν͔+β͔ι͕+β͕ι͔-δ͔η͕-δ͕η͔)e1253",
    ///     "+(+a͕η͔+b͕θ͔+c͕ι͔+d͕κ͔+e͕λ͔+f͕μ͔+g͕ν͔+h͕ξ͔+i͕ο͔+j͕π͔+k͕ρ͔+l͕σ͔+m͕τ͔+n͕υ͔+o͕φ͔)e234567",
    ///     "+(-a͕β͔-b͕γ͔-c͕δ͔-d͕ε͔-e͕ζ͔+p͕μ͔+q͕ν͔+r͕ξ͔+s͕ο͔+t͕π͔+u͕ρ͔+á͕σ͔+ä͕τ͔+å͕υ͔+æ͕φ͔)e134576",
    ///     "+(+a͕α͔-f͕γ͔-g͕δ͔-h͕ε͔-i͕ζ͔-p͕θ͔-q͕ι͔-r͕κ͔-s͕λ͔+ç͕π͔+é͕ρ͔+ë͕σ͔+í͕τ͔+ï͕υ͔+ñ͕φ͔)e124567",
    ///     "+(+b͕α͔+f͕β͔-j͕δ͔-k͕ε͔-l͕ζ͔+p͕η͔-t͕ι͔-u͕κ͔-á͕λ͔-ç͕ν͔-é͕ξ͔-ë͕ο͔+ó͕τ͔+ö͕υ͔+ú͕φ͔)e123576",
    ///     "+(+c͕α͔+g͕β͔+j͕γ͔-m͕ε͔-n͕ζ͔+q͕η͔+t͕θ͔-ä͕κ͔-å͕λ͔+ç͕μ͔-í͕ξ͔-ï͕ο͔-ó͕ρ͔-ö͕σ͔+ü͕φ͔)e123467",
    ///     "+(+d͕α͔+h͕β͔+k͕γ͔+m͕δ͔-o͕ζ͔+r͕η͔+u͕θ͔+ä͕ι͔-æ͕λ͔+é͕μ͔+í͕ν͔-ñ͕ο͔+ó͕π͔-ú͕σ͔-ü͕υ͔)e123475",
    ///     "+(+e͕α͔+i͕β͔+l͕γ͔+n͕δ͔+o͕ε͔+s͕η͔+á͕θ͔+å͕ι͔+æ͕κ͔+ë͕μ͔+ï͕ν͔+ñ͕ξ͔+ö͕π͔+ú͕ρ͔+ü͕τ͔)e123456",
    /// ]);
    ///
    /// let triple_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(triple_rotator.basis_blades(), Vee::triple_rotator().basis_blades());
    /// format_eq!(triple_rotator, [
    ///     "+a͔a͕+b͔b͕+c͔c͕+d͔d͕+e͔e͕+f͔f͕+g͔g͕+h͔h͕+i͔i͕+j͔j͕+k͔k͕+l͔l͕+m͔m͕+n͔n͕+o͔o͕+p͔p͕+q͔q͕\
    ///      +r͔r͕+s͔s͕+t͔t͕+u͔u͕+á͔á͕+ä͔ä͕+å͔å͕+æ͔æ͕+ç͔ç͕+é͔é͕+ë͔ë͕+í͔í͕+ï͔ï͕+ñ͔ñ͕+ó͔ó͕+ö͔ö͕+ú͔ú͕+ü͔ü͕",
    ///     "+(+f͔p͕-f͕p͔+g͔q͕-g͕q͔+h͔r͕-h͕r͔+i͔s͕-i͕s͔+j͔t͕-j͕t͔+k͔u͕-k͕u͔+l͔á͕-l͕á͔+m͔ä͕-m͕ä͔+n͔å͕-n͕å͔+o͔æ͕-o͕æ͔)e12",
    ///     "+(-b͔p͕+b͕p͔-c͔q͕+c͕q͔-d͔r͕+d͕r͔-e͔s͕+e͕s͔+j͔ç͕-j͕ç͔+k͔é͕-k͕é͔+l͔ë͕-l͕ë͔+m͔í͕-m͕í͔+n͔ï͕-n͕ï͔+o͔ñ͕-o͕ñ͔)e13",
    ///     "+(+a͔p͕-a͕p͔-c͔t͕+c͕t͔-d͔u͕+d͕u͔-e͔á͕+e͕á͔-g͔ç͕+g͕ç͔-h͔é͕+h͕é͔-i͔ë͕+i͕ë͔+m͔ó͕-m͕ó͔+n͔ö͕-n͕ö͔+o͔ú͕-o͕ú͔)e14",
    ///     "+(+a͔q͕-a͕q͔+b͔t͕-b͕t͔-d͔ä͕+d͕ä͔-e͔å͕+e͕å͔+f͔ç͕-f͕ç͔-h͔í͕+h͕í͔-i͔ï͕+i͕ï͔-k͔ó͕+k͕ó͔-l͔ö͕+l͕ö͔+o͔ü͕-o͕ü͔)e15",
    ///     "+(+a͔r͕-a͕r͔+b͔u͕-b͕u͔+c͔ä͕-c͕ä͔-e͔æ͕+e͕æ͔+f͔é͕-f͕é͔+g͔í͕-g͕í͔-i͔ñ͕+i͕ñ͔+j͔ó͕-j͕ó͔-l͔ú͕+l͕ú͔-n͔ü͕+n͕ü͔)e16",
    ///     "+(+a͔s͕-a͕s͔+b͔á͕-b͕á͔+c͔å͕-c͕å͔+d͔æ͕-d͕æ͔+f͔ë͕-f͕ë͔+g͔ï͕-g͕ï͔+h͔ñ͕-h͕ñ͔+j͔ö͕-j͕ö͔+k͔ú͕-k͕ú͔+m͔ü͕-m͕ü͔)e17",
    ///     "+(+b͔f͕-b͕f͔+c͔g͕-c͕g͔+d͔h͕-d͕h͔+e͔i͕-e͕i͔+t͔ç͕-t͕ç͔+u͔é͕-u͕é͔+á͔ë͕-á͕ë͔+ä͔í͕-ä͕í͔+å͔ï͕-å͕ï͔+æ͔ñ͕-æ͕ñ͔)e23",
    ///     "+(-a͔f͕+a͕f͔+c͔j͕-c͕j͔+d͔k͕-d͕k͔+e͔l͕-e͕l͔-q͔ç͕+q͕ç͔-r͔é͕+r͕é͔-s͔ë͕+s͕ë͔+ä͔ó͕-ä͕ó͔+å͔ö͕-å͕ö͔+æ͔ú͕-æ͕ú͔)e24",
    ///     "+(-a͔g͕+a͕g͔-b͔j͕+b͕j͔+d͔m͕-d͕m͔+e͔n͕-e͕n͔+p͔ç͕-p͕ç͔-r͔í͕+r͕í͔-s͔ï͕+s͕ï͔-u͔ó͕+u͕ó͔-á͔ö͕+á͕ö͔+æ͔ü͕-æ͕ü͔)e25",
    ///     "+(-a͔h͕+a͕h͔-b͔k͕+b͕k͔-c͔m͕+c͕m͔+e͔o͕-e͕o͔+p͔é͕-p͕é͔+q͔í͕-q͕í͔-s͔ñ͕+s͕ñ͔+t͔ó͕-t͕ó͔-á͔ú͕+á͕ú͔-å͔ü͕+å͕ü͔)e26",
    ///     "+(-a͔i͕+a͕i͔-b͔l͕+b͕l͔-c͔n͕+c͕n͔-d͔o͕+d͕o͔+p͔ë͕-p͕ë͔+q͔ï͕-q͕ï͔+r͔ñ͕-r͕ñ͔+t͔ö͕-t͕ö͔+u͔ú͕-u͕ú͔+ä͔ü͕-ä͕ü͔)e27",
    ///     "+(+a͔b͕-a͕b͔+g͔j͕-g͕j͔+h͔k͕-h͕k͔+i͔l͕-i͕l͔+q͔t͕-q͕t͔+r͔u͕-r͕u͔+s͔á͕-s͕á͔+í͔ó͕-í͕ó͔+ï͔ö͕-ï͕ö͔+ñ͔ú͕-ñ͕ú͔)e34",
    ///     "+(+a͔c͕-a͕c͔-f͔j͕+f͕j͔+h͔m͕-h͕m͔+i͔n͕-i͕n͔-p͔t͕+p͕t͔+r͔ä͕-r͕ä͔+s͔å͕-s͕å͔-é͔ó͕+é͕ó͔-ë͔ö͕+ë͕ö͔+ñ͔ü͕-ñ͕ü͔)e35",
    ///     "+(+a͔d͕-a͕d͔-f͔k͕+f͕k͔-g͔m͕+g͕m͔+i͔o͕-i͕o͔-p͔u͕+p͕u͔-q͔ä͕+q͕ä͔+s͔æ͕-s͕æ͔+ç͔ó͕-ç͕ó͔-ë͔ú͕+ë͕ú͔-ï͔ü͕+ï͕ü͔)e36",
    ///     "+(+a͔e͕-a͕e͔-f͔l͕+f͕l͔-g͔n͕+g͕n͔-h͔o͕+h͕o͔-p͔á͕+p͕á͔-q͔å͕+q͕å͔-r͔æ͕+r͕æ͔+ç͔ö͕-ç͕ö͔+é͔ú͕-é͕ú͔+í͔ü͕-í͕ü͔)e37",
    ///     "+(+b͔c͕-b͕c͔+f͔g͕-f͕g͔+k͔m͕-k͕m͔+l͔n͕-l͕n͔+p͔q͕-p͕q͔+u͔ä͕-u͕ä͔+á͔å͕-á͕å͔+é͔í͕-é͕í͔+ë͔ï͕-ë͕ï͔+ú͔ü͕-ú͕ü͔)e45",
    ///     "+(+b͔d͕-b͕d͔+f͔h͕-f͕h͔-j͔m͕+j͕m͔+l͔o͕-l͕o͔+p͔r͕-p͕r͔-t͔ä͕+t͕ä͔+á͔æ͕-á͕æ͔-ç͔í͕+ç͕í͔+ë͔ñ͕-ë͕ñ͔-ö͔ü͕+ö͕ü͔)e46",
    ///     "+(+b͔e͕-b͕e͔+f͔i͕-f͕i͔-j͔n͕+j͕n͔-k͔o͕+k͕o͔+p͔s͕-p͕s͔-t͔å͕+t͕å͔-u͔æ͕+u͕æ͔-ç͔ï͕+ç͕ï͔-é͔ñ͕+é͕ñ͔+ó͔ü͕-ó͕ü͔)e47",
    ///     "+(+c͔d͕-c͕d͔+g͔h͕-g͕h͔+j͔k͕-j͕k͔+n͔o͕-n͕o͔+q͔r͕-q͕r͔+t͔u͕-t͕u͔+å͔æ͕-å͕æ͔+ç͔é͕-ç͕é͔+ï͔ñ͕-ï͕ñ͔+ö͔ú͕-ö͕ú͔)e56",
    ///     "+(+c͔e͕-c͕e͔+g͔i͕-g͕i͔+j͔l͕-j͕l͔-m͔o͕+m͕o͔+q͔s͕-q͕s͔+t͔á͕-t͕á͔-ä͔æ͕+ä͕æ͔+ç͔ë͕-ç͕ë͔-í͔ñ͕+í͕ñ͔-ó͔ú͕+ó͕ú͔)e57",
    ///     "+(+d͔e͕-d͕e͔+h͔i͕-h͕i͔+k͔l͕-k͕l͔+m͔n͕-m͕n͔+r͔s͕-r͕s͔+u͔á͕-u͕á͔+ä͔å͕-ä͕å͔+é͔ë͕-é͕ë͔+í͔ï͕-í͕ï͔+ó͔ö͕-ó͕ö͔)e67",
    ///     "+(-c͔ç͕-c͕ç͔-d͔é͕-d͕é͔-e͔ë͕-e͕ë͔+g͔t͕+g͕t͔+h͔u͕+h͕u͔+i͔á͕+i͕á͔-j͔q͕-j͕q͔-k͔r͕-k͕r͔-l͔s͕-l͕s͔)e1234",
    ///     "+(-j͔o͕-j͕o͔+k͔n͕+k͕n͔-l͔m͕-l͕m͔-t͔æ͕-t͕æ͔+u͔å͕+u͕å͔-á͔ä͕-á͕ä͔-ç͔ñ͕-ç͕ñ͔+é͔ï͕+é͕ï͔-ë͔í͕-ë͕í͔)e4567",
    ///     "+(+g͔o͕+g͕o͔-h͔n͕-h͕n͔+i͔m͕+i͕m͔+q͔æ͕+q͕æ͔-r͔å͕-r͕å͔+s͔ä͕+s͕ä͔-ç͔ú͕-ç͕ú͔+é͔ö͕+é͕ö͔-ë͔ó͕-ë͕ó͔)e3576",
    ///     "+(-f͔o͕-f͕o͔+h͔l͕+h͕l͔-i͔k͕-i͕k͔-p͔æ͕-p͕æ͔+r͔á͕+r͕á͔-s͔u͕-s͕u͔-ç͔ü͕-ç͕ü͔+í͔ö͕+í͕ö͔-ï͔ó͕-ï͕ó͔)e3467",
    ///     "+(+f͔n͕+f͕n͔-g͔l͕-g͕l͔+i͔j͕+i͕j͔+p͔å͕+p͕å͔-q͔á͕-q͕á͔+s͔t͕+s͕t͔-é͔ü͕-é͕ü͔+í͔ú͕+í͕ú͔-ñ͔ó͕-ñ͕ó͔)e3475",
    ///     "+(-f͔m͕-f͕m͔+g͔k͕+g͕k͔-h͔j͕-h͕j͔-p͔ä͕-p͕ä͔+q͔u͕+q͕u͔-r͔t͕-r͕t͔-ë͔ü͕-ë͕ü͔+ï͔ú͕+ï͕ú͔-ñ͔ö͕-ñ͕ö͔)e3456",
    ///     "+(-c͔o͕-c͕o͔+d͔n͕+d͕n͔-e͔m͕-e͕m͔+q͔ñ͕+q͕ñ͔-r͔ï͕-r͕ï͔+s͔í͕+s͕í͔+t͔ú͕+t͕ú͔-u͔ö͕-u͕ö͔+á͔ó͕+á͕ó͔)e2567",
    ///     "+(+b͔o͕+b͕o͔-d͔l͕-d͕l͔+e͔k͕+e͕k͔-p͔ñ͕-p͕ñ͔+r͔ë͕+r͕ë͔-s͔é͕-s͕é͔+t͔ü͕+t͕ü͔-ä͔ö͕-ä͕ö͔+å͔ó͕+å͕ó͔)e2476",
    ///     "+(-b͔n͕-b͕n͔+c͔l͕+c͕l͔-e͔j͕-e͕j͔+p͔ï͕+p͕ï͔-q͔ë͕-q͕ë͔+s͔ç͕+s͕ç͔+u͔ü͕+u͕ü͔-ä͔ú͕-ä͕ú͔+æ͔ó͕+æ͕ó͔)e2457",
    ///     "+(+b͔m͕+b͕m͔-c͔k͕-c͕k͔+d͔j͕+d͕j͔-p͔í͕-p͕í͔+q͔é͕+q͕é͔-r͔ç͕-r͕ç͔+á͔ü͕+á͕ü͔-å͔ú͕-å͕ú͔+æ͔ö͕+æ͕ö͔)e2465",
    ///     "+(-a͔o͕-a͕o͔+d͔i͕+d͕i͔-e͔h͕-e͕h͔-p͔ú͕-p͕ú͔-q͔ü͕-q͕ü͔+u͔ë͕+u͕ë͔-á͔é͕-á͕é͔+ä͔ï͕+ä͕ï͔-å͔í͕-å͕í͔)e2367",
    ///     "+(+a͔n͕+a͕n͔-c͔i͕-c͕i͔+e͔g͕+e͕g͔+p͔ö͕+p͕ö͔-r͔ü͕-r͕ü͔-t͔ë͕-t͕ë͔+á͔ç͕+á͕ç͔+ä͔ñ͕+ä͕ñ͔-æ͔í͕-æ͕í͔)e2375",
    ///     "+(-a͔m͕-a͕m͔+c͔h͕+c͕h͔-d͔g͕-d͕g͔-p͔ó͕-p͕ó͔-s͔ü͕-s͕ü͔+t͔é͕+t͕é͔-u͔ç͕-u͕ç͔+å͔ñ͕+å͕ñ͔-æ͔ï͕-æ͕ï͔)e2356",
    ///     "+(-a͔l͕-a͕l͔+b͔i͕+b͕i͔-e͔f͕-e͕f͔+q͔ö͕+q͕ö͔+r͔ú͕+r͕ú͔-t͔ï͕-t͕ï͔-u͔ñ͕-u͕ñ͔+å͔ç͕+å͕ç͔+æ͔é͕+æ͕é͔)e2347",
    ///     "+(+a͔k͕+a͕k͔-b͔h͕-b͕h͔+d͔f͕+d͕f͔-q͔ó͕-q͕ó͔+s͔ú͕+s͕ú͔+t͔í͕+t͕í͔-á͔ñ͕-á͕ñ͔-ä͔ç͕-ä͕ç͔+æ͔ë͕+æ͕ë͔)e2364",
    ///     "+(-a͔j͕-a͕j͔+b͔g͕+b͕g͔-c͔f͕-c͕f͔-r͔ó͕-r͕ó͔-s͔ö͕-s͕ö͔+u͔í͕+u͕í͔+á͔ï͕+á͕ï͔-ä͔é͕-ä͕é͔-å͔ë͕-å͕ë͔)e2345",
    ///     "+(-c͔æ͕-c͕æ͔+d͔å͕+d͕å͔-e͔ä͕-e͕ä͔-g͔ñ͕-g͕ñ͔+h͔ï͕+h͕ï͔-i͔í͕-i͕í͔-j͔ú͕-j͕ú͔+k͔ö͕+k͕ö͔-l͔ó͕-l͕ó͔)e1576",
    ///     "+(+b͔æ͕+b͕æ͔-d͔á͕-d͕á͔+e͔u͕+e͕u͔+f͔ñ͕+f͕ñ͔-h͔ë͕-h͕ë͔+i͔é͕+i͕é͔-j͔ü͕-j͕ü͔+m͔ö͕+m͕ö͔-n͔ó͕-n͕ó͔)e1467",
    ///     "+(-b͔å͕-b͕å͔+c͔á͕+c͕á͔-e͔t͕-e͕t͔-f͔ï͕-f͕ï͔+g͔ë͕+g͕ë͔-i͔ç͕-i͕ç͔-k͔ü͕-k͕ü͔+m͔ú͕+m͕ú͔-o͔ó͕-o͕ó͔)e1475",
    ///     "+(+b͔ä͕+b͕ä͔-c͔u͕-c͕u͔+d͔t͕+d͕t͔+f͔í͕+f͕í͔-g͔é͕-g͕é͔+h͔ç͕+h͕ç͔-l͔ü͕-l͕ü͔+n͔ú͕+n͕ú͔-o͔ö͕-o͕ö͔)e1456",
    ///     "+(-a͔æ͕-a͕æ͔+d͔s͕+d͕s͔-e͔r͕-e͕r͔+f͔ú͕+f͕ú͔+g͔ü͕+g͕ü͔-k͔ë͕-k͕ë͔+l͔é͕+l͕é͔-m͔ï͕-m͕ï͔+n͔í͕+n͕í͔)e1376",
    ///     "+(+a͔å͕+a͕å͔-c͔s͕-c͕s͔+e͔q͕+e͕q͔-f͔ö͕-f͕ö͔+h͔ü͕+h͕ü͔+j͔ë͕+j͕ë͔-l͔ç͕-l͕ç͔-m͔ñ͕-m͕ñ͔+o͔í͕+o͕í͔)e1357",
    ///     "+(-a͔ä͕-a͕ä͔+c͔r͕+c͕r͔-d͔q͕-d͕q͔+f͔ó͕+f͕ó͔+i͔ü͕+i͕ü͔-j͔é͕-j͕é͔+k͔ç͕+k͕ç͔-n͔ñ͕-n͕ñ͔+o͔ï͕+o͕ï͔)e1365",
    ///     "+(-a͔á͕-a͕á͔+b͔s͕+b͕s͔-e͔p͕-e͕p͔-g͔ö͕-g͕ö͔-h͔ú͕-h͕ú͔+j͔ï͕+j͕ï͔+k͔ñ͕+k͕ñ͔-n͔ç͕-n͕ç͔-o͔é͕-o͕é͔)e1374",
    ///     "+(+a͔u͕+a͕u͔-b͔r͕-b͕r͔+d͔p͕+d͕p͔+g͔ó͕+g͕ó͔-i͔ú͕-i͕ú͔-j͔í͕-j͕í͔+l͔ñ͕+l͕ñ͔+m͔ç͕+m͕ç͔-o͔ë͕-o͕ë͔)e1346",
    ///     "+(-a͔t͕-a͕t͔+b͔q͕+b͕q͔-c͔p͕-c͕p͔+h͔ó͕+h͕ó͔+i͔ö͕+i͕ö͔-k͔í͕-k͕í͔-l͔ï͕-l͕ï͔+m͔é͕+m͕é͔+n͔ë͕+n͕ë͔)e1354",
    ///     "+(-a͔ñ͕-a͕ñ͔-b͔ú͕-b͕ú͔-c͔ü͕-c͕ü͔+h͔s͕+h͕s͔-i͔r͕-i͕r͔+k͔á͕+k͕á͔-l͔u͕-l͕u͔+m͔å͕+m͕å͔-n͔ä͕-n͕ä͔)e1267",
    ///     "+(+a͔ï͕+a͕ï͔+b͔ö͕+b͕ö͔-d͔ü͕-d͕ü͔-g͔s͕-g͕s͔+i͔q͕+i͕q͔-j͔á͕-j͕á͔+l͔t͕+l͕t͔+m͔æ͕+m͕æ͔-o͔ä͕-o͕ä͔)e1275",
    ///     "+(-a͔í͕-a͕í͔-b͔ó͕-b͕ó͔-e͔ü͕-e͕ü͔+g͔r͕+g͕r͔-h͔q͕-h͕q͔+j͔u͕+j͕u͔-k͔t͕-k͕t͔+n͔æ͕+n͕æ͔-o͔å͕-o͕å͔)e1256",
    ///     "+(-a͔ë͕-a͕ë͔+c͔ö͕+c͕ö͔+d͔ú͕+d͕ú͔+f͔s͕+f͕s͔-i͔p͕-i͕p͔-j͔å͕-j͕å͔-k͔æ͕-k͕æ͔+n͔t͕+n͕t͔+o͔u͕+o͕u͔)e1247",
    ///     "+(+a͔é͕+a͕é͔-c͔ó͕-c͕ó͔+e͔ú͕+e͕ú͔-f͔r͕-f͕r͔+h͔p͕+h͕p͔+j͔ä͕+j͕ä͔-l͔æ͕-l͕æ͔-m͔t͕-m͕t͔+o͔á͕+o͕á͔)e1264",
    ///     "+(-a͔ç͕-a͕ç͔-d͔ó͕-d͕ó͔-e͔ö͕-e͕ö͔+f͔q͕+f͕q͔-g͔p͕-g͕p͔+k͔ä͕+k͕ä͔+l͔å͕+l͕å͔-m͔u͕-m͕u͔-n͔á͕-n͕á͔)e1245",
    ///     "+(-b͔ë͕-b͕ë͔-c͔ï͕-c͕ï͔-d͔ñ͕-d͕ñ͔+f͔á͕+f͕á͔+g͔å͕+g͕å͔+h͔æ͕+h͕æ͔-l͔p͕-l͕p͔-n͔q͕-n͕q͔-o͔r͕-o͕r͔)e1273",
    ///     "+(+b͔é͕+b͕é͔+c͔í͕+c͕í͔-e͔ñ͕-e͕ñ͔-f͔u͕-f͕u͔-g͔ä͕-g͕ä͔+i͔æ͕+i͕æ͔+k͔p͕+k͕p͔+m͔q͕+m͕q͔-o͔s͕-o͕s͔)e1236",
    ///     "+(-b͔ç͕-b͕ç͔+d͔í͕+d͕í͔+e͔ï͕+e͕ï͔+f͔t͕+f͕t͔-h͔ä͕-h͕ä͔-i͔å͕-i͕å͔-j͔p͕-j͕p͔+m͔r͕+m͕r͔+n͔s͕+n͕s͔)e1253",
    ///     "+(-p͔ü͕+p͕ü͔+q͔ú͕-q͕ú͔-r͔ö͕+r͕ö͔+s͔ó͕-s͕ó͔-t͔ñ͕+t͕ñ͔+u͔ï͕-u͕ï͔-á͔í͕+á͕í͔-ä͔ë͕+ä͕ë͔+å͔é͕-å͕é͔-æ͔ç͕+æ͕ç͔)e234567",
    ///     "+(+f͔ü͕-f͕ü͔-g͔ú͕+g͕ú͔+h͔ö͕-h͕ö͔-i͔ó͕+i͕ó͔+j͔ñ͕-j͕ñ͔-k͔ï͕+k͕ï͔+l͔í͕-l͕í͔+m͔ë͕-m͕ë͔-n͔é͕+n͕é͔+o͔ç͕-o͕ç͔)e134576",
    ///     "+(-b͔ü͕+b͕ü͔+c͔ú͕-c͕ú͔-d͔ö͕+d͕ö͔+e͔ó͕-e͕ó͔-j͔æ͕+j͕æ͔+k͔å͕-k͕å͔-l͔ä͕+l͕ä͔-m͔á͕+m͕á͔+n͔u͕-n͕u͔-o͔t͕+o͕t͔)e124567",
    ///     "+(+a͔ü͕-a͕ü͔-c͔ñ͕+c͕ñ͔+d͔ï͕-d͕ï͔-e͔í͕+e͕í͔+g͔æ͕-g͕æ͔-h͔å͕+h͕å͔+i͔ä͕-i͕ä͔+m͔s͕-m͕s͔-n͔r͕+n͕r͔+o͔q͕-o͕q͔)e123576",
    ///     "+(-a͔ú͕+a͕ú͔+b͔ñ͕-b͕ñ͔-d͔ë͕+d͕ë͔+e͔é͕-e͕é͔-f͔æ͕+f͕æ͔+h͔á͕-h͕á͔-i͔u͕+i͕u͔-k͔s͕+k͕s͔+l͔r͕-l͕r͔-o͔p͕+o͕p͔)e123467",
    ///     "+(+a͔ö͕-a͕ö͔-b͔ï͕+b͕ï͔+c͔ë͕-c͕ë͔-e͔ç͕+e͕ç͔+f͔å͕-f͕å͔-g͔á͕+g͕á͔+i͔t͕-i͕t͔+j͔s͕-j͕s͔-l͔q͕+l͕q͔+n͔p͕-n͕p͔)e123475",
    ///     "+(-a͔ó͕+a͕ó͔+b͔í͕-b͕í͔-c͔é͕+c͕é͔+d͔ç͕-d͕ç͔-f͔ä͕+f͕ä͔+g͔u͕-g͕u͔-h͔t͕+h͕t͔-j͔r͕+j͕r͔+k͔q͕-k͕q͔-m͔p͕+m͕p͔)e123456",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_rotator() -> Self {
        Self::scalar()
            + Self::volume5_displacement()
            + Self::volume_displacement()
            + Self::line_displacement()
    }
    /// The multivector of translator $`t \equiv s + v^5_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// format_eq!(translator, [
    ///     "-w͔w͕",
    ///     "+(+X͔w͕-X͕w͔)e01",
    ///     "+(+Y͔w͕-Y͕w͔)e02",
    ///     "+(+Z͔w͕-Z͕w͔)e03",
    ///     "+(-w͔Ð͕+w͕Ð͔)e04",
    ///     "+(-w͔Ø͕+w͕Ø͔)e05",
    ///     "+(-w͔Þ͕+w͕Þ͔)e06",
    ///     "+(-w͔Œ͕+w͕Œ͔)e07",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::volume5_moment()
    }
    /// The multivector of simple single motor $`m_{s1} \equiv s + v^5`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let simple_single_motor = Vee::volume6().lhs() * Vee::volume6().rhs();
    ///
    /// assert_eq!(simple_single_motor.basis_blades(), Vee::simple_single_motor().basis_blades());
    /// format_eq!(simple_single_motor, [
    ///     "+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕+œ͔œ͕",
    ///     "+(+W͔x͕-W͕x͔)e01",
    ///     "+(+W͔y͕-W͕y͔)e02",
    ///     "+(+W͔z͕-W͕z͔)e03",
    ///     "+(+W͔ð͕-W͕ð͔)e04",
    ///     "+(+W͔ø͕-W͕ø͔)e05",
    ///     "+(+W͔þ͕-W͕þ͔)e06",
    ///     "+(+W͔œ͕-W͕œ͔)e07",
    ///     "+(+x͔y͕-x͕y͔)e12",
    ///     "+(+x͔z͕-x͕z͔)e13",
    ///     "+(+x͔ð͕-x͕ð͔)e14",
    ///     "+(+x͔ø͕-x͕ø͔)e15",
    ///     "+(+x͔þ͕-x͕þ͔)e16",
    ///     "+(+x͔œ͕-x͕œ͔)e17",
    ///     "+(+y͔z͕-y͕z͔)e23",
    ///     "+(+y͔ð͕-y͕ð͔)e24",
    ///     "+(+y͔ø͕-y͕ø͔)e25",
    ///     "+(+y͔þ͕-y͕þ͔)e26",
    ///     "+(+y͔œ͕-y͕œ͔)e27",
    ///     "+(+z͔ð͕-z͕ð͔)e34",
    ///     "+(+z͔ø͕-z͕ø͔)e35",
    ///     "+(+z͔þ͕-z͕þ͔)e36",
    ///     "+(+z͔œ͕-z͕œ͔)e37",
    ///     "+(+ð͔ø͕-ð͕ø͔)e45",
    ///     "+(+ð͔þ͕-ð͕þ͔)e46",
    ///     "+(+ð͔œ͕-ð͕œ͔)e47",
    ///     "+(+ø͔þ͕-ø͕þ͔)e56",
    ///     "+(+ø͔œ͕-ø͕œ͔)e57",
    ///     "+(+þ͔œ͕-þ͕œ͔)e67",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_motor() -> Self {
        Self::scalar() + Self::volume5()
    }
    /// The multivector of single motor $`m_1 \equiv s + v^5 + v_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔+Œ͕ζ͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕η͔+Ð͕θ͔+Ø͕ι͔+Þ͕κ͔+Œ͕λ͔)e02",
    ///     "+(-X͕β͔-Y͕η͔+Z͕v͔+Ð͕μ͔+Ø͕ν͔+Þ͕ξ͔+Œ͕ο͔)e03",
    ///     "+(-X͕γ͔-Y͕θ͔-Z͕μ͔+v͔Ð͕+Ø͕π͔+Þ͕ρ͔+Œ͕σ͔)e04",
    ///     "+(-X͕δ͔-Y͕ι͔-Z͕ν͔+v͔Ø͕-Ð͕π͔+Þ͕τ͔+Œ͕υ͔)e05",
    ///     "+(-X͕ε͔-Y͕κ͔-Z͕ξ͔+v͔Þ͕-Ð͕ρ͔-Ø͕τ͔+Œ͕φ͔)e06",
    ///     "+(-X͕ζ͔-Y͕λ͔-Z͕ο͔+v͔Œ͕-Ð͕σ͔-Ø͕υ͔-Þ͕φ͔)e07",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e17",
    ///     "+v͕η͔e23",
    ///     "+v͕θ͔e24",
    ///     "+v͕ι͔e25",
    ///     "+v͕κ͔e26",
    ///     "+v͕λ͔e27",
    ///     "+v͕μ͔e34",
    ///     "+v͕ν͔e35",
    ///     "+v͕ξ͔e36",
    ///     "+v͕ο͔e37",
    ///     "+v͕π͔e45",
    ///     "+v͕ρ͔e46",
    ///     "+v͕σ͔e47",
    ///     "+v͕τ͔e56",
    ///     "+v͕υ͔e57",
    ///     "+v͕φ͔e67",
    ///     "+(+X͕η͔-Y͕β͔+Z͕α͔)e0123",
    ///     "+(+X͕θ͔-Y͕γ͔+Ð͕α͔)e0124",
    ///     "+(+X͕ι͔-Y͕δ͔+Ø͕α͔)e0125",
    ///     "+(+X͕κ͔-Y͕ε͔+Þ͕α͔)e0126",
    ///     "+(+X͕λ͔-Y͕ζ͔+Œ͕α͔)e0127",
    ///     "+(+X͕μ͔-Z͕γ͔+Ð͕β͔)e0134",
    ///     "+(+X͕ν͔-Z͕δ͔+Ø͕β͔)e0135",
    ///     "+(+X͕ξ͔-Z͕ε͔+Þ͕β͔)e0136",
    ///     "+(+X͕ο͔-Z͕ζ͔+Œ͕β͔)e0137",
    ///     "+(+X͕π͔-Ð͕δ͔+Ø͕γ͔)e0145",
    ///     "+(+X͕ρ͔-Ð͕ε͔+Þ͕γ͔)e0146",
    ///     "+(+X͕σ͔-Ð͕ζ͔+Œ͕γ͔)e0147",
    ///     "+(+X͕τ͔-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(+X͕υ͔-Ø͕ζ͔+Œ͕δ͔)e0157",
    ///     "+(+X͕φ͔-Þ͕ζ͔+Œ͕ε͔)e0167",
    ///     "+(+Y͕μ͔-Z͕θ͔+Ð͕η͔)e0234",
    ///     "+(+Y͕ν͔-Z͕ι͔+Ø͕η͔)e0235",
    ///     "+(+Y͕ξ͔-Z͕κ͔+Þ͕η͔)e0236",
    ///     "+(+Y͕ο͔-Z͕λ͔+Œ͕η͔)e0237",
    ///     "+(+Y͕π͔-Ð͕ι͔+Ø͕θ͔)e0245",
    ///     "+(+Y͕ρ͔-Ð͕κ͔+Þ͕θ͔)e0246",
    ///     "+(+Y͕σ͔-Ð͕λ͔+Œ͕θ͔)e0247",
    ///     "+(+Y͕τ͔-Ø͕κ͔+Þ͕ι͔)e0256",
    ///     "+(+Y͕υ͔-Ø͕λ͔+Œ͕ι͔)e0257",
    ///     "+(+Y͕φ͔-Þ͕λ͔+Œ͕κ͔)e0267",
    ///     "+(+Z͕π͔-Ð͕ν͔+Ø͕μ͔)e0345",
    ///     "+(+Z͕ρ͔-Ð͕ξ͔+Þ͕μ͔)e0346",
    ///     "+(+Z͕σ͔-Ð͕ο͔+Œ͕μ͔)e0347",
    ///     "+(+Z͕τ͔-Ø͕ξ͔+Þ͕ν͔)e0356",
    ///     "+(+Z͕υ͔-Ø͕ο͔+Œ͕ν͔)e0357",
    ///     "+(+Z͕φ͔-Þ͕ο͔+Œ͕ξ͔)e0367",
    ///     "+(+Ð͕τ͔-Ø͕ρ͔+Þ͕π͔)e0456",
    ///     "+(+Ð͕υ͔-Ø͕σ͔+Œ͕π͔)e0457",
    ///     "+(+Ð͕φ͔-Þ͕σ͔+Œ͕ρ͔)e0467",
    ///     "+(+Ø͕φ͔-Þ͕υ͔+Œ͕τ͔)e0567",
    /// ]);
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// format_eq!(single_motor, [
    ///     "-x͔x͕-y͔y͕-z͔z͕-ð͔ð͕-ø͔ø͕-þ͔þ͕-œ͔œ͕",
    ///     "+(-y͔Α͕+y͕Α͔-z͔Β͕+z͕Β͔-ð͔Γ͕+ð͕Γ͔-ø͔Δ͕+ø͕Δ͔-þ͔Ε͕+þ͕Ε͔-œ͔Ζ͕+œ͕Ζ͔)e01",
    ///     "+(+x͔Α͕-x͕Α͔-z͔Η͕+z͕Η͔-ð͔Θ͕+ð͕Θ͔-ø͔Ι͕+ø͕Ι͔-þ͔Κ͕+þ͕Κ͔-œ͔Λ͕+œ͕Λ͔)e02",
    ///     "+(+x͔Β͕-x͕Β͔+y͔Η͕-y͕Η͔-ð͔Μ͕+ð͕Μ͔-ø͔Ν͕+ø͕Ν͔-þ͔Ξ͕+þ͕Ξ͔-œ͔Ο͕+œ͕Ο͔)e03",
    ///     "+(+x͔Γ͕-x͕Γ͔+y͔Θ͕-y͕Θ͔+z͔Μ͕-z͕Μ͔-ø͔Π͕+ø͕Π͔-þ͔Ρ͕+þ͕Ρ͔-œ͔Σ͕+œ͕Σ͔)e04",
    ///     "+(+x͔Δ͕-x͕Δ͔+y͔Ι͕-y͕Ι͔+z͔Ν͕-z͕Ν͔+ð͔Π͕-ð͕Π͔-þ͔Τ͕+þ͕Τ͔-œ͔Υ͕+œ͕Υ͔)e05",
    ///     "+(+x͔Ε͕-x͕Ε͔+y͔Κ͕-y͕Κ͔+z͔Ξ͕-z͕Ξ͔+ð͔Ρ͕-ð͕Ρ͔+ø͔Τ͕-ø͕Τ͔-œ͔Φ͕+œ͕Φ͔)e06",
    ///     "+(+x͔Ζ͕-x͕Ζ͔+y͔Λ͕-y͕Λ͔+z͔Ο͕-z͕Ο͔+ð͔Σ͕-ð͕Σ͔+ø͔Υ͕-ø͕Υ͔+þ͔Φ͕-þ͕Φ͔)e07",
    ///     "+(-x͔y͕+x͕y͔)e12",
    ///     "+(-x͔z͕+x͕z͔)e13",
    ///     "+(-x͔ð͕+x͕ð͔)e14",
    ///     "+(-x͔ø͕+x͕ø͔)e15",
    ///     "+(-x͔þ͕+x͕þ͔)e16",
    ///     "+(-x͔œ͕+x͕œ͔)e17",
    ///     "+(-y͔z͕+y͕z͔)e23",
    ///     "+(-y͔ð͕+y͕ð͔)e24",
    ///     "+(-y͔ø͕+y͕ø͔)e25",
    ///     "+(-y͔þ͕+y͕þ͔)e26",
    ///     "+(-y͔œ͕+y͕œ͔)e27",
    ///     "+(-z͔ð͕+z͕ð͔)e34",
    ///     "+(-z͔ø͕+z͕ø͔)e35",
    ///     "+(-z͔þ͕+z͕þ͔)e36",
    ///     "+(-z͔œ͕+z͕œ͔)e37",
    ///     "+(-ð͔ø͕+ð͕ø͔)e45",
    ///     "+(-ð͔þ͕+ð͕þ͔)e46",
    ///     "+(-ð͔œ͕+ð͕œ͔)e47",
    ///     "+(-ø͔þ͕+ø͕þ͔)e56",
    ///     "+(-ø͔œ͕+ø͕œ͔)e57",
    ///     "+(-þ͔œ͕+þ͕œ͔)e67",
    ///     "+(+x͔Η͕+x͕Η͔-y͔Β͕-y͕Β͔+z͔Α͕+z͕Α͔)e0123",
    ///     "+(+x͔Θ͕+x͕Θ͔-y͔Γ͕-y͕Γ͔+ð͔Α͕+ð͕Α͔)e0124",
    ///     "+(+x͔Ι͕+x͕Ι͔-y͔Δ͕-y͕Δ͔+ø͔Α͕+ø͕Α͔)e0125",
    ///     "+(+x͔Κ͕+x͕Κ͔-y͔Ε͕-y͕Ε͔+þ͔Α͕+þ͕Α͔)e0126",
    ///     "+(+x͔Λ͕+x͕Λ͔-y͔Ζ͕-y͕Ζ͔+œ͔Α͕+œ͕Α͔)e0127",
    ///     "+(+x͔Μ͕+x͕Μ͔-z͔Γ͕-z͕Γ͔+ð͔Β͕+ð͕Β͔)e0134",
    ///     "+(+x͔Ν͕+x͕Ν͔-z͔Δ͕-z͕Δ͔+ø͔Β͕+ø͕Β͔)e0135",
    ///     "+(+x͔Ξ͕+x͕Ξ͔-z͔Ε͕-z͕Ε͔+þ͔Β͕+þ͕Β͔)e0136",
    ///     "+(+x͔Ο͕+x͕Ο͔-z͔Ζ͕-z͕Ζ͔+œ͔Β͕+œ͕Β͔)e0137",
    ///     "+(+x͔Π͕+x͕Π͔-ð͔Δ͕-ð͕Δ͔+ø͔Γ͕+ø͕Γ͔)e0145",
    ///     "+(+x͔Ρ͕+x͕Ρ͔-ð͔Ε͕-ð͕Ε͔+þ͔Γ͕+þ͕Γ͔)e0146",
    ///     "+(+x͔Σ͕+x͕Σ͔-ð͔Ζ͕-ð͕Ζ͔+œ͔Γ͕+œ͕Γ͔)e0147",
    ///     "+(+x͔Τ͕+x͕Τ͔-ø͔Ε͕-ø͕Ε͔+þ͔Δ͕+þ͕Δ͔)e0156",
    ///     "+(+x͔Υ͕+x͕Υ͔-ø͔Ζ͕-ø͕Ζ͔+œ͔Δ͕+œ͕Δ͔)e0157",
    ///     "+(+x͔Φ͕+x͕Φ͔-þ͔Ζ͕-þ͕Ζ͔+œ͔Ε͕+œ͕Ε͔)e0167",
    ///     "+(+y͔Μ͕+y͕Μ͔-z͔Θ͕-z͕Θ͔+ð͔Η͕+ð͕Η͔)e0234",
    ///     "+(+y͔Ν͕+y͕Ν͔-z͔Ι͕-z͕Ι͔+ø͔Η͕+ø͕Η͔)e0235",
    ///     "+(+y͔Ξ͕+y͕Ξ͔-z͔Κ͕-z͕Κ͔+þ͔Η͕+þ͕Η͔)e0236",
    ///     "+(+y͔Ο͕+y͕Ο͔-z͔Λ͕-z͕Λ͔+œ͔Η͕+œ͕Η͔)e0237",
    ///     "+(+y͔Π͕+y͕Π͔-ð͔Ι͕-ð͕Ι͔+ø͔Θ͕+ø͕Θ͔)e0245",
    ///     "+(+y͔Ρ͕+y͕Ρ͔-ð͔Κ͕-ð͕Κ͔+þ͔Θ͕+þ͕Θ͔)e0246",
    ///     "+(+y͔Σ͕+y͕Σ͔-ð͔Λ͕-ð͕Λ͔+œ͔Θ͕+œ͕Θ͔)e0247",
    ///     "+(+y͔Τ͕+y͕Τ͔-ø͔Κ͕-ø͕Κ͔+þ͔Ι͕+þ͕Ι͔)e0256",
    ///     "+(+y͔Υ͕+y͕Υ͔-ø͔Λ͕-ø͕Λ͔+œ͔Ι͕+œ͕Ι͔)e0257",
    ///     "+(+y͔Φ͕+y͕Φ͔-þ͔Λ͕-þ͕Λ͔+œ͔Κ͕+œ͕Κ͔)e0267",
    ///     "+(+z͔Π͕+z͕Π͔-ð͔Ν͕-ð͕Ν͔+ø͔Μ͕+ø͕Μ͔)e0345",
    ///     "+(+z͔Ρ͕+z͕Ρ͔-ð͔Ξ͕-ð͕Ξ͔+þ͔Μ͕+þ͕Μ͔)e0346",
    ///     "+(+z͔Σ͕+z͕Σ͔-ð͔Ο͕-ð͕Ο͔+œ͔Μ͕+œ͕Μ͔)e0347",
    ///     "+(+z͔Τ͕+z͕Τ͔-ø͔Ξ͕-ø͕Ξ͔+þ͔Ν͕+þ͕Ν͔)e0356",
    ///     "+(+z͔Υ͕+z͕Υ͔-ø͔Ο͕-ø͕Ο͔+œ͔Ν͕+œ͕Ν͔)e0357",
    ///     "+(+z͔Φ͕+z͕Φ͔-þ͔Ο͕-þ͕Ο͔+œ͔Ξ͕+œ͕Ξ͔)e0367",
    ///     "+(+ð͔Τ͕+ð͕Τ͔-ø͔Ρ͕-ø͕Ρ͔+þ͔Π͕+þ͕Π͔)e0456",
    ///     "+(+ð͔Υ͕+ð͕Υ͔-ø͔Σ͕-ø͕Σ͔+œ͔Π͕+œ͕Π͔)e0457",
    ///     "+(+ð͔Φ͕+ð͕Φ͔-þ͔Σ͕-þ͕Σ͔+œ͔Ρ͕+œ͕Ρ͔)e0467",
    ///     "+(+ø͔Φ͕+ø͕Φ͔-þ͔Υ͕-þ͕Υ͔+œ͔Τ͕+œ͕Τ͔)e0567",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::volume5() + Self::volume_moment()
    }
    /// The multivector of double motor $`m_2 \equiv s + v^5 + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔+Œ͕ζ͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕η͔+Ð͕θ͔+Ø͕ι͔+Þ͕κ͔+Œ͕λ͔)e02",
    ///     "+(-X͕β͔-Y͕η͔+Z͕v͔+Ð͕μ͔+Ø͕ν͔+Þ͕ξ͔+Œ͕ο͔)e03",
    ///     "+(-X͕γ͔-Y͕θ͔-Z͕μ͔+v͔Ð͕+Ø͕π͔+Þ͕ρ͔+Œ͕σ͔)e04",
    ///     "+(-X͕δ͔-Y͕ι͔-Z͕ν͔+v͔Ø͕-Ð͕π͔+Þ͕τ͔+Œ͕υ͔)e05",
    ///     "+(-X͕ε͔-Y͕κ͔-Z͕ξ͔+v͔Þ͕-Ð͕ρ͔-Ø͕τ͔+Œ͕φ͔)e06",
    ///     "+(-X͕ζ͔-Y͕λ͔-Z͕ο͔+v͔Œ͕-Ð͕σ͔-Ø͕υ͔-Þ͕φ͔)e07",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e17",
    ///     "+v͕η͔e23",
    ///     "+v͕θ͔e24",
    ///     "+v͕ι͔e25",
    ///     "+v͕κ͔e26",
    ///     "+v͕λ͔e27",
    ///     "+v͕μ͔e34",
    ///     "+v͕ν͔e35",
    ///     "+v͕ξ͔e36",
    ///     "+v͕ο͔e37",
    ///     "+v͕π͔e45",
    ///     "+v͕ρ͔e46",
    ///     "+v͕σ͔e47",
    ///     "+v͕τ͔e56",
    ///     "+v͕υ͔e57",
    ///     "+v͕φ͔e67",
    ///     "+v͕ü͔e1234",
    ///     "+(+X͕η͔-Y͕β͔+Z͕α͔+Ð͕ü͔-Ø͕ú͔+Þ͕ö͔-ó͔Œ͕)e0123",
    ///     "+(+X͕θ͔-Y͕γ͔-Z͕ü͔+Ð͕α͔+Ø͕ñ͔-Þ͕ï͔+í͔Œ͕)e0124",
    ///     "+(+X͕ι͔-Y͕δ͔+Z͕ú͔-Ð͕ñ͔+Ø͕α͔+Þ͕ë͔-é͔Œ͕)e0125",
    ///     "+(+X͕κ͔-Y͕ε͔-Z͕ö͔+Ð͕ï͔-Ø͕ë͔+Þ͕α͔+ç͔Œ͕)e0126",
    ///     "+(+X͕λ͔-Y͕ζ͔+Z͕ó͔-Ð͕í͔+Ø͕é͔-Þ͕ç͔+Œ͕α͔)e0127",
    ///     "+(+X͕μ͔+Y͕ü͔-Z͕γ͔+Ð͕β͔-Ø͕æ͔+Þ͕å͔-ä͔Œ͕)e0134",
    ///     "+(+X͕ν͔-Y͕ú͔-Z͕δ͔+u͔Œ͕+Ð͕æ͔+Ø͕β͔-Þ͕á͔)e0135",
    ///     "+(+X͕ξ͔+Y͕ö͔-Z͕ε͔-t͔Œ͕-Ð͕å͔+Ø͕á͔+Þ͕β͔)e0136",
    ///     "+(+X͕ο͔-Y͕ó͔-Z͕ζ͔+t͔Þ͕-u͔Ø͕+Ð͕ä͔+Œ͕β͔)e0137",
    ///     "+(+X͕π͔+Y͕ñ͔-Z͕æ͔-r͔Œ͕+s͔Þ͕-Ð͕δ͔+Ø͕γ͔)e0145",
    ///     "+(+X͕ρ͔-Y͕ï͔+Z͕å͔+q͔Œ͕-s͔Ø͕-Ð͕ε͔+Þ͕γ͔)e0146",
    ///     "+(+X͕σ͔+Y͕í͔-Z͕ä͔-q͔Þ͕+r͔Ø͕-Ð͕ζ͔+Œ͕γ͔)e0147",
    ///     "+(+X͕τ͔+Y͕ë͔-Z͕á͔-p͔Œ͕+s͔Ð͕-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(+X͕υ͔-Y͕é͔+Z͕u͔+p͔Þ͕-r͔Ð͕-Ø͕ζ͔+Œ͕δ͔)e0157",
    ///     "+(+X͕φ͔+Y͕ç͔-Z͕t͔-p͔Ø͕+q͔Ð͕-Þ͕ζ͔+Œ͕ε͔)e0167",
    ///     "+(-X͕ü͔+Y͕μ͔-Z͕θ͔+m͔Œ͕-n͔Þ͕+o͔Ø͕+Ð͕η͔)e0234",
    ///     "+(+X͕ú͔+Y͕ν͔-Z͕ι͔-k͔Œ͕+l͔Þ͕-o͔Ð͕+Ø͕η͔)e0235",
    ///     "+(-X͕ö͔+Y͕ξ͔-Z͕κ͔+j͔Œ͕-l͔Ø͕+n͔Ð͕+Þ͕η͔)e0236",
    ///     "+(+X͕ó͔+Y͕ο͔-Z͕λ͔-j͔Þ͕+k͔Ø͕-m͔Ð͕+Œ͕η͔)e0237",
    ///     "+(-X͕ñ͔+Y͕π͔+Z͕o͔+h͔Œ͕-i͔Þ͕-Ð͕ι͔+Ø͕θ͔)e0245",
    ///     "+(+X͕ï͔+Y͕ρ͔-Z͕n͔-g͔Œ͕+i͔Ø͕-Ð͕κ͔+Þ͕θ͔)e0246",
    ///     "+(-X͕í͔+Y͕σ͔+Z͕m͔+g͔Þ͕-h͔Ø͕-Ð͕λ͔+Œ͕θ͔)e0247",
    ///     "+(-X͕ë͔+Y͕τ͔+Z͕l͔+f͔Œ͕-i͔Ð͕-Ø͕κ͔+Þ͕ι͔)e0256",
    ///     "+(+X͕é͔+Y͕υ͔-Z͕k͔-f͔Þ͕+h͔Ð͕-Ø͕λ͔+Œ͕ι͔)e0257",
    ///     "+(-X͕ç͔+Y͕φ͔+Z͕j͔+f͔Ø͕-g͔Ð͕-Þ͕λ͔+Œ͕κ͔)e0267",
    ///     "+(+X͕æ͔-Y͕o͔+Z͕π͔-d͔Œ͕+e͔Þ͕-Ð͕ν͔+Ø͕μ͔)e0345",
    ///     "+(-X͕å͔+Y͕n͔+Z͕ρ͔+c͔Œ͕-e͔Ø͕-Ð͕ξ͔+Þ͕μ͔)e0346",
    ///     "+(+X͕ä͔-Y͕m͔+Z͕σ͔-c͔Þ͕+d͔Ø͕-Ð͕ο͔+Œ͕μ͔)e0347",
    ///     "+(+X͕á͔-Y͕l͔+Z͕τ͔-b͔Œ͕+e͔Ð͕-Ø͕ξ͔+Þ͕ν͔)e0356",
    ///     "+(-X͕u͔+Y͕k͔+Z͕υ͔+b͔Þ͕-d͔Ð͕-Ø͕ο͔+Œ͕ν͔)e0357",
    ///     "+(+X͕t͔-Y͕j͔+Z͕φ͔-b͔Ø͕+c͔Ð͕-Þ͕ο͔+Œ͕ξ͔)e0367",
    ///     "+(-X͕s͔+Y͕i͔-Z͕e͔+a͔Œ͕+Ð͕τ͔-Ø͕ρ͔+Þ͕π͔)e0456",
    ///     "+(+X͕r͔-Y͕h͔+Z͕d͔-a͔Þ͕+Ð͕υ͔-Ø͕σ͔+Œ͕π͔)e0457",
    ///     "+(-X͕q͔+Y͕g͔-Z͕c͔+a͔Ø͕+Ð͕φ͔-Þ͕σ͔+Œ͕ρ͔)e0467",
    ///     "+(+X͕p͔-Y͕f͔+Z͕b͔-a͔Ð͕+Ø͕φ͔-Þ͕υ͔+Œ͕τ͔)e0567",
    ///     "+a͔v͕e4567",
    ///     "+b͔v͕e3576",
    ///     "+c͔v͕e3467",
    ///     "+d͔v͕e3475",
    ///     "+e͔v͕e3456",
    ///     "+f͔v͕e2567",
    ///     "+g͔v͕e2476",
    ///     "+h͔v͕e2457",
    ///     "+i͔v͕e2465",
    ///     "+j͔v͕e2367",
    ///     "+k͔v͕e2375",
    ///     "+l͔v͕e2356",
    ///     "+m͔v͕e2347",
    ///     "+n͔v͕e2364",
    ///     "+o͔v͕e2345",
    ///     "+p͔v͕e1576",
    ///     "+q͔v͕e1467",
    ///     "+r͔v͕e1475",
    ///     "+s͔v͕e1456",
    ///     "+t͔v͕e1376",
    ///     "+u͔v͕e1357",
    ///     "+v͕á͔e1365",
    ///     "+v͕ä͔e1374",
    ///     "+v͕å͔e1346",
    ///     "+v͕æ͔e1354",
    ///     "+v͕ç͔e1267",
    ///     "+v͕é͔e1275",
    ///     "+v͕ë͔e1256",
    ///     "+v͕í͔e1247",
    ///     "+v͕ï͔e1264",
    ///     "+v͕ñ͔e1245",
    ///     "+v͕ó͔e1273",
    ///     "+v͕ö͔e1236",
    ///     "+v͕ú͔e1253",
    ///     "+(+Z͕a͔+b͔Ð͕+c͔Ø͕+d͔Þ͕+e͔Œ͕)e034567",
    ///     "+(-Y͕a͔+f͔Ð͕+g͔Ø͕+h͔Þ͕+i͔Œ͕)e024576",
    ///     "+(-Y͕b͔-Z͕f͔+j͔Ø͕+k͔Þ͕+l͔Œ͕)e023567",
    ///     "+(-Y͕c͔-Z͕g͔-j͔Ð͕+m͔Þ͕+n͔Œ͕)e023476",
    ///     "+(-Y͕d͔-Z͕h͔-k͔Ð͕-m͔Ø͕+o͔Œ͕)e023457",
    ///     "+(-Y͕e͔-Z͕i͔-l͔Ð͕-n͔Ø͕-o͔Þ͕)e023465",
    ///     "+(+X͕a͔+p͔Ð͕+q͔Ø͕+r͔Þ͕+s͔Œ͕)e014567",
    ///     "+(+X͕b͔-Z͕p͔+t͔Ø͕+u͔Þ͕+á͔Œ͕)e013576",
    ///     "+(+X͕c͔-Z͕q͔-t͔Ð͕+Þ͕ä͔+å͔Œ͕)e013467",
    ///     "+(+X͕d͔-Z͕r͔-u͔Ð͕-Ø͕ä͔+æ͔Œ͕)e013475",
    ///     "+(+X͕e͔-Z͕s͔-Ð͕á͔-Ø͕å͔-Þ͕æ͔)e013456",
    ///     "+(+X͕f͔+Y͕p͔+Ø͕ç͔+Þ͕é͔+ë͔Œ͕)e012567",
    ///     "+(+X͕g͔+Y͕q͔-Ð͕ç͔+Þ͕í͔+ï͔Œ͕)e012476",
    ///     "+(+X͕h͔+Y͕r͔-Ð͕é͔-Ø͕í͔+ñ͔Œ͕)e012457",
    ///     "+(+X͕i͔+Y͕s͔-Ð͕ë͔-Ø͕ï͔-Þ͕ñ͔)e012465",
    ///     "+(+X͕j͔+Y͕t͔+Z͕ç͔+Þ͕ó͔+ö͔Œ͕)e012367",
    ///     "+(+X͕k͔+Y͕u͔+Z͕é͔-Ø͕ó͔+ú͔Œ͕)e012375",
    ///     "+(+X͕l͔+Y͕á͔+Z͕ë͔-Ø͕ö͔-Þ͕ú͔)e012356",
    ///     "+(+X͕m͔+Y͕ä͔+Z͕í͔+Ð͕ó͔+ü͔Œ͕)e012347",
    ///     "+(+X͕n͔+Y͕å͔+Z͕ï͔+Ð͕ö͔-Þ͕ü͔)e012364",
    ///     "+(+X͕o͔+Y͕æ͔+Z͕ñ͔+Ð͕ú͔+Ø͕ü͔)e012345",
    /// ]);
    ///
    /// let simple_double_motor = Vee::volume5().lhs() * Vee::volume5().rhs();
    ///
    /// assert_eq!(simple_double_motor.basis_blades(), Vee::simple_double_motor().basis_blades());
    /// format_eq!(simple_double_motor, [
    ///     "-α͔α͕-β͔β͕-γ͔γ͕-δ͔δ͕-ε͔ε͕-ζ͔ζ͕-η͔η͕-θ͔θ͕-ι͔ι͕-κ͔κ͕-λ͔λ͕-μ͔μ͕-ν͔ν͕-ξ͔ξ͕-ο͔ο͕-π͔π͕-ρ͔ρ͕-σ͔σ͕-τ͔τ͕-υ͔υ͕-φ͔φ͕",
    ///     "+(-Y͔α͕+Y͕α͔-Z͔β͕+Z͕β͔-Ð͔γ͕+Ð͕γ͔-Ø͔δ͕+Ø͕δ͔-Þ͔ε͕+Þ͕ε͔-Œ͔ζ͕+Œ͕ζ͔)e01",
    ///     "+(+X͔α͕-X͕α͔-Z͔η͕+Z͕η͔-Ð͔θ͕+Ð͕θ͔-Ø͔ι͕+Ø͕ι͔-Þ͔κ͕+Þ͕κ͔-Œ͔λ͕+Œ͕λ͔)e02",
    ///     "+(+X͔β͕-X͕β͔+Y͔η͕-Y͕η͔-Ð͔μ͕+Ð͕μ͔-Ø͔ν͕+Ø͕ν͔-Þ͔ξ͕+Þ͕ξ͔-Œ͔ο͕+Œ͕ο͔)e03",
    ///     "+(+X͔γ͕-X͕γ͔+Y͔θ͕-Y͕θ͔+Z͔μ͕-Z͕μ͔-Ø͔π͕+Ø͕π͔-Þ͔ρ͕+Þ͕ρ͔-Œ͔σ͕+Œ͕σ͔)e04",
    ///     "+(+X͔δ͕-X͕δ͔+Y͔ι͕-Y͕ι͔+Z͔ν͕-Z͕ν͔+Ð͔π͕-Ð͕π͔-Þ͔τ͕+Þ͕τ͔-Œ͔υ͕+Œ͕υ͔)e05",
    ///     "+(+X͔ε͕-X͕ε͔+Y͔κ͕-Y͕κ͔+Z͔ξ͕-Z͕ξ͔+Ð͔ρ͕-Ð͕ρ͔+Ø͔τ͕-Ø͕τ͔-Œ͔φ͕+Œ͕φ͔)e06",
    ///     "+(+X͔ζ͕-X͕ζ͔+Y͔λ͕-Y͕λ͔+Z͔ο͕-Z͕ο͔+Ð͔σ͕-Ð͕σ͔+Ø͔υ͕-Ø͕υ͔+Þ͔φ͕-Þ͕φ͔)e07",
    ///     "+(-β͔η͕+β͕η͔-γ͔θ͕+γ͕θ͔-δ͔ι͕+δ͕ι͔-ε͔κ͕+ε͕κ͔-ζ͔λ͕+ζ͕λ͔)e12",
    ///     "+(+α͔η͕-α͕η͔-γ͔μ͕+γ͕μ͔-δ͔ν͕+δ͕ν͔-ε͔ξ͕+ε͕ξ͔-ζ͔ο͕+ζ͕ο͔)e13",
    ///     "+(+α͔θ͕-α͕θ͔+β͔μ͕-β͕μ͔-δ͔π͕+δ͕π͔-ε͔ρ͕+ε͕ρ͔-ζ͔σ͕+ζ͕σ͔)e14",
    ///     "+(+α͔ι͕-α͕ι͔+β͔ν͕-β͕ν͔+γ͔π͕-γ͕π͔-ε͔τ͕+ε͕τ͔-ζ͔υ͕+ζ͕υ͔)e15",
    ///     "+(+α͔κ͕-α͕κ͔+β͔ξ͕-β͕ξ͔+γ͔ρ͕-γ͕ρ͔+δ͔τ͕-δ͕τ͔-ζ͔φ͕+ζ͕φ͔)e16",
    ///     "+(+α͔λ͕-α͕λ͔+β͔ο͕-β͕ο͔+γ͔σ͕-γ͕σ͔+δ͔υ͕-δ͕υ͔+ε͔φ͕-ε͕φ͔)e17",
    ///     "+(-α͔β͕+α͕β͔-θ͔μ͕+θ͕μ͔-ι͔ν͕+ι͕ν͔-κ͔ξ͕+κ͕ξ͔-λ͔ο͕+λ͕ο͔)e23",
    ///     "+(-α͔γ͕+α͕γ͔+η͔μ͕-η͕μ͔-ι͔π͕+ι͕π͔-κ͔ρ͕+κ͕ρ͔-λ͔σ͕+λ͕σ͔)e24",
    ///     "+(-α͔δ͕+α͕δ͔+η͔ν͕-η͕ν͔+θ͔π͕-θ͕π͔-κ͔τ͕+κ͕τ͔-λ͔υ͕+λ͕υ͔)e25",
    ///     "+(-α͔ε͕+α͕ε͔+η͔ξ͕-η͕ξ͔+θ͔ρ͕-θ͕ρ͔+ι͔τ͕-ι͕τ͔-λ͔φ͕+λ͕φ͔)e26",
    ///     "+(-α͔ζ͕+α͕ζ͔+η͔ο͕-η͕ο͔+θ͔σ͕-θ͕σ͔+ι͔υ͕-ι͕υ͔+κ͔φ͕-κ͕φ͔)e27",
    ///     "+(-β͔γ͕+β͕γ͔-η͔θ͕+η͕θ͔-ν͔π͕+ν͕π͔-ξ͔ρ͕+ξ͕ρ͔-ο͔σ͕+ο͕σ͔)e34",
    ///     "+(-β͔δ͕+β͕δ͔-η͔ι͕+η͕ι͔+μ͔π͕-μ͕π͔-ξ͔τ͕+ξ͕τ͔-ο͔υ͕+ο͕υ͔)e35",
    ///     "+(-β͔ε͕+β͕ε͔-η͔κ͕+η͕κ͔+μ͔ρ͕-μ͕ρ͔+ν͔τ͕-ν͕τ͔-ο͔φ͕+ο͕φ͔)e36",
    ///     "+(-β͔ζ͕+β͕ζ͔-η͔λ͕+η͕λ͔+μ͔σ͕-μ͕σ͔+ν͔υ͕-ν͕υ͔+ξ͔φ͕-ξ͕φ͔)e37",
    ///     "+(-γ͔δ͕+γ͕δ͔-θ͔ι͕+θ͕ι͔-μ͔ν͕+μ͕ν͔-ρ͔τ͕+ρ͕τ͔-σ͔υ͕+σ͕υ͔)e45",
    ///     "+(-γ͔ε͕+γ͕ε͔-θ͔κ͕+θ͕κ͔-μ͔ξ͕+μ͕ξ͔+π͔τ͕-π͕τ͔-σ͔φ͕+σ͕φ͔)e46",
    ///     "+(-γ͔ζ͕+γ͕ζ͔-θ͔λ͕+θ͕λ͔-μ͔ο͕+μ͕ο͔+π͔υ͕-π͕υ͔+ρ͔φ͕-ρ͕φ͔)e47",
    ///     "+(-δ͔ε͕+δ͕ε͔-ι͔κ͕+ι͕κ͔-ν͔ξ͕+ν͕ξ͔-π͔ρ͕+π͕ρ͔-υ͔φ͕+υ͕φ͔)e56",
    ///     "+(-δ͔ζ͕+δ͕ζ͔-ι͔λ͕+ι͕λ͔-ν͔ο͕+ν͕ο͔-π͔σ͕+π͕σ͔+τ͔φ͕-τ͕φ͔)e57",
    ///     "+(-ε͔ζ͕+ε͕ζ͔-κ͔λ͕+κ͕λ͔-ξ͔ο͕+ξ͕ο͔-ρ͔σ͕+ρ͕σ͔-τ͔υ͕+τ͕υ͔)e67",
    ///     "+(+α͔μ͕+α͕μ͔-β͔θ͕-β͕θ͔+γ͔η͕+γ͕η͔)e1234",
    ///     "+(+X͔η͕+X͕η͔-Y͔β͕-Y͕β͔+Z͔α͕+Z͕α͔)e0123",
    ///     "+(+X͔θ͕+X͕θ͔-Y͔γ͕-Y͕γ͔+Ð͔α͕+Ð͕α͔)e0124",
    ///     "+(+X͔ι͕+X͕ι͔-Y͔δ͕-Y͕δ͔+Ø͔α͕+Ø͕α͔)e0125",
    ///     "+(+X͔κ͕+X͕κ͔-Y͔ε͕-Y͕ε͔+Þ͔α͕+Þ͕α͔)e0126",
    ///     "+(+X͔λ͕+X͕λ͔-Y͔ζ͕-Y͕ζ͔+Œ͔α͕+Œ͕α͔)e0127",
    ///     "+(+X͔μ͕+X͕μ͔-Z͔γ͕-Z͕γ͔+Ð͔β͕+Ð͕β͔)e0134",
    ///     "+(+X͔ν͕+X͕ν͔-Z͔δ͕-Z͕δ͔+Ø͔β͕+Ø͕β͔)e0135",
    ///     "+(+X͔ξ͕+X͕ξ͔-Z͔ε͕-Z͕ε͔+Þ͔β͕+Þ͕β͔)e0136",
    ///     "+(+X͔ο͕+X͕ο͔-Z͔ζ͕-Z͕ζ͔+Œ͔β͕+Œ͕β͔)e0137",
    ///     "+(+X͔π͕+X͕π͔-Ð͔δ͕-Ð͕δ͔+Ø͔γ͕+Ø͕γ͔)e0145",
    ///     "+(+X͔ρ͕+X͕ρ͔-Ð͔ε͕-Ð͕ε͔+Þ͔γ͕+Þ͕γ͔)e0146",
    ///     "+(+X͔σ͕+X͕σ͔-Ð͔ζ͕-Ð͕ζ͔+Œ͔γ͕+Œ͕γ͔)e0147",
    ///     "+(+X͔τ͕+X͕τ͔-Ø͔ε͕-Ø͕ε͔+Þ͔δ͕+Þ͕δ͔)e0156",
    ///     "+(+X͔υ͕+X͕υ͔-Ø͔ζ͕-Ø͕ζ͔+Œ͔δ͕+Œ͕δ͔)e0157",
    ///     "+(+X͔φ͕+X͕φ͔-Þ͔ζ͕-Þ͕ζ͔+Œ͔ε͕+Œ͕ε͔)e0167",
    ///     "+(+Y͔μ͕+Y͕μ͔-Z͔θ͕-Z͕θ͔+Ð͔η͕+Ð͕η͔)e0234",
    ///     "+(+Y͔ν͕+Y͕ν͔-Z͔ι͕-Z͕ι͔+Ø͔η͕+Ø͕η͔)e0235",
    ///     "+(+Y͔ξ͕+Y͕ξ͔-Z͔κ͕-Z͕κ͔+Þ͔η͕+Þ͕η͔)e0236",
    ///     "+(+Y͔ο͕+Y͕ο͔-Z͔λ͕-Z͕λ͔+Œ͔η͕+Œ͕η͔)e0237",
    ///     "+(+Y͔π͕+Y͕π͔-Ð͔ι͕-Ð͕ι͔+Ø͔θ͕+Ø͕θ͔)e0245",
    ///     "+(+Y͔ρ͕+Y͕ρ͔-Ð͔κ͕-Ð͕κ͔+Þ͔θ͕+Þ͕θ͔)e0246",
    ///     "+(+Y͔σ͕+Y͕σ͔-Ð͔λ͕-Ð͕λ͔+Œ͔θ͕+Œ͕θ͔)e0247",
    ///     "+(+Y͔τ͕+Y͕τ͔-Ø͔κ͕-Ø͕κ͔+Þ͔ι͕+Þ͕ι͔)e0256",
    ///     "+(+Y͔υ͕+Y͕υ͔-Ø͔λ͕-Ø͕λ͔+Œ͔ι͕+Œ͕ι͔)e0257",
    ///     "+(+Y͔φ͕+Y͕φ͔-Þ͔λ͕-Þ͕λ͔+Œ͔κ͕+Œ͕κ͔)e0267",
    ///     "+(+Z͔π͕+Z͕π͔-Ð͔ν͕-Ð͕ν͔+Ø͔μ͕+Ø͕μ͔)e0345",
    ///     "+(+Z͔ρ͕+Z͕ρ͔-Ð͔ξ͕-Ð͕ξ͔+Þ͔μ͕+Þ͕μ͔)e0346",
    ///     "+(+Z͔σ͕+Z͕σ͔-Ð͔ο͕-Ð͕ο͔+Œ͔μ͕+Œ͕μ͔)e0347",
    ///     "+(+Z͔τ͕+Z͕τ͔-Ø͔ξ͕-Ø͕ξ͔+Þ͔ν͕+Þ͕ν͔)e0356",
    ///     "+(+Z͔υ͕+Z͕υ͔-Ø͔ο͕-Ø͕ο͔+Œ͔ν͕+Œ͕ν͔)e0357",
    ///     "+(+Z͔φ͕+Z͕φ͔-Þ͔ο͕-Þ͕ο͔+Œ͔ξ͕+Œ͕ξ͔)e0367",
    ///     "+(+Ð͔τ͕+Ð͕τ͔-Ø͔ρ͕-Ø͕ρ͔+Þ͔π͕+Þ͕π͔)e0456",
    ///     "+(+Ð͔υ͕+Ð͕υ͔-Ø͔σ͕-Ø͕σ͔+Œ͔π͕+Œ͕π͔)e0457",
    ///     "+(+Ð͔φ͕+Ð͕φ͔-Þ͔σ͕-Þ͕σ͔+Œ͔ρ͕+Œ͕ρ͔)e0467",
    ///     "+(+Ø͔φ͕+Ø͕φ͔-Þ͔υ͕-Þ͕υ͔+Œ͔τ͕+Œ͕τ͔)e0567",
    ///     "+(+π͔φ͕+π͕φ͔-ρ͔υ͕-ρ͕υ͔+σ͔τ͕+σ͕τ͔)e4567",
    ///     "+(-ν͔φ͕-ν͕φ͔+ξ͔υ͕+ξ͕υ͔-ο͔τ͕-ο͕τ͔)e3576",
    ///     "+(+μ͔φ͕+μ͕φ͔-ξ͔σ͕-ξ͕σ͔+ο͔ρ͕+ο͕ρ͔)e3467",
    ///     "+(-μ͔υ͕-μ͕υ͔+ν͔σ͕+ν͕σ͔-ο͔π͕-ο͕π͔)e3475",
    ///     "+(+μ͔τ͕+μ͕τ͔-ν͔ρ͕-ν͕ρ͔+ξ͔π͕+ξ͕π͔)e3456",
    ///     "+(+ι͔φ͕+ι͕φ͔-κ͔υ͕-κ͕υ͔+λ͔τ͕+λ͕τ͔)e2567",
    ///     "+(-θ͔φ͕-θ͕φ͔+κ͔σ͕+κ͕σ͔-λ͔ρ͕-λ͕ρ͔)e2476",
    ///     "+(+θ͔υ͕+θ͕υ͔-ι͔σ͕-ι͕σ͔+λ͔π͕+λ͕π͔)e2457",
    ///     "+(-θ͔τ͕-θ͕τ͔+ι͔ρ͕+ι͕ρ͔-κ͔π͕-κ͕π͔)e2465",
    ///     "+(+η͔φ͕+η͕φ͔-κ͔ο͕-κ͕ο͔+λ͔ξ͕+λ͕ξ͔)e2367",
    ///     "+(-η͔υ͕-η͕υ͔+ι͔ο͕+ι͕ο͔-λ͔ν͕-λ͕ν͔)e2375",
    ///     "+(+η͔τ͕+η͕τ͔-ι͔ξ͕-ι͕ξ͔+κ͔ν͕+κ͕ν͔)e2356",
    ///     "+(+η͔σ͕+η͕σ͔-θ͔ο͕-θ͕ο͔+λ͔μ͕+λ͕μ͔)e2347",
    ///     "+(-η͔ρ͕-η͕ρ͔+θ͔ξ͕+θ͕ξ͔-κ͔μ͕-κ͕μ͔)e2364",
    ///     "+(+η͔π͕+η͕π͔-θ͔ν͕-θ͕ν͔+ι͔μ͕+ι͕μ͔)e2345",
    ///     "+(-δ͔φ͕-δ͕φ͔+ε͔υ͕+ε͕υ͔-ζ͔τ͕-ζ͕τ͔)e1576",
    ///     "+(+γ͔φ͕+γ͕φ͔-ε͔σ͕-ε͕σ͔+ζ͔ρ͕+ζ͕ρ͔)e1467",
    ///     "+(-γ͔υ͕-γ͕υ͔+δ͔σ͕+δ͕σ͔-ζ͔π͕-ζ͕π͔)e1475",
    ///     "+(+γ͔τ͕+γ͕τ͔-δ͔ρ͕-δ͕ρ͔+ε͔π͕+ε͕π͔)e1456",
    ///     "+(-β͔φ͕-β͕φ͔+ε͔ο͕+ε͕ο͔-ζ͔ξ͕-ζ͕ξ͔)e1376",
    ///     "+(+β͔υ͕+β͕υ͔-δ͔ο͕-δ͕ο͔+ζ͔ν͕+ζ͕ν͔)e1357",
    ///     "+(-β͔τ͕-β͕τ͔+δ͔ξ͕+δ͕ξ͔-ε͔ν͕-ε͕ν͔)e1365",
    ///     "+(-β͔σ͕-β͕σ͔+γ͔ο͕+γ͕ο͔-ζ͔μ͕-ζ͕μ͔)e1374",
    ///     "+(+β͔ρ͕+β͕ρ͔-γ͔ξ͕-γ͕ξ͔+ε͔μ͕+ε͕μ͔)e1346",
    ///     "+(-β͔π͕-β͕π͔+γ͔ν͕+γ͕ν͔-δ͔μ͕-δ͕μ͔)e1354",
    ///     "+(+α͔φ͕+α͕φ͔-ε͔λ͕-ε͕λ͔+ζ͔κ͕+ζ͕κ͔)e1267",
    ///     "+(-α͔υ͕-α͕υ͔+δ͔λ͕+δ͕λ͔-ζ͔ι͕-ζ͕ι͔)e1275",
    ///     "+(+α͔τ͕+α͕τ͔-δ͔κ͕-δ͕κ͔+ε͔ι͕+ε͕ι͔)e1256",
    ///     "+(+α͔σ͕+α͕σ͔-γ͔λ͕-γ͕λ͔+ζ͔θ͕+ζ͕θ͔)e1247",
    ///     "+(-α͔ρ͕-α͕ρ͔+γ͔κ͕+γ͕κ͔-ε͔θ͕-ε͕θ͔)e1264",
    ///     "+(+α͔π͕+α͕π͔-γ͔ι͕-γ͕ι͔+δ͔θ͕+δ͕θ͔)e1245",
    ///     "+(-α͔ο͕-α͕ο͔+β͔λ͕+β͕λ͔-ζ͔η͕-ζ͕η͔)e1273",
    ///     "+(+α͔ξ͕+α͕ξ͔-β͔κ͕-β͕κ͔+ε͔η͕+ε͕η͔)e1236",
    ///     "+(-α͔ν͕-α͕ν͔+β͔ι͕+β͕ι͔-δ͔η͕-δ͕η͔)e1253",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_motor() -> Self {
        Self::scalar() + Self::volume5() + Self::volume()
    }
    /// The multivector of double motor $`m_2 \equiv s + v^5 + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// format_eq!(double_motor, [
    ///     "+α͔α͕+β͔β͕+γ͔γ͕+δ͔δ͕+ε͔ε͕+ζ͔ζ͕+η͔η͕+θ͔θ͕+ι͔ι͕+κ͔κ͕+λ͔λ͕+μ͔μ͕+ν͔ν͕+ξ͔ξ͕+ο͔ο͕+π͔π͕+ρ͔ρ͕+σ͔σ͕+τ͔τ͕+υ͔υ͕+φ͔φ͕",
    ///     "+(-A͔η͕+A͕η͔-B͔θ͕+B͕θ͔-C͔ι͕+C͕ι͔-D͔κ͕+D͕κ͔-E͔λ͕+E͕λ͔-F͔μ͕+F͕μ͔-G͔ν͕+G͕ν͔-H͔ξ͕\
    ///        +H͕ξ͔-I͔ο͕+I͕ο͔-J͔π͕+J͕π͔-K͔ρ͕+K͕ρ͔-L͔σ͕+L͕σ͔-M͔τ͕+M͕τ͔-N͔υ͕+N͕υ͔-O͔φ͕+O͕φ͔)e01",
    ///     "+(+A͔β͕-A͕β͔+B͔γ͕-B͕γ͔+C͔δ͕-C͕δ͔+D͔ε͕-D͕ε͔+E͔ζ͕-E͕ζ͔-P͔μ͕+P͕μ͔-Q͔ν͕+Q͕ν͔-R͔ξ͕\
    ///        +R͕ξ͔-S͔ο͕+S͕ο͔-T͔π͕+T͕π͔-U͔ρ͕+U͕ρ͔-Á͔σ͕+Á͕σ͔-Ä͔τ͕+Ä͕τ͔-Å͔υ͕+Å͕υ͔-Æ͔φ͕+Æ͕φ͔)e02",
    ///     "+(-A͔α͕+A͕α͔+F͔γ͕-F͕γ͔+G͔δ͕-G͕δ͔+H͔ε͕-H͕ε͔+I͔ζ͕-I͕ζ͔+P͔θ͕-P͕θ͔+Q͔ι͕-Q͕ι͔+R͔κ͕\
    ///        -R͕κ͔+S͔λ͕-S͕λ͔-Ç͔π͕+Ç͕π͔-É͔ρ͕+É͕ρ͔-Ë͔σ͕+Ë͕σ͔-Í͔τ͕+Í͕τ͔-Ï͔υ͕+Ï͕υ͔-Ñ͔φ͕+Ñ͕φ͔)e03",
    ///     "+(-B͔α͕+B͕α͔-F͔β͕+F͕β͔+J͔δ͕-J͕δ͔+K͔ε͕-K͕ε͔+L͔ζ͕-L͕ζ͔-P͔η͕+P͕η͔+T͔ι͕-T͕ι͔+U͔κ͕\
    ///        -U͕κ͔+Á͔λ͕-Á͕λ͔+Ç͔ν͕-Ç͕ν͔+É͔ξ͕-É͕ξ͔+Ë͔ο͕-Ë͕ο͔-Ó͔τ͕+Ó͕τ͔-Ö͔υ͕+Ö͕υ͔-Ú͔φ͕+Ú͕φ͔)e04",
    ///     "+(-C͔α͕+C͕α͔-G͔β͕+G͕β͔-J͔γ͕+J͕γ͔+M͔ε͕-M͕ε͔+N͔ζ͕-N͕ζ͔-Q͔η͕+Q͕η͔-T͔θ͕+T͕θ͔+Ä͔κ͕\
    ///        -Ä͕κ͔+Å͔λ͕-Å͕λ͔-Ç͔μ͕+Ç͕μ͔+Í͔ξ͕-Í͕ξ͔+Ï͔ο͕-Ï͕ο͔+Ó͔ρ͕-Ó͕ρ͔+Ö͔σ͕-Ö͕σ͔-Ü͔φ͕+Ü͕φ͔)e05",
    ///     "+(-D͔α͕+D͕α͔-H͔β͕+H͕β͔-K͔γ͕+K͕γ͔-M͔δ͕+M͕δ͔+O͔ζ͕-O͕ζ͔-R͔η͕+R͕η͔-U͔θ͕+U͕θ͔-Ä͔ι͕\
    ///        +Ä͕ι͔+Æ͔λ͕-Æ͕λ͔-É͔μ͕+É͕μ͔-Í͔ν͕+Í͕ν͔+Ñ͔ο͕-Ñ͕ο͔-Ó͔π͕+Ó͕π͔+Ú͔σ͕-Ú͕σ͔+Ü͔υ͕-Ü͕υ͔)e06",
    ///     "+(-E͔α͕+E͕α͔-I͔β͕+I͕β͔-L͔γ͕+L͕γ͔-N͔δ͕+N͕δ͔-O͔ε͕+O͕ε͔-S͔η͕+S͕η͔-Á͔θ͕+Á͕θ͔-Å͔ι͕\
    ///        +Å͕ι͔-Æ͔κ͕+Æ͕κ͔-Ë͔μ͕+Ë͕μ͔-Ï͔ν͕+Ï͕ν͔-Ñ͔ξ͕+Ñ͕ξ͔-Ö͔π͕+Ö͕π͔-Ú͔ρ͕+Ú͕ρ͔-Ü͔τ͕+Ü͕τ͔)e07",
    ///     "+(+β͔η͕-β͕η͔+γ͔θ͕-γ͕θ͔+δ͔ι͕-δ͕ι͔+ε͔κ͕-ε͕κ͔+ζ͔λ͕-ζ͕λ͔)e12",
    ///     "+(-α͔η͕+α͕η͔+γ͔μ͕-γ͕μ͔+δ͔ν͕-δ͕ν͔+ε͔ξ͕-ε͕ξ͔+ζ͔ο͕-ζ͕ο͔)e13",
    ///     "+(-α͔θ͕+α͕θ͔-β͔μ͕+β͕μ͔+δ͔π͕-δ͕π͔+ε͔ρ͕-ε͕ρ͔+ζ͔σ͕-ζ͕σ͔)e14",
    ///     "+(-α͔ι͕+α͕ι͔-β͔ν͕+β͕ν͔-γ͔π͕+γ͕π͔+ε͔τ͕-ε͕τ͔+ζ͔υ͕-ζ͕υ͔)e15",
    ///     "+(-α͔κ͕+α͕κ͔-β͔ξ͕+β͕ξ͔-γ͔ρ͕+γ͕ρ͔-δ͔τ͕+δ͕τ͔+ζ͔φ͕-ζ͕φ͔)e16",
    ///     "+(-α͔λ͕+α͕λ͔-β͔ο͕+β͕ο͔-γ͔σ͕+γ͕σ͔-δ͔υ͕+δ͕υ͔-ε͔φ͕+ε͕φ͔)e17",
    ///     "+(+α͔β͕-α͕β͔+θ͔μ͕-θ͕μ͔+ι͔ν͕-ι͕ν͔+κ͔ξ͕-κ͕ξ͔+λ͔ο͕-λ͕ο͔)e23",
    ///     "+(+α͔γ͕-α͕γ͔-η͔μ͕+η͕μ͔+ι͔π͕-ι͕π͔+κ͔ρ͕-κ͕ρ͔+λ͔σ͕-λ͕σ͔)e24",
    ///     "+(+α͔δ͕-α͕δ͔-η͔ν͕+η͕ν͔-θ͔π͕+θ͕π͔+κ͔τ͕-κ͕τ͔+λ͔υ͕-λ͕υ͔)e25",
    ///     "+(+α͔ε͕-α͕ε͔-η͔ξ͕+η͕ξ͔-θ͔ρ͕+θ͕ρ͔-ι͔τ͕+ι͕τ͔+λ͔φ͕-λ͕φ͔)e26",
    ///     "+(+α͔ζ͕-α͕ζ͔-η͔ο͕+η͕ο͔-θ͔σ͕+θ͕σ͔-ι͔υ͕+ι͕υ͔-κ͔φ͕+κ͕φ͔)e27",
    ///     "+(+β͔γ͕-β͕γ͔+η͔θ͕-η͕θ͔+ν͔π͕-ν͕π͔+ξ͔ρ͕-ξ͕ρ͔+ο͔σ͕-ο͕σ͔)e34",
    ///     "+(+β͔δ͕-β͕δ͔+η͔ι͕-η͕ι͔-μ͔π͕+μ͕π͔+ξ͔τ͕-ξ͕τ͔+ο͔υ͕-ο͕υ͔)e35",
    ///     "+(+β͔ε͕-β͕ε͔+η͔κ͕-η͕κ͔-μ͔ρ͕+μ͕ρ͔-ν͔τ͕+ν͕τ͔+ο͔φ͕-ο͕φ͔)e36",
    ///     "+(+β͔ζ͕-β͕ζ͔+η͔λ͕-η͕λ͔-μ͔σ͕+μ͕σ͔-ν͔υ͕+ν͕υ͔-ξ͔φ͕+ξ͕φ͔)e37",
    ///     "+(+γ͔δ͕-γ͕δ͔+θ͔ι͕-θ͕ι͔+μ͔ν͕-μ͕ν͔+ρ͔τ͕-ρ͕τ͔+σ͔υ͕-σ͕υ͔)e45",
    ///     "+(+γ͔ε͕-γ͕ε͔+θ͔κ͕-θ͕κ͔+μ͔ξ͕-μ͕ξ͔-π͔τ͕+π͕τ͔+σ͔φ͕-σ͕φ͔)e46",
    ///     "+(+γ͔ζ͕-γ͕ζ͔+θ͔λ͕-θ͕λ͔+μ͔ο͕-μ͕ο͔-π͔υ͕+π͕υ͔-ρ͔φ͕+ρ͕φ͔)e47",
    ///     "+(+δ͔ε͕-δ͕ε͔+ι͔κ͕-ι͕κ͔+ν͔ξ͕-ν͕ξ͔+π͔ρ͕-π͕ρ͔+υ͔φ͕-υ͕φ͔)e56",
    ///     "+(+δ͔ζ͕-δ͕ζ͔+ι͔λ͕-ι͕λ͔+ν͔ο͕-ν͕ο͔+π͔σ͕-π͕σ͔-τ͔φ͕+τ͕φ͔)e57",
    ///     "+(+ε͔ζ͕-ε͕ζ͔+κ͔λ͕-κ͕λ͔+ξ͔ο͕-ξ͕ο͔+ρ͔σ͕-ρ͕σ͔+τ͔υ͕-τ͕υ͔)e67",
    ///     "+(-α͔μ͕-α͕μ͔+β͔θ͕+β͕θ͔-γ͔η͕-γ͕η͔)e1234",
    ///     "+(-B͔μ͕-B͕μ͔-C͔ν͕-C͕ν͔-D͔ξ͕-D͕ξ͔-E͔ο͕-E͕ο͔+F͔θ͕+F͕θ͔+G͔ι͕+G͕ι͔+H͔κ͕+H͕κ͔+I͔λ͕+I͕λ͔-P͔γ͕-P͕γ͔-Q͔δ͕-Q͕δ͔-R͔ε͕-R͕ε͔-S͔ζ͕-S͕ζ͔)e0123",
    ///     "+(+A͔μ͕+A͕μ͔-C͔π͕-C͕π͔-D͔ρ͕-D͕ρ͔-E͔σ͕-E͕σ͔-F͔η͕-F͕η͔+J͔ι͕+J͕ι͔+K͔κ͕+K͕κ͔+L͔λ͕+L͕λ͔+P͔β͕+P͕β͔-T͔δ͕-T͕δ͔-U͔ε͕-U͕ε͔-Á͔ζ͕-Á͕ζ͔)e0124",
    ///     "+(+A͔ν͕+A͕ν͔+B͔π͕+B͕π͔-D͔τ͕-D͕τ͔-E͔υ͕-E͕υ͔-G͔η͕-G͕η͔-J͔θ͕-J͕θ͔+M͔κ͕+M͕κ͔+N͔λ͕+N͕λ͔+Q͔β͕+Q͕β͔+T͔γ͕+T͕γ͔-Ä͔ε͕-Ä͕ε͔-Å͔ζ͕-Å͕ζ͔)e0125",
    ///     "+(+A͔ξ͕+A͕ξ͔+B͔ρ͕+B͕ρ͔+C͔τ͕+C͕τ͔-E͔φ͕-E͕φ͔-H͔η͕-H͕η͔-K͔θ͕-K͕θ͔-M͔ι͕-M͕ι͔+O͔λ͕+O͕λ͔+R͔β͕+R͕β͔+U͔γ͕+U͕γ͔+Ä͔δ͕+Ä͕δ͔-Æ͔ζ͕-Æ͕ζ͔)e0126",
    ///     "+(+A͔ο͕+A͕ο͔+B͔σ͕+B͕σ͔+C͔υ͕+C͕υ͔+D͔φ͕+D͕φ͔-I͔η͕-I͕η͔-L͔θ͕-L͕θ͔-N͔ι͕-N͕ι͔-O͔κ͕-O͕κ͔+S͔β͕+S͕β͔+Á͔γ͕+Á͕γ͔+Å͔δ͕+Å͕δ͔+Æ͔ε͕+Æ͕ε͔)e0127",
    ///     "+(-A͔θ͕-A͕θ͔+B͔η͕+B͕η͔-G͔π͕-G͕π͔-H͔ρ͕-H͕ρ͔-I͔σ͕-I͕σ͔+J͔ν͕+J͕ν͔+K͔ξ͕+K͕ξ͔+L͔ο͕+L͕ο͔-P͔α͕-P͕α͔-Ç͔δ͕-Ç͕δ͔-É͔ε͕-É͕ε͔-Ë͔ζ͕-Ë͕ζ͔)e0134",
    ///     "+(-A͔ι͕-A͕ι͔+C͔η͕+C͕η͔+F͔π͕+F͕π͔-H͔τ͕-H͕τ͔-I͔υ͕-I͕υ͔-J͔μ͕-J͕μ͔+M͔ξ͕+M͕ξ͔+N͔ο͕+N͕ο͔-Q͔α͕-Q͕α͔+Ç͔γ͕+Ç͕γ͔-Í͔ε͕-Í͕ε͔-Ï͔ζ͕-Ï͕ζ͔)e0135",
    ///     "+(-A͔κ͕-A͕κ͔+D͔η͕+D͕η͔+F͔ρ͕+F͕ρ͔+G͔τ͕+G͕τ͔-I͔φ͕-I͕φ͔-K͔μ͕-K͕μ͔-M͔ν͕-M͕ν͔+O͔ο͕+O͕ο͔-R͔α͕-R͕α͔+É͔γ͕+É͕γ͔+Í͔δ͕+Í͕δ͔-Ñ͔ζ͕-Ñ͕ζ͔)e0136",
    ///     "+(-A͔λ͕-A͕λ͔+E͔η͕+E͕η͔+F͔σ͕+F͕σ͔+G͔υ͕+G͕υ͔+H͔φ͕+H͕φ͔-L͔μ͕-L͕μ͔-N͔ν͕-N͕ν͔-O͔ξ͕-O͕ξ͔-S͔α͕-S͕α͔+Ë͔γ͕+Ë͕γ͔+Ï͔δ͕+Ï͕δ͔+Ñ͔ε͕+Ñ͕ε͔)e0137",
    ///     "+(-B͔ι͕-B͕ι͔+C͔θ͕+C͕θ͔-F͔ν͕-F͕ν͔+G͔μ͕+G͕μ͔-K͔τ͕-K͕τ͔-L͔υ͕-L͕υ͔+M͔ρ͕+M͕ρ͔+N͔σ͕+N͕σ͔-T͔α͕-T͕α͔-Ç͔β͕-Ç͕β͔-Ó͔ε͕-Ó͕ε͔-Ö͔ζ͕-Ö͕ζ͔)e0145",
    ///     "+(-B͔κ͕-B͕κ͔+D͔θ͕+D͕θ͔-F͔ξ͕-F͕ξ͔+H͔μ͕+H͕μ͔+J͔τ͕+J͕τ͔-L͔φ͕-L͕φ͔-M͔π͕-M͕π͔+O͔σ͕+O͕σ͔-U͔α͕-U͕α͔-É͔β͕-É͕β͔+Ó͔δ͕+Ó͕δ͔-Ú͔ζ͕-Ú͕ζ͔)e0146",
    ///     "+(-B͔λ͕-B͕λ͔+E͔θ͕+E͕θ͔-F͔ο͕-F͕ο͔+I͔μ͕+I͕μ͔+J͔υ͕+J͕υ͔+K͔φ͕+K͕φ͔-N͔π͕-N͕π͔-O͔ρ͕-O͕ρ͔-Á͔α͕-Á͕α͔-Ë͔β͕-Ë͕β͔+Ö͔δ͕+Ö͕δ͔+Ú͔ε͕+Ú͕ε͔)e0147",
    ///     "+(-C͔κ͕-C͕κ͔+D͔ι͕+D͕ι͔-G͔ξ͕-G͕ξ͔+H͔ν͕+H͕ν͔-J͔ρ͕-J͕ρ͔+K͔π͕+K͕π͔-N͔φ͕-N͕φ͔+O͔υ͕+O͕υ͔-Ä͔α͕-Ä͕α͔-Í͔β͕-Í͕β͔-Ó͔γ͕-Ó͕γ͔-Ü͔ζ͕-Ü͕ζ͔)e0156",
    ///     "+(-C͔λ͕-C͕λ͔+E͔ι͕+E͕ι͔-G͔ο͕-G͕ο͔+I͔ν͕+I͕ν͔-J͔σ͕-J͕σ͔+L͔π͕+L͕π͔+M͔φ͕+M͕φ͔-O͔τ͕-O͕τ͔-Å͔α͕-Å͕α͔-Ï͔β͕-Ï͕β͔-Ö͔γ͕-Ö͕γ͔+Ü͔ε͕+Ü͕ε͔)e0157",
    ///     "+(-D͔λ͕-D͕λ͔+E͔κ͕+E͕κ͔-H͔ο͕-H͕ο͔+I͔ξ͕+I͕ξ͔-K͔σ͕-K͕σ͔+L͔ρ͕+L͕ρ͔-M͔υ͕-M͕υ͔+N͔τ͕+N͕τ͔-Æ͔α͕-Æ͕α͔-Ñ͔β͕-Ñ͕β͔-Ú͔γ͕-Ú͕γ͔-Ü͔δ͕-Ü͕δ͔)e0167",
    ///     "+(+A͔γ͕+A͕γ͔-B͔β͕-B͕β͔+F͔α͕+F͕α͔-Q͔π͕-Q͕π͔-R͔ρ͕-R͕ρ͔-S͔σ͕-S͕σ͔+T͔ν͕+T͕ν͔+U͔ξ͕+U͕ξ͔+Á͔ο͕+Á͕ο͔-Ç͔ι͕-Ç͕ι͔-É͔κ͕-É͕κ͔-Ë͔λ͕-Ë͕λ͔)e0234",
    ///     "+(+A͔δ͕+A͕δ͔-C͔β͕-C͕β͔+G͔α͕+G͕α͔+P͔π͕+P͕π͔-R͔τ͕-R͕τ͔-S͔υ͕-S͕υ͔-T͔μ͕-T͕μ͔+Ä͔ξ͕+Ä͕ξ͔+Å͔ο͕+Å͕ο͔+Ç͔θ͕+Ç͕θ͔-Í͔κ͕-Í͕κ͔-Ï͔λ͕-Ï͕λ͔)e0235",
    ///     "+(+A͔ε͕+A͕ε͔-D͔β͕-D͕β͔+H͔α͕+H͕α͔+P͔ρ͕+P͕ρ͔+Q͔τ͕+Q͕τ͔-S͔φ͕-S͕φ͔-U͔μ͕-U͕μ͔-Ä͔ν͕-Ä͕ν͔+Æ͔ο͕+Æ͕ο͔+É͔θ͕+É͕θ͔+Í͔ι͕+Í͕ι͔-Ñ͔λ͕-Ñ͕λ͔)e0236",
    ///     "+(+A͔ζ͕+A͕ζ͔-E͔β͕-E͕β͔+I͔α͕+I͕α͔+P͔σ͕+P͕σ͔+Q͔υ͕+Q͕υ͔+R͔φ͕+R͕φ͔-Á͔μ͕-Á͕μ͔-Å͔ν͕-Å͕ν͔-Æ͔ξ͕-Æ͕ξ͔+Ë͔θ͕+Ë͕θ͔+Ï͔ι͕+Ï͕ι͔+Ñ͔κ͕+Ñ͕κ͔)e0237",
    ///     "+(+B͔δ͕+B͕δ͔-C͔γ͕-C͕γ͔+J͔α͕+J͕α͔-P͔ν͕-P͕ν͔+Q͔μ͕+Q͕μ͔-U͔τ͕-U͕τ͔-Á͔υ͕-Á͕υ͔+Ä͔ρ͕+Ä͕ρ͔+Å͔σ͕+Å͕σ͔-Ç͔η͕-Ç͕η͔-Ó͔κ͕-Ó͕κ͔-Ö͔λ͕-Ö͕λ͔)e0245",
    ///     "+(+B͔ε͕+B͕ε͔-D͔γ͕-D͕γ͔+K͔α͕+K͕α͔-P͔ξ͕-P͕ξ͔+R͔μ͕+R͕μ͔+T͔τ͕+T͕τ͔-Á͔φ͕-Á͕φ͔-Ä͔π͕-Ä͕π͔+Æ͔σ͕+Æ͕σ͔-É͔η͕-É͕η͔+Ó͔ι͕+Ó͕ι͔-Ú͔λ͕-Ú͕λ͔)e0246",
    ///     "+(+B͔ζ͕+B͕ζ͔-E͔γ͕-E͕γ͔+L͔α͕+L͕α͔-P͔ο͕-P͕ο͔+S͔μ͕+S͕μ͔+T͔υ͕+T͕υ͔+U͔φ͕+U͕φ͔-Å͔π͕-Å͕π͔-Æ͔ρ͕-Æ͕ρ͔-Ë͔η͕-Ë͕η͔+Ö͔ι͕+Ö͕ι͔+Ú͔κ͕+Ú͕κ͔)e0247",
    ///     "+(+C͔ε͕+C͕ε͔-D͔δ͕-D͕δ͔+M͔α͕+M͕α͔-Q͔ξ͕-Q͕ξ͔+R͔ν͕+R͕ν͔-T͔ρ͕-T͕ρ͔+U͔π͕+U͕π͔-Å͔φ͕-Å͕φ͔+Æ͔υ͕+Æ͕υ͔-Í͔η͕-Í͕η͔-Ó͔θ͕-Ó͕θ͔-Ü͔λ͕-Ü͕λ͔)e0256",
    ///     "+(+C͔ζ͕+C͕ζ͔-E͔δ͕-E͕δ͔+N͔α͕+N͕α͔-Q͔ο͕-Q͕ο͔+S͔ν͕+S͕ν͔-T͔σ͕-T͕σ͔+Á͔π͕+Á͕π͔+Ä͔φ͕+Ä͕φ͔-Æ͔τ͕-Æ͕τ͔-Ï͔η͕-Ï͕η͔-Ö͔θ͕-Ö͕θ͔+Ü͔κ͕+Ü͕κ͔)e0257",
    ///     "+(+D͔ζ͕+D͕ζ͔-E͔ε͕-E͕ε͔+O͔α͕+O͕α͔-R͔ο͕-R͕ο͔+S͔ξ͕+S͕ξ͔-U͔σ͕-U͕σ͔+Á͔ρ͕+Á͕ρ͔-Ä͔υ͕-Ä͕υ͔+Å͔τ͕+Å͕τ͔-Ñ͔η͕-Ñ͕η͔-Ú͔θ͕-Ú͕θ͔-Ü͔ι͕-Ü͕ι͔)e0267",
    ///     "+(+F͔δ͕+F͕δ͔-G͔γ͕-G͕γ͔+J͔β͕+J͕β͔+P͔ι͕+P͕ι͔-Q͔θ͕-Q͕θ͔+T͔η͕+T͕η͔-É͔τ͕-É͕τ͔-Ë͔υ͕-Ë͕υ͔+Í͔ρ͕+Í͕ρ͔+Ï͔σ͕+Ï͕σ͔-Ó͔ξ͕-Ó͕ξ͔-Ö͔ο͕-Ö͕ο͔)e0345",
    ///     "+(+F͔ε͕+F͕ε͔-H͔γ͕-H͕γ͔+K͔β͕+K͕β͔+P͔κ͕+P͕κ͔-R͔θ͕-R͕θ͔+U͔η͕+U͕η͔+Ç͔τ͕+Ç͕τ͔-Ë͔φ͕-Ë͕φ͔-Í͔π͕-Í͕π͔+Ñ͔σ͕+Ñ͕σ͔+Ó͔ν͕+Ó͕ν͔-Ú͔ο͕-Ú͕ο͔)e0346",
    ///     "+(+F͔ζ͕+F͕ζ͔-I͔γ͕-I͕γ͔+L͔β͕+L͕β͔+P͔λ͕+P͕λ͔-S͔θ͕-S͕θ͔+Á͔η͕+Á͕η͔+Ç͔υ͕+Ç͕υ͔+É͔φ͕+É͕φ͔-Ï͔π͕-Ï͕π͔-Ñ͔ρ͕-Ñ͕ρ͔+Ö͔ν͕+Ö͕ν͔+Ú͔ξ͕+Ú͕ξ͔)e0347",
    ///     "+(+G͔ε͕+G͕ε͔-H͔δ͕-H͕δ͔+M͔β͕+M͕β͔+Q͔κ͕+Q͕κ͔-R͔ι͕-R͕ι͔+Ä͔η͕+Ä͕η͔-Ç͔ρ͕-Ç͕ρ͔+É͔π͕+É͕π͔-Ï͔φ͕-Ï͕φ͔+Ñ͔υ͕+Ñ͕υ͔-Ó͔μ͕-Ó͕μ͔-Ü͔ο͕-Ü͕ο͔)e0356",
    ///     "+(+G͔ζ͕+G͕ζ͔-I͔δ͕-I͕δ͔+N͔β͕+N͕β͔+Q͔λ͕+Q͕λ͔-S͔ι͕-S͕ι͔+Å͔η͕+Å͕η͔-Ç͔σ͕-Ç͕σ͔+Ë͔π͕+Ë͕π͔+Í͔φ͕+Í͕φ͔-Ñ͔τ͕-Ñ͕τ͔-Ö͔μ͕-Ö͕μ͔+Ü͔ξ͕+Ü͕ξ͔)e0357",
    ///     "+(+H͔ζ͕+H͕ζ͔-I͔ε͕-I͕ε͔+O͔β͕+O͕β͔+R͔λ͕+R͕λ͔-S͔κ͕-S͕κ͔+Æ͔η͕+Æ͕η͔-É͔σ͕-É͕σ͔+Ë͔ρ͕+Ë͕ρ͔-Í͔υ͕-Í͕υ͔+Ï͔τ͕+Ï͕τ͔-Ú͔μ͕-Ú͕μ͔-Ü͔ν͕-Ü͕ν͔)e0367",
    ///     "+(+J͔ε͕+J͕ε͔-K͔δ͕-K͕δ͔+M͔γ͕+M͕γ͔+T͔κ͕+T͕κ͔-U͔ι͕-U͕ι͔+Ä͔θ͕+Ä͕θ͔+Ç͔ξ͕+Ç͕ξ͔-É͔ν͕-É͕ν͔+Í͔μ͕+Í͕μ͔-Ö͔φ͕-Ö͕φ͔+Ú͔υ͕+Ú͕υ͔-Ü͔σ͕-Ü͕σ͔)e0456",
    ///     "+(+J͔ζ͕+J͕ζ͔-L͔δ͕-L͕δ͔+N͔γ͕+N͕γ͔+T͔λ͕+T͕λ͔-Á͔ι͕-Á͕ι͔+Å͔θ͕+Å͕θ͔+Ç͔ο͕+Ç͕ο͔-Ë͔ν͕-Ë͕ν͔+Ï͔μ͕+Ï͕μ͔+Ó͔φ͕+Ó͕φ͔-Ú͔τ͕-Ú͕τ͔+Ü͔ρ͕+Ü͕ρ͔)e0457",
    ///     "+(+K͔ζ͕+K͕ζ͔-L͔ε͕-L͕ε͔+O͔γ͕+O͕γ͔+U͔λ͕+U͕λ͔-Á͔κ͕-Á͕κ͔+Æ͔θ͕+Æ͕θ͔+É͔ο͕+É͕ο͔-Ë͔ξ͕-Ë͕ξ͔+Ñ͔μ͕+Ñ͕μ͔-Ó͔υ͕-Ó͕υ͔+Ö͔τ͕+Ö͕τ͔-Ü͔π͕-Ü͕π͔)e0467",
    ///     "+(+M͔ζ͕+M͕ζ͔-N͔ε͕-N͕ε͔+O͔δ͕+O͕δ͔+Ä͔λ͕+Ä͕λ͔-Å͔κ͕-Å͕κ͔+Æ͔ι͕+Æ͕ι͔+Í͔ο͕+Í͕ο͔-Ï͔ξ͕-Ï͕ξ͔+Ñ͔ν͕+Ñ͕ν͔+Ó͔σ͕+Ó͕σ͔-Ö͔ρ͕-Ö͕ρ͔+Ú͔π͕+Ú͕π͔)e0567",
    ///     "+(-π͔φ͕-π͕φ͔+ρ͔υ͕+ρ͕υ͔-σ͔τ͕-σ͕τ͔)e4567",
    ///     "+(+ν͔φ͕+ν͕φ͔-ξ͔υ͕-ξ͕υ͔+ο͔τ͕+ο͕τ͔)e3576",
    ///     "+(-μ͔φ͕-μ͕φ͔+ξ͔σ͕+ξ͕σ͔-ο͔ρ͕-ο͕ρ͔)e3467",
    ///     "+(+μ͔υ͕+μ͕υ͔-ν͔σ͕-ν͕σ͔+ο͔π͕+ο͕π͔)e3475",
    ///     "+(-μ͔τ͕-μ͕τ͔+ν͔ρ͕+ν͕ρ͔-ξ͔π͕-ξ͕π͔)e3456",
    ///     "+(-ι͔φ͕-ι͕φ͔+κ͔υ͕+κ͕υ͔-λ͔τ͕-λ͕τ͔)e2567",
    ///     "+(+θ͔φ͕+θ͕φ͔-κ͔σ͕-κ͕σ͔+λ͔ρ͕+λ͕ρ͔)e2476",
    ///     "+(-θ͔υ͕-θ͕υ͔+ι͔σ͕+ι͕σ͔-λ͔π͕-λ͕π͔)e2457",
    ///     "+(+θ͔τ͕+θ͕τ͔-ι͔ρ͕-ι͕ρ͔+κ͔π͕+κ͕π͔)e2465",
    ///     "+(-η͔φ͕-η͕φ͔+κ͔ο͕+κ͕ο͔-λ͔ξ͕-λ͕ξ͔)e2367",
    ///     "+(+η͔υ͕+η͕υ͔-ι͔ο͕-ι͕ο͔+λ͔ν͕+λ͕ν͔)e2375",
    ///     "+(-η͔τ͕-η͕τ͔+ι͔ξ͕+ι͕ξ͔-κ͔ν͕-κ͕ν͔)e2356",
    ///     "+(-η͔σ͕-η͕σ͔+θ͔ο͕+θ͕ο͔-λ͔μ͕-λ͕μ͔)e2347",
    ///     "+(+η͔ρ͕+η͕ρ͔-θ͔ξ͕-θ͕ξ͔+κ͔μ͕+κ͕μ͔)e2364",
    ///     "+(-η͔π͕-η͕π͔+θ͔ν͕+θ͕ν͔-ι͔μ͕-ι͕μ͔)e2345",
    ///     "+(+δ͔φ͕+δ͕φ͔-ε͔υ͕-ε͕υ͔+ζ͔τ͕+ζ͕τ͔)e1576",
    ///     "+(-γ͔φ͕-γ͕φ͔+ε͔σ͕+ε͕σ͔-ζ͔ρ͕-ζ͕ρ͔)e1467",
    ///     "+(+γ͔υ͕+γ͕υ͔-δ͔σ͕-δ͕σ͔+ζ͔π͕+ζ͕π͔)e1475",
    ///     "+(-γ͔τ͕-γ͕τ͔+δ͔ρ͕+δ͕ρ͔-ε͔π͕-ε͕π͔)e1456",
    ///     "+(+β͔φ͕+β͕φ͔-ε͔ο͕-ε͕ο͔+ζ͔ξ͕+ζ͕ξ͔)e1376",
    ///     "+(-β͔υ͕-β͕υ͔+δ͔ο͕+δ͕ο͔-ζ͔ν͕-ζ͕ν͔)e1357",
    ///     "+(+β͔τ͕+β͕τ͔-δ͔ξ͕-δ͕ξ͔+ε͔ν͕+ε͕ν͔)e1365",
    ///     "+(+β͔σ͕+β͕σ͔-γ͔ο͕-γ͕ο͔+ζ͔μ͕+ζ͕μ͔)e1374",
    ///     "+(-β͔ρ͕-β͕ρ͔+γ͔ξ͕+γ͕ξ͔-ε͔μ͕-ε͕μ͔)e1346",
    ///     "+(+β͔π͕+β͕π͔-γ͔ν͕-γ͕ν͔+δ͔μ͕+δ͕μ͔)e1354",
    ///     "+(-α͔φ͕-α͕φ͔+ε͔λ͕+ε͕λ͔-ζ͔κ͕-ζ͕κ͔)e1267",
    ///     "+(+α͔υ͕+α͕υ͔-δ͔λ͕-δ͕λ͔+ζ͔ι͕+ζ͕ι͔)e1275",
    ///     "+(-α͔τ͕-α͕τ͔+δ͔κ͕+δ͕κ͔-ε͔ι͕-ε͕ι͔)e1256",
    ///     "+(-α͔σ͕-α͕σ͔+γ͔λ͕+γ͕λ͔-ζ͔θ͕-ζ͕θ͔)e1247",
    ///     "+(+α͔ρ͕+α͕ρ͔-γ͔κ͕-γ͕κ͔+ε͔θ͕+ε͕θ͔)e1264",
    ///     "+(-α͔π͕-α͕π͔+γ͔ι͕+γ͕ι͔-δ͔θ͕-δ͕θ͔)e1245",
    ///     "+(+α͔ο͕+α͕ο͔-β͔λ͕-β͕λ͔+ζ͔η͕+ζ͕η͔)e1273",
    ///     "+(-α͔ξ͕-α͕ξ͔+β͔κ͕+β͕κ͔-ε͔η͕-ε͕η͔)e1236",
    ///     "+(+α͔ν͕+α͕ν͔-β͔ι͕-β͕ι͔+δ͔η͕+δ͕η͔)e1253",
    ///     "+(+Ç͔φ͕-Ç͕φ͔-É͔υ͕+É͕υ͔+Ë͔τ͕-Ë͕τ͔+Í͔σ͕-Í͕σ͔-Ï͔ρ͕+Ï͕ρ͔+Ñ͔π͕-Ñ͕π͔-Ó͔ο͕+Ó͕ο͔+Ö͔ξ͕-Ö͕ξ͔-Ú͔ν͕+Ú͕ν͔+Ü͔μ͕-Ü͕μ͔)e034567",
    ///     "+(-T͔φ͕+T͕φ͔+U͔υ͕-U͕υ͔-Á͔τ͕+Á͕τ͔-Ä͔σ͕+Ä͕σ͔+Å͔ρ͕-Å͕ρ͔-Æ͔π͕+Æ͕π͔+Ó͔λ͕-Ó͕λ͔-Ö͔κ͕+Ö͕κ͔+Ú͔ι͕-Ú͕ι͔-Ü͔θ͕+Ü͕θ͔)e024576",
    ///     "+(+Q͔φ͕-Q͕φ͔-R͔υ͕+R͕υ͔+S͔τ͕-S͕τ͔+Ä͔ο͕-Ä͕ο͔-Å͔ξ͕+Å͕ξ͔+Æ͔ν͕-Æ͕ν͔-Í͔λ͕+Í͕λ͔+Ï͔κ͕-Ï͕κ͔-Ñ͔ι͕+Ñ͕ι͔+Ü͔η͕-Ü͕η͔)e023567",
    ///     "+(-P͔φ͕+P͕φ͔+R͔σ͕-R͕σ͔-S͔ρ͕+S͕ρ͔-U͔ο͕+U͕ο͔+Á͔ξ͕-Á͕ξ͔-Æ͔μ͕+Æ͕μ͔+É͔λ͕-É͕λ͔-Ë͔κ͕+Ë͕κ͔+Ñ͔θ͕-Ñ͕θ͔-Ú͔η͕+Ú͕η͔)e023476",
    ///     "+(+P͔υ͕-P͕υ͔-Q͔σ͕+Q͕σ͔+S͔π͕-S͕π͔+T͔ο͕-T͕ο͔-Á͔ν͕+Á͕ν͔+Å͔μ͕-Å͕μ͔-Ç͔λ͕+Ç͕λ͔+Ë͔ι͕-Ë͕ι͔-Ï͔θ͕+Ï͕θ͔+Ö͔η͕-Ö͕η͔)e023457",
    ///     "+(-P͔τ͕+P͕τ͔+Q͔ρ͕-Q͕ρ͔-R͔π͕+R͕π͔-T͔ξ͕+T͕ξ͔+U͔ν͕-U͕ν͔-Ä͔μ͕+Ä͕μ͔+Ç͔κ͕-Ç͕κ͔-É͔ι͕+É͕ι͔+Í͔θ͕-Í͕θ͔-Ó͔η͕+Ó͕η͔)e023465",
    ///     "+(+J͔φ͕-J͕φ͔-K͔υ͕+K͕υ͔+L͔τ͕-L͕τ͔+M͔σ͕-M͕σ͔-N͔ρ͕+N͕ρ͔+O͔π͕-O͕π͔-Ó͔ζ͕+Ó͕ζ͔+Ö͔ε͕-Ö͕ε͔-Ú͔δ͕+Ú͕δ͔+Ü͔γ͕-Ü͕γ͔)e014567",
    ///     "+(-G͔φ͕+G͕φ͔+H͔υ͕-H͕υ͔-I͔τ͕+I͕τ͔-M͔ο͕+M͕ο͔+N͔ξ͕-N͕ξ͔-O͔ν͕+O͕ν͔+Í͔ζ͕-Í͕ζ͔-Ï͔ε͕+Ï͕ε͔+Ñ͔δ͕-Ñ͕δ͔-Ü͔β͕+Ü͕β͔)e013576",
    ///     "+(+F͔φ͕-F͕φ͔-H͔σ͕+H͕σ͔+I͔ρ͕-I͕ρ͔+K͔ο͕-K͕ο͔-L͔ξ͕+L͕ξ͔+O͔μ͕-O͕μ͔-É͔ζ͕+É͕ζ͔+Ë͔ε͕-Ë͕ε͔-Ñ͔γ͕+Ñ͕γ͔+Ú͔β͕-Ú͕β͔)e013467",
    ///     "+(-F͔υ͕+F͕υ͔+G͔σ͕-G͕σ͔-I͔π͕+I͕π͔-J͔ο͕+J͕ο͔+L͔ν͕-L͕ν͔-N͔μ͕+N͕μ͔+Ç͔ζ͕-Ç͕ζ͔-Ë͔δ͕+Ë͕δ͔+Ï͔γ͕-Ï͕γ͔-Ö͔β͕+Ö͕β͔)e013475",
    ///     "+(+F͔τ͕-F͕τ͔-G͔ρ͕+G͕ρ͔+H͔π͕-H͕π͔+J͔ξ͕-J͕ξ͔-K͔ν͕+K͕ν͔+M͔μ͕-M͕μ͔-Ç͔ε͕+Ç͕ε͔+É͔δ͕-É͕δ͔-Í͔γ͕+Í͕γ͔+Ó͔β͕-Ó͕β͔)e013456",
    ///     "+(+C͔φ͕-C͕φ͔-D͔υ͕+D͕υ͔+E͔τ͕-E͕τ͔+M͔λ͕-M͕λ͔-N͔κ͕+N͕κ͔+O͔ι͕-O͕ι͔-Ä͔ζ͕+Ä͕ζ͔+Å͔ε͕-Å͕ε͔-Æ͔δ͕+Æ͕δ͔+Ü͔α͕-Ü͕α͔)e012567",
    ///     "+(-B͔φ͕+B͕φ͔+D͔σ͕-D͕σ͔-E͔ρ͕+E͕ρ͔-K͔λ͕+K͕λ͔+L͔κ͕-L͕κ͔-O͔θ͕+O͕θ͔+U͔ζ͕-U͕ζ͔-Á͔ε͕+Á͕ε͔+Æ͔γ͕-Æ͕γ͔-Ú͔α͕+Ú͕α͔)e012476",
    ///     "+(+B͔υ͕-B͕υ͔-C͔σ͕+C͕σ͔+E͔π͕-E͕π͔+J͔λ͕-J͕λ͔-L͔ι͕+L͕ι͔+N͔θ͕-N͕θ͔-T͔ζ͕+T͕ζ͔+Á͔δ͕-Á͕δ͔-Å͔γ͕+Å͕γ͔+Ö͔α͕-Ö͕α͔)e012457",
    ///     "+(-B͔τ͕+B͕τ͔+C͔ρ͕-C͕ρ͔-D͔π͕+D͕π͔-J͔κ͕+J͕κ͔+K͔ι͕-K͕ι͔-M͔θ͕+M͕θ͔+T͔ε͕-T͕ε͔-U͔δ͕+U͕δ͔+Ä͔γ͕-Ä͕γ͔-Ó͔α͕+Ó͕α͔)e012465",
    ///     "+(+A͔φ͕-A͕φ͔-D͔ο͕+D͕ο͔+E͔ξ͕-E͕ξ͔+H͔λ͕-H͕λ͔-I͔κ͕+I͕κ͔+O͔η͕-O͕η͔-R͔ζ͕+R͕ζ͔+S͔ε͕-S͕ε͔-Æ͔β͕+Æ͕β͔+Ñ͔α͕-Ñ͕α͔)e012367",
    ///     "+(-A͔υ͕+A͕υ͔+C͔ο͕-C͕ο͔-E͔ν͕+E͕ν͔-G͔λ͕+G͕λ͔+I͔ι͕-I͕ι͔-N͔η͕+N͕η͔+Q͔ζ͕-Q͕ζ͔-S͔δ͕+S͕δ͔+Å͔β͕-Å͕β͔-Ï͔α͕+Ï͕α͔)e012375",
    ///     "+(+A͔τ͕-A͕τ͔-C͔ξ͕+C͕ξ͔+D͔ν͕-D͕ν͔+G͔κ͕-G͕κ͔-H͔ι͕+H͕ι͔+M͔η͕-M͕η͔-Q͔ε͕+Q͕ε͔+R͔δ͕-R͕δ͔-Ä͔β͕+Ä͕β͔+Í͔α͕-Í͕α͔)e012356",
    ///     "+(+A͔σ͕-A͕σ͔-B͔ο͕+B͕ο͔+E͔μ͕-E͕μ͔+F͔λ͕-F͕λ͔-I͔θ͕+I͕θ͔+L͔η͕-L͕η͔-P͔ζ͕+P͕ζ͔+S͔γ͕-S͕γ͔-Á͔β͕+Á͕β͔+Ë͔α͕-Ë͕α͔)e012347",
    ///     "+(-A͔ρ͕+A͕ρ͔+B͔ξ͕-B͕ξ͔-D͔μ͕+D͕μ͔-F͔κ͕+F͕κ͔+H͔θ͕-H͕θ͔-K͔η͕+K͕η͔+P͔ε͕-P͕ε͔-R͔γ͕+R͕γ͔+U͔β͕-U͕β͔-É͔α͕+É͕α͔)e012364",
    ///     "+(+A͔π͕-A͕π͔-B͔ν͕+B͕ν͔+C͔μ͕-C͕μ͔+F͔ι͕-F͕ι͔-G͔θ͕+G͕θ͔+J͔η͕-J͕η͔-P͔δ͕+P͕δ͔+Q͔γ͕-Q͕γ͔-T͔β͕+T͕β͔+Ç͔α͕-Ç͕α͔)e012345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume5() + Self::volume() + Self::line_moment()
    }
    /// The multivector of double motor $`m_2 \equiv s + v^5 + v + \ell`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let simple_triple_motor = Vee::volume4().lhs() * Vee::volume4().rhs();
    ///
    /// assert_eq!(simple_triple_motor.basis_blades(), Vee::simple_triple_motor().basis_blades());
    /// format_eq!(simple_triple_motor, [
    ///     "-a͔a͕-b͔b͕-c͔c͕-d͔d͕-e͔e͕-f͔f͕-g͔g͕-h͔h͕-i͔i͕-j͔j͕-k͔k͕-l͔l͕-m͔m͕-n͔n͕-o͔o͕-p͔p͕-q͔q͕\
    ///      -r͔r͕-s͔s͕-t͔t͕-u͔u͕-á͔á͕-ä͔ä͕-å͔å͕-æ͔æ͕-ç͔ç͕-é͔é͕-ë͔ë͕-í͔í͕-ï͔ï͕-ñ͔ñ͕-ó͔ó͕-ö͔ö͕-ú͔ú͕-ü͔ü͕",
    ///     "+(+a͔Η͕-a͕Η͔+b͔Θ͕-b͕Θ͔+c͔Ι͕-c͕Ι͔+d͔Κ͕-d͕Κ͔+e͔Λ͕-e͕Λ͔+f͔Μ͕-f͕Μ͔+g͔Ν͕-g͕Ν͔+h͔Ξ͕\
    ///        -h͕Ξ͔+i͔Ο͕-i͕Ο͔+j͔Π͕-j͕Π͔+k͔Ρ͕-k͕Ρ͔+l͔Σ͕-l͕Σ͔+m͔Τ͕-m͕Τ͔+n͔Υ͕-n͕Υ͔+o͔Φ͕-o͕Φ͔)e01",
    ///     "+(-a͔Β͕+a͕Β͔-b͔Γ͕+b͕Γ͔-c͔Δ͕+c͕Δ͔-d͔Ε͕+d͕Ε͔-e͔Ζ͕+e͕Ζ͔+p͔Μ͕-p͕Μ͔+q͔Ν͕-q͕Ν͔+r͔Ξ͕\
    ///        -r͕Ξ͔+s͔Ο͕-s͕Ο͔+t͔Π͕-t͕Π͔+u͔Ρ͕-u͕Ρ͔+á͔Σ͕-á͕Σ͔+ä͔Τ͕-ä͕Τ͔+å͔Υ͕-å͕Υ͔+æ͔Φ͕-æ͕Φ͔)e02",
    ///     "+(+a͔Α͕-a͕Α͔-f͔Γ͕+f͕Γ͔-g͔Δ͕+g͕Δ͔-h͔Ε͕+h͕Ε͔-i͔Ζ͕+i͕Ζ͔-p͔Θ͕+p͕Θ͔-q͔Ι͕+q͕Ι͔-r͔Κ͕\
    ///        +r͕Κ͔-s͔Λ͕+s͕Λ͔+ç͔Π͕-ç͕Π͔+é͔Ρ͕-é͕Ρ͔+ë͔Σ͕-ë͕Σ͔+í͔Τ͕-í͕Τ͔+ï͔Υ͕-ï͕Υ͔+ñ͔Φ͕-ñ͕Φ͔)e03",
    ///     "+(+b͔Α͕-b͕Α͔+f͔Β͕-f͕Β͔-j͔Δ͕+j͕Δ͔-k͔Ε͕+k͕Ε͔-l͔Ζ͕+l͕Ζ͔+p͔Η͕-p͕Η͔-t͔Ι͕+t͕Ι͔-u͔Κ͕\
    ///        +u͕Κ͔-á͔Λ͕+á͕Λ͔-ç͔Ν͕+ç͕Ν͔-é͔Ξ͕+é͕Ξ͔-ë͔Ο͕+ë͕Ο͔+ó͔Τ͕-ó͕Τ͔+ö͔Υ͕-ö͕Υ͔+ú͔Φ͕-ú͕Φ͔)e04",
    ///     "+(+c͔Α͕-c͕Α͔+g͔Β͕-g͕Β͔+j͔Γ͕-j͕Γ͔-m͔Ε͕+m͕Ε͔-n͔Ζ͕+n͕Ζ͔+q͔Η͕-q͕Η͔+t͔Θ͕-t͕Θ͔-ä͔Κ͕\
    ///        +ä͕Κ͔-å͔Λ͕+å͕Λ͔+ç͔Μ͕-ç͕Μ͔-í͔Ξ͕+í͕Ξ͔-ï͔Ο͕+ï͕Ο͔-ó͔Ρ͕+ó͕Ρ͔-ö͔Σ͕+ö͕Σ͔+ü͔Φ͕-ü͕Φ͔)e05",
    ///     "+(+d͔Α͕-d͕Α͔+h͔Β͕-h͕Β͔+k͔Γ͕-k͕Γ͔+m͔Δ͕-m͕Δ͔-o͔Ζ͕+o͕Ζ͔+r͔Η͕-r͕Η͔+u͔Θ͕-u͕Θ͔+ä͔Ι͕\
    ///        -ä͕Ι͔-æ͔Λ͕+æ͕Λ͔+é͔Μ͕-é͕Μ͔+í͔Ν͕-í͕Ν͔-ñ͔Ο͕+ñ͕Ο͔+ó͔Π͕-ó͕Π͔-ú͔Σ͕+ú͕Σ͔-ü͔Υ͕+ü͕Υ͔)e06",
    ///     "+(+e͔Α͕-e͕Α͔+i͔Β͕-i͕Β͔+l͔Γ͕-l͕Γ͔+n͔Δ͕-n͕Δ͔+o͔Ε͕-o͕Ε͔+s͔Η͕-s͕Η͔+á͔Θ͕-á͕Θ͔+å͔Ι͕\
    ///        -å͕Ι͔+æ͔Κ͕-æ͕Κ͔+ë͔Μ͕-ë͕Μ͔+ï͔Ν͕-ï͕Ν͔+ñ͔Ξ͕-ñ͕Ξ͔+ö͔Π͕-ö͕Π͔+ú͔Ρ͕-ú͕Ρ͔+ü͔Τ͕-ü͕Τ͔)e07",
    ///     "+(-f͔p͕+f͕p͔-g͔q͕+g͕q͔-h͔r͕+h͕r͔-i͔s͕+i͕s͔-j͔t͕+j͕t͔-k͔u͕+k͕u͔-l͔á͕+l͕á͔-m͔ä͕+m͕ä͔-n͔å͕+n͕å͔-o͔æ͕+o͕æ͔)e12",
    ///     "+(+b͔p͕-b͕p͔+c͔q͕-c͕q͔+d͔r͕-d͕r͔+e͔s͕-e͕s͔-j͔ç͕+j͕ç͔-k͔é͕+k͕é͔-l͔ë͕+l͕ë͔-m͔í͕+m͕í͔-n͔ï͕+n͕ï͔-o͔ñ͕+o͕ñ͔)e13",
    ///     "+(-a͔p͕+a͕p͔+c͔t͕-c͕t͔+d͔u͕-d͕u͔+e͔á͕-e͕á͔+g͔ç͕-g͕ç͔+h͔é͕-h͕é͔+i͔ë͕-i͕ë͔-m͔ó͕+m͕ó͔-n͔ö͕+n͕ö͔-o͔ú͕+o͕ú͔)e14",
    ///     "+(-a͔q͕+a͕q͔-b͔t͕+b͕t͔+d͔ä͕-d͕ä͔+e͔å͕-e͕å͔-f͔ç͕+f͕ç͔+h͔í͕-h͕í͔+i͔ï͕-i͕ï͔+k͔ó͕-k͕ó͔+l͔ö͕-l͕ö͔-o͔ü͕+o͕ü͔)e15",
    ///     "+(-a͔r͕+a͕r͔-b͔u͕+b͕u͔-c͔ä͕+c͕ä͔+e͔æ͕-e͕æ͔-f͔é͕+f͕é͔-g͔í͕+g͕í͔+i͔ñ͕-i͕ñ͔-j͔ó͕+j͕ó͔+l͔ú͕-l͕ú͔+n͔ü͕-n͕ü͔)e16",
    ///     "+(-a͔s͕+a͕s͔-b͔á͕+b͕á͔-c͔å͕+c͕å͔-d͔æ͕+d͕æ͔-f͔ë͕+f͕ë͔-g͔ï͕+g͕ï͔-h͔ñ͕+h͕ñ͔-j͔ö͕+j͕ö͔-k͔ú͕+k͕ú͔-m͔ü͕+m͕ü͔)e17",
    ///     "+(-b͔f͕+b͕f͔-c͔g͕+c͕g͔-d͔h͕+d͕h͔-e͔i͕+e͕i͔-t͔ç͕+t͕ç͔-u͔é͕+u͕é͔-á͔ë͕+á͕ë͔-ä͔í͕+ä͕í͔-å͔ï͕+å͕ï͔-æ͔ñ͕+æ͕ñ͔)e23",
    ///     "+(+a͔f͕-a͕f͔-c͔j͕+c͕j͔-d͔k͕+d͕k͔-e͔l͕+e͕l͔+q͔ç͕-q͕ç͔+r͔é͕-r͕é͔+s͔ë͕-s͕ë͔-ä͔ó͕+ä͕ó͔-å͔ö͕+å͕ö͔-æ͔ú͕+æ͕ú͔)e24",
    ///     "+(+a͔g͕-a͕g͔+b͔j͕-b͕j͔-d͔m͕+d͕m͔-e͔n͕+e͕n͔-p͔ç͕+p͕ç͔+r͔í͕-r͕í͔+s͔ï͕-s͕ï͔+u͔ó͕-u͕ó͔+á͔ö͕-á͕ö͔-æ͔ü͕+æ͕ü͔)e25",
    ///     "+(+a͔h͕-a͕h͔+b͔k͕-b͕k͔+c͔m͕-c͕m͔-e͔o͕+e͕o͔-p͔é͕+p͕é͔-q͔í͕+q͕í͔+s͔ñ͕-s͕ñ͔-t͔ó͕+t͕ó͔+á͔ú͕-á͕ú͔+å͔ü͕-å͕ü͔)e26",
    ///     "+(+a͔i͕-a͕i͔+b͔l͕-b͕l͔+c͔n͕-c͕n͔+d͔o͕-d͕o͔-p͔ë͕+p͕ë͔-q͔ï͕+q͕ï͔-r͔ñ͕+r͕ñ͔-t͔ö͕+t͕ö͔-u͔ú͕+u͕ú͔-ä͔ü͕+ä͕ü͔)e27",
    ///     "+(-a͔b͕+a͕b͔-g͔j͕+g͕j͔-h͔k͕+h͕k͔-i͔l͕+i͕l͔-q͔t͕+q͕t͔-r͔u͕+r͕u͔-s͔á͕+s͕á͔-í͔ó͕+í͕ó͔-ï͔ö͕+ï͕ö͔-ñ͔ú͕+ñ͕ú͔)e34",
    ///     "+(-a͔c͕+a͕c͔+f͔j͕-f͕j͔-h͔m͕+h͕m͔-i͔n͕+i͕n͔+p͔t͕-p͕t͔-r͔ä͕+r͕ä͔-s͔å͕+s͕å͔+é͔ó͕-é͕ó͔+ë͔ö͕-ë͕ö͔-ñ͔ü͕+ñ͕ü͔)e35",
    ///     "+(-a͔d͕+a͕d͔+f͔k͕-f͕k͔+g͔m͕-g͕m͔-i͔o͕+i͕o͔+p͔u͕-p͕u͔+q͔ä͕-q͕ä͔-s͔æ͕+s͕æ͔-ç͔ó͕+ç͕ó͔+ë͔ú͕-ë͕ú͔+ï͔ü͕-ï͕ü͔)e36",
    ///     "+(-a͔e͕+a͕e͔+f͔l͕-f͕l͔+g͔n͕-g͕n͔+h͔o͕-h͕o͔+p͔á͕-p͕á͔+q͔å͕-q͕å͔+r͔æ͕-r͕æ͔-ç͔ö͕+ç͕ö͔-é͔ú͕+é͕ú͔-í͔ü͕+í͕ü͔)e37",
    ///     "+(-b͔c͕+b͕c͔-f͔g͕+f͕g͔-k͔m͕+k͕m͔-l͔n͕+l͕n͔-p͔q͕+p͕q͔-u͔ä͕+u͕ä͔-á͔å͕+á͕å͔-é͔í͕+é͕í͔-ë͔ï͕+ë͕ï͔-ú͔ü͕+ú͕ü͔)e45",
    ///     "+(-b͔d͕+b͕d͔-f͔h͕+f͕h͔+j͔m͕-j͕m͔-l͔o͕+l͕o͔-p͔r͕+p͕r͔+t͔ä͕-t͕ä͔-á͔æ͕+á͕æ͔+ç͔í͕-ç͕í͔-ë͔ñ͕+ë͕ñ͔+ö͔ü͕-ö͕ü͔)e46",
    ///     "+(-b͔e͕+b͕e͔-f͔i͕+f͕i͔+j͔n͕-j͕n͔+k͔o͕-k͕o͔-p͔s͕+p͕s͔+t͔å͕-t͕å͔+u͔æ͕-u͕æ͔+ç͔ï͕-ç͕ï͔+é͔ñ͕-é͕ñ͔-ó͔ü͕+ó͕ü͔)e47",
    ///     "+(-c͔d͕+c͕d͔-g͔h͕+g͕h͔-j͔k͕+j͕k͔-n͔o͕+n͕o͔-q͔r͕+q͕r͔-t͔u͕+t͕u͔-å͔æ͕+å͕æ͔-ç͔é͕+ç͕é͔-ï͔ñ͕+ï͕ñ͔-ö͔ú͕+ö͕ú͔)e56",
    ///     "+(-c͔e͕+c͕e͔-g͔i͕+g͕i͔-j͔l͕+j͕l͔+m͔o͕-m͕o͔-q͔s͕+q͕s͔-t͔á͕+t͕á͔+ä͔æ͕-ä͕æ͔-ç͔ë͕+ç͕ë͔+í͔ñ͕-í͕ñ͔+ó͔ú͕-ó͕ú͔)e57",
    ///     "+(-d͔e͕+d͕e͔-h͔i͕+h͕i͔-k͔l͕+k͕l͔-m͔n͕+m͕n͔-r͔s͕+r͕s͔-u͔á͕+u͕á͔-ä͔å͕+ä͕å͔-é͔ë͕+é͕ë͔-í͔ï͕+í͕ï͔-ó͔ö͕+ó͕ö͔)e67",
    ///     "+(+c͔ç͕+c͕ç͔+d͔é͕+d͕é͔+e͔ë͕+e͕ë͔-g͔t͕-g͕t͔-h͔u͕-h͕u͔-i͔á͕-i͕á͔+j͔q͕+j͕q͔+k͔r͕+k͕r͔+l͔s͕+l͕s͔)e1234",
    ///     "+(+b͔Μ͕+b͕Μ͔+c͔Ν͕+c͕Ν͔+d͔Ξ͕+d͕Ξ͔+e͔Ο͕+e͕Ο͔-f͔Θ͕-f͕Θ͔-g͔Ι͕-g͕Ι͔-h͔Κ͕-h͕Κ͔-i͔Λ͕-i͕Λ͔+p͔Γ͕+p͕Γ͔+q͔Δ͕+q͕Δ͔+r͔Ε͕+r͕Ε͔+s͔Ζ͕+s͕Ζ͔)e0123",
    ///     "+(-a͔Μ͕-a͕Μ͔+c͔Π͕+c͕Π͔+d͔Ρ͕+d͕Ρ͔+e͔Σ͕+e͕Σ͔+f͔Η͕+f͕Η͔-j͔Ι͕-j͕Ι͔-k͔Κ͕-k͕Κ͔-l͔Λ͕-l͕Λ͔-p͔Β͕-p͕Β͔+t͔Δ͕+t͕Δ͔+u͔Ε͕+u͕Ε͔+á͔Ζ͕+á͕Ζ͔)e0124",
    ///     "+(-a͔Ν͕-a͕Ν͔-b͔Π͕-b͕Π͔+d͔Τ͕+d͕Τ͔+e͔Υ͕+e͕Υ͔+g͔Η͕+g͕Η͔+j͔Θ͕+j͕Θ͔-m͔Κ͕-m͕Κ͔-n͔Λ͕-n͕Λ͔-q͔Β͕-q͕Β͔-t͔Γ͕-t͕Γ͔+ä͔Ε͕+ä͕Ε͔+å͔Ζ͕+å͕Ζ͔)e0125",
    ///     "+(-a͔Ξ͕-a͕Ξ͔-b͔Ρ͕-b͕Ρ͔-c͔Τ͕-c͕Τ͔+e͔Φ͕+e͕Φ͔+h͔Η͕+h͕Η͔+k͔Θ͕+k͕Θ͔+m͔Ι͕+m͕Ι͔-o͔Λ͕-o͕Λ͔-r͔Β͕-r͕Β͔-u͔Γ͕-u͕Γ͔-ä͔Δ͕-ä͕Δ͔+æ͔Ζ͕+æ͕Ζ͔)e0126",
    ///     "+(-a͔Ο͕-a͕Ο͔-b͔Σ͕-b͕Σ͔-c͔Υ͕-c͕Υ͔-d͔Φ͕-d͕Φ͔+i͔Η͕+i͕Η͔+l͔Θ͕+l͕Θ͔+n͔Ι͕+n͕Ι͔+o͔Κ͕+o͕Κ͔-s͔Β͕-s͕Β͔-á͔Γ͕-á͕Γ͔-å͔Δ͕-å͕Δ͔-æ͔Ε͕-æ͕Ε͔)e0127",
    ///     "+(+a͔Θ͕+a͕Θ͔-b͔Η͕-b͕Η͔+g͔Π͕+g͕Π͔+h͔Ρ͕+h͕Ρ͔+i͔Σ͕+i͕Σ͔-j͔Ν͕-j͕Ν͔-k͔Ξ͕-k͕Ξ͔-l͔Ο͕-l͕Ο͔+p͔Α͕+p͕Α͔+ç͔Δ͕+ç͕Δ͔+é͔Ε͕+é͕Ε͔+ë͔Ζ͕+ë͕Ζ͔)e0134",
    ///     "+(+a͔Ι͕+a͕Ι͔-c͔Η͕-c͕Η͔-f͔Π͕-f͕Π͔+h͔Τ͕+h͕Τ͔+i͔Υ͕+i͕Υ͔+j͔Μ͕+j͕Μ͔-m͔Ξ͕-m͕Ξ͔-n͔Ο͕-n͕Ο͔+q͔Α͕+q͕Α͔-ç͔Γ͕-ç͕Γ͔+í͔Ε͕+í͕Ε͔+ï͔Ζ͕+ï͕Ζ͔)e0135",
    ///     "+(+a͔Κ͕+a͕Κ͔-d͔Η͕-d͕Η͔-f͔Ρ͕-f͕Ρ͔-g͔Τ͕-g͕Τ͔+i͔Φ͕+i͕Φ͔+k͔Μ͕+k͕Μ͔+m͔Ν͕+m͕Ν͔-o͔Ο͕-o͕Ο͔+r͔Α͕+r͕Α͔-é͔Γ͕-é͕Γ͔-í͔Δ͕-í͕Δ͔+ñ͔Ζ͕+ñ͕Ζ͔)e0136",
    ///     "+(+a͔Λ͕+a͕Λ͔-e͔Η͕-e͕Η͔-f͔Σ͕-f͕Σ͔-g͔Υ͕-g͕Υ͔-h͔Φ͕-h͕Φ͔+l͔Μ͕+l͕Μ͔+n͔Ν͕+n͕Ν͔+o͔Ξ͕+o͕Ξ͔+s͔Α͕+s͕Α͔-ë͔Γ͕-ë͕Γ͔-ï͔Δ͕-ï͕Δ͔-ñ͔Ε͕-ñ͕Ε͔)e0137",
    ///     "+(+b͔Ι͕+b͕Ι͔-c͔Θ͕-c͕Θ͔+f͔Ν͕+f͕Ν͔-g͔Μ͕-g͕Μ͔+k͔Τ͕+k͕Τ͔+l͔Υ͕+l͕Υ͔-m͔Ρ͕-m͕Ρ͔-n͔Σ͕-n͕Σ͔+t͔Α͕+t͕Α͔+ç͔Β͕+ç͕Β͔+ó͔Ε͕+ó͕Ε͔+ö͔Ζ͕+ö͕Ζ͔)e0145",
    ///     "+(+b͔Κ͕+b͕Κ͔-d͔Θ͕-d͕Θ͔+f͔Ξ͕+f͕Ξ͔-h͔Μ͕-h͕Μ͔-j͔Τ͕-j͕Τ͔+l͔Φ͕+l͕Φ͔+m͔Π͕+m͕Π͔-o͔Σ͕-o͕Σ͔+u͔Α͕+u͕Α͔+é͔Β͕+é͕Β͔-ó͔Δ͕-ó͕Δ͔+ú͔Ζ͕+ú͕Ζ͔)e0146",
    ///     "+(+b͔Λ͕+b͕Λ͔-e͔Θ͕-e͕Θ͔+f͔Ο͕+f͕Ο͔-i͔Μ͕-i͕Μ͔-j͔Υ͕-j͕Υ͔-k͔Φ͕-k͕Φ͔+n͔Π͕+n͕Π͔+o͔Ρ͕+o͕Ρ͔+á͔Α͕+á͕Α͔+ë͔Β͕+ë͕Β͔-ö͔Δ͕-ö͕Δ͔-ú͔Ε͕-ú͕Ε͔)e0147",
    ///     "+(+c͔Κ͕+c͕Κ͔-d͔Ι͕-d͕Ι͔+g͔Ξ͕+g͕Ξ͔-h͔Ν͕-h͕Ν͔+j͔Ρ͕+j͕Ρ͔-k͔Π͕-k͕Π͔+n͔Φ͕+n͕Φ͔-o͔Υ͕-o͕Υ͔+ä͔Α͕+ä͕Α͔+í͔Β͕+í͕Β͔+ó͔Γ͕+ó͕Γ͔+ü͔Ζ͕+ü͕Ζ͔)e0156",
    ///     "+(+c͔Λ͕+c͕Λ͔-e͔Ι͕-e͕Ι͔+g͔Ο͕+g͕Ο͔-i͔Ν͕-i͕Ν͔+j͔Σ͕+j͕Σ͔-l͔Π͕-l͕Π͔-m͔Φ͕-m͕Φ͔+o͔Τ͕+o͕Τ͔+å͔Α͕+å͕Α͔+ï͔Β͕+ï͕Β͔+ö͔Γ͕+ö͕Γ͔-ü͔Ε͕-ü͕Ε͔)e0157",
    ///     "+(+d͔Λ͕+d͕Λ͔-e͔Κ͕-e͕Κ͔+h͔Ο͕+h͕Ο͔-i͔Ξ͕-i͕Ξ͔+k͔Σ͕+k͕Σ͔-l͔Ρ͕-l͕Ρ͔+m͔Υ͕+m͕Υ͔-n͔Τ͕-n͕Τ͔+æ͔Α͕+æ͕Α͔+ñ͔Β͕+ñ͕Β͔+ú͔Γ͕+ú͕Γ͔+ü͔Δ͕+ü͕Δ͔)e0167",
    ///     "+(-a͔Γ͕-a͕Γ͔+b͔Β͕+b͕Β͔-f͔Α͕-f͕Α͔+q͔Π͕+q͕Π͔+r͔Ρ͕+r͕Ρ͔+s͔Σ͕+s͕Σ͔-t͔Ν͕-t͕Ν͔-u͔Ξ͕-u͕Ξ͔-á͔Ο͕-á͕Ο͔+ç͔Ι͕+ç͕Ι͔+é͔Κ͕+é͕Κ͔+ë͔Λ͕+ë͕Λ͔)e0234",
    ///     "+(-a͔Δ͕-a͕Δ͔+c͔Β͕+c͕Β͔-g͔Α͕-g͕Α͔-p͔Π͕-p͕Π͔+r͔Τ͕+r͕Τ͔+s͔Υ͕+s͕Υ͔+t͔Μ͕+t͕Μ͔-ä͔Ξ͕-ä͕Ξ͔-å͔Ο͕-å͕Ο͔-ç͔Θ͕-ç͕Θ͔+í͔Κ͕+í͕Κ͔+ï͔Λ͕+ï͕Λ͔)e0235",
    ///     "+(-a͔Ε͕-a͕Ε͔+d͔Β͕+d͕Β͔-h͔Α͕-h͕Α͔-p͔Ρ͕-p͕Ρ͔-q͔Τ͕-q͕Τ͔+s͔Φ͕+s͕Φ͔+u͔Μ͕+u͕Μ͔+ä͔Ν͕+ä͕Ν͔-æ͔Ο͕-æ͕Ο͔-é͔Θ͕-é͕Θ͔-í͔Ι͕-í͕Ι͔+ñ͔Λ͕+ñ͕Λ͔)e0236",
    ///     "+(-a͔Ζ͕-a͕Ζ͔+e͔Β͕+e͕Β͔-i͔Α͕-i͕Α͔-p͔Σ͕-p͕Σ͔-q͔Υ͕-q͕Υ͔-r͔Φ͕-r͕Φ͔+á͔Μ͕+á͕Μ͔+å͔Ν͕+å͕Ν͔+æ͔Ξ͕+æ͕Ξ͔-ë͔Θ͕-ë͕Θ͔-ï͔Ι͕-ï͕Ι͔-ñ͔Κ͕-ñ͕Κ͔)e0237",
    ///     "+(-b͔Δ͕-b͕Δ͔+c͔Γ͕+c͕Γ͔-j͔Α͕-j͕Α͔+p͔Ν͕+p͕Ν͔-q͔Μ͕-q͕Μ͔+u͔Τ͕+u͕Τ͔+á͔Υ͕+á͕Υ͔-ä͔Ρ͕-ä͕Ρ͔-å͔Σ͕-å͕Σ͔+ç͔Η͕+ç͕Η͔+ó͔Κ͕+ó͕Κ͔+ö͔Λ͕+ö͕Λ͔)e0245",
    ///     "+(-b͔Ε͕-b͕Ε͔+d͔Γ͕+d͕Γ͔-k͔Α͕-k͕Α͔+p͔Ξ͕+p͕Ξ͔-r͔Μ͕-r͕Μ͔-t͔Τ͕-t͕Τ͔+á͔Φ͕+á͕Φ͔+ä͔Π͕+ä͕Π͔-æ͔Σ͕-æ͕Σ͔+é͔Η͕+é͕Η͔-ó͔Ι͕-ó͕Ι͔+ú͔Λ͕+ú͕Λ͔)e0246",
    ///     "+(-b͔Ζ͕-b͕Ζ͔+e͔Γ͕+e͕Γ͔-l͔Α͕-l͕Α͔+p͔Ο͕+p͕Ο͔-s͔Μ͕-s͕Μ͔-t͔Υ͕-t͕Υ͔-u͔Φ͕-u͕Φ͔+å͔Π͕+å͕Π͔+æ͔Ρ͕+æ͕Ρ͔+ë͔Η͕+ë͕Η͔-ö͔Ι͕-ö͕Ι͔-ú͔Κ͕-ú͕Κ͔)e0247",
    ///     "+(-c͔Ε͕-c͕Ε͔+d͔Δ͕+d͕Δ͔-m͔Α͕-m͕Α͔+q͔Ξ͕+q͕Ξ͔-r͔Ν͕-r͕Ν͔+t͔Ρ͕+t͕Ρ͔-u͔Π͕-u͕Π͔+å͔Φ͕+å͕Φ͔-æ͔Υ͕-æ͕Υ͔+í͔Η͕+í͕Η͔+ó͔Θ͕+ó͕Θ͔+ü͔Λ͕+ü͕Λ͔)e0256",
    ///     "+(-c͔Ζ͕-c͕Ζ͔+e͔Δ͕+e͕Δ͔-n͔Α͕-n͕Α͔+q͔Ο͕+q͕Ο͔-s͔Ν͕-s͕Ν͔+t͔Σ͕+t͕Σ͔-á͔Π͕-á͕Π͔-ä͔Φ͕-ä͕Φ͔+æ͔Τ͕+æ͕Τ͔+ï͔Η͕+ï͕Η͔+ö͔Θ͕+ö͕Θ͔-ü͔Κ͕-ü͕Κ͔)e0257",
    ///     "+(-d͔Ζ͕-d͕Ζ͔+e͔Ε͕+e͕Ε͔-o͔Α͕-o͕Α͔+r͔Ο͕+r͕Ο͔-s͔Ξ͕-s͕Ξ͔+u͔Σ͕+u͕Σ͔-á͔Ρ͕-á͕Ρ͔+ä͔Υ͕+ä͕Υ͔-å͔Τ͕-å͕Τ͔+ñ͔Η͕+ñ͕Η͔+ú͔Θ͕+ú͕Θ͔+ü͔Ι͕+ü͕Ι͔)e0267",
    ///     "+(-f͔Δ͕-f͕Δ͔+g͔Γ͕+g͕Γ͔-j͔Β͕-j͕Β͔-p͔Ι͕-p͕Ι͔+q͔Θ͕+q͕Θ͔-t͔Η͕-t͕Η͔+é͔Τ͕+é͕Τ͔+ë͔Υ͕+ë͕Υ͔-í͔Ρ͕-í͕Ρ͔-ï͔Σ͕-ï͕Σ͔+ó͔Ξ͕+ó͕Ξ͔+ö͔Ο͕+ö͕Ο͔)e0345",
    ///     "+(-f͔Ε͕-f͕Ε͔+h͔Γ͕+h͕Γ͔-k͔Β͕-k͕Β͔-p͔Κ͕-p͕Κ͔+r͔Θ͕+r͕Θ͔-u͔Η͕-u͕Η͔-ç͔Τ͕-ç͕Τ͔+ë͔Φ͕+ë͕Φ͔+í͔Π͕+í͕Π͔-ñ͔Σ͕-ñ͕Σ͔-ó͔Ν͕-ó͕Ν͔+ú͔Ο͕+ú͕Ο͔)e0346",
    ///     "+(-f͔Ζ͕-f͕Ζ͔+i͔Γ͕+i͕Γ͔-l͔Β͕-l͕Β͔-p͔Λ͕-p͕Λ͔+s͔Θ͕+s͕Θ͔-á͔Η͕-á͕Η͔-ç͔Υ͕-ç͕Υ͔-é͔Φ͕-é͕Φ͔+ï͔Π͕+ï͕Π͔+ñ͔Ρ͕+ñ͕Ρ͔-ö͔Ν͕-ö͕Ν͔-ú͔Ξ͕-ú͕Ξ͔)e0347",
    ///     "+(-g͔Ε͕-g͕Ε͔+h͔Δ͕+h͕Δ͔-m͔Β͕-m͕Β͔-q͔Κ͕-q͕Κ͔+r͔Ι͕+r͕Ι͔-ä͔Η͕-ä͕Η͔+ç͔Ρ͕+ç͕Ρ͔-é͔Π͕-é͕Π͔+ï͔Φ͕+ï͕Φ͔-ñ͔Υ͕-ñ͕Υ͔+ó͔Μ͕+ó͕Μ͔+ü͔Ο͕+ü͕Ο͔)e0356",
    ///     "+(-g͔Ζ͕-g͕Ζ͔+i͔Δ͕+i͕Δ͔-n͔Β͕-n͕Β͔-q͔Λ͕-q͕Λ͔+s͔Ι͕+s͕Ι͔-å͔Η͕-å͕Η͔+ç͔Σ͕+ç͕Σ͔-ë͔Π͕-ë͕Π͔-í͔Φ͕-í͕Φ͔+ñ͔Τ͕+ñ͕Τ͔+ö͔Μ͕+ö͕Μ͔-ü͔Ξ͕-ü͕Ξ͔)e0357",
    ///     "+(-h͔Ζ͕-h͕Ζ͔+i͔Ε͕+i͕Ε͔-o͔Β͕-o͕Β͔-r͔Λ͕-r͕Λ͔+s͔Κ͕+s͕Κ͔-æ͔Η͕-æ͕Η͔+é͔Σ͕+é͕Σ͔-ë͔Ρ͕-ë͕Ρ͔+í͔Υ͕+í͕Υ͔-ï͔Τ͕-ï͕Τ͔+ú͔Μ͕+ú͕Μ͔+ü͔Ν͕+ü͕Ν͔)e0367",
    ///     "+(-j͔Ε͕-j͕Ε͔+k͔Δ͕+k͕Δ͔-m͔Γ͕-m͕Γ͔-t͔Κ͕-t͕Κ͔+u͔Ι͕+u͕Ι͔-ä͔Θ͕-ä͕Θ͔-ç͔Ξ͕-ç͕Ξ͔+é͔Ν͕+é͕Ν͔-í͔Μ͕-í͕Μ͔+ö͔Φ͕+ö͕Φ͔-ú͔Υ͕-ú͕Υ͔+ü͔Σ͕+ü͕Σ͔)e0456",
    ///     "+(-j͔Ζ͕-j͕Ζ͔+l͔Δ͕+l͕Δ͔-n͔Γ͕-n͕Γ͔-t͔Λ͕-t͕Λ͔+á͔Ι͕+á͕Ι͔-å͔Θ͕-å͕Θ͔-ç͔Ο͕-ç͕Ο͔+ë͔Ν͕+ë͕Ν͔-ï͔Μ͕-ï͕Μ͔-ó͔Φ͕-ó͕Φ͔+ú͔Τ͕+ú͕Τ͔-ü͔Ρ͕-ü͕Ρ͔)e0457",
    ///     "+(-k͔Ζ͕-k͕Ζ͔+l͔Ε͕+l͕Ε͔-o͔Γ͕-o͕Γ͔-u͔Λ͕-u͕Λ͔+á͔Κ͕+á͕Κ͔-æ͔Θ͕-æ͕Θ͔-é͔Ο͕-é͕Ο͔+ë͔Ξ͕+ë͕Ξ͔-ñ͔Μ͕-ñ͕Μ͔+ó͔Υ͕+ó͕Υ͔-ö͔Τ͕-ö͕Τ͔+ü͔Π͕+ü͕Π͔)e0467",
    ///     "+(-m͔Ζ͕-m͕Ζ͔+n͔Ε͕+n͕Ε͔-o͔Δ͕-o͕Δ͔-ä͔Λ͕-ä͕Λ͔+å͔Κ͕+å͕Κ͔-æ͔Ι͕-æ͕Ι͔-í͔Ο͕-í͕Ο͔+ï͔Ξ͕+ï͕Ξ͔-ñ͔Ν͕-ñ͕Ν͔-ó͔Σ͕-ó͕Σ͔+ö͔Ρ͕+ö͕Ρ͔-ú͔Π͕-ú͕Π͔)e0567",
    ///     "+(+j͔o͕+j͕o͔-k͔n͕-k͕n͔+l͔m͕+l͕m͔+t͔æ͕+t͕æ͔-u͔å͕-u͕å͔+á͔ä͕+á͕ä͔+ç͔ñ͕+ç͕ñ͔-é͔ï͕-é͕ï͔+ë͔í͕+ë͕í͔)e4567",
    ///     "+(-g͔o͕-g͕o͔+h͔n͕+h͕n͔-i͔m͕-i͕m͔-q͔æ͕-q͕æ͔+r͔å͕+r͕å͔-s͔ä͕-s͕ä͔+ç͔ú͕+ç͕ú͔-é͔ö͕-é͕ö͔+ë͔ó͕+ë͕ó͔)e3576",
    ///     "+(+f͔o͕+f͕o͔-h͔l͕-h͕l͔+i͔k͕+i͕k͔+p͔æ͕+p͕æ͔-r͔á͕-r͕á͔+s͔u͕+s͕u͔+ç͔ü͕+ç͕ü͔-í͔ö͕-í͕ö͔+ï͔ó͕+ï͕ó͔)e3467",
    ///     "+(-f͔n͕-f͕n͔+g͔l͕+g͕l͔-i͔j͕-i͕j͔-p͔å͕-p͕å͔+q͔á͕+q͕á͔-s͔t͕-s͕t͔+é͔ü͕+é͕ü͔-í͔ú͕-í͕ú͔+ñ͔ó͕+ñ͕ó͔)e3475",
    ///     "+(+f͔m͕+f͕m͔-g͔k͕-g͕k͔+h͔j͕+h͕j͔+p͔ä͕+p͕ä͔-q͔u͕-q͕u͔+r͔t͕+r͕t͔+ë͔ü͕+ë͕ü͔-ï͔ú͕-ï͕ú͔+ñ͔ö͕+ñ͕ö͔)e3456",
    ///     "+(+c͔o͕+c͕o͔-d͔n͕-d͕n͔+e͔m͕+e͕m͔-q͔ñ͕-q͕ñ͔+r͔ï͕+r͕ï͔-s͔í͕-s͕í͔-t͔ú͕-t͕ú͔+u͔ö͕+u͕ö͔-á͔ó͕-á͕ó͔)e2567",
    ///     "+(-b͔o͕-b͕o͔+d͔l͕+d͕l͔-e͔k͕-e͕k͔+p͔ñ͕+p͕ñ͔-r͔ë͕-r͕ë͔+s͔é͕+s͕é͔-t͔ü͕-t͕ü͔+ä͔ö͕+ä͕ö͔-å͔ó͕-å͕ó͔)e2476",
    ///     "+(+b͔n͕+b͕n͔-c͔l͕-c͕l͔+e͔j͕+e͕j͔-p͔ï͕-p͕ï͔+q͔ë͕+q͕ë͔-s͔ç͕-s͕ç͔-u͔ü͕-u͕ü͔+ä͔ú͕+ä͕ú͔-æ͔ó͕-æ͕ó͔)e2457",
    ///     "+(-b͔m͕-b͕m͔+c͔k͕+c͕k͔-d͔j͕-d͕j͔+p͔í͕+p͕í͔-q͔é͕-q͕é͔+r͔ç͕+r͕ç͔-á͔ü͕-á͕ü͔+å͔ú͕+å͕ú͔-æ͔ö͕-æ͕ö͔)e2465",
    ///     "+(+a͔o͕+a͕o͔-d͔i͕-d͕i͔+e͔h͕+e͕h͔+p͔ú͕+p͕ú͔+q͔ü͕+q͕ü͔-u͔ë͕-u͕ë͔+á͔é͕+á͕é͔-ä͔ï͕-ä͕ï͔+å͔í͕+å͕í͔)e2367",
    ///     "+(-a͔n͕-a͕n͔+c͔i͕+c͕i͔-e͔g͕-e͕g͔-p͔ö͕-p͕ö͔+r͔ü͕+r͕ü͔+t͔ë͕+t͕ë͔-á͔ç͕-á͕ç͔-ä͔ñ͕-ä͕ñ͔+æ͔í͕+æ͕í͔)e2375",
    ///     "+(+a͔m͕+a͕m͔-c͔h͕-c͕h͔+d͔g͕+d͕g͔+p͔ó͕+p͕ó͔+s͔ü͕+s͕ü͔-t͔é͕-t͕é͔+u͔ç͕+u͕ç͔-å͔ñ͕-å͕ñ͔+æ͔ï͕+æ͕ï͔)e2356",
    ///     "+(+a͔l͕+a͕l͔-b͔i͕-b͕i͔+e͔f͕+e͕f͔-q͔ö͕-q͕ö͔-r͔ú͕-r͕ú͔+t͔ï͕+t͕ï͔+u͔ñ͕+u͕ñ͔-å͔ç͕-å͕ç͔-æ͔é͕-æ͕é͔)e2347",
    ///     "+(-a͔k͕-a͕k͔+b͔h͕+b͕h͔-d͔f͕-d͕f͔+q͔ó͕+q͕ó͔-s͔ú͕-s͕ú͔-t͔í͕-t͕í͔+á͔ñ͕+á͕ñ͔+ä͔ç͕+ä͕ç͔-æ͔ë͕-æ͕ë͔)e2364",
    ///     "+(+a͔j͕+a͕j͔-b͔g͕-b͕g͔+c͔f͕+c͕f͔+r͔ó͕+r͕ó͔+s͔ö͕+s͕ö͔-u͔í͕-u͕í͔-á͔ï͕-á͕ï͔+ä͔é͕+ä͕é͔+å͔ë͕+å͕ë͔)e2345",
    ///     "+(+c͔æ͕+c͕æ͔-d͔å͕-d͕å͔+e͔ä͕+e͕ä͔+g͔ñ͕+g͕ñ͔-h͔ï͕-h͕ï͔+i͔í͕+i͕í͔+j͔ú͕+j͕ú͔-k͔ö͕-k͕ö͔+l͔ó͕+l͕ó͔)e1576",
    ///     "+(-b͔æ͕-b͕æ͔+d͔á͕+d͕á͔-e͔u͕-e͕u͔-f͔ñ͕-f͕ñ͔+h͔ë͕+h͕ë͔-i͔é͕-i͕é͔+j͔ü͕+j͕ü͔-m͔ö͕-m͕ö͔+n͔ó͕+n͕ó͔)e1467",
    ///     "+(+b͔å͕+b͕å͔-c͔á͕-c͕á͔+e͔t͕+e͕t͔+f͔ï͕+f͕ï͔-g͔ë͕-g͕ë͔+i͔ç͕+i͕ç͔+k͔ü͕+k͕ü͔-m͔ú͕-m͕ú͔+o͔ó͕+o͕ó͔)e1475",
    ///     "+(-b͔ä͕-b͕ä͔+c͔u͕+c͕u͔-d͔t͕-d͕t͔-f͔í͕-f͕í͔+g͔é͕+g͕é͔-h͔ç͕-h͕ç͔+l͔ü͕+l͕ü͔-n͔ú͕-n͕ú͔+o͔ö͕+o͕ö͔)e1456",
    ///     "+(+a͔æ͕+a͕æ͔-d͔s͕-d͕s͔+e͔r͕+e͕r͔-f͔ú͕-f͕ú͔-g͔ü͕-g͕ü͔+k͔ë͕+k͕ë͔-l͔é͕-l͕é͔+m͔ï͕+m͕ï͔-n͔í͕-n͕í͔)e1376",
    ///     "+(-a͔å͕-a͕å͔+c͔s͕+c͕s͔-e͔q͕-e͕q͔+f͔ö͕+f͕ö͔-h͔ü͕-h͕ü͔-j͔ë͕-j͕ë͔+l͔ç͕+l͕ç͔+m͔ñ͕+m͕ñ͔-o͔í͕-o͕í͔)e1357",
    ///     "+(+a͔ä͕+a͕ä͔-c͔r͕-c͕r͔+d͔q͕+d͕q͔-f͔ó͕-f͕ó͔-i͔ü͕-i͕ü͔+j͔é͕+j͕é͔-k͔ç͕-k͕ç͔+n͔ñ͕+n͕ñ͔-o͔ï͕-o͕ï͔)e1365",
    ///     "+(+a͔á͕+a͕á͔-b͔s͕-b͕s͔+e͔p͕+e͕p͔+g͔ö͕+g͕ö͔+h͔ú͕+h͕ú͔-j͔ï͕-j͕ï͔-k͔ñ͕-k͕ñ͔+n͔ç͕+n͕ç͔+o͔é͕+o͕é͔)e1374",
    ///     "+(-a͔u͕-a͕u͔+b͔r͕+b͕r͔-d͔p͕-d͕p͔-g͔ó͕-g͕ó͔+i͔ú͕+i͕ú͔+j͔í͕+j͕í͔-l͔ñ͕-l͕ñ͔-m͔ç͕-m͕ç͔+o͔ë͕+o͕ë͔)e1346",
    ///     "+(+a͔t͕+a͕t͔-b͔q͕-b͕q͔+c͔p͕+c͕p͔-h͔ó͕-h͕ó͔-i͔ö͕-i͕ö͔+k͔í͕+k͕í͔+l͔ï͕+l͕ï͔-m͔é͕-m͕é͔-n͔ë͕-n͕ë͔)e1354",
    ///     "+(+a͔ñ͕+a͕ñ͔+b͔ú͕+b͕ú͔+c͔ü͕+c͕ü͔-h͔s͕-h͕s͔+i͔r͕+i͕r͔-k͔á͕-k͕á͔+l͔u͕+l͕u͔-m͔å͕-m͕å͔+n͔ä͕+n͕ä͔)e1267",
    ///     "+(-a͔ï͕-a͕ï͔-b͔ö͕-b͕ö͔+d͔ü͕+d͕ü͔+g͔s͕+g͕s͔-i͔q͕-i͕q͔+j͔á͕+j͕á͔-l͔t͕-l͕t͔-m͔æ͕-m͕æ͔+o͔ä͕+o͕ä͔)e1275",
    ///     "+(+a͔í͕+a͕í͔+b͔ó͕+b͕ó͔+e͔ü͕+e͕ü͔-g͔r͕-g͕r͔+h͔q͕+h͕q͔-j͔u͕-j͕u͔+k͔t͕+k͕t͔-n͔æ͕-n͕æ͔+o͔å͕+o͕å͔)e1256",
    ///     "+(+a͔ë͕+a͕ë͔-c͔ö͕-c͕ö͔-d͔ú͕-d͕ú͔-f͔s͕-f͕s͔+i͔p͕+i͕p͔+j͔å͕+j͕å͔+k͔æ͕+k͕æ͔-n͔t͕-n͕t͔-o͔u͕-o͕u͔)e1247",
    ///     "+(-a͔é͕-a͕é͔+c͔ó͕+c͕ó͔-e͔ú͕-e͕ú͔+f͔r͕+f͕r͔-h͔p͕-h͕p͔-j͔ä͕-j͕ä͔+l͔æ͕+l͕æ͔+m͔t͕+m͕t͔-o͔á͕-o͕á͔)e1264",
    ///     "+(+a͔ç͕+a͕ç͔+d͔ó͕+d͕ó͔+e͔ö͕+e͕ö͔-f͔q͕-f͕q͔+g͔p͕+g͕p͔-k͔ä͕-k͕ä͔-l͔å͕-l͕å͔+m͔u͕+m͕u͔+n͔á͕+n͕á͔)e1245",
    ///     "+(+b͔ë͕+b͕ë͔+c͔ï͕+c͕ï͔+d͔ñ͕+d͕ñ͔-f͔á͕-f͕á͔-g͔å͕-g͕å͔-h͔æ͕-h͕æ͔+l͔p͕+l͕p͔+n͔q͕+n͕q͔+o͔r͕+o͕r͔)e1273",
    ///     "+(-b͔é͕-b͕é͔-c͔í͕-c͕í͔+e͔ñ͕+e͕ñ͔+f͔u͕+f͕u͔+g͔ä͕+g͕ä͔-i͔æ͕-i͕æ͔-k͔p͕-k͕p͔-m͔q͕-m͕q͔+o͔s͕+o͕s͔)e1236",
    ///     "+(+b͔ç͕+b͕ç͔-d͔í͕-d͕í͔-e͔ï͕-e͕ï͔-f͔t͕-f͕t͔+h͔ä͕+h͕ä͔+i͔å͕+i͕å͔+j͔p͕+j͕p͔-m͔r͕-m͕r͔-n͔s͕-n͕s͔)e1253",
    ///     "+(+p͔ü͕-p͕ü͔-q͔ú͕+q͕ú͔+r͔ö͕-r͕ö͔-s͔ó͕+s͕ó͔+t͔ñ͕-t͕ñ͔-u͔ï͕+u͕ï͔+á͔í͕-á͕í͔+ä͔ë͕-ä͕ë͔-å͔é͕+å͕é͔+æ͔ç͕-æ͕ç͔)e234567",
    ///     "+(-f͔ü͕+f͕ü͔+g͔ú͕-g͕ú͔-h͔ö͕+h͕ö͔+i͔ó͕-i͕ó͔-j͔ñ͕+j͕ñ͔+k͔ï͕-k͕ï͔-l͔í͕+l͕í͔-m͔ë͕+m͕ë͔+n͔é͕-n͕é͔-o͔ç͕+o͕ç͔)e134576",
    ///     "+(+b͔ü͕-b͕ü͔-c͔ú͕+c͕ú͔+d͔ö͕-d͕ö͔-e͔ó͕+e͕ó͔+j͔æ͕-j͕æ͔-k͔å͕+k͕å͔+l͔ä͕-l͕ä͔+m͔á͕-m͕á͔-n͔u͕+n͕u͔+o͔t͕-o͕t͔)e124567",
    ///     "+(-a͔ü͕+a͕ü͔+c͔ñ͕-c͕ñ͔-d͔ï͕+d͕ï͔+e͔í͕-e͕í͔-g͔æ͕+g͕æ͔+h͔å͕-h͕å͔-i͔ä͕+i͕ä͔-m͔s͕+m͕s͔+n͔r͕-n͕r͔-o͔q͕+o͕q͔)e123576",
    ///     "+(+a͔ú͕-a͕ú͔-b͔ñ͕+b͕ñ͔+d͔ë͕-d͕ë͔-e͔é͕+e͕é͔+f͔æ͕-f͕æ͔-h͔á͕+h͕á͔+i͔u͕-i͕u͔+k͔s͕-k͕s͔-l͔r͕+l͕r͔+o͔p͕-o͕p͔)e123467",
    ///     "+(-a͔ö͕+a͕ö͔+b͔ï͕-b͕ï͔-c͔ë͕+c͕ë͔+e͔ç͕-e͕ç͔-f͔å͕+f͕å͔+g͔á͕-g͕á͔-i͔t͕+i͕t͔-j͔s͕+j͕s͔+l͔q͕-l͕q͔-n͔p͕+n͕p͔)e123475",
    ///     "+(+a͔ó͕-a͕ó͔-b͔í͕+b͕í͔+c͔é͕-c͕é͔-d͔ç͕+d͕ç͔+f͔ä͕-f͕ä͔-g͔u͕+g͕u͔+h͔t͕-h͕t͔+j͔r͕-j͕r͔-k͔q͕+k͕q͔+m͔p͕-m͕p͔)e123456",
    ///     "+(-ç͔Φ͕+ç͕Φ͔+é͔Υ͕-é͕Υ͔-ë͔Τ͕+ë͕Τ͔-í͔Σ͕+í͕Σ͔+ï͔Ρ͕-ï͕Ρ͔-ñ͔Π͕+ñ͕Π͔+ó͔Ο͕-ó͕Ο͔-ö͔Ξ͕+ö͕Ξ͔+ú͔Ν͕-ú͕Ν͔-ü͔Μ͕+ü͕Μ͔)e034567",
    ///     "+(+t͔Φ͕-t͕Φ͔-u͔Υ͕+u͕Υ͔+á͔Τ͕-á͕Τ͔+ä͔Σ͕-ä͕Σ͔-å͔Ρ͕+å͕Ρ͔+æ͔Π͕-æ͕Π͔-ó͔Λ͕+ó͕Λ͔+ö͔Κ͕-ö͕Κ͔-ú͔Ι͕+ú͕Ι͔+ü͔Θ͕-ü͕Θ͔)e024576",
    ///     "+(-q͔Φ͕+q͕Φ͔+r͔Υ͕-r͕Υ͔-s͔Τ͕+s͕Τ͔-ä͔Ο͕+ä͕Ο͔+å͔Ξ͕-å͕Ξ͔-æ͔Ν͕+æ͕Ν͔+í͔Λ͕-í͕Λ͔-ï͔Κ͕+ï͕Κ͔+ñ͔Ι͕-ñ͕Ι͔-ü͔Η͕+ü͕Η͔)e023567",
    ///     "+(+p͔Φ͕-p͕Φ͔-r͔Σ͕+r͕Σ͔+s͔Ρ͕-s͕Ρ͔+u͔Ο͕-u͕Ο͔-á͔Ξ͕+á͕Ξ͔+æ͔Μ͕-æ͕Μ͔-é͔Λ͕+é͕Λ͔+ë͔Κ͕-ë͕Κ͔-ñ͔Θ͕+ñ͕Θ͔+ú͔Η͕-ú͕Η͔)e023476",
    ///     "+(-p͔Υ͕+p͕Υ͔+q͔Σ͕-q͕Σ͔-s͔Π͕+s͕Π͔-t͔Ο͕+t͕Ο͔+á͔Ν͕-á͕Ν͔-å͔Μ͕+å͕Μ͔+ç͔Λ͕-ç͕Λ͔-ë͔Ι͕+ë͕Ι͔+ï͔Θ͕-ï͕Θ͔-ö͔Η͕+ö͕Η͔)e023457",
    ///     "+(+p͔Τ͕-p͕Τ͔-q͔Ρ͕+q͕Ρ͔+r͔Π͕-r͕Π͔+t͔Ξ͕-t͕Ξ͔-u͔Ν͕+u͕Ν͔+ä͔Μ͕-ä͕Μ͔-ç͔Κ͕+ç͕Κ͔+é͔Ι͕-é͕Ι͔-í͔Θ͕+í͕Θ͔+ó͔Η͕-ó͕Η͔)e023465",
    ///     "+(-j͔Φ͕+j͕Φ͔+k͔Υ͕-k͕Υ͔-l͔Τ͕+l͕Τ͔-m͔Σ͕+m͕Σ͔+n͔Ρ͕-n͕Ρ͔-o͔Π͕+o͕Π͔+ó͔Ζ͕-ó͕Ζ͔-ö͔Ε͕+ö͕Ε͔+ú͔Δ͕-ú͕Δ͔-ü͔Γ͕+ü͕Γ͔)e014567",
    ///     "+(+g͔Φ͕-g͕Φ͔-h͔Υ͕+h͕Υ͔+i͔Τ͕-i͕Τ͔+m͔Ο͕-m͕Ο͔-n͔Ξ͕+n͕Ξ͔+o͔Ν͕-o͕Ν͔-í͔Ζ͕+í͕Ζ͔+ï͔Ε͕-ï͕Ε͔-ñ͔Δ͕+ñ͕Δ͔+ü͔Β͕-ü͕Β͔)e013576",
    ///     "+(-f͔Φ͕+f͕Φ͔+h͔Σ͕-h͕Σ͔-i͔Ρ͕+i͕Ρ͔-k͔Ο͕+k͕Ο͔+l͔Ξ͕-l͕Ξ͔-o͔Μ͕+o͕Μ͔+é͔Ζ͕-é͕Ζ͔-ë͔Ε͕+ë͕Ε͔+ñ͔Γ͕-ñ͕Γ͔-ú͔Β͕+ú͕Β͔)e013467",
    ///     "+(+f͔Υ͕-f͕Υ͔-g͔Σ͕+g͕Σ͔+i͔Π͕-i͕Π͔+j͔Ο͕-j͕Ο͔-l͔Ν͕+l͕Ν͔+n͔Μ͕-n͕Μ͔-ç͔Ζ͕+ç͕Ζ͔+ë͔Δ͕-ë͕Δ͔-ï͔Γ͕+ï͕Γ͔+ö͔Β͕-ö͕Β͔)e013475",
    ///     "+(-f͔Τ͕+f͕Τ͔+g͔Ρ͕-g͕Ρ͔-h͔Π͕+h͕Π͔-j͔Ξ͕+j͕Ξ͔+k͔Ν͕-k͕Ν͔-m͔Μ͕+m͕Μ͔+ç͔Ε͕-ç͕Ε͔-é͔Δ͕+é͕Δ͔+í͔Γ͕-í͕Γ͔-ó͔Β͕+ó͕Β͔)e013456",
    ///     "+(-c͔Φ͕+c͕Φ͔+d͔Υ͕-d͕Υ͔-e͔Τ͕+e͕Τ͔-m͔Λ͕+m͕Λ͔+n͔Κ͕-n͕Κ͔-o͔Ι͕+o͕Ι͔+ä͔Ζ͕-ä͕Ζ͔-å͔Ε͕+å͕Ε͔+æ͔Δ͕-æ͕Δ͔-ü͔Α͕+ü͕Α͔)e012567",
    ///     "+(+b͔Φ͕-b͕Φ͔-d͔Σ͕+d͕Σ͔+e͔Ρ͕-e͕Ρ͔+k͔Λ͕-k͕Λ͔-l͔Κ͕+l͕Κ͔+o͔Θ͕-o͕Θ͔-u͔Ζ͕+u͕Ζ͔+á͔Ε͕-á͕Ε͔-æ͔Γ͕+æ͕Γ͔+ú͔Α͕-ú͕Α͔)e012476",
    ///     "+(-b͔Υ͕+b͕Υ͔+c͔Σ͕-c͕Σ͔-e͔Π͕+e͕Π͔-j͔Λ͕+j͕Λ͔+l͔Ι͕-l͕Ι͔-n͔Θ͕+n͕Θ͔+t͔Ζ͕-t͕Ζ͔-á͔Δ͕+á͕Δ͔+å͔Γ͕-å͕Γ͔-ö͔Α͕+ö͕Α͔)e012457",
    ///     "+(+b͔Τ͕-b͕Τ͔-c͔Ρ͕+c͕Ρ͔+d͔Π͕-d͕Π͔+j͔Κ͕-j͕Κ͔-k͔Ι͕+k͕Ι͔+m͔Θ͕-m͕Θ͔-t͔Ε͕+t͕Ε͔+u͔Δ͕-u͕Δ͔-ä͔Γ͕+ä͕Γ͔+ó͔Α͕-ó͕Α͔)e012465",
    ///     "+(-a͔Φ͕+a͕Φ͔+d͔Ο͕-d͕Ο͔-e͔Ξ͕+e͕Ξ͔-h͔Λ͕+h͕Λ͔+i͔Κ͕-i͕Κ͔-o͔Η͕+o͕Η͔+r͔Ζ͕-r͕Ζ͔-s͔Ε͕+s͕Ε͔+æ͔Β͕-æ͕Β͔-ñ͔Α͕+ñ͕Α͔)e012367",
    ///     "+(+a͔Υ͕-a͕Υ͔-c͔Ο͕+c͕Ο͔+e͔Ν͕-e͕Ν͔+g͔Λ͕-g͕Λ͔-i͔Ι͕+i͕Ι͔+n͔Η͕-n͕Η͔-q͔Ζ͕+q͕Ζ͔+s͔Δ͕-s͕Δ͔-å͔Β͕+å͕Β͔+ï͔Α͕-ï͕Α͔)e012375",
    ///     "+(-a͔Τ͕+a͕Τ͔+c͔Ξ͕-c͕Ξ͔-d͔Ν͕+d͕Ν͔-g͔Κ͕+g͕Κ͔+h͔Ι͕-h͕Ι͔-m͔Η͕+m͕Η͔+q͔Ε͕-q͕Ε͔-r͔Δ͕+r͕Δ͔+ä͔Β͕-ä͕Β͔-í͔Α͕+í͕Α͔)e012356",
    ///     "+(-a͔Σ͕+a͕Σ͔+b͔Ο͕-b͕Ο͔-e͔Μ͕+e͕Μ͔-f͔Λ͕+f͕Λ͔+i͔Θ͕-i͕Θ͔-l͔Η͕+l͕Η͔+p͔Ζ͕-p͕Ζ͔-s͔Γ͕+s͕Γ͔+á͔Β͕-á͕Β͔-ë͔Α͕+ë͕Α͔)e012347",
    ///     "+(+a͔Ρ͕-a͕Ρ͔-b͔Ξ͕+b͕Ξ͔+d͔Μ͕-d͕Μ͔+f͔Κ͕-f͕Κ͔-h͔Θ͕+h͕Θ͔+k͔Η͕-k͕Η͔-p͔Ε͕+p͕Ε͔+r͔Γ͕-r͕Γ͔-u͔Β͕+u͕Β͔+é͔Α͕-é͕Α͔)e012364",
    ///     "+(-a͔Π͕+a͕Π͔+b͔Ν͕-b͕Ν͔-c͔Μ͕+c͕Μ͔-f͔Ι͕+f͕Ι͔+g͔Θ͕-g͕Θ͔-j͔Η͕+j͕Η͔+p͔Δ͕-p͕Δ͔-q͔Γ͕+q͕Γ͔+t͔Β͕-t͕Β͔-ç͔Α͕+ç͕Α͔)e012345",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_triple_motor() -> Self {
        Self::scalar() + Self::volume5() + Self::volume() + Self::line()
    }
    /// The multivector of triple motor $`m_3 \equiv s + v^5 + v + \ell + S`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let triple_motor = Vee::triple_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(triple_motor.basis_blades(), Vee::triple_motor().basis_blades());
    /// format_eq!(triple_motor, [
    ///     "+v͔v͕",
    ///     "+(+X͕v͔+Y͕α͔+Z͕β͔+Ð͕γ͔+Ø͕δ͔+Þ͕ε͔+Œ͕ζ͔)e01",
    ///     "+(-X͕α͔+Y͕v͔+Z͕η͔+Ð͕θ͔+Ø͕ι͔+Þ͕κ͔+Œ͕λ͔)e02",
    ///     "+(-X͕β͔-Y͕η͔+Z͕v͔+Ð͕μ͔+Ø͕ν͔+Þ͕ξ͔+Œ͕ο͔)e03",
    ///     "+(-X͕γ͔-Y͕θ͔-Z͕μ͔+v͔Ð͕+Ø͕π͔+Þ͕ρ͔+Œ͕σ͔)e04",
    ///     "+(-X͕δ͔-Y͕ι͔-Z͕ν͔+v͔Ø͕-Ð͕π͔+Þ͕τ͔+Œ͕υ͔)e05",
    ///     "+(-X͕ε͔-Y͕κ͔-Z͕ξ͔+v͔Þ͕-Ð͕ρ͔-Ø͕τ͔+Œ͕φ͔)e06",
    ///     "+(-X͕ζ͔-Y͕λ͔-Z͕ο͔+v͔Œ͕-Ð͕σ͔-Ø͕υ͔-Þ͕φ͔)e07",
    ///     "+v͕α͔e12",
    ///     "+v͕β͔e13",
    ///     "+v͕γ͔e14",
    ///     "+v͕δ͔e15",
    ///     "+v͕ε͔e16",
    ///     "+v͕ζ͔e17",
    ///     "+v͕η͔e23",
    ///     "+v͕θ͔e24",
    ///     "+v͕ι͔e25",
    ///     "+v͕κ͔e26",
    ///     "+v͕λ͔e27",
    ///     "+v͕μ͔e34",
    ///     "+v͕ν͔e35",
    ///     "+v͕ξ͔e36",
    ///     "+v͕ο͔e37",
    ///     "+v͕π͔e45",
    ///     "+v͕ρ͔e46",
    ///     "+v͕σ͔e47",
    ///     "+v͕τ͔e56",
    ///     "+v͕υ͔e57",
    ///     "+v͕φ͔e67",
    ///     "+v͕ü͔e1234",
    ///     "+(+X͕η͔-Y͕β͔+Z͕α͔+Ð͕ü͔-Ø͕ú͔+Þ͕ö͔-ó͔Œ͕)e0123",
    ///     "+(+X͕θ͔-Y͕γ͔-Z͕ü͔+Ð͕α͔+Ø͕ñ͔-Þ͕ï͔+í͔Œ͕)e0124",
    ///     "+(+X͕ι͔-Y͕δ͔+Z͕ú͔-Ð͕ñ͔+Ø͕α͔+Þ͕ë͔-é͔Œ͕)e0125",
    ///     "+(+X͕κ͔-Y͕ε͔-Z͕ö͔+Ð͕ï͔-Ø͕ë͔+Þ͕α͔+ç͔Œ͕)e0126",
    ///     "+(+X͕λ͔-Y͕ζ͔+Z͕ó͔-Ð͕í͔+Ø͕é͔-Þ͕ç͔+Œ͕α͔)e0127",
    ///     "+(+X͕μ͔+Y͕ü͔-Z͕γ͔+Ð͕β͔-Ø͕æ͔+Þ͕å͔-ä͔Œ͕)e0134",
    ///     "+(+X͕ν͔-Y͕ú͔-Z͕δ͔+u͔Œ͕+Ð͕æ͔+Ø͕β͔-Þ͕á͔)e0135",
    ///     "+(+X͕ξ͔+Y͕ö͔-Z͕ε͔-t͔Œ͕-Ð͕å͔+Ø͕á͔+Þ͕β͔)e0136",
    ///     "+(+X͕ο͔-Y͕ó͔-Z͕ζ͔+t͔Þ͕-u͔Ø͕+Ð͕ä͔+Œ͕β͔)e0137",
    ///     "+(+X͕π͔+Y͕ñ͔-Z͕æ͔-r͔Œ͕+s͔Þ͕-Ð͕δ͔+Ø͕γ͔)e0145",
    ///     "+(+X͕ρ͔-Y͕ï͔+Z͕å͔+q͔Œ͕-s͔Ø͕-Ð͕ε͔+Þ͕γ͔)e0146",
    ///     "+(+X͕σ͔+Y͕í͔-Z͕ä͔-q͔Þ͕+r͔Ø͕-Ð͕ζ͔+Œ͕γ͔)e0147",
    ///     "+(+X͕τ͔+Y͕ë͔-Z͕á͔-p͔Œ͕+s͔Ð͕-Ø͕ε͔+Þ͕δ͔)e0156",
    ///     "+(+X͕υ͔-Y͕é͔+Z͕u͔+p͔Þ͕-r͔Ð͕-Ø͕ζ͔+Œ͕δ͔)e0157",
    ///     "+(+X͕φ͔+Y͕ç͔-Z͕t͔-p͔Ø͕+q͔Ð͕-Þ͕ζ͔+Œ͕ε͔)e0167",
    ///     "+(-X͕ü͔+Y͕μ͔-Z͕θ͔+m͔Œ͕-n͔Þ͕+o͔Ø͕+Ð͕η͔)e0234",
    ///     "+(+X͕ú͔+Y͕ν͔-Z͕ι͔-k͔Œ͕+l͔Þ͕-o͔Ð͕+Ø͕η͔)e0235",
    ///     "+(-X͕ö͔+Y͕ξ͔-Z͕κ͔+j͔Œ͕-l͔Ø͕+n͔Ð͕+Þ͕η͔)e0236",
    ///     "+(+X͕ó͔+Y͕ο͔-Z͕λ͔-j͔Þ͕+k͔Ø͕-m͔Ð͕+Œ͕η͔)e0237",
    ///     "+(-X͕ñ͔+Y͕π͔+Z͕o͔+h͔Œ͕-i͔Þ͕-Ð͕ι͔+Ø͕θ͔)e0245",
    ///     "+(+X͕ï͔+Y͕ρ͔-Z͕n͔-g͔Œ͕+i͔Ø͕-Ð͕κ͔+Þ͕θ͔)e0246",
    ///     "+(-X͕í͔+Y͕σ͔+Z͕m͔+g͔Þ͕-h͔Ø͕-Ð͕λ͔+Œ͕θ͔)e0247",
    ///     "+(-X͕ë͔+Y͕τ͔+Z͕l͔+f͔Œ͕-i͔Ð͕-Ø͕κ͔+Þ͕ι͔)e0256",
    ///     "+(+X͕é͔+Y͕υ͔-Z͕k͔-f͔Þ͕+h͔Ð͕-Ø͕λ͔+Œ͕ι͔)e0257",
    ///     "+(-X͕ç͔+Y͕φ͔+Z͕j͔+f͔Ø͕-g͔Ð͕-Þ͕λ͔+Œ͕κ͔)e0267",
    ///     "+(+X͕æ͔-Y͕o͔+Z͕π͔-d͔Œ͕+e͔Þ͕-Ð͕ν͔+Ø͕μ͔)e0345",
    ///     "+(-X͕å͔+Y͕n͔+Z͕ρ͔+c͔Œ͕-e͔Ø͕-Ð͕ξ͔+Þ͕μ͔)e0346",
    ///     "+(+X͕ä͔-Y͕m͔+Z͕σ͔-c͔Þ͕+d͔Ø͕-Ð͕ο͔+Œ͕μ͔)e0347",
    ///     "+(+X͕á͔-Y͕l͔+Z͕τ͔-b͔Œ͕+e͔Ð͕-Ø͕ξ͔+Þ͕ν͔)e0356",
    ///     "+(-X͕u͔+Y͕k͔+Z͕υ͔+b͔Þ͕-d͔Ð͕-Ø͕ο͔+Œ͕ν͔)e0357",
    ///     "+(+X͕t͔-Y͕j͔+Z͕φ͔-b͔Ø͕+c͔Ð͕-Þ͕ο͔+Œ͕ξ͔)e0367",
    ///     "+(-X͕s͔+Y͕i͔-Z͕e͔+a͔Œ͕+Ð͕τ͔-Ø͕ρ͔+Þ͕π͔)e0456",
    ///     "+(+X͕r͔-Y͕h͔+Z͕d͔-a͔Þ͕+Ð͕υ͔-Ø͕σ͔+Œ͕π͔)e0457",
    ///     "+(-X͕q͔+Y͕g͔-Z͕c͔+a͔Ø͕+Ð͕φ͔-Þ͕σ͔+Œ͕ρ͔)e0467",
    ///     "+(+X͕p͔-Y͕f͔+Z͕b͔-a͔Ð͕+Ø͕φ͔-Þ͕υ͔+Œ͕τ͔)e0567",
    ///     "+a͔v͕e4567",
    ///     "+b͔v͕e3576",
    ///     "+c͔v͕e3467",
    ///     "+d͔v͕e3475",
    ///     "+e͔v͕e3456",
    ///     "+f͔v͕e2567",
    ///     "+g͔v͕e2476",
    ///     "+h͔v͕e2457",
    ///     "+i͔v͕e2465",
    ///     "+j͔v͕e2367",
    ///     "+k͔v͕e2375",
    ///     "+l͔v͕e2356",
    ///     "+m͔v͕e2347",
    ///     "+n͔v͕e2364",
    ///     "+o͔v͕e2345",
    ///     "+p͔v͕e1576",
    ///     "+q͔v͕e1467",
    ///     "+r͔v͕e1475",
    ///     "+s͔v͕e1456",
    ///     "+t͔v͕e1376",
    ///     "+u͔v͕e1357",
    ///     "+v͕á͔e1365",
    ///     "+v͕ä͔e1374",
    ///     "+v͕å͔e1346",
    ///     "+v͕æ͔e1354",
    ///     "+v͕ç͔e1267",
    ///     "+v͕é͔e1275",
    ///     "+v͕ë͔e1256",
    ///     "+v͕í͔e1247",
    ///     "+v͕ï͔e1264",
    ///     "+v͕ñ͔e1245",
    ///     "+v͕ó͔e1273",
    ///     "+v͕ö͔e1236",
    ///     "+v͕ú͔e1253",
    ///     "+v͕x͔e234567",
    ///     "+v͕y͔e134576",
    ///     "+v͕z͔e124567",
    ///     "+v͕ð͔e123576",
    ///     "+v͕ø͔e123467",
    ///     "+v͕þ͔e123475",
    ///     "+v͕œ͔e123456",
    ///     "+(+X͕y͔-Y͕x͔+Z͕a͔+b͔Ð͕+c͔Ø͕+d͔Þ͕+e͔Œ͕)e034567",
    ///     "+(+X͕z͔-Y͕a͔-Z͕x͔+f͔Ð͕+g͔Ø͕+h͔Þ͕+i͔Œ͕)e024576",
    ///     "+(+X͕ð͔-Y͕b͔-Z͕f͔+j͔Ø͕+k͔Þ͕+l͔Œ͕-x͔Ð͕)e023567",
    ///     "+(+X͕ø͔-Y͕c͔-Z͕g͔-j͔Ð͕+m͔Þ͕+n͔Œ͕-x͔Ø͕)e023476",
    ///     "+(+X͕þ͔-Y͕d͔-Z͕h͔-k͔Ð͕-m͔Ø͕+o͔Œ͕-x͔Þ͕)e023457",
    ///     "+(+X͕œ͔-Y͕e͔-Z͕i͔-l͔Ð͕-n͔Ø͕-o͔Þ͕-x͔Œ͕)e023465",
    ///     "+(+X͕a͔+Y͕z͔-Z͕y͔+p͔Ð͕+q͔Ø͕+r͔Þ͕+s͔Œ͕)e014567",
    ///     "+(+X͕b͔+Y͕ð͔-Z͕p͔+t͔Ø͕+u͔Þ͕-y͔Ð͕+á͔Œ͕)e013576",
    ///     "+(+X͕c͔+Y͕ø͔-Z͕q͔-t͔Ð͕-y͔Ø͕+Þ͕ä͔+å͔Œ͕)e013467",
    ///     "+(+X͕d͔+Y͕þ͔-Z͕r͔-u͔Ð͕-y͔Þ͕-Ø͕ä͔+æ͔Œ͕)e013475",
    ///     "+(+X͕e͔+Y͕œ͔-Z͕s͔-y͔Œ͕-Ð͕á͔-Ø͕å͔-Þ͕æ͔)e013456",
    ///     "+(+X͕f͔+Y͕p͔+Z͕ð͔-z͔Ð͕+Ø͕ç͔+Þ͕é͔+ë͔Œ͕)e012567",
    ///     "+(+X͕g͔+Y͕q͔+Z͕ø͔-z͔Ø͕-Ð͕ç͔+Þ͕í͔+ï͔Œ͕)e012476",
    ///     "+(+X͕h͔+Y͕r͔+Z͕þ͔-z͔Þ͕-Ð͕é͔-Ø͕í͔+ñ͔Œ͕)e012457",
    ///     "+(+X͕i͔+Y͕s͔+Z͕œ͔-z͔Œ͕-Ð͕ë͔-Ø͕ï͔-Þ͕ñ͔)e012465",
    ///     "+(+X͕j͔+Y͕t͔+Z͕ç͔+Ð͕ø͔-Ø͕ð͔+Þ͕ó͔+ö͔Œ͕)e012367",
    ///     "+(+X͕k͔+Y͕u͔+Z͕é͔+Ð͕þ͔-Ø͕ó͔-Þ͕ð͔+ú͔Œ͕)e012375",
    ///     "+(+X͕l͔+Y͕á͔+Z͕ë͔+Ð͕œ͔-Ø͕ö͔-Þ͕ú͔-ð͔Œ͕)e012356",
    ///     "+(+X͕m͔+Y͕ä͔+Z͕í͔+Ð͕ó͔+Ø͕þ͔-Þ͕ø͔+ü͔Œ͕)e012347",
    ///     "+(+X͕n͔+Y͕å͔+Z͕ï͔+Ð͕ö͔+Ø͕œ͔-Þ͕ü͔-ø͔Œ͕)e012364",
    ///     "+(+X͕o͔+Y͕æ͔+Z͕ñ͔+Ð͕ú͔+Ø͕ü͔+Þ͕œ͔-þ͔Œ͕)e012345",
    ///     "+(+X͕x͔+Y͕y͔+Z͕z͔+Ð͕ð͔+Ø͕ø͔+Þ͕þ͔+Œ͕œ͔)I",
    /// ]);
    ///
    /// let norm_squared = Vee::triple_motor().norm_squared();
    ///
    /// assert_eq!(norm_squared.basis_blades(), Vee::norm().basis_blades());
    /// format_eq!(norm_squared, [
    ///     "+aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk+ll+mm+nn+oo+pp+qq+rr+ss+tt+uu+vv\
    ///      +xx+yy+zz+áá+ää+åå+ææ+çç+éé+ëë+íí+ïï+ðð+ññ+óó+öö+øø+úú+üü+þþ+œœ\
    ///      +αα+ββ+γγ+δδ+εε+ζζ+ηη+θθ+ιι+κκ+λλ+μμ+νν+ξξ+οο+ππ+ρρ+σσ+ττ+υυ+φφ",
    ///     "+2(-að+bz-cç-dé-eë-fy+gt+hu+iá-jq-kr-ls-mζ+nε-oδ+px\
    ///        +vü-äλ+åκ-æι-íο+ïξ-ñν-óσ+öρ+øφ-úπ-þυ+œτ-αμ+βθ-γη)e1234",
    ///     "+2(+Av+Bμ+Cν+Dξ+Eο-Fθ-Gι-Hκ-Iλ-Jo+Kn-Lm-Ml+Nk-Oj+Pγ\
    ///         +Qδ+Rε+Sζ-Tæ+Uå+Va-Xη+Yβ-Zα-bΜ-cΝ-dΞ-eΟ+fΘ+gΙ+hΚ\
    ///         +iΛ-pΓ-qΔ-rΕ-sΖ-tÆ+uÅ-xΗ+yΒ-zΑ-Áä-Äá-Çñ+Éï-Ëí-Íë\
    ///         +Ïé-Ðü-Ñç-Óœ+Öþ+Øú-Úø+Üð-Þö+óŒ+Πφ-Ρυ+Στ+Τσ-Υρ+Φπ)e0123",
    ///     "+2(-Aμ+Bv+Cπ+Dρ+Eσ+Fη+Go-Hn+Im-Jι-Kκ-Lλ+Mi-Nh+Og-Pβ\
    ///         +Qæ-Rå+Sä+Tδ+Uε+Vb-Xθ+Yγ+Zü+aΜ-cΠ-dΡ-eΣ-fΗ+jΙ+kΚ\
    ///         +lΛ+pΒ+qÆ-rÅ+sÄ-tΔ-uΕ-xΘ+yΓ-zÜ+Áζ-Çú+Éö-Ëó+Íœ-Ïþ\
    ///         -Ðα+Ñø-Óë+Öé-Øñ-Úç+Þï-áΖ-íŒ-ðΑ-Νφ+Ξυ-Οτ-Το+Υξ-Φν)e0124",
    ///     "+2(-Aν-Bπ+Cv+Dτ+Eυ-Fo+Gη+Hl-Ik+Jθ-Ki+Lh-Mκ-Nλ-Of-Pæ\
    ///         -Qβ+Rá-Su-Tγ-Us+Vc-Xι+Yδ-Zú+aΝ+bΠ-dΤ-eΥ-gΗ-jΘ+mΚ\
    ///         +nΛ-pÆ+qΒ+rÁ+tΓ-xΙ+yΔ+zÚ+Äε+Åζ-Çü-Éœ+Ëþ+Íö-Ïó+Ðñ\
    ///         -Ñð-Óï+Öí-Øα-Üç-Þë-äΕ-åΖ+éŒ-øΑ+Μφ-Ξσ+Ορ+Ρο-Σξ+Φμ)e0125",
    ///     "+2(-Aξ-Bρ-Cτ+Dv+Eφ+Fn-Gl+Hη+Ij+Ji+Kθ-Lg+Mι+Nf-Oλ+På\
    ///         -Qá-Rβ+St+Ts-Uγ+Vd-Xκ+Yε+Zö+aΞ+bΡ+cΤ-eΦ-hΗ-kΘ-mΙ\
    ///         +oΛ+pÅ-qÁ+rΒ+uΓ-xΚ+yΕ-zÖ-Äδ+Æζ+Çœ-Éü-Ëø+Íú+Ïð-Ðï\
    ///         -Ñó-Óñ+Øë+Úí-Üé-Þα+äΔ-æΖ-çŒ-þΑ-Μυ+Νσ-Οπ-Πο+Σν-Υμ)e0126",
    ///     "+2(-Aο-Bσ-Cυ-Dφ+Ev-Fm+Gk-Hj+Iη-Jh+Kg+Lθ-Mf+Nι+Oκ-Pä\
    ///         +Qu-Rt-Sβ-Tr+Uq+Ve-Xλ+Yζ-Zó+aΟ+bΣ+cΥ+dΦ-iΗ-lΘ-nΙ\
    ///         -oΚ-pÄ+sΒ-xΛ+yΖ+zÓ-Áγ-Åδ-Æε-Çþ+Éø-Ëü-Íð+Ïú+Ðí-Ñö\
    ///         -Öñ-Øé+Úï-Üë+Þç+áΓ+åΔ+æΕ-Œα-œΑ+Μτ-Νρ+Ξπ+Πξ-Ρν+Τμ)e0127",
    ///     "+2(+Aθ-Bη-Co+Dn-Em+Fv+Gπ+Hρ+Iσ-Jν-Kξ-Lο-Me+Nd-Oc+Pα\
    ///         +Qñ-Rï+Sí+Tú-Uö+Vf-Xμ-Yü+Zγ-aΘ+bΗ-gΠ-hΡ-iΣ+jΝ+kΞ\
    ///         +lΟ-pΑ+qÑ-rÏ+sÍ+tÚ-uÖ-xΜ+yÜ+zΓ+Áó-Äœ+Åþ-Æø+Çδ+Éε\
    ///         +Ëζ-Ðβ+Óá+Øæ-Þå+äŒ-çΔ-éΕ-ëΖ-ðΒ+Ιφ-Κυ+Λτ+Τλ-Υκ+Φι)e0134",
    ///     "+2(+Aι+Bo-Cη-Dl+Ek-Fπ+Gv+Hτ+Iυ+Jμ+Ke-Ld-Mξ-Nο+Ob-Pñ\
    ///         +Qα+Rë-Sé+Tü+Uœ+Vg-Xν+Yú+Zδ-aΙ+cΗ+fΠ-hΤ-iΥ-jΜ+mΞ\
    ///         +nΟ-pÑ-qΑ+rË-sÉ+tÜ-uŒ-xΝ-yÚ+zΔ-Áþ-Äö+Åó+Æð-Çγ+Íε\
    ///         +Ïζ-Ðæ+Óå-Öä-Øβ+Þá+çΓ-íΕ-ïΖ-øΒ-Θφ+Κσ-Λρ-Ρλ+Σκ-Φθ)e0135",
    ///     "+2(+Aκ-Bn+Cl-Dη-Ej-Fρ-Gτ+Hv+Iφ-Je+Kμ+Lc+Mν-Nb-Oο+Pï\
    ///         -Që+Rα+Sç-Tœ+Uü+Vh-Xξ-Yö+Zε-aΚ+dΗ+fΡ+gΤ-iΦ-kΜ-mΝ\
    ///         +oΟ+pÏ-qË-rΑ+sÇ+tŒ+uÜ-xΞ+yÖ+zΕ+Áø-Äú-Åð+Æó-Éγ-Íδ\
    ///         +Ðå+Ñζ+Óæ-Øá-Úä-Þβ+éΓ+íΔ-ñΖ-þΒ+Θυ-Ισ+Λπ+Πλ-Σι+Υθ)e0136",
    ///     "+2(+Aλ+Bm-Ck+Dj-Eη-Fσ-Gυ-Hφ+Iv+Jd-Kc+Lμ+Mb+Nν+Oξ-Pí\
    ///         +Qé-Rç+Sα+Tþ-Uø+Vi-Xο+Yó+Zζ-aΛ+eΗ+fΣ+gΥ+hΦ-lΜ-nΝ\
    ///         -oΞ-pÍ+qÉ-rÇ-sΑ-tÞ+uØ-xΟ-yÓ+zΖ+Áü+Äð-Åú+Æö-Ëγ-Ïδ\
    ///         -Ðä-Ñε+Öæ-Úå+Üá+ëΓ+ïΔ+ñΕ-Œβ-œΒ-Θτ+Ιρ-Κπ-Πκ+Ρι-Τθ)e0137",
    ///     "+2(-Ao+Bι-Cθ+Di-Eh+Fν-Gμ-He+Id+Jv+Kτ+Lυ-Mρ-Nσ-Oa-Pú\
    ///         -Qü-Rœ+Sþ+Tα+Uë+Vj-Xπ-Yñ+Zæ-bΙ+cΘ-fΝ+gΜ-kΤ-lΥ+mΡ\
    ///         +nΣ-pÚ-qÜ+rŒ-sÞ-tΑ+uË-xΠ+yÑ-zÆ-Áé+Äï-Åí+Çβ-Éá-Íå\
    ///         +Ïä+Ðδ+Óε+Öζ-Øγ-çΒ+ðΔ-óΕ-öΖ-øΓ+Ηφ-Κο+Λξ+Ξλ-Οκ+Φη)e0145",
    ///     "+2(+An+Bκ-Ci-Dθ+Eg+Fξ+Ge-Hμ-Ic-Jτ+Kv+Lφ+Mπ+Na-Oσ+Pö\
    ///         +Qœ-Rü-Sø-Të+Uα+Vk-Xρ+Yï-Zå-bΚ+dΘ-fΞ+hΜ+jΤ-lΦ-mΠ\
    ///         +oΣ+pÖ-qŒ-rÜ+sØ-tË-uΑ-xΡ-yÏ+zÅ+Áç+Äñ-Æí+Çá+Éβ-Íæ\
    ///         +Ðε+Ñä-Óδ+Úζ-Þγ-éΒ+ðΕ+óΔ-úΖ-þΓ-Ηυ+Ιο-Λν-Νλ+Οι-Υη)e0146",
    ///     "+2(-Am+Bλ+Ch-Dg-Eθ+Fο-Gd+Hc-Iμ-Jυ-Kφ+Lv-Ma+Nπ+Oρ-Pó\
    ///         -Qþ+Rø-Sü+Té-Uç+Vl-Xσ-Yí+Zä-bΛ+eΘ-fΟ+iΜ+jΥ+kΦ-nΠ\
    ///         -oΡ-pÓ+qÞ-rØ-sÜ+tÉ-uÇ-xΣ+yÍ-zÄ+Áα+Åñ-Æï+Ëβ-Ïæ+Ðζ\
    ///         +Ñå-Öδ-Úε-áΑ-ëΒ+ðΖ+öΔ+úΕ-Œγ-œΓ+Ητ-Ιξ+Κν+Νκ-Ξι+Τη)e0147",
    ///     "+2(-Al+Bi+Cκ-Dι-Ef-Fe+Gξ-Hν+Ib+Jρ-Kπ-La+Mv+Nφ-Oυ-Pœ\
    ///         +Qö+Rú+Sð-Tï-Uñ+Vm-Xτ-Yë+Zá-cΚ+dΙ-gΞ+hΝ-jΡ+kΠ-nΦ\
    ///         +oΥ+pŒ+qÖ+rÚ-sÐ-tÏ-uÑ-xΤ+yË-zÁ+Äα+Åç+Æé+Çå+Éæ+Íβ\
    ///         +Óγ+Øε+Üζ-Þδ-äΑ-íΒ-óΓ+øΕ-üΖ-þΔ+Ησ-Θο+Λμ+Μλ-Οθ+Ση)e0156",
    ///     "+2(+Ak-Bh+Cλ+Df-Eι+Fd+Gο-Hb-Iν+Jσ+Ka-Lπ-Mφ+Nv+Oτ+Pþ\
    ///         -Qó-Rð+Sú+Tí+Uz+Vn-Xυ+Yé-Zu-cΛ+eΙ-gΟ+iΝ-jΣ+lΠ+mΦ\
    ///         -oΤ-pÞ-qÓ+rÐ+sÚ+tÍ-xΥ-yÉ-Áñ-Äç+Åα+Æë-Çä+Ëæ+Ïβ-Ñá\
    ///         +Öγ+Øζ-Üε-åΑ-ïΒ-öΓ+øΖ+üΕ-Œδ-œΔ-Ηρ+Θξ-Κμ-Μκ+Ξθ-Ρη)e0157",
    ///     "+2(-Aj+Bg-Cf+Dλ-Eκ-Fc+Gb+Hο-Iξ-Ja+Kσ-Lρ+Mυ-Nτ+Ov-Pø\
    ///         +Qð-Ró-Sö-Tz+Uí+Vo-Xφ-Yç+Zt-dΛ+eΚ-hΟ+iΞ-kΣ+lΡ-mΥ\
    ///         +nΤ+pØ-qÐ-rÓ-sÖ+uÍ-xΦ+yÇ+Áï-Äé-Åë+Æα-Éä-Ëå+Ïá+Ñβ\
    ///         +Úγ+Üδ+Þζ-æΑ-ñΒ-úΓ-üΔ+þΖ-Œε-œΕ+Ηπ-Θν+Ιμ+Μι-Νθ+Πη)e0167",
    ///     "+2(-Aγ+Bβ-Cæ+Då-Eä-Fα-Gñ+Hï-Ií-Jú+Kö-Ló+Mœ-Nþ+Oø+Pv\
    ///         +Qπ+Rρ+Sσ-Tν-Uξ+Vp+Xü-Yμ+Zθ+aΓ-bΒ-cÆ+dÅ-eÄ+fΑ-gÑ\
    ///         +hÏ-iÍ-jÚ+kÖ-lÓ-mŒ+nÞ-oØ-qΠ-rΡ-sΣ+tΝ+uΞ-xÜ-yΜ+zΘ\
    ///         -Áο+Çι+Éκ+Ëλ-Ðη+áΟ-çΙ-éΚ-ëΛ-ðΗ-Δφ+Ευ-Ζτ-Τζ+Υε-Φδ)e0234",
    ///     "+2(-Aδ+Bæ+Cβ-Dá+Eu+Fñ-Gα-Hë+Ié-Jü-Kœ+Lþ+Mö-Nó-Oð-Pπ\
    ///         +Qv+Rτ+Sυ+Tμ+Ue+Vq-Xú-Yν+Zι+aΔ+bÆ-cΒ-dÁ+fÑ+gΑ-hË\
    ///         +iÉ-jÜ+kŒ-lÞ+mÖ-nÓ+oÐ+pΠ-rΤ-sΥ-tΜ+xÚ-yΝ+zΙ-Äξ-Åο\
    ///         -Çθ+Íκ+Ïλ-Øη+äΞ+åΟ+çΘ-íΚ-ïΛ-øΗ+Γφ-Εσ+Ζρ+Ρζ-Σε+Φγ)e0235",
    ///     "+2(-Aε-Bå+Cá+Dβ-Et-Fï+Gë-Hα-Iç+Jœ-Kü-Lø+Mú+Nð-Oó-Pρ\
    ///         -Qτ+Rv+Sφ-Te+Uμ+Vr+Xö-Yξ+Zκ+aΕ-bÅ+cÁ-dΒ-fÏ+gË+hΑ\
    ///         -iÇ-jŒ-kÜ+lØ+mÚ-nÐ-oÓ+pΡ+qΤ-sΦ-uΜ-xÖ-yΞ+zΚ+Äν-Æο\
    ///         -Éθ-Íι+Ñλ-Þη-äΝ+æΟ+éΘ+íΙ-ñΛ-þΗ-Γυ+Δσ-Ζπ-Πζ+Σδ-Υγ)e0236",
    ///     "+2(-Aζ+Bä-Cu+Dt+Eβ+Fí-Gé+Hç-Iα-Jþ+Kø-Lü-Mð+Nú-Oö-Pσ\
    ///         -Qυ-Rφ+Sv+Td-Uc+Vs-Xó-Yο+Zλ+aΖ+bÄ-eΒ+fÍ-gÉ+hÇ+iΑ\
    ///         +jÞ-kØ-lÜ+mÐ+nÚ-oÖ+pΣ+qΥ+rΦ+xÓ-yΟ+zΛ+Áμ+Åν+Æξ-Ëθ\
    ///         -Ïι-Ñκ-áΜ-åΝ-æΞ+ëΘ+ïΙ+ñΚ-Œη-œΗ+Γτ-Δρ+Επ+Πε-Ρδ+Τγ)e0237",
    ///     "+2(-Aæ-Bδ+Cγ+Ds-Er+Fú+Gü+Hœ-Iþ-Jα-Kë+Lé-Mï+Ní+Oz+Pν\
    ///         -Qμ-Re+Sd+Tv+Uτ+Vt+Xñ-Yπ-Zo-aÆ+bΔ-cΓ+fÚ+gÜ-hŒ+iÞ\
    ///         +jΑ-kË+lÉ-mÏ+nÍ-pΝ+qΜ-uΤ-xÑ-yΠ+Áυ-Äρ-Åσ+Çη+Ðι+Óκ\
    ///         +Öλ-Øθ-áΥ+äΡ+åΣ-çΗ+ðΙ-óΚ-öΛ-øΘ-Βφ+Εο-Ζξ-Ξζ+Οε-Φβ)e0245",
    ///     "+2(+Aå-Bε-Cs+Dγ+Eq-Fö-Gœ+Hü+Iø+Jë-Kα-Lç-Mñ-Nz+Oí+Pξ\
    ///         +Qe-Rμ-Sc-Tτ+Uv+Vu-Xï-Yρ+Zn+aÅ+bΕ-dΓ-fÖ+gŒ+hÜ-iØ\
    ///         +jË+kΑ-lÇ-mÑ+oÍ-pΞ+rΜ+tΤ+xÏ-yΡ+Áφ+Äπ-Æσ+Éη+Ðκ-Óι\
    ///         +Úλ-Þθ-áΦ-äΠ+æΣ-éΗ+ðΚ+óΙ-úΛ-þΘ+Βυ-Δο+Ζν+Νζ-Οδ+Υβ)e0246",
    ///     "+2(-Aä-Bζ+Cr-Dq+Eγ+Fó+Gþ-Hø+Iü-Jé+Kç-Lα+Mz-Nñ+Oï+Pο\
    ///         -Qd+Rc-Sμ-Tυ-Uφ+Vá+Xí-Yσ-Zm-aÄ+bΖ-eΓ+fÓ-gÞ+hØ+iÜ\
    ///         -jÉ+kÇ+lΑ-nÑ+oÏ-pΟ+sΜ+tΥ+uΦ+vÁ-xÍ-yΣ+Åπ+Æρ+Ëη+Ðλ\
    ///         -Öι-Úκ-åΠ-æΡ-ëΗ+ðΛ+öΙ+úΚ-Œθ-œΘ-Βτ+Δξ-Εν-Νε+Ξδ-Τβ)e0247",
    ///     "+2(-Aá+Bs-Cε+Dδ-Ep+Fœ-Gö-Hú-Ið+Jï+Kñ+Lz-Mα-Nç-Oé-Pe\
    ///         +Qξ-Rν+Sb+Tρ-Uπ+Vä+Xë-Yτ-Zl-aÁ+cΕ-dΔ-fŒ-gÖ-hÚ+iÐ\
    ///         +jÏ+kÑ+mΑ-nÇ-oÉ-qΞ+rΝ-tΡ+uΠ+vÄ-xË-yΤ+Åφ-Æυ+Íη+Óθ\
    ///         +Øκ+Üλ-Þι-åΦ+æΥ-íΗ-óΘ+øΚ-üΛ-þΙ-Βσ+Γο-Ζμ-Μζ+Ογ-Σβ)e0256",
    ///     "+2(+Au-Br-Cζ+Dp+Eδ-Fþ+Gó+Hð-Iú-Jí-Kz+Lñ+Mç-Nα-Oë+Pd\
    ///         +Qο-Rb-Sν+Tσ+Ua+Vå-Xé-Yυ+Zk+cΖ-eΔ+fÞ+gÓ-hÐ-iÚ-jÍ\
    ///         +lÑ+mÇ+nΑ-oË-qΟ+sΝ-tΣ+vÅ+xÉ-yΥ-Áπ-Äφ+Æτ+Ïη+Öθ+Øλ\
    ///         -Üκ+áΠ+äΦ-æΤ-ïΗ-öΘ+øΛ+üΚ-Œι-œΙ+Βρ-Γξ+Εμ+Με-Ξγ+Ρβ)e0257",
    ///     "+2(-At+Bq-Cp-Dζ+Eε+Fø-Gð+Hó+Iö+Jz-Kí-Lï+Mé+Në-Oα-Pc\
    ///         +Qb+Rο-Sξ-Ta+Uσ+Væ+Xç-Yφ-Zj+dΖ-eΕ-fØ+gÐ+hÓ+iÖ-kÍ\
    ///         -lÏ+mÉ+nË+oΑ-rΟ+sΞ-uΣ+vÆ-xÇ-yΦ-Áρ+Äυ-Åτ+Ñη+Úθ+Üι\
    ///         +Þλ+áΡ-äΥ+åΤ-ñΗ-úΘ-üΙ+þΛ-Œκ-œΚ-Βπ+Γν-Δμ-Μδ+Νγ-Πβ)e0267",
    ///     "+2(-Añ-Bú-Cü-Dœ+Eþ-Fδ+Gγ+Hs-Ir-Jβ+Ká-Lu+Må-Nä-Oy-Pι\
    ///         +Qθ-Ri+Sh-Tη-Ul+Vç-Xæ+Yo-Zπ-aÑ-bÚ-cÜ+dŒ-eÞ+fΔ-gΓ\
    ///         +jΒ+kÁ+mÅ-nÄ+pΙ-qΘ+tΗ+vÇ+xÆ-zΠ+Éτ+Ëυ-Íρ-Ïσ+Ðν+Óξ\
    ///         +Öο-Øμ-éΤ-ëΥ+íΡ+ïΣ+ðΝ-óΞ-öΟ-øΜ+Αφ-Ελ+Ζκ+Κζ-Λε+Φα)e0345",
    ///     "+2(+Aï+Bö+Cœ-Dü-Eø-Fε-Gs+Hγ+Iq-Já-Kβ+Lt+Mæ+Ny-Oä-Pκ\
    ///         +Qi+Rθ-Sg+Tl-Uη+Vé+Xå-Yn-Zρ+aÏ+bÖ-cŒ-dÜ+eØ+fΕ-hΓ\
    ///         -jÁ+kΒ+mÆ-oÄ+pΚ-rΘ+uΗ+vÉ-xÅ-zΡ-Çτ+Ëφ+Íπ+Ðξ-Ñσ-Óν\
    ///         +Úο-Þμ+çΤ-ëΦ-íΠ+ðΞ+ñΣ+óΝ-úΟ-þΜ-Αυ+Δλ-Ζι-Ιζ+Λδ-Υα)e0346",
    ///     "+2(-Aí-Bó-Cþ+Dø-Eü-Fζ+Gr-Hq+Iγ+Ju-Kt-Lβ-My+Næ-Oå-Pλ\
    ///         -Qh+Rg+Sθ-Tk+Uj+Vë-Xä+Ym-Zσ-aÍ-bÓ+cÞ-dØ-eÜ+fΖ-iΓ\
    ///         +lΒ+nÆ-oÅ+pΛ-sΘ+vË+xÄ-zΣ-Áη-Çυ-Éφ+Ïπ+Ðο+Ñρ-Öν-Úξ\
    ///         +áΗ+çΥ+éΦ-ïΠ+ðΟ-ñΡ+öΝ+úΞ-Œμ-œΜ+Ατ-Δκ+Ει+Ιε-Κδ+Τα)e0347",
    ///     "+2(-Aë-Bœ+Cö+Dú+Eð+Fs-Gε+Hδ-Ip-Jå-Kæ-Ly-Mβ+Nt+Ou-Pi\
    ///         -Qκ+Rι+Sf+Tn+Uo+Ví-Xá+Yl-Zτ-aË+bŒ+cÖ+dÚ-eÐ+gΕ-hΔ\
    ///         -jÅ-kÆ+mΒ+qΚ-rΙ+vÍ+xÁ-zΤ-Äη+Çρ-Éπ+Ïφ-Ñυ+Óμ+Øξ+Üο\
    ///         -Þν+äΗ-çΡ+éΠ-ïΦ+ñΥ-óΜ+øΞ-üΟ-þΝ+Ασ-Γλ+Ζθ+Θζ-Λγ+Σα)e0356",
    ///     "+2(+Aé+Bþ-Có-Dð+Eú-Fr-Gζ+Hp+Iδ+Jä+Ky-Læ-Mt-Nβ+Oá+Ph\
    ///         -Qλ-Rf+Sι-Tm-Ux+Vï+Xu-Yk-Zυ+aÉ-bÞ-cÓ+dÐ+eÚ+gΖ-iΔ\
    ///         +jÄ-lÆ+nΒ+oÁ+qΛ-sΙ+vÏ-zΥ-Åη+Çσ-Ëπ-Íφ+Ñτ+Öμ+Øο-Üξ\
    ///         +åΗ-çΣ+ëΠ+íΦ-ñΤ-öΜ+øΟ+üΞ-Œν-œΝ-Αρ+Γκ-Εθ-Θε+Κγ-Ρα)e0357",
    ///     "+2(-Aç-Bø+Cð-Dó-Eö+Fq-Gp-Hζ+Iε-Jy+Kä+Lå-Mu-Ná-Oβ-Pg\
    ///         +Qf-Rλ+Sκ+Tx-Um+Vñ-Xt+Yj-Zφ-aÇ+bØ-cÐ-dÓ-eÖ+hΖ-iΕ\
    ///         +kÄ+lÅ-nÁ+oΒ+rΛ-sΚ+vÑ-zΦ-Æη+Éσ-Ëρ+Íυ-Ïτ+Úμ+Üν+Þο\
    ///         +æΗ-éΣ+ëΡ-íΥ+ïΤ-úΜ-üΝ+þΟ-Œξ-œΞ+Απ-Γι+Δθ+Θδ-Ιγ+Πα)e0367",
    ///     "+2(+Aœ-Bë-Cï-Dñ-Ez+Fá+Gå+Hæ+Iy-Jε+Kδ-Lp-Mγ-Nq-Or-Pl\
    ///         -Qn-Ro-Sx-Tκ+Uι+Vó+Xs-Yi+Ze-aŒ-bË-cÏ-dÑ+fÁ+gÅ+hÆ\
    ///         +jΕ-kΔ+mΓ+tΚ-uΙ+vÓ-Äθ-Çξ+Éν-Íμ-Ðτ+Öφ+Øρ-Úυ+Üσ-Þπ\
    ///         +äΘ+çΞ-éΝ+íΜ-ðΤ-öΦ+øΡ+úΥ-üΣ-þΠ-Αο+Βλ-Ζη-Ηζ+Λβ-Οα)e0456",
    ///     "+2(-Aþ+Bé+Cí+Dz-Eñ-Fu-Gä-Hy+Iæ-Jζ+Kp+Lδ+Mq-Nγ-Os+Pk\
    ///         +Qm+Rx-So-Tλ-Uf+Vö-Xr+Yh-Zd+aÞ+bÉ+cÍ-eÑ-gÄ+iÆ+jΖ\
    ///         -lΔ+nΓ+tΛ+vÖ+Áι-Åθ-Çο+Ëν-Ïμ-Ðυ-Óφ+Øσ+Úτ-Üρ-áΙ+åΘ\
    ///         +çΟ-ëΝ+ïΜ-ðΥ+óΦ+øΣ-úΤ+üΡ-Œπ-œΠ+Αξ-Βκ+Εη+Ηε-Κβ+Ξα)e0457",
    ///     "+2(+Aø-Bç-Cz+Dí+Eï+Ft+Gy-Hä-Iå-Jp-Kζ+Lε+Mr+Ns-Oγ-Pj\
    ///         -Qx+Rm+Sn+Tf-Uλ+Vú+Xq-Yg+Zc-aØ-bÇ+dÍ+eÏ-hÄ-iÅ+kΖ\
    ///         -lΕ+oΓ+uΛ+vÚ+Áκ-Æθ-Éο+Ëξ-Ðφ-Ñμ+Óυ-Öτ+Üπ+Þσ-áΚ+æΘ\
    ///         +éΟ-ëΞ-ðΦ+ñΜ-óΥ+öΤ-üΠ+þΣ-Œρ-œΡ-Αν+Βι-Δη-Ηδ+Ιβ-Να)e0467",
    ///     "+2(-Að+Bz-Cç-Dé-Eë-Fy+Gt+Hu+Iá-Jq-Kr-Ls-Mζ+Nε-Oδ+Px\
    ///         -Qj-Rk-Sl+Tg+Uh+Vü-Xp+Yf-Zb+aÐ-cÇ-dÉ-eË+iÁ+mΖ-nΕ\
    ///         +oΔ+vÜ-Äλ+Åκ-Æι-Íο+Ïξ-Ñν-Óσ+Öρ-Øφ-Úπ+Þυ+äΛ-åΚ+æΙ\
    ///         +íΟ-ïΞ+ñΝ+óΣ-öΡ-øΦ+úΠ+þΥ-Œτ-œΤ+Αμ-Βθ+Γη+Ηγ-Θβ+Μα)e0567",
    ///     "+2(+av+bμ+cν+dξ+eο-fθ-gι-hκ-iλ-jo+kn-lm+pγ+qδ+rε+sζ\
    ///         -tæ+uå+xη-yβ+zα-áä-çñ+éï-ëí+ðü-óœ+öþ-øú-πφ+ρυ-στ)e4567",
    ///     "+2(-aμ+bv+cπ+dρ+eσ+fη+go-hn+im-jι-kκ-lλ-pβ+qæ-rå+sä\
    ///         +tδ+uε+xθ-yγ-zü+áζ-çú+éö-ëó+íœ-ïþ+ðα+ñø+νφ-ξυ+οτ)e3576",
    ///     "+2(-aν-bπ+cv+dτ+eυ-fo+gη+hl-ik+jθ-mκ-nλ-pæ-qβ+rá-su\
    ///         -tγ+xι-yδ+zú+äε+åζ-çü-éœ+ëþ+íö-ïó-ðñ+øα-μφ+ξσ-ορ)e3467",
    ///     "+2(-aξ-bρ-cτ+dv+eφ+fn-gl+hη+ij+kθ+mι-oλ+på-qá-rβ+st\
    ///         -uγ+xκ-yε-zö-äδ+æζ+çœ-éü-ëø+íú+ïð-ñó+þα+μυ-νσ+οπ)e3475",
    ///     "+2(-aο-bσ-cυ-dφ+ev-fm+gk-hj+iη+lθ+nι+oκ-pä+qu-rt-sβ\
    ///         +xλ-yζ+zó-áγ-åδ-æε-çþ+éø-ëü-íð+ïú-ñö+œα-μτ+νρ-ξπ)e3456",
    ///     "+2(+aθ-bη-co+dn-em+fv+gπ+hρ+iσ-jν-kξ-lο+pα+qñ-rï+sí\
    ///         +tú-uö+xμ+yü-zγ+áó-äœ+åþ-æø+çδ+éε+ëζ+ðβ-ιφ+κυ-λτ)e2567",
    ///     "+2(+aι+bo-cη-dl+ek-fπ+gv+hτ+iυ+jμ-mξ-nο-pñ+qα+rë-sé\
    ///         +tü+uœ+xν-yú-zδ-áþ-äö+åó+æð-çγ+íε+ïζ+øβ+θφ-κσ+λρ)e2476",
    ///     "+2(+aκ-bn+cl-dη-ej-fρ-gτ+hv+iφ+kμ+mν-oο+pï-që+rα+sç\
    ///         -tœ+uü+xξ+yö-zε+áø-äú-åð+æó-éγ-íδ+ñζ+þβ-θυ+ισ-λπ)e2457",
    ///     "+2(+aλ+bm-ck+dj-eη-fσ-gυ-hφ+iv+lμ+nν+oξ-pí+qé-rç+sα\
    ///         +tþ-uø+xο-yó-zζ+áü+äð-åú+æö-ëγ-ïδ-ñε+œβ+θτ-ιρ+κπ)e2465",
    ///     "+2(-ao+bι-cθ+di-eh+fν-gμ+jv+kτ+lυ-mρ-nσ-pú-qü-rœ+sþ\
    ///         +tα+uë+xπ+yñ-zæ-áé+äï-åí+çβ-ðδ+óε+öζ+øγ-ηφ+κο-λξ)e2367",
    ///     "+2(+an+bκ-ci-dθ+eg+fξ-hμ-jτ+kv+lφ+mπ-oσ+pö+qœ-rü-sø\
    ///         -të+uα+xρ-yï+zå+áç+äñ-æí+éβ-ðε-óδ+úζ+þγ+ηυ-ιο+λν)e2375",
    ///     "+2(-am+bλ+ch-dg-eθ+fο-iμ-jυ-kφ+lv+nπ+oρ-pó-qþ+rø-sü\
    ///         +té-uç+xσ+yí-zä+áα+åñ-æï+ëβ-ðζ-öδ-úε+œγ-ητ+ιξ-κν)e2356",
    ///     "+2(-al+bi+cκ-dι-ef+gξ-hν+jρ-kπ+mv+nφ-oυ-pœ+qö+rú+sð\
    ///         -tï-uñ+xτ+yë-zá+äα+åç+æé+íβ+óγ-øε+üζ+þδ-ησ+θο-λμ)e2347",
    ///     "+2(+ak-bh+cλ+df-eι+gο-iν+jσ-lπ-mφ+nv+oτ+pþ-qó-rð+sú\
    ///         +tí+uz+xυ-yé-áñ-äç+åα+æë+ïβ+öγ-øζ-üε+œδ+ηρ-θξ+κμ)e2364",
    ///     "+2(-aj+bg-cf+dλ-eκ+hο-iξ+kσ-lρ+mυ-nτ+ov-pø+qð-ró-sö\
    ///         -tz+uí+xφ+yç+áï-äé-åë+æα+ñβ+úγ+üδ-þζ+œε-ηπ+θν-ιμ)e2345",
    ///     "+2(-aγ+bβ-cæ+då-eä-fα-gñ+hï-ií-jú+kö-ló+mœ-nþ+oø+pv\
    ///         +qπ+rρ+sσ-tν-uξ-xü+yμ-zθ-áο+çι+éκ+ëλ+ðη+δφ-ευ+ζτ)e1576",
    ///     "+2(-aδ+bæ+cβ-dá+eu+fñ-gα-hë+ié-jü-kœ+lþ+mö-nó-oð-pπ\
    ///         +qv+rτ+sυ+tμ+xú+yν-zι-äξ-åο-çθ+íκ+ïλ+øη-γφ+εσ-ζρ)e1467",
    ///     "+2(-aε-bå+cá+dβ-et-fï+gë-hα-iç+jœ-kü-lø+mú+nð-oó-pρ\
    ///         -qτ+rv+sφ+uμ-xö+yξ-zκ+äν-æο-éθ-íι+ñλ+þη+γυ-δσ+ζπ)e1475",
    ///     "+2(-aζ+bä-cu+dt+eβ+fí-gé+hç-iα-jþ+kø-lü-mð+nú-oö-pσ\
    ///         -qυ-rφ+sv+xó+yο-zλ+áμ+åν+æξ-ëθ-ïι-ñκ+œη-γτ+δρ-επ)e1456",
    ///     "+2(-aæ-bδ+cγ+ds-er+fú+gü+hœ-iþ-jα-kë+lé-mï+ní+oz+pν\
    ///         -qμ+tv+uτ-xñ+yπ+áυ-äρ-åσ+çη-ðι+óκ+öλ+øθ+βφ-εο+ζξ)e1376",
    ///     "+2(+aå-bε-cs+dγ+eq-fö-gœ+hü+iø+jë-kα-lç-mñ-nz+oí+pξ\
    ///         -rμ-tτ+uv+xï+yρ+áφ+äπ-æσ+éη-ðκ-óι+úλ+þθ-βυ+δο-ζν)e1357",
    ///     "+2(-aä-bζ+cr-dq+eγ+fó+gþ-hø+iü-jé+kç-lα+mz-nñ+oï+pο\
    ///         -sμ-tυ-uφ+vá-xí+yσ+åπ+æρ+ëη-ðλ-öι-úκ+œθ+βτ-δξ+εν)e1365",
    ///     "+2(-aá+bs-cε+dδ-ep+fœ-gö-hú-ið+jï+kñ+lz-mα-nç-oé+qξ\
    ///         -rν+tρ-uπ+vä-xë+yτ+åφ-æυ+íη+óθ-øκ+üλ+þι+βσ-γο+ζμ)e1374",
    ///     "+2(+au-br-cζ+dp+eδ-fþ+gó+hð-iú-jí-kz+lñ+mç-nα-oë+qο\
    ///         -sν+tσ+vå+xé+yυ-áπ-äφ+æτ+ïη+öθ-øλ-üκ+œι-βρ+γξ-εμ)e1346",
    ///     "+2(-at+bq-cp-dζ+eε+fø-gð+hó+iö+jz-kí-lï+mé+në-oα+rο\
    ///         -sξ+uσ+væ-xç+yφ-áρ+äυ-åτ+ñη+úθ+üι-þλ+œκ+βπ-γν+δμ)e1354",
    ///     "+2(-añ-bú-cü-dœ+eþ-fδ+gγ+hs-ir-jβ+ká-lu+må-nä-oy-pι\
    ///         +qθ-tη+vç+xæ+zπ+éτ+ëυ-íρ-ïσ-ðν+óξ+öο+øμ-αφ+ελ-ζκ)e1267",
    ///     "+2(+aï+bö+cœ-dü-eø-fε-gs+hγ+iq-já-kβ+lt+mæ+ny-oä-pκ\
    ///         +rθ-uη+vé-xå+zρ-çτ+ëφ+íπ-ðξ-ñσ-óν+úο+þμ+αυ-δλ+ζι)e1275",
    ///     "+2(-aí-bó-cþ+dø-eü-fζ+gr-hq+iγ+ju-kt-lβ-my+næ-oå-pλ\
    ///         +sθ+vë+xä+zσ-áη-çυ-éφ+ïπ-ðο+ñρ-öν-úξ+œμ-ατ+δκ-ει)e1256",
    ///     "+2(-aë-bœ+cö+dú+eð+fs-gε+hδ-ip-jå-kæ-ly-mβ+nt+ou-qκ\
    ///         +rι+ví+xá+zτ-äη+çρ-éπ+ïφ-ñυ+óμ-øξ+üο+þν-ασ+γλ-ζθ)e1247",
    ///     "+2(+aé+bþ-có-dð+eú-fr-gζ+hp+iδ+jä+ky-læ-mt-nβ+oá-qλ\
    ///         +sι-ux+vï+zυ-åη+çσ-ëπ-íφ+ñτ+öμ-øο-üξ+œν+αρ-γκ+εθ)e1264",
    ///     "+2(-aç-bø+cð-dó-eö+fq-gp-hζ+iε-jy+kä+lå-mu-ná-oβ-rλ\
    ///         +sκ+tx+vñ+zφ-æη+éσ-ëρ+íυ-ïτ+úμ+üν-þο+œξ-απ+γι-δθ)e1245",
    ///     "+2(+aœ-bë-cï-dñ-ez+fá+gå+hæ+iy-jε+kδ-lp-mγ-nq-or-sx\
    ///         -tκ+uι+vó-äθ-çξ+éν-íμ+ðτ+öφ-øρ-úυ+üσ+þπ+αο-βλ+ζη)e1273",
    ///     "+2(-aþ+bé+cí+dz-eñ-fu-gä-hy+iæ-jζ+kp+lδ+mq-nγ-os+rx\
    ///         -tλ+vö+áι-åθ-çο+ëν-ïμ+ðυ-óφ-øσ+úτ-üρ+œπ-αξ+βκ-εη)e1236",
    ///     "+2(+aø-bç-cz+dí+eï+ft+gy-hä-iå-jp-kζ+lε+mr+ns-oγ-qx\
    ///         -uλ+vú+áκ-æθ-éο+ëξ+ðφ-ñμ+óυ-öτ+üπ-þσ+œρ+αν-βι+δη)e1253",
    ///     "+2(+Aa+Bb+Cc+Dd+Ee+Ff+Gg+Hh+Ii+Jj+Kk+Ll+Mm+Nn+Oo+Pp\
    ///         +Qq+Rr+Ss+Tt+Uu+Vv-Xx-Yy-Zz+Áá+Ää+Åå+Ææ+Çç+Éé+Ëë\
    ///         +Íí+Ïï-Ðð+Ññ+Óó+Öö-Øø+Úú+Üü-Þþ-Œœ-Αα-Ββ-Γγ-Δδ-Εε\
    ///         -Ζζ-Ηη-Θθ-Ιι-Κκ-Λλ-Μμ-Νν-Ξξ-Οο-Ππ-Ρρ-Σσ-Ττ-Υυ-Φφ)I",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_motor() -> Self {
        Self::scalar() + Self::volume5() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
    /// The multivector of single rotoreflector $`r_1 \equiv v^6 + v^4`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let single_rotoreflector = Vee::volume6().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(single_rotoreflector.basis_blades(), Vee::single_rotoreflector().basis_blades());
    /// format_eq!(single_rotoreflector, [
    ///     "+W͔v͕e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+W͔α͕e012",
    ///     "+W͔β͕e013",
    ///     "+W͔γ͕e014",
    ///     "+W͔δ͕e015",
    ///     "+W͔ε͕e016",
    ///     "+W͔ζ͕e017",
    ///     "+W͔η͕e023",
    ///     "+W͔θ͕e024",
    ///     "+W͔ι͕e025",
    ///     "+W͔κ͕e026",
    ///     "+W͔λ͕e027",
    ///     "+W͔μ͕e034",
    ///     "+W͔ν͕e035",
    ///     "+W͔ξ͕e036",
    ///     "+W͔ο͕e037",
    ///     "+W͔π͕e045",
    ///     "+W͔ρ͕e046",
    ///     "+W͔σ͕e047",
    ///     "+W͔τ͕e056",
    ///     "+W͔υ͕e057",
    ///     "+W͔φ͕e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔ι͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕+œ͔α͕)e127",
    ///     "+(+x͔μ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔ν͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔ξ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ο͕-z͔ζ͕+œ͔β͕)e137",
    ///     "+(+x͔π͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ρ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔σ͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+x͔τ͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔υ͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+x͔φ͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::volume6() + Self::volume4()
    }
    /// The multivector of double rotoreflector $`r_2 \equiv v^6 + v^4 + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let double_rotoreflector = Vee::volume6().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(double_rotoreflector.basis_blades(), Vee::double_rotoreflector().basis_blades());
    /// format_eq!(double_rotoreflector, [
    ///     "+W͔v͕e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+W͔α͕e012",
    ///     "+W͔β͕e013",
    ///     "+W͔γ͕e014",
    ///     "+W͔δ͕e015",
    ///     "+W͔ε͕e016",
    ///     "+W͔ζ͕e017",
    ///     "+W͔η͕e023",
    ///     "+W͔θ͕e024",
    ///     "+W͔ι͕e025",
    ///     "+W͔κ͕e026",
    ///     "+W͔λ͕e027",
    ///     "+W͔μ͕e034",
    ///     "+W͔ν͕e035",
    ///     "+W͔ξ͕e036",
    ///     "+W͔ο͕e037",
    ///     "+W͔π͕e045",
    ///     "+W͔ρ͕e046",
    ///     "+W͔σ͕e047",
    ///     "+W͔τ͕e056",
    ///     "+W͔υ͕e057",
    ///     "+W͔φ͕e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕-ð͔ü͕+ó͕œ͔-ö͕þ͔+ø͔ú͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+z͔ü͕-í͕œ͔+ï͕þ͔+ð͔α͕-ñ͕ø͔)e124",
    ///     "+(+x͔ι͕-y͔δ͕-z͔ú͕+é͕œ͔-ë͕þ͔+ð͔ñ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+z͔ö͕-ç͕œ͔+ë͕ø͔-ï͕ð͔+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕-z͔ó͕+ç͕þ͔-é͕ø͔+í͕ð͔+œ͔α͕)e127",
    ///     "+(+x͔μ͕-y͔ü͕-z͔γ͕+ä͕œ͔-å͕þ͔+æ͕ø͔+ð͔β͕)e134",
    ///     "+(-u͕œ͔+x͔ν͕+y͔ú͕-z͔δ͕+á͕þ͔-æ͕ð͔+ø͔β͕)e135",
    ///     "+(+t͕œ͔+x͔ξ͕-y͔ö͕-z͔ε͕-á͕ø͔+å͕ð͔+þ͔β͕)e136",
    ///     "+(-t͕þ͔+u͕ø͔+x͔ο͕+y͔ó͕-z͔ζ͕-ä͕ð͔+œ͔β͕)e137",
    ///     "+(+r͕œ͔-s͕þ͔+x͔π͕-y͔ñ͕+z͔æ͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(-q͕œ͔+s͕ø͔+x͔ρ͕+y͔ï͕-z͔å͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+q͕þ͔-r͕ø͔+x͔σ͕-y͔í͕+z͔ä͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+p͕œ͔-s͕ð͔+x͔τ͕-y͔ë͕+z͔á͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(-p͕þ͔+r͕ð͔-u͕z͔+x͔υ͕+y͔é͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+p͕ø͔-q͕ð͔+t͕z͔+x͔φ͕-y͔ç͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(-m͕œ͔+n͕þ͔-o͕ø͔+x͔ü͕+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+k͕œ͔-l͕þ͔+o͕ð͔-x͔ú͕+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(-j͕œ͔+l͕ø͔-n͕ð͔+x͔ö͕+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+j͕þ͔-k͕ø͔+m͕ð͔-x͔ó͕+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(-h͕œ͔+i͕þ͔-o͕z͔+x͔ñ͕+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+g͕œ͔-i͕ø͔+n͕z͔-x͔ï͕+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(-g͕þ͔+h͕ø͔-m͕z͔+x͔í͕+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(-f͕œ͔+i͕ð͔-l͕z͔+x͔ë͕+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+f͕þ͔-h͕ð͔+k͕z͔-x͔é͕+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(-f͕ø͔+g͕ð͔-j͕z͔+x͔ç͕+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+d͕œ͔-e͕þ͔+o͕y͔-x͔æ͕+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(-c͕œ͔+e͕ø͔-n͕y͔+x͔å͕+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+c͕þ͔-d͕ø͔+m͕y͔-x͔ä͕+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+b͕œ͔-e͕ð͔+l͕y͔-x͔á͕+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(-b͕þ͔+d͕ð͔-k͕y͔+u͕x͔+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+b͕ø͔-c͕ð͔+j͕y͔-t͕x͔+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(-a͕œ͔+e͕z͔-i͕y͔+s͕x͔+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+a͕þ͔-d͕z͔+h͕y͔-r͕x͔+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(-a͕ø͔+c͕z͔-g͕y͔+q͕x͔+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+a͕ð͔-b͕z͔+f͕y͔-p͕x͔+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(+a͕z͔+b͕ð͔+c͕ø͔+d͕þ͔+e͕œ͔)e34567",
    ///     "+(-a͕y͔+f͕ð͔+g͕ø͔+h͕þ͔+i͕œ͔)e24576",
    ///     "+(-b͕y͔-f͕z͔+j͕ø͔+k͕þ͔+l͕œ͔)e23567",
    ///     "+(-c͕y͔-g͕z͔-j͕ð͔+m͕þ͔+n͕œ͔)e23476",
    ///     "+(-d͕y͔-h͕z͔-k͕ð͔-m͕ø͔+o͕œ͔)e23457",
    ///     "+(-e͕y͔-i͕z͔-l͕ð͔-n͕ø͔-o͕þ͔)e23465",
    ///     "+(+a͕x͔+p͕ð͔+q͕ø͔+r͕þ͔+s͕œ͔)e14567",
    ///     "+(+b͕x͔-p͕z͔+t͕ø͔+u͕þ͔+á͕œ͔)e13576",
    ///     "+(+c͕x͔-q͕z͔-t͕ð͔+ä͕þ͔+å͕œ͔)e13467",
    ///     "+(+d͕x͔-r͕z͔-u͕ð͔-ä͕ø͔+æ͕œ͔)e13475",
    ///     "+(+e͕x͔-s͕z͔-á͕ð͔-å͕ø͔-æ͕þ͔)e13456",
    ///     "+(+f͕x͔+p͕y͔+ç͕ø͔+é͕þ͔+ë͕œ͔)e12567",
    ///     "+(+g͕x͔+q͕y͔-ç͕ð͔+í͕þ͔+ï͕œ͔)e12476",
    ///     "+(+h͕x͔+r͕y͔-é͕ð͔-í͕ø͔+ñ͕œ͔)e12457",
    ///     "+(+i͕x͔+s͕y͔-ë͕ð͔-ï͕ø͔-ñ͕þ͔)e12465",
    ///     "+(+j͕x͔+t͕y͔+z͔ç͕+ó͕þ͔+ö͕œ͔)e12367",
    ///     "+(+k͕x͔+u͕y͔+z͔é͕-ó͕ø͔+ú͕œ͔)e12375",
    ///     "+(+l͕x͔+y͔á͕+z͔ë͕-ö͕ø͔-ú͕þ͔)e12356",
    ///     "+(+m͕x͔+y͔ä͕+z͔í͕+ð͔ó͕+ü͕œ͔)e12347",
    ///     "+(+n͕x͔+y͔å͕+z͔ï͕+ð͔ö͕-ü͕þ͔)e12364",
    ///     "+(+o͕x͔+y͔æ͕+z͔ñ͕+ð͔ú͕+ø͔ü͕)e12345",
    ///     "-W͔a͕e04576",
    ///     "-W͔b͕e03567",
    ///     "-W͔c͕e03476",
    ///     "-W͔d͕e03457",
    ///     "-W͔e͕e03465",
    ///     "-W͔f͕e02576",
    ///     "-W͔g͕e02467",
    ///     "-W͔h͕e02475",
    ///     "-W͔i͕e02456",
    ///     "-W͔j͕e02376",
    ///     "-W͔k͕e02357",
    ///     "-W͔l͕e02365",
    ///     "-W͔m͕e02374",
    ///     "-W͔n͕e02346",
    ///     "-W͔o͕e02354",
    ///     "-W͔p͕e01567",
    ///     "-W͔q͕e01476",
    ///     "-W͔r͕e01457",
    ///     "-W͔s͕e01465",
    ///     "-W͔t͕e01367",
    ///     "-W͔u͕e01375",
    ///     "-W͔á͕e01356",
    ///     "-W͔ä͕e01347",
    ///     "-W͔å͕e01364",
    ///     "-W͔æ͕e01345",
    ///     "-W͔ç͕e01276",
    ///     "-W͔é͕e01257",
    ///     "-W͔ë͕e01265",
    ///     "-W͔í͕e01274",
    ///     "-W͔ï͕e01246",
    ///     "-W͔ñ͕e01254",
    ///     "-W͔ó͕e01237",
    ///     "-W͔ö͕e01263",
    ///     "-W͔ú͕e01235",
    ///     "-W͔ü͕e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane()
    }
    /// The multivector of triple rotoreflector $`r_3 \equiv v^6 + v^4 + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let triple_rotoreflector = Vee::volume6().lhs() * Vee::triple_rotator().rhs();
    ///
    /// assert_eq!(triple_rotoreflector.basis_blades(), Vee::triple_rotoreflector().basis_blades());
    /// format_eq!(triple_rotoreflector, [
    ///     "+W͔v͕e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+W͔α͕e012",
    ///     "+W͔β͕e013",
    ///     "+W͔γ͕e014",
    ///     "+W͔δ͕e015",
    ///     "+W͔ε͕e016",
    ///     "+W͔ζ͕e017",
    ///     "+W͔η͕e023",
    ///     "+W͔θ͕e024",
    ///     "+W͔ι͕e025",
    ///     "+W͔κ͕e026",
    ///     "+W͔λ͕e027",
    ///     "+W͔μ͕e034",
    ///     "+W͔ν͕e035",
    ///     "+W͔ξ͕e036",
    ///     "+W͔ο͕e037",
    ///     "+W͔π͕e045",
    ///     "+W͔ρ͕e046",
    ///     "+W͔σ͕e047",
    ///     "+W͔τ͕e056",
    ///     "+W͔υ͕e057",
    ///     "+W͔φ͕e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕-ð͔ü͕+ó͕œ͔-ö͕þ͔+ø͔ú͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+z͔ü͕-í͕œ͔+ï͕þ͔+ð͔α͕-ñ͕ø͔)e124",
    ///     "+(+x͔ι͕-y͔δ͕-z͔ú͕+é͕œ͔-ë͕þ͔+ð͔ñ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+z͔ö͕-ç͕œ͔+ë͕ø͔-ï͕ð͔+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕-z͔ó͕+ç͕þ͔-é͕ø͔+í͕ð͔+œ͔α͕)e127",
    ///     "+(+x͔μ͕-y͔ü͕-z͔γ͕+ä͕œ͔-å͕þ͔+æ͕ø͔+ð͔β͕)e134",
    ///     "+(-u͕œ͔+x͔ν͕+y͔ú͕-z͔δ͕+á͕þ͔-æ͕ð͔+ø͔β͕)e135",
    ///     "+(+t͕œ͔+x͔ξ͕-y͔ö͕-z͔ε͕-á͕ø͔+å͕ð͔+þ͔β͕)e136",
    ///     "+(-t͕þ͔+u͕ø͔+x͔ο͕+y͔ó͕-z͔ζ͕-ä͕ð͔+œ͔β͕)e137",
    ///     "+(+r͕œ͔-s͕þ͔+x͔π͕-y͔ñ͕+z͔æ͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(-q͕œ͔+s͕ø͔+x͔ρ͕+y͔ï͕-z͔å͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+q͕þ͔-r͕ø͔+x͔σ͕-y͔í͕+z͔ä͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+p͕œ͔-s͕ð͔+x͔τ͕-y͔ë͕+z͔á͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(-p͕þ͔+r͕ð͔-u͕z͔+x͔υ͕+y͔é͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+p͕ø͔-q͕ð͔+t͕z͔+x͔φ͕-y͔ç͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(-m͕œ͔+n͕þ͔-o͕ø͔+x͔ü͕+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+k͕œ͔-l͕þ͔+o͕ð͔-x͔ú͕+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(-j͕œ͔+l͕ø͔-n͕ð͔+x͔ö͕+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+j͕þ͔-k͕ø͔+m͕ð͔-x͔ó͕+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(-h͕œ͔+i͕þ͔-o͕z͔+x͔ñ͕+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+g͕œ͔-i͕ø͔+n͕z͔-x͔ï͕+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(-g͕þ͔+h͕ø͔-m͕z͔+x͔í͕+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(-f͕œ͔+i͕ð͔-l͕z͔+x͔ë͕+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+f͕þ͔-h͕ð͔+k͕z͔-x͔é͕+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(-f͕ø͔+g͕ð͔-j͕z͔+x͔ç͕+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+d͕œ͔-e͕þ͔+o͕y͔-x͔æ͕+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(-c͕œ͔+e͕ø͔-n͕y͔+x͔å͕+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+c͕þ͔-d͕ø͔+m͕y͔-x͔ä͕+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+b͕œ͔-e͕ð͔+l͕y͔-x͔á͕+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(-b͕þ͔+d͕ð͔-k͕y͔+u͕x͔+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+b͕ø͔-c͕ð͔+j͕y͔-t͕x͔+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(-a͕œ͔+e͕z͔-i͕y͔+s͕x͔+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+a͕þ͔-d͕z͔+h͕y͔-r͕x͔+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(-a͕ø͔+c͕z͔-g͕y͔+q͕x͔+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+a͕ð͔-b͕z͔+f͕y͔-p͕x͔+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(+a͕z͔+b͕ð͔+c͕ø͔+d͕þ͔+e͕œ͔-x͔y͕+x͕y͔)e34567",
    ///     "+(-a͕y͔+f͕ð͔+g͕ø͔+h͕þ͔+i͕œ͔-x͔z͕+x͕z͔)e24576",
    ///     "+(-b͕y͔-f͕z͔+j͕ø͔+k͕þ͔+l͕œ͔-x͔ð͕+x͕ð͔)e23567",
    ///     "+(-c͕y͔-g͕z͔-j͕ð͔+m͕þ͔+n͕œ͔-x͔ø͕+x͕ø͔)e23476",
    ///     "+(-d͕y͔-h͕z͔-k͕ð͔-m͕ø͔+o͕œ͔-x͔þ͕+x͕þ͔)e23457",
    ///     "+(-e͕y͔-i͕z͔-l͕ð͔-n͕ø͔-o͕þ͔-x͔œ͕+x͕œ͔)e23465",
    ///     "+(+a͕x͔+p͕ð͔+q͕ø͔+r͕þ͔+s͕œ͔-y͔z͕+y͕z͔)e14567",
    ///     "+(+b͕x͔-p͕z͔+t͕ø͔+u͕þ͔-y͔ð͕+y͕ð͔+á͕œ͔)e13576",
    ///     "+(+c͕x͔-q͕z͔-t͕ð͔-y͔ø͕+y͕ø͔+ä͕þ͔+å͕œ͔)e13467",
    ///     "+(+d͕x͔-r͕z͔-u͕ð͔-y͔þ͕+y͕þ͔-ä͕ø͔+æ͕œ͔)e13475",
    ///     "+(+e͕x͔-s͕z͔-y͔œ͕+y͕œ͔-á͕ð͔-å͕ø͔-æ͕þ͔)e13456",
    ///     "+(+f͕x͔+p͕y͔-z͔ð͕+z͕ð͔+ç͕ø͔+é͕þ͔+ë͕œ͔)e12567",
    ///     "+(+g͕x͔+q͕y͔-z͔ø͕+z͕ø͔-ç͕ð͔+í͕þ͔+ï͕œ͔)e12476",
    ///     "+(+h͕x͔+r͕y͔-z͔þ͕+z͕þ͔-é͕ð͔-í͕ø͔+ñ͕œ͔)e12457",
    ///     "+(+i͕x͔+s͕y͔-z͔œ͕+z͕œ͔-ë͕ð͔-ï͕ø͔-ñ͕þ͔)e12465",
    ///     "+(+j͕x͔+t͕y͔+z͔ç͕-ð͔ø͕+ð͕ø͔+ó͕þ͔+ö͕œ͔)e12367",
    ///     "+(+k͕x͔+u͕y͔+z͔é͕-ð͔þ͕+ð͕þ͔-ó͕ø͔+ú͕œ͔)e12375",
    ///     "+(+l͕x͔+y͔á͕+z͔ë͕-ð͔œ͕+ð͕œ͔-ö͕ø͔-ú͕þ͔)e12356",
    ///     "+(+m͕x͔+y͔ä͕+z͔í͕+ð͔ó͕-ø͔þ͕+ø͕þ͔+ü͕œ͔)e12347",
    ///     "+(+n͕x͔+y͔å͕+z͔ï͕+ð͔ö͕-ø͔œ͕+ø͕œ͔-ü͕þ͔)e12364",
    ///     "+(+o͕x͔+y͔æ͕+z͔ñ͕+ð͔ú͕+ø͔ü͕-þ͔œ͕+þ͕œ͔)e12345",
    ///     "-W͔a͕e04576",
    ///     "-W͔b͕e03567",
    ///     "-W͔c͕e03476",
    ///     "-W͔d͕e03457",
    ///     "-W͔e͕e03465",
    ///     "-W͔f͕e02576",
    ///     "-W͔g͕e02467",
    ///     "-W͔h͕e02475",
    ///     "-W͔i͕e02456",
    ///     "-W͔j͕e02376",
    ///     "-W͔k͕e02357",
    ///     "-W͔l͕e02365",
    ///     "-W͔m͕e02374",
    ///     "-W͔n͕e02346",
    ///     "-W͔o͕e02354",
    ///     "-W͔p͕e01567",
    ///     "-W͔q͕e01476",
    ///     "-W͔r͕e01457",
    ///     "-W͔s͕e01465",
    ///     "-W͔t͕e01367",
    ///     "-W͔u͕e01375",
    ///     "-W͔á͕e01356",
    ///     "-W͔ä͕e01347",
    ///     "-W͔å͕e01364",
    ///     "-W͔æ͕e01345",
    ///     "-W͔ç͕e01276",
    ///     "-W͔é͕e01257",
    ///     "-W͔ë͕e01265",
    ///     "-W͔í͕e01274",
    ///     "-W͔ï͕e01246",
    ///     "-W͔ñ͕e01254",
    ///     "-W͔ó͕e01237",
    ///     "-W͔ö͕e01263",
    ///     "-W͔ú͕e01235",
    ///     "-W͔ü͕e01243",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕+œ͔œ͕)e1234567",
    ///     "-W͔x͕e0234576",
    ///     "-W͔y͕e0134567",
    ///     "-W͔z͕e0124576",
    ///     "-W͔ð͕e0123567",
    ///     "-W͔ø͕e0123476",
    ///     "-W͔þ͕e0123457",
    ///     "-W͔œ͕e0123465",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn triple_rotoreflector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane() + Self::point()
    }
    /// The multivector of transflector $`f_t \equiv s + v^6 + v^4_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let transflector = Vee::volume6().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// format_eq!(transflector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+v͕x͔e1",
    ///     "+v͕y͔e2",
    ///     "+v͕z͔e3",
    ///     "+v͕ð͔e4",
    ///     "+v͕ø͔e5",
    ///     "+v͕þ͔e6",
    ///     "+v͕œ͔e7",
    ///     "+(+X͕y͔-Y͕x͔)e012",
    ///     "+(+X͕z͔-Z͕x͔)e013",
    ///     "+(+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+X͕œ͔-x͔Œ͕)e017",
    ///     "+(+Y͕z͔-Z͕y͔)e023",
    ///     "+(+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+Y͕œ͔-y͔Œ͕)e027",
    ///     "+(+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+Z͕œ͔-z͔Œ͕)e037",
    ///     "+(+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+Ð͕œ͔-ð͔Œ͕)e047",
    ///     "+(+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+Ø͕œ͔-ø͔Œ͕)e057",
    ///     "+(+Þ͕œ͔-þ͔Œ͕)e067",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume6() + Self::volume4_moment()
    }
    /// The multivector of simple single flector $`f_{s1} \equiv v^6 + v^4`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let simple_single_flector = Vee::volume6().lhs() * Vee::simple_single_motor().rhs();
    ///
    /// assert_eq!(simple_single_flector.basis_blades(),
    ///     Vee::simple_single_flector().basis_blades());
    /// format_eq!(simple_single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+(+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+W͔ζ͕+X͕œ͔-x͔Œ͕)e017",
    ///     "+(+W͔η͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(+W͔θ͕+Y͕ð͔-y͔Ð͕)e024",
    ///     "+(+W͔ι͕+Y͕ø͔-y͔Ø͕)e025",
    ///     "+(+W͔κ͕+Y͕þ͔-y͔Þ͕)e026",
    ///     "+(+W͔λ͕+Y͕œ͔-y͔Œ͕)e027",
    ///     "+(+W͔μ͕+Z͕ð͔-z͔Ð͕)e034",
    ///     "+(+W͔ν͕+Z͕ø͔-z͔Ø͕)e035",
    ///     "+(+W͔ξ͕+Z͕þ͔-z͔Þ͕)e036",
    ///     "+(+W͔ο͕+Z͕œ͔-z͔Œ͕)e037",
    ///     "+(+W͔π͕+Ð͕ø͔-Ø͕ð͔)e045",
    ///     "+(+W͔ρ͕+Ð͕þ͔-Þ͕ð͔)e046",
    ///     "+(+W͔σ͕+Ð͕œ͔-ð͔Œ͕)e047",
    ///     "+(+W͔τ͕+Ø͕þ͔-Þ͕ø͔)e056",
    ///     "+(+W͔υ͕+Ø͕œ͔-ø͔Œ͕)e057",
    ///     "+(+W͔φ͕+Þ͕œ͔-þ͔Œ͕)e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔ι͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕+œ͔α͕)e127",
    ///     "+(+x͔μ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔ν͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔ξ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ο͕-z͔ζ͕+œ͔β͕)e137",
    ///     "+(+x͔π͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ρ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔σ͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+x͔τ͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔υ͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+x͔φ͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_flector() -> Self {
        Self::volume6() + Self::volume4()
    }
    /// The multivector of single flector $`f_1 \equiv v^6 + v^4 + p_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let single_flector = Vee::volume6().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(single_flector.basis_blades(), Vee::single_flector().basis_blades());
    /// format_eq!(single_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+(-A͕z͔-B͕ð͔-C͕ø͔-D͕þ͔-E͕œ͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+A͕y͔-F͕ð͔-G͕ø͔-H͕þ͔-I͕œ͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+B͕y͔+F͕z͔-J͕ø͔-K͕þ͔-L͕œ͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+C͕y͔+G͕z͔+J͕ð͔-M͕þ͔-N͕œ͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+D͕y͔+H͕z͔+K͕ð͔+M͕ø͔-O͕œ͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+E͕y͔+I͕z͔+L͕ð͔+N͕ø͔+O͕þ͔+W͔ζ͕+X͕œ͔-x͔Œ͕)e017",
    ///     "+(-A͕x͔-P͕ð͔-Q͕ø͔-R͕þ͔-S͕œ͔+W͔η͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-B͕x͔+P͕z͔-T͕ø͔-U͕þ͔+W͔θ͕+Y͕ð͔-y͔Ð͕-Á͕œ͔)e024",
    ///     "+(-C͕x͔+Q͕z͔+T͕ð͔+W͔ι͕+Y͕ø͔-y͔Ø͕-Ä͕þ͔-Å͕œ͔)e025",
    ///     "+(-D͕x͔+R͕z͔+U͕ð͔+W͔κ͕+Y͕þ͔-y͔Þ͕+Ä͕ø͔-Æ͕œ͔)e026",
    ///     "+(-E͕x͔+S͕z͔+W͔λ͕+Y͕œ͔-y͔Œ͕+Á͕ð͔+Å͕ø͔+Æ͕þ͔)e027",
    ///     "+(-F͕x͔-P͕y͔+W͔μ͕+Z͕ð͔-z͔Ð͕-Ç͕ø͔-É͕þ͔-Ë͕œ͔)e034",
    ///     "+(-G͕x͔-Q͕y͔+W͔ν͕+Z͕ø͔-z͔Ø͕+Ç͕ð͔-Í͕þ͔-Ï͕œ͔)e035",
    ///     "+(-H͕x͔-R͕y͔+W͔ξ͕+Z͕þ͔-z͔Þ͕+É͕ð͔+Í͕ø͔-Ñ͕œ͔)e036",
    ///     "+(-I͕x͔-S͕y͔+W͔ο͕+Z͕œ͔-z͔Œ͕+Ë͕ð͔+Ï͕ø͔+Ñ͕þ͔)e037",
    ///     "+(-J͕x͔-T͕y͔+W͔π͕-z͔Ç͕+Ð͕ø͔-Ó͕þ͔-Ö͕œ͔-Ø͕ð͔)e045",
    ///     "+(-K͕x͔-U͕y͔+W͔ρ͕-z͔É͕+Ð͕þ͔+Ó͕ø͔-Ú͕œ͔-Þ͕ð͔)e046",
    ///     "+(-L͕x͔+W͔σ͕-y͔Á͕-z͔Ë͕+Ð͕œ͔+Ö͕ø͔+Ú͕þ͔-ð͔Œ͕)e047",
    ///     "+(-M͕x͔+W͔τ͕-y͔Ä͕-z͔Í͕-Ó͕ð͔+Ø͕þ͔-Ü͕œ͔-Þ͕ø͔)e056",
    ///     "+(-N͕x͔+W͔υ͕-y͔Å͕-z͔Ï͕-Ö͕ð͔+Ø͕œ͔+Ü͕þ͔-ø͔Œ͕)e057",
    ///     "+(-O͕x͔+W͔φ͕-y͔Æ͕-z͔Ñ͕-Ú͕ð͔-Ü͕ø͔+Þ͕œ͔-þ͔Œ͕)e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+ð͔α͕)e124",
    ///     "+(+x͔ι͕-y͔δ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕+œ͔α͕)e127",
    ///     "+(+x͔μ͕-z͔γ͕+ð͔β͕)e134",
    ///     "+(+x͔ν͕-z͔δ͕+ø͔β͕)e135",
    ///     "+(+x͔ξ͕-z͔ε͕+þ͔β͕)e136",
    ///     "+(+x͔ο͕-z͔ζ͕+œ͔β͕)e137",
    ///     "+(+x͔π͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(+x͔ρ͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+x͔σ͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+x͔τ͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(+x͔υ͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+x͔φ͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(-Ó͕œ͔+Ö͕þ͔-Ú͕ø͔+Ü͕ð͔)e04576",
    ///     "+(-z͔Ü͕+Í͕œ͔-Ï͕þ͔+Ñ͕ø͔)e03567",
    ///     "+(+z͔Ú͕-É͕œ͔+Ë͕þ͔-Ñ͕ð͔)e03476",
    ///     "+(-z͔Ö͕+Ç͕œ͔-Ë͕ø͔+Ï͕ð͔)e03457",
    ///     "+(+z͔Ó͕-Ç͕þ͔+É͕ø͔-Í͕ð͔)e03465",
    ///     "+(+y͔Ü͕-Ä͕œ͔+Å͕þ͔-Æ͕ø͔)e02576",
    ///     "+(+U͕œ͔-y͔Ú͕-Á͕þ͔+Æ͕ð͔)e02467",
    ///     "+(-T͕œ͔+y͔Ö͕+Á͕ø͔-Å͕ð͔)e02475",
    ///     "+(+T͕þ͔-U͕ø͔-y͔Ó͕+Ä͕ð͔)e02456",
    ///     "+(-R͕œ͔+S͕þ͔+y͔Ñ͕-z͔Æ͕)e02376",
    ///     "+(+Q͕œ͔-S͕ø͔-y͔Ï͕+z͔Å͕)e02357",
    ///     "+(-Q͕þ͔+R͕ø͔+y͔Í͕-z͔Ä͕)e02365",
    ///     "+(-P͕œ͔+S͕ð͔+y͔Ë͕-z͔Á͕)e02374",
    ///     "+(+P͕þ͔-R͕ð͔+U͕z͔-y͔É͕)e02346",
    ///     "+(-P͕ø͔+Q͕ð͔-T͕z͔+y͔Ç͕)e02354",
    ///     "+(+M͕œ͔-N͕þ͔+O͕ø͔-x͔Ü͕)e01567",
    ///     "+(-K͕œ͔+L͕þ͔-O͕ð͔+x͔Ú͕)e01476",
    ///     "+(+J͕œ͔-L͕ø͔+N͕ð͔-x͔Ö͕)e01457",
    ///     "+(-J͕þ͔+K͕ø͔-M͕ð͔+x͔Ó͕)e01465",
    ///     "+(+H͕œ͔-I͕þ͔+O͕z͔-x͔Ñ͕)e01367",
    ///     "+(-G͕œ͔+I͕ø͔-N͕z͔+x͔Ï͕)e01375",
    ///     "+(+G͕þ͔-H͕ø͔+M͕z͔-x͔Í͕)e01356",
    ///     "+(+F͕œ͔-I͕ð͔+L͕z͔-x͔Ë͕)e01347",
    ///     "+(-F͕þ͔+H͕ð͔-K͕z͔+x͔É͕)e01364",
    ///     "+(+F͕ø͔-G͕ð͔+J͕z͔-x͔Ç͕)e01345",
    ///     "+(-D͕œ͔+E͕þ͔-O͕y͔+x͔Æ͕)e01276",
    ///     "+(+C͕œ͔-E͕ø͔+N͕y͔-x͔Å͕)e01257",
    ///     "+(-C͕þ͔+D͕ø͔-M͕y͔+x͔Ä͕)e01265",
    ///     "+(-B͕œ͔+E͕ð͔-L͕y͔+x͔Á͕)e01274",
    ///     "+(+B͕þ͔-D͕ð͔+K͕y͔-U͕x͔)e01246",
    ///     "+(-B͕ø͔+C͕ð͔-J͕y͔+T͕x͔)e01254",
    ///     "+(+A͕œ͔-E͕z͔+I͕y͔-S͕x͔)e01237",
    ///     "+(-A͕þ͔+D͕z͔-H͕y͔+R͕x͔)e01263",
    ///     "+(+A͕ø͔-C͕z͔+G͕y͔-Q͕x͔)e01235",
    ///     "+(-A͕ð͔+B͕z͔-F͕y͔+P͕x͔)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane_moment()
    }
    /// The multivector of simple double flector $`f_{s2} \equiv v^6 + v^4 + p`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let simple_double_flector = Vee::volume6().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(simple_double_flector.basis_blades(),
    ///     Vee::simple_double_flector().basis_blades());
    /// format_eq!(simple_double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+(-A͕z͔-B͕ð͔-C͕ø͔-D͕þ͔-E͕œ͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+A͕y͔-F͕ð͔-G͕ø͔-H͕þ͔-I͕œ͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+B͕y͔+F͕z͔-J͕ø͔-K͕þ͔-L͕œ͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+C͕y͔+G͕z͔+J͕ð͔-M͕þ͔-N͕œ͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+D͕y͔+H͕z͔+K͕ð͔+M͕ø͔-O͕œ͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+E͕y͔+I͕z͔+L͕ð͔+N͕ø͔+O͕þ͔+W͔ζ͕+X͕œ͔-x͔Œ͕)e017",
    ///     "+(-A͕x͔-P͕ð͔-Q͕ø͔-R͕þ͔-S͕œ͔+W͔η͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-B͕x͔+P͕z͔-T͕ø͔-U͕þ͔+W͔θ͕+Y͕ð͔-y͔Ð͕-Á͕œ͔)e024",
    ///     "+(-C͕x͔+Q͕z͔+T͕ð͔+W͔ι͕+Y͕ø͔-y͔Ø͕-Ä͕þ͔-Å͕œ͔)e025",
    ///     "+(-D͕x͔+R͕z͔+U͕ð͔+W͔κ͕+Y͕þ͔-y͔Þ͕+Ä͕ø͔-Æ͕œ͔)e026",
    ///     "+(-E͕x͔+S͕z͔+W͔λ͕+Y͕œ͔-y͔Œ͕+Á͕ð͔+Å͕ø͔+Æ͕þ͔)e027",
    ///     "+(-F͕x͔-P͕y͔+W͔μ͕+Z͕ð͔-z͔Ð͕-Ç͕ø͔-É͕þ͔-Ë͕œ͔)e034",
    ///     "+(-G͕x͔-Q͕y͔+W͔ν͕+Z͕ø͔-z͔Ø͕+Ç͕ð͔-Í͕þ͔-Ï͕œ͔)e035",
    ///     "+(-H͕x͔-R͕y͔+W͔ξ͕+Z͕þ͔-z͔Þ͕+É͕ð͔+Í͕ø͔-Ñ͕œ͔)e036",
    ///     "+(-I͕x͔-S͕y͔+W͔ο͕+Z͕œ͔-z͔Œ͕+Ë͕ð͔+Ï͕ø͔+Ñ͕þ͔)e037",
    ///     "+(-J͕x͔-T͕y͔+W͔π͕-z͔Ç͕+Ð͕ø͔-Ó͕þ͔-Ö͕œ͔-Ø͕ð͔)e045",
    ///     "+(-K͕x͔-U͕y͔+W͔ρ͕-z͔É͕+Ð͕þ͔+Ó͕ø͔-Ú͕œ͔-Þ͕ð͔)e046",
    ///     "+(-L͕x͔+W͔σ͕-y͔Á͕-z͔Ë͕+Ð͕œ͔+Ö͕ø͔+Ú͕þ͔-ð͔Œ͕)e047",
    ///     "+(-M͕x͔+W͔τ͕-y͔Ä͕-z͔Í͕-Ó͕ð͔+Ø͕þ͔-Ü͕œ͔-Þ͕ø͔)e056",
    ///     "+(-N͕x͔+W͔υ͕-y͔Å͕-z͔Ï͕-Ö͕ð͔+Ø͕œ͔+Ü͕þ͔-ø͔Œ͕)e057",
    ///     "+(-O͕x͔+W͔φ͕-y͔Æ͕-z͔Ñ͕-Ú͕ð͔-Ü͕ø͔+Þ͕œ͔-þ͔Œ͕)e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕-ð͔ü͕+ó͕œ͔-ö͕þ͔+ø͔ú͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+z͔ü͕-í͕œ͔+ï͕þ͔+ð͔α͕-ñ͕ø͔)e124",
    ///     "+(+x͔ι͕-y͔δ͕-z͔ú͕+é͕œ͔-ë͕þ͔+ð͔ñ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+z͔ö͕-ç͕œ͔+ë͕ø͔-ï͕ð͔+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕-z͔ó͕+ç͕þ͔-é͕ø͔+í͕ð͔+œ͔α͕)e127",
    ///     "+(+x͔μ͕-y͔ü͕-z͔γ͕+ä͕œ͔-å͕þ͔+æ͕ø͔+ð͔β͕)e134",
    ///     "+(-u͕œ͔+x͔ν͕+y͔ú͕-z͔δ͕+á͕þ͔-æ͕ð͔+ø͔β͕)e135",
    ///     "+(+t͕œ͔+x͔ξ͕-y͔ö͕-z͔ε͕-á͕ø͔+å͕ð͔+þ͔β͕)e136",
    ///     "+(-t͕þ͔+u͕ø͔+x͔ο͕+y͔ó͕-z͔ζ͕-ä͕ð͔+œ͔β͕)e137",
    ///     "+(+r͕œ͔-s͕þ͔+x͔π͕-y͔ñ͕+z͔æ͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(-q͕œ͔+s͕ø͔+x͔ρ͕+y͔ï͕-z͔å͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+q͕þ͔-r͕ø͔+x͔σ͕-y͔í͕+z͔ä͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+p͕œ͔-s͕ð͔+x͔τ͕-y͔ë͕+z͔á͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(-p͕þ͔+r͕ð͔-u͕z͔+x͔υ͕+y͔é͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+p͕ø͔-q͕ð͔+t͕z͔+x͔φ͕-y͔ç͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(-m͕œ͔+n͕þ͔-o͕ø͔+x͔ü͕+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+k͕œ͔-l͕þ͔+o͕ð͔-x͔ú͕+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(-j͕œ͔+l͕ø͔-n͕ð͔+x͔ö͕+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+j͕þ͔-k͕ø͔+m͕ð͔-x͔ó͕+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(-h͕œ͔+i͕þ͔-o͕z͔+x͔ñ͕+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+g͕œ͔-i͕ø͔+n͕z͔-x͔ï͕+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(-g͕þ͔+h͕ø͔-m͕z͔+x͔í͕+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(-f͕œ͔+i͕ð͔-l͕z͔+x͔ë͕+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+f͕þ͔-h͕ð͔+k͕z͔-x͔é͕+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(-f͕ø͔+g͕ð͔-j͕z͔+x͔ç͕+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+d͕œ͔-e͕þ͔+o͕y͔-x͔æ͕+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(-c͕œ͔+e͕ø͔-n͕y͔+x͔å͕+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+c͕þ͔-d͕ø͔+m͕y͔-x͔ä͕+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+b͕œ͔-e͕ð͔+l͕y͔-x͔á͕+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(-b͕þ͔+d͕ð͔-k͕y͔+u͕x͔+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+b͕ø͔-c͕ð͔+j͕y͔-t͕x͔+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(-a͕œ͔+e͕z͔-i͕y͔+s͕x͔+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+a͕þ͔-d͕z͔+h͕y͔-r͕x͔+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(-a͕ø͔+c͕z͔-g͕y͔+q͕x͔+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+a͕ð͔-b͕z͔+f͕y͔-p͕x͔+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(+a͕z͔+b͕ð͔+c͕ø͔+d͕þ͔+e͕œ͔)e34567",
    ///     "+(-a͕y͔+f͕ð͔+g͕ø͔+h͕þ͔+i͕œ͔)e24576",
    ///     "+(-b͕y͔-f͕z͔+j͕ø͔+k͕þ͔+l͕œ͔)e23567",
    ///     "+(-c͕y͔-g͕z͔-j͕ð͔+m͕þ͔+n͕œ͔)e23476",
    ///     "+(-d͕y͔-h͕z͔-k͕ð͔-m͕ø͔+o͕œ͔)e23457",
    ///     "+(-e͕y͔-i͕z͔-l͕ð͔-n͕ø͔-o͕þ͔)e23465",
    ///     "+(+a͕x͔+p͕ð͔+q͕ø͔+r͕þ͔+s͕œ͔)e14567",
    ///     "+(+b͕x͔-p͕z͔+t͕ø͔+u͕þ͔+á͕œ͔)e13576",
    ///     "+(+c͕x͔-q͕z͔-t͕ð͔+ä͕þ͔+å͕œ͔)e13467",
    ///     "+(+d͕x͔-r͕z͔-u͕ð͔-ä͕ø͔+æ͕œ͔)e13475",
    ///     "+(+e͕x͔-s͕z͔-á͕ð͔-å͕ø͔-æ͕þ͔)e13456",
    ///     "+(+f͕x͔+p͕y͔+ç͕ø͔+é͕þ͔+ë͕œ͔)e12567",
    ///     "+(+g͕x͔+q͕y͔-ç͕ð͔+í͕þ͔+ï͕œ͔)e12476",
    ///     "+(+h͕x͔+r͕y͔-é͕ð͔-í͕ø͔+ñ͕œ͔)e12457",
    ///     "+(+i͕x͔+s͕y͔-ë͕ð͔-ï͕ø͔-ñ͕þ͔)e12465",
    ///     "+(+j͕x͔+t͕y͔+z͔ç͕+ó͕þ͔+ö͕œ͔)e12367",
    ///     "+(+k͕x͔+u͕y͔+z͔é͕-ó͕ø͔+ú͕œ͔)e12375",
    ///     "+(+l͕x͔+y͔á͕+z͔ë͕-ö͕ø͔-ú͕þ͔)e12356",
    ///     "+(+m͕x͔+y͔ä͕+z͔í͕+ð͔ó͕+ü͕œ͔)e12347",
    ///     "+(+n͕x͔+y͔å͕+z͔ï͕+ð͔ö͕-ü͕þ͔)e12364",
    ///     "+(+o͕x͔+y͔æ͕+z͔ñ͕+ð͔ú͕+ø͔ü͕)e12345",
    ///     "+(-W͔a͕-Ó͕œ͔+Ö͕þ͔-Ú͕ø͔+Ü͕ð͔)e04576",
    ///     "+(-W͔b͕-z͔Ü͕+Í͕œ͔-Ï͕þ͔+Ñ͕ø͔)e03567",
    ///     "+(-W͔c͕+z͔Ú͕-É͕œ͔+Ë͕þ͔-Ñ͕ð͔)e03476",
    ///     "+(-W͔d͕-z͔Ö͕+Ç͕œ͔-Ë͕ø͔+Ï͕ð͔)e03457",
    ///     "+(-W͔e͕+z͔Ó͕-Ç͕þ͔+É͕ø͔-Í͕ð͔)e03465",
    ///     "+(-W͔f͕+y͔Ü͕-Ä͕œ͔+Å͕þ͔-Æ͕ø͔)e02576",
    ///     "+(+U͕œ͔-W͔g͕-y͔Ú͕-Á͕þ͔+Æ͕ð͔)e02467",
    ///     "+(-T͕œ͔-W͔h͕+y͔Ö͕+Á͕ø͔-Å͕ð͔)e02475",
    ///     "+(+T͕þ͔-U͕ø͔-W͔i͕-y͔Ó͕+Ä͕ð͔)e02456",
    ///     "+(-R͕œ͔+S͕þ͔-W͔j͕+y͔Ñ͕-z͔Æ͕)e02376",
    ///     "+(+Q͕œ͔-S͕ø͔-W͔k͕-y͔Ï͕+z͔Å͕)e02357",
    ///     "+(-Q͕þ͔+R͕ø͔-W͔l͕+y͔Í͕-z͔Ä͕)e02365",
    ///     "+(-P͕œ͔+S͕ð͔-W͔m͕+y͔Ë͕-z͔Á͕)e02374",
    ///     "+(+P͕þ͔-R͕ð͔+U͕z͔-W͔n͕-y͔É͕)e02346",
    ///     "+(-P͕ø͔+Q͕ð͔-T͕z͔-W͔o͕+y͔Ç͕)e02354",
    ///     "+(+M͕œ͔-N͕þ͔+O͕ø͔-W͔p͕-x͔Ü͕)e01567",
    ///     "+(-K͕œ͔+L͕þ͔-O͕ð͔-W͔q͕+x͔Ú͕)e01476",
    ///     "+(+J͕œ͔-L͕ø͔+N͕ð͔-W͔r͕-x͔Ö͕)e01457",
    ///     "+(-J͕þ͔+K͕ø͔-M͕ð͔-W͔s͕+x͔Ó͕)e01465",
    ///     "+(+H͕œ͔-I͕þ͔+O͕z͔-W͔t͕-x͔Ñ͕)e01367",
    ///     "+(-G͕œ͔+I͕ø͔-N͕z͔-W͔u͕+x͔Ï͕)e01375",
    ///     "+(+G͕þ͔-H͕ø͔+M͕z͔-W͔á͕-x͔Í͕)e01356",
    ///     "+(+F͕œ͔-I͕ð͔+L͕z͔-W͔ä͕-x͔Ë͕)e01347",
    ///     "+(-F͕þ͔+H͕ð͔-K͕z͔-W͔å͕+x͔É͕)e01364",
    ///     "+(+F͕ø͔-G͕ð͔+J͕z͔-W͔æ͕-x͔Ç͕)e01345",
    ///     "+(-D͕œ͔+E͕þ͔-O͕y͔-W͔ç͕+x͔Æ͕)e01276",
    ///     "+(+C͕œ͔-E͕ø͔+N͕y͔-W͔é͕-x͔Å͕)e01257",
    ///     "+(-C͕þ͔+D͕ø͔-M͕y͔-W͔ë͕+x͔Ä͕)e01265",
    ///     "+(-B͕œ͔+E͕ð͔-L͕y͔-W͔í͕+x͔Á͕)e01274",
    ///     "+(+B͕þ͔-D͕ð͔+K͕y͔-U͕x͔-W͔ï͕)e01246",
    ///     "+(-B͕ø͔+C͕ð͔-J͕y͔+T͕x͔-W͔ñ͕)e01254",
    ///     "+(+A͕œ͔-E͕z͔+I͕y͔-S͕x͔-W͔ó͕)e01237",
    ///     "+(-A͕þ͔+D͕z͔-H͕y͔+R͕x͔-W͔ö͕)e01263",
    ///     "+(+A͕ø͔-C͕z͔+G͕y͔-Q͕x͔-W͔ú͕)e01235",
    ///     "+(-A͕ð͔+B͕z͔-F͕y͔+P͕x͔-W͔ü͕)e01243",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_flector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane()
    }
    /// The multivector of double flector $`f_2 \equiv v^6 + v^4 + p + P_\infty`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let double_flector = Vee::volume6().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// format_eq!(double_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+(-A͕z͔-B͕ð͔-C͕ø͔-D͕þ͔-E͕œ͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+A͕y͔-F͕ð͔-G͕ø͔-H͕þ͔-I͕œ͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+B͕y͔+F͕z͔-J͕ø͔-K͕þ͔-L͕œ͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+C͕y͔+G͕z͔+J͕ð͔-M͕þ͔-N͕œ͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+D͕y͔+H͕z͔+K͕ð͔+M͕ø͔-O͕œ͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+E͕y͔+I͕z͔+L͕ð͔+N͕ø͔+O͕þ͔+W͔ζ͕+X͕œ͔-x͔Œ͕)e017",
    ///     "+(-A͕x͔-P͕ð͔-Q͕ø͔-R͕þ͔-S͕œ͔+W͔η͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-B͕x͔+P͕z͔-T͕ø͔-U͕þ͔+W͔θ͕+Y͕ð͔-y͔Ð͕-Á͕œ͔)e024",
    ///     "+(-C͕x͔+Q͕z͔+T͕ð͔+W͔ι͕+Y͕ø͔-y͔Ø͕-Ä͕þ͔-Å͕œ͔)e025",
    ///     "+(-D͕x͔+R͕z͔+U͕ð͔+W͔κ͕+Y͕þ͔-y͔Þ͕+Ä͕ø͔-Æ͕œ͔)e026",
    ///     "+(-E͕x͔+S͕z͔+W͔λ͕+Y͕œ͔-y͔Œ͕+Á͕ð͔+Å͕ø͔+Æ͕þ͔)e027",
    ///     "+(-F͕x͔-P͕y͔+W͔μ͕+Z͕ð͔-z͔Ð͕-Ç͕ø͔-É͕þ͔-Ë͕œ͔)e034",
    ///     "+(-G͕x͔-Q͕y͔+W͔ν͕+Z͕ø͔-z͔Ø͕+Ç͕ð͔-Í͕þ͔-Ï͕œ͔)e035",
    ///     "+(-H͕x͔-R͕y͔+W͔ξ͕+Z͕þ͔-z͔Þ͕+É͕ð͔+Í͕ø͔-Ñ͕œ͔)e036",
    ///     "+(-I͕x͔-S͕y͔+W͔ο͕+Z͕œ͔-z͔Œ͕+Ë͕ð͔+Ï͕ø͔+Ñ͕þ͔)e037",
    ///     "+(-J͕x͔-T͕y͔+W͔π͕-z͔Ç͕+Ð͕ø͔-Ó͕þ͔-Ö͕œ͔-Ø͕ð͔)e045",
    ///     "+(-K͕x͔-U͕y͔+W͔ρ͕-z͔É͕+Ð͕þ͔+Ó͕ø͔-Ú͕œ͔-Þ͕ð͔)e046",
    ///     "+(-L͕x͔+W͔σ͕-y͔Á͕-z͔Ë͕+Ð͕œ͔+Ö͕ø͔+Ú͕þ͔-ð͔Œ͕)e047",
    ///     "+(-M͕x͔+W͔τ͕-y͔Ä͕-z͔Í͕-Ó͕ð͔+Ø͕þ͔-Ü͕œ͔-Þ͕ø͔)e056",
    ///     "+(-N͕x͔+W͔υ͕-y͔Å͕-z͔Ï͕-Ö͕ð͔+Ø͕œ͔+Ü͕þ͔-ø͔Œ͕)e057",
    ///     "+(-O͕x͔+W͔φ͕-y͔Æ͕-z͔Ñ͕-Ú͕ð͔-Ü͕ø͔+Þ͕œ͔-þ͔Œ͕)e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕-ð͔ü͕+ó͕œ͔-ö͕þ͔+ø͔ú͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+z͔ü͕-í͕œ͔+ï͕þ͔+ð͔α͕-ñ͕ø͔)e124",
    ///     "+(+x͔ι͕-y͔δ͕-z͔ú͕+é͕œ͔-ë͕þ͔+ð͔ñ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+z͔ö͕-ç͕œ͔+ë͕ø͔-ï͕ð͔+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕-z͔ó͕+ç͕þ͔-é͕ø͔+í͕ð͔+œ͔α͕)e127",
    ///     "+(+x͔μ͕-y͔ü͕-z͔γ͕+ä͕œ͔-å͕þ͔+æ͕ø͔+ð͔β͕)e134",
    ///     "+(-u͕œ͔+x͔ν͕+y͔ú͕-z͔δ͕+á͕þ͔-æ͕ð͔+ø͔β͕)e135",
    ///     "+(+t͕œ͔+x͔ξ͕-y͔ö͕-z͔ε͕-á͕ø͔+å͕ð͔+þ͔β͕)e136",
    ///     "+(-t͕þ͔+u͕ø͔+x͔ο͕+y͔ó͕-z͔ζ͕-ä͕ð͔+œ͔β͕)e137",
    ///     "+(+r͕œ͔-s͕þ͔+x͔π͕-y͔ñ͕+z͔æ͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(-q͕œ͔+s͕ø͔+x͔ρ͕+y͔ï͕-z͔å͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+q͕þ͔-r͕ø͔+x͔σ͕-y͔í͕+z͔ä͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+p͕œ͔-s͕ð͔+x͔τ͕-y͔ë͕+z͔á͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(-p͕þ͔+r͕ð͔-u͕z͔+x͔υ͕+y͔é͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+p͕ø͔-q͕ð͔+t͕z͔+x͔φ͕-y͔ç͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(-m͕œ͔+n͕þ͔-o͕ø͔+x͔ü͕+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+k͕œ͔-l͕þ͔+o͕ð͔-x͔ú͕+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(-j͕œ͔+l͕ø͔-n͕ð͔+x͔ö͕+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+j͕þ͔-k͕ø͔+m͕ð͔-x͔ó͕+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(-h͕œ͔+i͕þ͔-o͕z͔+x͔ñ͕+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+g͕œ͔-i͕ø͔+n͕z͔-x͔ï͕+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(-g͕þ͔+h͕ø͔-m͕z͔+x͔í͕+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(-f͕œ͔+i͕ð͔-l͕z͔+x͔ë͕+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+f͕þ͔-h͕ð͔+k͕z͔-x͔é͕+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(-f͕ø͔+g͕ð͔-j͕z͔+x͔ç͕+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+d͕œ͔-e͕þ͔+o͕y͔-x͔æ͕+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(-c͕œ͔+e͕ø͔-n͕y͔+x͔å͕+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+c͕þ͔-d͕ø͔+m͕y͔-x͔ä͕+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+b͕œ͔-e͕ð͔+l͕y͔-x͔á͕+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(-b͕þ͔+d͕ð͔-k͕y͔+u͕x͔+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+b͕ø͔-c͕ð͔+j͕y͔-t͕x͔+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(-a͕œ͔+e͕z͔-i͕y͔+s͕x͔+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+a͕þ͔-d͕z͔+h͕y͔-r͕x͔+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(-a͕ø͔+c͕z͔-g͕y͔+q͕x͔+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+a͕ð͔-b͕z͔+f͕y͔-p͕x͔+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(+a͕z͔+b͕ð͔+c͕ø͔+d͕þ͔+e͕œ͔)e34567",
    ///     "+(-a͕y͔+f͕ð͔+g͕ø͔+h͕þ͔+i͕œ͔)e24576",
    ///     "+(-b͕y͔-f͕z͔+j͕ø͔+k͕þ͔+l͕œ͔)e23567",
    ///     "+(-c͕y͔-g͕z͔-j͕ð͔+m͕þ͔+n͕œ͔)e23476",
    ///     "+(-d͕y͔-h͕z͔-k͕ð͔-m͕ø͔+o͕œ͔)e23457",
    ///     "+(-e͕y͔-i͕z͔-l͕ð͔-n͕ø͔-o͕þ͔)e23465",
    ///     "+(+a͕x͔+p͕ð͔+q͕ø͔+r͕þ͔+s͕œ͔)e14567",
    ///     "+(+b͕x͔-p͕z͔+t͕ø͔+u͕þ͔+á͕œ͔)e13576",
    ///     "+(+c͕x͔-q͕z͔-t͕ð͔+ä͕þ͔+å͕œ͔)e13467",
    ///     "+(+d͕x͔-r͕z͔-u͕ð͔-ä͕ø͔+æ͕œ͔)e13475",
    ///     "+(+e͕x͔-s͕z͔-á͕ð͔-å͕ø͔-æ͕þ͔)e13456",
    ///     "+(+f͕x͔+p͕y͔+ç͕ø͔+é͕þ͔+ë͕œ͔)e12567",
    ///     "+(+g͕x͔+q͕y͔-ç͕ð͔+í͕þ͔+ï͕œ͔)e12476",
    ///     "+(+h͕x͔+r͕y͔-é͕ð͔-í͕ø͔+ñ͕œ͔)e12457",
    ///     "+(+i͕x͔+s͕y͔-ë͕ð͔-ï͕ø͔-ñ͕þ͔)e12465",
    ///     "+(+j͕x͔+t͕y͔+z͔ç͕+ó͕þ͔+ö͕œ͔)e12367",
    ///     "+(+k͕x͔+u͕y͔+z͔é͕-ó͕ø͔+ú͕œ͔)e12375",
    ///     "+(+l͕x͔+y͔á͕+z͔ë͕-ö͕ø͔-ú͕þ͔)e12356",
    ///     "+(+m͕x͔+y͔ä͕+z͔í͕+ð͔ó͕+ü͕œ͔)e12347",
    ///     "+(+n͕x͔+y͔å͕+z͔ï͕+ð͔ö͕-ü͕þ͔)e12364",
    ///     "+(+o͕x͔+y͔æ͕+z͔ñ͕+ð͔ú͕+ø͔ü͕)e12345",
    ///     "+(-W͔a͕+x͔Η͕-y͔Β͕+z͔Α͕-Ó͕œ͔+Ö͕þ͔-Ú͕ø͔+Ü͕ð͔)e04576",
    ///     "+(-W͔b͕+x͔Θ͕-y͔Γ͕-z͔Ü͕+Í͕œ͔-Ï͕þ͔+Ñ͕ø͔+ð͔Α͕)e03567",
    ///     "+(-W͔c͕+x͔Ι͕-y͔Δ͕+z͔Ú͕-É͕œ͔+Ë͕þ͔-Ñ͕ð͔+ø͔Α͕)e03476",
    ///     "+(-W͔d͕+x͔Κ͕-y͔Ε͕-z͔Ö͕+Ç͕œ͔-Ë͕ø͔+Ï͕ð͔+þ͔Α͕)e03457",
    ///     "+(-W͔e͕+x͔Λ͕-y͔Ζ͕+z͔Ó͕-Ç͕þ͔+É͕ø͔-Í͕ð͔+œ͔Α͕)e03465",
    ///     "+(-W͔f͕+x͔Μ͕+y͔Ü͕-z͔Γ͕-Ä͕œ͔+Å͕þ͔-Æ͕ø͔+ð͔Β͕)e02576",
    ///     "+(+U͕œ͔-W͔g͕+x͔Ν͕-y͔Ú͕-z͔Δ͕-Á͕þ͔+Æ͕ð͔+ø͔Β͕)e02467",
    ///     "+(-T͕œ͔-W͔h͕+x͔Ξ͕+y͔Ö͕-z͔Ε͕+Á͕ø͔-Å͕ð͔+þ͔Β͕)e02475",
    ///     "+(+T͕þ͔-U͕ø͔-W͔i͕+x͔Ο͕-y͔Ó͕-z͔Ζ͕+Ä͕ð͔+œ͔Β͕)e02456",
    ///     "+(-R͕œ͔+S͕þ͔-W͔j͕+x͔Π͕+y͔Ñ͕-z͔Æ͕-ð͔Δ͕+ø͔Γ͕)e02376",
    ///     "+(+Q͕œ͔-S͕ø͔-W͔k͕+x͔Ρ͕-y͔Ï͕+z͔Å͕-ð͔Ε͕+þ͔Γ͕)e02357",
    ///     "+(-Q͕þ͔+R͕ø͔-W͔l͕+x͔Σ͕+y͔Í͕-z͔Ä͕-ð͔Ζ͕+œ͔Γ͕)e02365",
    ///     "+(-P͕œ͔+S͕ð͔-W͔m͕+x͔Τ͕+y͔Ë͕-z͔Á͕-ø͔Ε͕+þ͔Δ͕)e02374",
    ///     "+(+P͕þ͔-R͕ð͔+U͕z͔-W͔n͕+x͔Υ͕-y͔É͕-ø͔Ζ͕+œ͔Δ͕)e02346",
    ///     "+(-P͕ø͔+Q͕ð͔-T͕z͔-W͔o͕+x͔Φ͕+y͔Ç͕-þ͔Ζ͕+œ͔Ε͕)e02354",
    ///     "+(+M͕œ͔-N͕þ͔+O͕ø͔-W͔p͕-x͔Ü͕+y͔Μ͕-z͔Θ͕+ð͔Η͕)e01567",
    ///     "+(-K͕œ͔+L͕þ͔-O͕ð͔-W͔q͕+x͔Ú͕+y͔Ν͕-z͔Ι͕+ø͔Η͕)e01476",
    ///     "+(+J͕œ͔-L͕ø͔+N͕ð͔-W͔r͕-x͔Ö͕+y͔Ξ͕-z͔Κ͕+þ͔Η͕)e01457",
    ///     "+(-J͕þ͔+K͕ø͔-M͕ð͔-W͔s͕+x͔Ó͕+y͔Ο͕-z͔Λ͕+œ͔Η͕)e01465",
    ///     "+(+H͕œ͔-I͕þ͔+O͕z͔-W͔t͕-x͔Ñ͕+y͔Π͕-ð͔Ι͕+ø͔Θ͕)e01367",
    ///     "+(-G͕œ͔+I͕ø͔-N͕z͔-W͔u͕+x͔Ï͕+y͔Ρ͕-ð͔Κ͕+þ͔Θ͕)e01375",
    ///     "+(+G͕þ͔-H͕ø͔+M͕z͔-W͔á͕-x͔Í͕+y͔Σ͕-ð͔Λ͕+œ͔Θ͕)e01356",
    ///     "+(+F͕œ͔-I͕ð͔+L͕z͔-W͔ä͕-x͔Ë͕+y͔Τ͕-ø͔Κ͕+þ͔Ι͕)e01347",
    ///     "+(-F͕þ͔+H͕ð͔-K͕z͔-W͔å͕+x͔É͕+y͔Υ͕-ø͔Λ͕+œ͔Ι͕)e01364",
    ///     "+(+F͕ø͔-G͕ð͔+J͕z͔-W͔æ͕-x͔Ç͕+y͔Φ͕-þ͔Λ͕+œ͔Κ͕)e01345",
    ///     "+(-D͕œ͔+E͕þ͔-O͕y͔-W͔ç͕+x͔Æ͕+z͔Π͕-ð͔Ν͕+ø͔Μ͕)e01276",
    ///     "+(+C͕œ͔-E͕ø͔+N͕y͔-W͔é͕-x͔Å͕+z͔Ρ͕-ð͔Ξ͕+þ͔Μ͕)e01257",
    ///     "+(-C͕þ͔+D͕ø͔-M͕y͔-W͔ë͕+x͔Ä͕+z͔Σ͕-ð͔Ο͕+œ͔Μ͕)e01265",
    ///     "+(-B͕œ͔+E͕ð͔-L͕y͔-W͔í͕+x͔Á͕+z͔Τ͕-ø͔Ξ͕+þ͔Ν͕)e01274",
    ///     "+(+B͕þ͔-D͕ð͔+K͕y͔-U͕x͔-W͔ï͕+z͔Υ͕-ø͔Ο͕+œ͔Ν͕)e01246",
    ///     "+(-B͕ø͔+C͕ð͔-J͕y͔+T͕x͔-W͔ñ͕+z͔Φ͕-þ͔Ο͕+œ͔Ξ͕)e01254",
    ///     "+(+A͕œ͔-E͕z͔+I͕y͔-S͕x͔-W͔ó͕+ð͔Τ͕-ø͔Ρ͕+þ͔Π͕)e01237",
    ///     "+(-A͕þ͔+D͕z͔-H͕y͔+R͕x͔-W͔ö͕+ð͔Υ͕-ø͔Σ͕+œ͔Π͕)e01263",
    ///     "+(+A͕ø͔-C͕z͔+G͕y͔-Q͕x͔-W͔ú͕+ð͔Φ͕-þ͔Σ͕+œ͔Ρ͕)e01235",
    ///     "+(-A͕ð͔+B͕z͔-F͕y͔+P͕x͔-W͔ü͕+ø͔Φ͕-þ͔Υ͕+œ͔Τ͕)e01243",
    ///     "+(+y͔Α͕+z͔Β͕+ð͔Γ͕+ø͔Δ͕+þ͔Ε͕+œ͔Ζ͕)e0234576",
    ///     "+(-x͔Α͕+z͔Η͕+ð͔Θ͕+ø͔Ι͕+þ͔Κ͕+œ͔Λ͕)e0134567",
    ///     "+(-x͔Β͕-y͔Η͕+ð͔Μ͕+ø͔Ν͕+þ͔Ξ͕+œ͔Ο͕)e0124576",
    ///     "+(-x͔Γ͕-y͔Θ͕-z͔Μ͕+ø͔Π͕+þ͔Ρ͕+œ͔Σ͕)e0123567",
    ///     "+(-x͔Δ͕-y͔Ι͕-z͔Ν͕-ð͔Π͕+þ͔Τ͕+œ͔Υ͕)e0123476",
    ///     "+(-x͔Ε͕-y͔Κ͕-z͔Ξ͕-ð͔Ρ͕-ø͔Τ͕+œ͔Φ͕)e0123457",
    ///     "+(-x͔Ζ͕-y͔Λ͕-z͔Ο͕-ð͔Σ͕-ø͔Υ͕-þ͔Φ͕)e0123465",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn double_flector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane() + Self::direction()
    }
    /// The multivector of simple triple flector $`f_{s3} \equiv v^6 + v^4 + p + P`$.
    ///
    /// ```
    /// use vee::{format_eq, PgaP7 as Vee};
    ///
    /// let simple_triple_flector = Vee::volume6().lhs() * Vee::simple_triple_motor().rhs();
    ///
    /// assert_eq!(simple_triple_flector.basis_blades(),
    ///     Vee::simple_triple_flector().basis_blades());
    /// format_eq!(simple_triple_flector, [
    ///     "+(+W͔v͕-X͕x͔-Y͕y͔-Z͕z͔-Ð͕ð͔-Ø͕ø͔-Þ͕þ͔-Œ͕œ͔)e0",
    ///     "+(+v͕x͔-y͔α͕-z͔β͕-ð͔γ͕-ø͔δ͕-þ͔ε͕-œ͔ζ͕)e1",
    ///     "+(+v͕y͔+x͔α͕-z͔η͕-ð͔θ͕-ø͔ι͕-þ͔κ͕-œ͔λ͕)e2",
    ///     "+(+v͕z͔+x͔β͕+y͔η͕-ð͔μ͕-ø͔ν͕-þ͔ξ͕-œ͔ο͕)e3",
    ///     "+(+v͕ð͔+x͔γ͕+y͔θ͕+z͔μ͕-ø͔π͕-þ͔ρ͕-œ͔σ͕)e4",
    ///     "+(+v͕ø͔+x͔δ͕+y͔ι͕+z͔ν͕+ð͔π͕-þ͔τ͕-œ͔υ͕)e5",
    ///     "+(+v͕þ͔+x͔ε͕+y͔κ͕+z͔ξ͕+ð͔ρ͕+ø͔τ͕-œ͔φ͕)e6",
    ///     "+(+v͕œ͔+x͔ζ͕+y͔λ͕+z͔ο͕+ð͔σ͕+ø͔υ͕+þ͔φ͕)e7",
    ///     "+(-A͕z͔-B͕ð͔-C͕ø͔-D͕þ͔-E͕œ͔+W͔α͕+X͕y͔-Y͕x͔)e012",
    ///     "+(+A͕y͔-F͕ð͔-G͕ø͔-H͕þ͔-I͕œ͔+W͔β͕+X͕z͔-Z͕x͔)e013",
    ///     "+(+B͕y͔+F͕z͔-J͕ø͔-K͕þ͔-L͕œ͔+W͔γ͕+X͕ð͔-x͔Ð͕)e014",
    ///     "+(+C͕y͔+G͕z͔+J͕ð͔-M͕þ͔-N͕œ͔+W͔δ͕+X͕ø͔-x͔Ø͕)e015",
    ///     "+(+D͕y͔+H͕z͔+K͕ð͔+M͕ø͔-O͕œ͔+W͔ε͕+X͕þ͔-x͔Þ͕)e016",
    ///     "+(+E͕y͔+I͕z͔+L͕ð͔+N͕ø͔+O͕þ͔+W͔ζ͕+X͕œ͔-x͔Œ͕)e017",
    ///     "+(-A͕x͔-P͕ð͔-Q͕ø͔-R͕þ͔-S͕œ͔+W͔η͕+Y͕z͔-Z͕y͔)e023",
    ///     "+(-B͕x͔+P͕z͔-T͕ø͔-U͕þ͔+W͔θ͕+Y͕ð͔-y͔Ð͕-Á͕œ͔)e024",
    ///     "+(-C͕x͔+Q͕z͔+T͕ð͔+W͔ι͕+Y͕ø͔-y͔Ø͕-Ä͕þ͔-Å͕œ͔)e025",
    ///     "+(-D͕x͔+R͕z͔+U͕ð͔+W͔κ͕+Y͕þ͔-y͔Þ͕+Ä͕ø͔-Æ͕œ͔)e026",
    ///     "+(-E͕x͔+S͕z͔+W͔λ͕+Y͕œ͔-y͔Œ͕+Á͕ð͔+Å͕ø͔+Æ͕þ͔)e027",
    ///     "+(-F͕x͔-P͕y͔+W͔μ͕+Z͕ð͔-z͔Ð͕-Ç͕ø͔-É͕þ͔-Ë͕œ͔)e034",
    ///     "+(-G͕x͔-Q͕y͔+W͔ν͕+Z͕ø͔-z͔Ø͕+Ç͕ð͔-Í͕þ͔-Ï͕œ͔)e035",
    ///     "+(-H͕x͔-R͕y͔+W͔ξ͕+Z͕þ͔-z͔Þ͕+É͕ð͔+Í͕ø͔-Ñ͕œ͔)e036",
    ///     "+(-I͕x͔-S͕y͔+W͔ο͕+Z͕œ͔-z͔Œ͕+Ë͕ð͔+Ï͕ø͔+Ñ͕þ͔)e037",
    ///     "+(-J͕x͔-T͕y͔+W͔π͕-z͔Ç͕+Ð͕ø͔-Ó͕þ͔-Ö͕œ͔-Ø͕ð͔)e045",
    ///     "+(-K͕x͔-U͕y͔+W͔ρ͕-z͔É͕+Ð͕þ͔+Ó͕ø͔-Ú͕œ͔-Þ͕ð͔)e046",
    ///     "+(-L͕x͔+W͔σ͕-y͔Á͕-z͔Ë͕+Ð͕œ͔+Ö͕ø͔+Ú͕þ͔-ð͔Œ͕)e047",
    ///     "+(-M͕x͔+W͔τ͕-y͔Ä͕-z͔Í͕-Ó͕ð͔+Ø͕þ͔-Ü͕œ͔-Þ͕ø͔)e056",
    ///     "+(-N͕x͔+W͔υ͕-y͔Å͕-z͔Ï͕-Ö͕ð͔+Ø͕œ͔+Ü͕þ͔-ø͔Œ͕)e057",
    ///     "+(-O͕x͔+W͔φ͕-y͔Æ͕-z͔Ñ͕-Ú͕ð͔-Ü͕ø͔+Þ͕œ͔-þ͔Œ͕)e067",
    ///     "+(+x͔η͕-y͔β͕+z͔α͕-ð͔ü͕+ó͕œ͔-ö͕þ͔+ø͔ú͕)e123",
    ///     "+(+x͔θ͕-y͔γ͕+z͔ü͕-í͕œ͔+ï͕þ͔+ð͔α͕-ñ͕ø͔)e124",
    ///     "+(+x͔ι͕-y͔δ͕-z͔ú͕+é͕œ͔-ë͕þ͔+ð͔ñ͕+ø͔α͕)e125",
    ///     "+(+x͔κ͕-y͔ε͕+z͔ö͕-ç͕œ͔+ë͕ø͔-ï͕ð͔+þ͔α͕)e126",
    ///     "+(+x͔λ͕-y͔ζ͕-z͔ó͕+ç͕þ͔-é͕ø͔+í͕ð͔+œ͔α͕)e127",
    ///     "+(+x͔μ͕-y͔ü͕-z͔γ͕+ä͕œ͔-å͕þ͔+æ͕ø͔+ð͔β͕)e134",
    ///     "+(-u͕œ͔+x͔ν͕+y͔ú͕-z͔δ͕+á͕þ͔-æ͕ð͔+ø͔β͕)e135",
    ///     "+(+t͕œ͔+x͔ξ͕-y͔ö͕-z͔ε͕-á͕ø͔+å͕ð͔+þ͔β͕)e136",
    ///     "+(-t͕þ͔+u͕ø͔+x͔ο͕+y͔ó͕-z͔ζ͕-ä͕ð͔+œ͔β͕)e137",
    ///     "+(+r͕œ͔-s͕þ͔+x͔π͕-y͔ñ͕+z͔æ͕-ð͔δ͕+ø͔γ͕)e145",
    ///     "+(-q͕œ͔+s͕ø͔+x͔ρ͕+y͔ï͕-z͔å͕-ð͔ε͕+þ͔γ͕)e146",
    ///     "+(+q͕þ͔-r͕ø͔+x͔σ͕-y͔í͕+z͔ä͕-ð͔ζ͕+œ͔γ͕)e147",
    ///     "+(+p͕œ͔-s͕ð͔+x͔τ͕-y͔ë͕+z͔á͕-ø͔ε͕+þ͔δ͕)e156",
    ///     "+(-p͕þ͔+r͕ð͔-u͕z͔+x͔υ͕+y͔é͕-ø͔ζ͕+œ͔δ͕)e157",
    ///     "+(+p͕ø͔-q͕ð͔+t͕z͔+x͔φ͕-y͔ç͕-þ͔ζ͕+œ͔ε͕)e167",
    ///     "+(-m͕œ͔+n͕þ͔-o͕ø͔+x͔ü͕+y͔μ͕-z͔θ͕+ð͔η͕)e234",
    ///     "+(+k͕œ͔-l͕þ͔+o͕ð͔-x͔ú͕+y͔ν͕-z͔ι͕+ø͔η͕)e235",
    ///     "+(-j͕œ͔+l͕ø͔-n͕ð͔+x͔ö͕+y͔ξ͕-z͔κ͕+þ͔η͕)e236",
    ///     "+(+j͕þ͔-k͕ø͔+m͕ð͔-x͔ó͕+y͔ο͕-z͔λ͕+œ͔η͕)e237",
    ///     "+(-h͕œ͔+i͕þ͔-o͕z͔+x͔ñ͕+y͔π͕-ð͔ι͕+ø͔θ͕)e245",
    ///     "+(+g͕œ͔-i͕ø͔+n͕z͔-x͔ï͕+y͔ρ͕-ð͔κ͕+þ͔θ͕)e246",
    ///     "+(-g͕þ͔+h͕ø͔-m͕z͔+x͔í͕+y͔σ͕-ð͔λ͕+œ͔θ͕)e247",
    ///     "+(-f͕œ͔+i͕ð͔-l͕z͔+x͔ë͕+y͔τ͕-ø͔κ͕+þ͔ι͕)e256",
    ///     "+(+f͕þ͔-h͕ð͔+k͕z͔-x͔é͕+y͔υ͕-ø͔λ͕+œ͔ι͕)e257",
    ///     "+(-f͕ø͔+g͕ð͔-j͕z͔+x͔ç͕+y͔φ͕-þ͔λ͕+œ͔κ͕)e267",
    ///     "+(+d͕œ͔-e͕þ͔+o͕y͔-x͔æ͕+z͔π͕-ð͔ν͕+ø͔μ͕)e345",
    ///     "+(-c͕œ͔+e͕ø͔-n͕y͔+x͔å͕+z͔ρ͕-ð͔ξ͕+þ͔μ͕)e346",
    ///     "+(+c͕þ͔-d͕ø͔+m͕y͔-x͔ä͕+z͔σ͕-ð͔ο͕+œ͔μ͕)e347",
    ///     "+(+b͕œ͔-e͕ð͔+l͕y͔-x͔á͕+z͔τ͕-ø͔ξ͕+þ͔ν͕)e356",
    ///     "+(-b͕þ͔+d͕ð͔-k͕y͔+u͕x͔+z͔υ͕-ø͔ο͕+œ͔ν͕)e357",
    ///     "+(+b͕ø͔-c͕ð͔+j͕y͔-t͕x͔+z͔φ͕-þ͔ο͕+œ͔ξ͕)e367",
    ///     "+(-a͕œ͔+e͕z͔-i͕y͔+s͕x͔+ð͔τ͕-ø͔ρ͕+þ͔π͕)e456",
    ///     "+(+a͕þ͔-d͕z͔+h͕y͔-r͕x͔+ð͔υ͕-ø͔σ͕+œ͔π͕)e457",
    ///     "+(-a͕ø͔+c͕z͔-g͕y͔+q͕x͔+ð͔φ͕-þ͔σ͕+œ͔ρ͕)e467",
    ///     "+(+a͕ð͔-b͕z͔+f͕y͔-p͕x͔+ø͔φ͕-þ͔υ͕+œ͔τ͕)e567",
    ///     "+(+a͕z͔+b͕ð͔+c͕ø͔+d͕þ͔+e͕œ͔-x͔y͕+x͕y͔)e34567",
    ///     "+(-a͕y͔+f͕ð͔+g͕ø͔+h͕þ͔+i͕œ͔-x͔z͕+x͕z͔)e24576",
    ///     "+(-b͕y͔-f͕z͔+j͕ø͔+k͕þ͔+l͕œ͔-x͔ð͕+x͕ð͔)e23567",
    ///     "+(-c͕y͔-g͕z͔-j͕ð͔+m͕þ͔+n͕œ͔-x͔ø͕+x͕ø͔)e23476",
    ///     "+(-d͕y͔-h͕z͔-k͕ð͔-m͕ø͔+o͕œ͔-x͔þ͕+x͕þ͔)e23457",
    ///     "+(-e͕y͔-i͕z͔-l͕ð͔-n͕ø͔-o͕þ͔-x͔œ͕+x͕œ͔)e23465",
    ///     "+(+a͕x͔+p͕ð͔+q͕ø͔+r͕þ͔+s͕œ͔-y͔z͕+y͕z͔)e14567",
    ///     "+(+b͕x͔-p͕z͔+t͕ø͔+u͕þ͔-y͔ð͕+y͕ð͔+á͕œ͔)e13576",
    ///     "+(+c͕x͔-q͕z͔-t͕ð͔-y͔ø͕+y͕ø͔+ä͕þ͔+å͕œ͔)e13467",
    ///     "+(+d͕x͔-r͕z͔-u͕ð͔-y͔þ͕+y͕þ͔-ä͕ø͔+æ͕œ͔)e13475",
    ///     "+(+e͕x͔-s͕z͔-y͔œ͕+y͕œ͔-á͕ð͔-å͕ø͔-æ͕þ͔)e13456",
    ///     "+(+f͕x͔+p͕y͔-z͔ð͕+z͕ð͔+ç͕ø͔+é͕þ͔+ë͕œ͔)e12567",
    ///     "+(+g͕x͔+q͕y͔-z͔ø͕+z͕ø͔-ç͕ð͔+í͕þ͔+ï͕œ͔)e12476",
    ///     "+(+h͕x͔+r͕y͔-z͔þ͕+z͕þ͔-é͕ð͔-í͕ø͔+ñ͕œ͔)e12457",
    ///     "+(+i͕x͔+s͕y͔-z͔œ͕+z͕œ͔-ë͕ð͔-ï͕ø͔-ñ͕þ͔)e12465",
    ///     "+(+j͕x͔+t͕y͔+z͔ç͕-ð͔ø͕+ð͕ø͔+ó͕þ͔+ö͕œ͔)e12367",
    ///     "+(+k͕x͔+u͕y͔+z͔é͕-ð͔þ͕+ð͕þ͔-ó͕ø͔+ú͕œ͔)e12375",
    ///     "+(+l͕x͔+y͔á͕+z͔ë͕-ð͔œ͕+ð͕œ͔-ö͕ø͔-ú͕þ͔)e12356",
    ///     "+(+m͕x͔+y͔ä͕+z͔í͕+ð͔ó͕-ø͔þ͕+ø͕þ͔+ü͕œ͔)e12347",
    ///     "+(+n͕x͔+y͔å͕+z͔ï͕+ð͔ö͕-ø͔œ͕+ø͕œ͔-ü͕þ͔)e12364",
    ///     "+(+o͕x͔+y͔æ͕+z͔ñ͕+ð͔ú͕+ø͔ü͕-þ͔œ͕+þ͕œ͔)e12345",
    ///     "+(-W͔a͕+x͔Η͕-y͔Β͕+z͔Α͕-Ó͕œ͔+Ö͕þ͔-Ú͕ø͔+Ü͕ð͔)e04576",
    ///     "+(-W͔b͕+x͔Θ͕-y͔Γ͕-z͔Ü͕+Í͕œ͔-Ï͕þ͔+Ñ͕ø͔+ð͔Α͕)e03567",
    ///     "+(-W͔c͕+x͔Ι͕-y͔Δ͕+z͔Ú͕-É͕œ͔+Ë͕þ͔-Ñ͕ð͔+ø͔Α͕)e03476",
    ///     "+(-W͔d͕+x͔Κ͕-y͔Ε͕-z͔Ö͕+Ç͕œ͔-Ë͕ø͔+Ï͕ð͔+þ͔Α͕)e03457",
    ///     "+(-W͔e͕+x͔Λ͕-y͔Ζ͕+z͔Ó͕-Ç͕þ͔+É͕ø͔-Í͕ð͔+œ͔Α͕)e03465",
    ///     "+(-W͔f͕+x͔Μ͕+y͔Ü͕-z͔Γ͕-Ä͕œ͔+Å͕þ͔-Æ͕ø͔+ð͔Β͕)e02576",
    ///     "+(+U͕œ͔-W͔g͕+x͔Ν͕-y͔Ú͕-z͔Δ͕-Á͕þ͔+Æ͕ð͔+ø͔Β͕)e02467",
    ///     "+(-T͕œ͔-W͔h͕+x͔Ξ͕+y͔Ö͕-z͔Ε͕+Á͕ø͔-Å͕ð͔+þ͔Β͕)e02475",
    ///     "+(+T͕þ͔-U͕ø͔-W͔i͕+x͔Ο͕-y͔Ó͕-z͔Ζ͕+Ä͕ð͔+œ͔Β͕)e02456",
    ///     "+(-R͕œ͔+S͕þ͔-W͔j͕+x͔Π͕+y͔Ñ͕-z͔Æ͕-ð͔Δ͕+ø͔Γ͕)e02376",
    ///     "+(+Q͕œ͔-S͕ø͔-W͔k͕+x͔Ρ͕-y͔Ï͕+z͔Å͕-ð͔Ε͕+þ͔Γ͕)e02357",
    ///     "+(-Q͕þ͔+R͕ø͔-W͔l͕+x͔Σ͕+y͔Í͕-z͔Ä͕-ð͔Ζ͕+œ͔Γ͕)e02365",
    ///     "+(-P͕œ͔+S͕ð͔-W͔m͕+x͔Τ͕+y͔Ë͕-z͔Á͕-ø͔Ε͕+þ͔Δ͕)e02374",
    ///     "+(+P͕þ͔-R͕ð͔+U͕z͔-W͔n͕+x͔Υ͕-y͔É͕-ø͔Ζ͕+œ͔Δ͕)e02346",
    ///     "+(-P͕ø͔+Q͕ð͔-T͕z͔-W͔o͕+x͔Φ͕+y͔Ç͕-þ͔Ζ͕+œ͔Ε͕)e02354",
    ///     "+(+M͕œ͔-N͕þ͔+O͕ø͔-W͔p͕-x͔Ü͕+y͔Μ͕-z͔Θ͕+ð͔Η͕)e01567",
    ///     "+(-K͕œ͔+L͕þ͔-O͕ð͔-W͔q͕+x͔Ú͕+y͔Ν͕-z͔Ι͕+ø͔Η͕)e01476",
    ///     "+(+J͕œ͔-L͕ø͔+N͕ð͔-W͔r͕-x͔Ö͕+y͔Ξ͕-z͔Κ͕+þ͔Η͕)e01457",
    ///     "+(-J͕þ͔+K͕ø͔-M͕ð͔-W͔s͕+x͔Ó͕+y͔Ο͕-z͔Λ͕+œ͔Η͕)e01465",
    ///     "+(+H͕œ͔-I͕þ͔+O͕z͔-W͔t͕-x͔Ñ͕+y͔Π͕-ð͔Ι͕+ø͔Θ͕)e01367",
    ///     "+(-G͕œ͔+I͕ø͔-N͕z͔-W͔u͕+x͔Ï͕+y͔Ρ͕-ð͔Κ͕+þ͔Θ͕)e01375",
    ///     "+(+G͕þ͔-H͕ø͔+M͕z͔-W͔á͕-x͔Í͕+y͔Σ͕-ð͔Λ͕+œ͔Θ͕)e01356",
    ///     "+(+F͕œ͔-I͕ð͔+L͕z͔-W͔ä͕-x͔Ë͕+y͔Τ͕-ø͔Κ͕+þ͔Ι͕)e01347",
    ///     "+(-F͕þ͔+H͕ð͔-K͕z͔-W͔å͕+x͔É͕+y͔Υ͕-ø͔Λ͕+œ͔Ι͕)e01364",
    ///     "+(+F͕ø͔-G͕ð͔+J͕z͔-W͔æ͕-x͔Ç͕+y͔Φ͕-þ͔Λ͕+œ͔Κ͕)e01345",
    ///     "+(-D͕œ͔+E͕þ͔-O͕y͔-W͔ç͕+x͔Æ͕+z͔Π͕-ð͔Ν͕+ø͔Μ͕)e01276",
    ///     "+(+C͕œ͔-E͕ø͔+N͕y͔-W͔é͕-x͔Å͕+z͔Ρ͕-ð͔Ξ͕+þ͔Μ͕)e01257",
    ///     "+(-C͕þ͔+D͕ø͔-M͕y͔-W͔ë͕+x͔Ä͕+z͔Σ͕-ð͔Ο͕+œ͔Μ͕)e01265",
    ///     "+(-B͕œ͔+E͕ð͔-L͕y͔-W͔í͕+x͔Á͕+z͔Τ͕-ø͔Ξ͕+þ͔Ν͕)e01274",
    ///     "+(+B͕þ͔-D͕ð͔+K͕y͔-U͕x͔-W͔ï͕+z͔Υ͕-ø͔Ο͕+œ͔Ν͕)e01246",
    ///     "+(-B͕ø͔+C͕ð͔-J͕y͔+T͕x͔-W͔ñ͕+z͔Φ͕-þ͔Ο͕+œ͔Ξ͕)e01254",
    ///     "+(+A͕œ͔-E͕z͔+I͕y͔-S͕x͔-W͔ó͕+ð͔Τ͕-ø͔Ρ͕+þ͔Π͕)e01237",
    ///     "+(-A͕þ͔+D͕z͔-H͕y͔+R͕x͔-W͔ö͕+ð͔Υ͕-ø͔Σ͕+œ͔Π͕)e01263",
    ///     "+(+A͕ø͔-C͕z͔+G͕y͔-Q͕x͔-W͔ú͕+ð͔Φ͕-þ͔Σ͕+œ͔Ρ͕)e01235",
    ///     "+(-A͕ð͔+B͕z͔-F͕y͔+P͕x͔-W͔ü͕+ø͔Φ͕-þ͔Υ͕+œ͔Τ͕)e01243",
    ///     "+(+x͔x͕+y͔y͕+z͔z͕+ð͔ð͕+ø͔ø͕+þ͔þ͕+œ͔œ͕)e1234567",
    ///     "+(-W͔x͕+y͔Α͕+z͔Β͕+ð͔Γ͕+ø͔Δ͕+þ͔Ε͕+œ͔Ζ͕)e0234576",
    ///     "+(-W͔y͕-x͔Α͕+z͔Η͕+ð͔Θ͕+ø͔Ι͕+þ͔Κ͕+œ͔Λ͕)e0134567",
    ///     "+(-W͔z͕-x͔Β͕-y͔Η͕+ð͔Μ͕+ø͔Ν͕+þ͔Ξ͕+œ͔Ο͕)e0124576",
    ///     "+(-W͔ð͕-x͔Γ͕-y͔Θ͕-z͔Μ͕+ø͔Π͕+þ͔Ρ͕+œ͔Σ͕)e0123567",
    ///     "+(-W͔ø͕-x͔Δ͕-y͔Ι͕-z͔Ν͕-ð͔Π͕+þ͔Τ͕+œ͔Υ͕)e0123476",
    ///     "+(-W͔þ͕-x͔Ε͕-y͔Κ͕-z͔Ξ͕-ð͔Ρ͕-ø͔Τ͕+œ͔Φ͕)e0123457",
    ///     "+(-W͔œ͕-x͔Ζ͕-y͔Λ͕-z͔Ο͕-ð͔Σ͕-ø͔Υ͕-þ͔Φ͕)e0123465",
    /// ]);
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_triple_flector() -> Self {
        Self::volume6() + Self::volume4() + Self::plane() + Self::point()
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
fn cnv() {
    use super::{PgaP4 as Vee, Tree};

    let zero = Vee::zero();
    assert_eq!(Vee::try_from(Tree::from(zero.clone())), Ok(zero));

    let one = Vee::one();
    assert_eq!(Vee::try_from(Tree::from(one.clone())), Ok(one));

    let point = Vee::point().pin() << Vee::double_motor().unit();
    assert_eq!(
        Vee::try_from(Tree::with_factorization(point.clone(), true)),
        Ok(point.clone())
    );
    assert_eq!(
        Vee::try_from(Tree::with_factorization(point.clone(), false)),
        Ok(point)
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn ops() {
    use super::{PgaP0, PgaP1, PgaP2, PgaP3, PgaP4, PgaP5, PgaP6, PgaP7};

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

#[test]
#[allow(clippy::iter_on_single_items, clippy::too_many_lines)]
fn sym() {
    use super::{PgaP0, PgaP1, PgaP2, PgaP3, PgaP4, PgaP5, PgaP6, PgaP7};
    #[rustfmt::skip]
    assert!(
        [
            PgaP0::norm(),
        ]
        .iter()
        .all(PgaP0::is_entity)
    );
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

#[test]
fn tab() {
    fn pga<const N: u32>() {
        use core::str::from_utf8;
        let basis_len = 1 << (N + 1);
        let mut basis = Vec::with_capacity(basis_len);
        for idx in 0..basis_len {
            let mut fmt = Vec::with_capacity(idx.count_ones() as usize + 1);
            fmt.push(b'e');
            for nth in 0..=N {
                if (idx & (1 << nth)) != 0 {
                    fmt.push(b'0' + u8::try_from(nth).unwrap());
                }
            }
            basis.push((
                Pga::<0, N> {
                    idx: u8::try_from(idx).unwrap(),
                },
                fmt,
            ));
        }
        basis.sort_by(|(_, a), (_, b)| a.cmp(b));
        basis.sort_by_key(|(_, b)| b.len());
        for (i, (b, sym)) in basis.iter_mut().enumerate() {
            let lut = TAB[N as usize][i];
            if i >= basis_len / 2 {
                let (_sig, not) = !*b;
                let cnt = if not.cnt(*b) & 1 != 0 && sym.len() > 2 {
                    let len = sym.len();
                    sym.swap(len - 2, len - 1);
                    1
                } else {
                    0
                };
                assert_eq!(
                    lut.cnt,
                    cnt,
                    "`TAB[{N}][{i}].cnt`: expected {}, found {}",
                    from_utf8(sym).unwrap(),
                    lut.sym,
                );
            }
            assert_eq!(
                lut.sym.as_bytes(),
                sym,
                "`TAB[{N}][{i}].sym`: expected {}, found {}",
                from_utf8(sym).unwrap(),
                lut.sym
            );
        }
    }
    pga::<0>();
    pga::<1>();
    pga::<2>();
    pga::<3>();
    pga::<4>();
    pga::<5>();
    pga::<6>();
    pga::<7>();
}

#[test]
#[allow(clippy::cognitive_complexity)]
fn not() {
    use super::{PgaP0, PgaP1, PgaP2, PgaP3, PgaP4, PgaP5, PgaP6, PgaP7};

    assert_eq!(!PgaP0::norm(), PgaP0::norm().swp());

    assert_eq!(
        !PgaP1::point(),
        (PgaP1::weight() - PgaP1::direction()).swp()
    );
    assert_eq!(!!PgaP1::point(), -PgaP1::point());

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
        (PgaP5::plane_displacement() - PgaP5::plane_moment()).swp()
    );
    assert_eq!(!!PgaP5::plane(), -PgaP5::plane());
    assert_eq!(!PgaP5::line(), PgaP5::volume().swp());
    assert_eq!(!PgaP5::point(), -PgaP5::volume4().swp());

    assert_eq!(!PgaP6::volume5().alt(), PgaP6::point().swp());
    assert_eq!(!PgaP6::volume4().alt(), PgaP6::line().swp());
    assert_eq!(!PgaP6::volume().alt(), PgaP6::plane().swp());
    assert_eq!(!PgaP6::plane(), PgaP6::volume().alt().swp());
    assert_eq!(!PgaP6::line(), PgaP6::volume4().alt().swp());
    assert_eq!(!PgaP6::point(), PgaP6::volume5().alt().swp());

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
