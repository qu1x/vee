// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Plane-Based Pistachio Flavor -- Projective Geometric Algebra (PGA)

use super::{Algebra, Multivector, Symbol, choose};
use core::{
    cmp::Ordering,
    fmt::{self, Debug, Display, Error, Write},
    ops::{Mul, Not},
};
use std::fmt::Alignment;

/// Basis blade of Elliptic 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaE0 = Pga<1, 0>;
/// Basis blade of Elliptic 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaE1 = Pga<1, 1>;
/// Basis blade of Elliptic 2D PGA.
pub type PgaE2 = Pga<1, 2>;
/// Basis blade of Elliptic 3D PGA.
pub type PgaE3 = Pga<1, 3>;
/// Basis blade of Elliptic 4D PGA.
pub type PgaE4 = Pga<1, 4>;
/// Basis blade of Elliptic 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE5 = Pga<1, 5>;
/// Basis blade of Elliptic 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE6 = Pga<1, 6>;
/// Basis blade of Elliptic 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaE7 = Pga<1, 7>;

/// Basis blade of Hyperbolic 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaH0 = Pga<-1, 0>;
/// Basis blade of Hyperbolic 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaH1 = Pga<-1, 1>;
/// Basis blade of Hyperbolic 2D PGA.
pub type PgaH2 = Pga<-1, 2>;
/// Basis blade of Hyperbolic 3D PGA.
pub type PgaH3 = Pga<-1, 3>;
/// Basis blade of Hyperbolic 4D PGA.
pub type PgaH4 = Pga<-1, 4>;
/// Basis blade of Hyperbolic 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH5 = Pga<-1, 5>;
/// Basis blade of Hyperbolic 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH6 = Pga<-1, 6>;
/// Basis blade of Hyperbolic 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaH7 = Pga<-1, 7>;

/// Basis blade of Parabolic (Euclidean) 0D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaP0 = Pga<0, 0>;
/// Basis blade of Parabolic (Euclidean) 1D PGA.
#[cfg(feature = "rudimentary")]
pub type PgaP1 = Pga<0, 1>;
/// Basis blade of Parabolic (Euclidean) 2D PGA.
pub type PgaP2 = Pga<0, 2>;
/// Basis blade of Parabolic (Euclidean) 3D PGA.
pub type PgaP3 = Pga<0, 3>;
/// Basis blade of Parabolic (Euclidean) 4D PGA.
pub type PgaP4 = Pga<0, 4>;
/// Basis blade of Parabolic (Euclidean) 5D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP5 = Pga<0, 5>;
/// Basis blade of Parabolic (Euclidean) 6D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP6 = Pga<0, 6>;
/// Basis blade of Parabolic (Euclidean) 7D PGA (exploratory).
#[cfg(feature = "exploratory")]
pub type PgaP7 = Pga<0, 7>;

/// Basis blade of PGA with metric $`M\in\{\pm 1,0\}`$ and embedded dimension $`N\in[0, 7]`$.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Pga<const M: i8, const N: u32> {
    idx: u8,
}

/// Flavor-specific methods.
impl<const M: i8, const N: u32> Pga<M, N> {
    #[inline]
    const fn const_assert() {
        assert!(M.abs() <= 1 && N <= 7);
    }
    /// Creates basis blade from $`\e`$-notation.
    ///
    /// ```
    /// use vee::pga::PgaP3 as Bee;
    ///
    /// let e = Bee::new("e").unwrap();
    /// let e0 = Bee::new("e0").unwrap();
    /// let e1 = Bee::new("e1").unwrap();
    /// let e2 = Bee::new("e2").unwrap();
    /// let e12 = Bee::new("e12").unwrap();
    ///
    /// assert_eq!(e1 * e2, (1, e12));
    /// assert_eq!(e2 * e1, (-1, e12));
    /// assert_eq!(e0 * e0, (0, e));
    ///
    /// // The basis blade must be part of the flavor.
    ///
    /// assert!(Bee::new("e3").is_some());
    /// assert!(Bee::new("e4").is_none());
    ///
    /// // The basis blade must be postively permuted as defined by the flavor.
    ///
    /// assert!(Bee::new("e123").is_some());
    /// assert!(Bee::new("e321").is_some());
    /// assert!(Bee::new("e213").is_none());
    /// ```
    #[must_use]
    #[inline]
    pub const fn new(sym: &str) -> Option<Self> {
        const { Self::const_assert() };
        let lut = LUT[N as usize];
        let BasisBlade { cnt, idx, .. } = BasisBlade::new(sym);
        if (idx as usize) < lut.len() && (lut[idx as usize].cnt + cnt) & 1 == 0 {
            Some(Self { idx })
        } else {
            None
        }
    }
    /// Creates basis blade for scalar.
    #[must_use]
    #[inline]
    pub const fn scalar() -> Self {
        const { Self::const_assert() };
        Self { idx: 0 }
    }
    /// Creates basis blade for pseudoscalar.
    #[must_use]
    #[inline]
    pub const fn pss() -> Self {
        const { Self::const_assert() };
        Self {
            idx: u8::MAX >> (u8::BITS - (N + 1)),
        }
    }
    /// Counts swaps until `self * other` is ordered.
    #[must_use]
    #[inline]
    pub fn cnt(self, other: Self) -> u32 {
        ((1..=N).fold(0, |p, n| p ^ (self.idx >> n)) & other.idx).count_ones()
    }
    /// Constructs Cayley table.
    ///
    /// # Errors
    ///
    /// Fails in case of formatting errors.
    pub fn table() -> Result<String, Error> {
        const { Self::const_assert() };
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

impl<const M: i8, const N: u32> From<(Self, Pga<M, N>)> for Symbol {
    #[inline]
    fn from((s, _b): (Self, Pga<M, N>)) -> Self {
        s
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
        Self::new(s.lab).ok_or(s)
    }
}

impl<const M: i8, const N: u32> Algebra for Pga<M, N> {
    const N: u32 = N;

    #[inline]
    fn scalar() -> Self {
        Self::scalar()
    }
    #[inline]
    fn pseudoscalar() -> Self {
        Self::pss()
    }
    #[inline]
    fn basis() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator<Item = Self> {
        const { Self::const_assert() };
        TAB[N as usize].iter().map(|b| Self { idx: b.idx })
    }
    #[inline]
    fn blade_len(&self) -> usize {
        choose(N + 1, self.grade()) as usize
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
        let sig = if cnt & 1 != 0 { -1 } else { 1 };
        let sig = if lhs & rhs & 1 != 0 { sig * M } else { sig };
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
        let code = fmt.alternate() || fmt.align() == Some(Alignment::Center);
        match self.idx {
            0 if !code => write!(fmt, "1"),
            idx if !code && idx == Self::pss().idx => write!(fmt, "I"),
            _ => {
                let deref = if fmt.align() == Some(Alignment::Center) {
                    "o."
                } else {
                    ""
                };
                write!(fmt, "{deref}{}", LUT[N as usize][self.idx as usize].sym)
            }
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
            assert!(b'0' <= s[i] && s[i] <= b'7');
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

macro_rules! basis {
    ($t:ident, $u:ident, $n:tt, [$(($s:tt, $b:tt),)*]) => {
        #[doc = concat!(
            "The symbols and basis blades of the PGA with embedded dimension $`N = ", $n, "`$."
        )]
        ///
        /// ```gdef
        /// \gdef\var#1#2#3{#2}
        /// \gdef\sym#1{\var #1}
        /// \gdef\idx#1{\expandafter\sub#1\relax}
        /// \gdef\sub#1#2\relax{#2}
        /// \gdef\fmt#1{\e_{\idx{#1}}}
        /// ```
        impl<const M: i8> Pga<M, $n> {
            $(
                #[doc = concat!(
                    "The symbol and basis blade $`\\sym{",
                    stringify!($s),
                    "}\\fmt{",
                    stringify!($b),
                    "}`$.",
                )]
                #[must_use]
                #[inline]
                pub const fn $b() -> (Symbol, Self) {
                    (Symbol::new($s, stringify!($b)), Pga::new(stringify!($b)).unwrap())
                }
            )*
        }
        #[doc = concat!("The basis blades of the PGA with embedded dimension $`N = ", $n, "`$.")]
        ///
        /// ```gdef
        /// \gdef\var#1#2#3{#2}
        /// \gdef\sym#1{\var #1}
        /// \gdef\idx#1{\expandafter\sub#1\relax}
        /// \gdef\sub#1#2\relax{#2}
        /// \gdef\fmt#1{\e_{\idx{#1}}}
        /// ```
        impl<const M: i8> Multivector<Pga<M, $n>> {
            $(
                #[doc = concat!(
                    "The multivector of basis blade $`\\sym{",
                    stringify!($s),
                    "}\\fmt{",
                    stringify!($b),
                    "}`$.",
                )]
                #[must_use]
                #[inline]
                pub fn $b() -> Self {
                    Self::new([const { Pga::<M, $n>::$b() }])
                }
            )*
        }
        pub const $t: [BasisBlade; 1 << ($n + 1)] = BasisBlade::tab([$(stringify!($b),)*]);
        pub const $u: [BasisBlade; 1 << ($n + 1)] = BasisBlade::lut($t);
    };
}

#[cfg(feature = "rudimentary")]
mod n0;
#[cfg(feature = "rudimentary")]
mod n1;
mod n2;
mod n3;
mod n4;
#[cfg(feature = "exploratory")]
mod n5;
#[cfg(feature = "exploratory")]
mod n6;
#[cfg(feature = "exploratory")]
mod n7;

#[cfg(feature = "rudimentary")]
use n0::LUT0;
#[cfg(not(feature = "rudimentary"))]
const LUT0: [BasisBlade; 0] = [];
#[cfg(feature = "rudimentary")]
use n1::LUT1;
#[cfg(not(feature = "rudimentary"))]
const LUT1: [BasisBlade; 0] = [];
use n2::LUT2;
use n3::LUT3;
use n4::LUT4;
#[cfg(feature = "exploratory")]
use n5::LUT5;
#[cfg(not(feature = "exploratory"))]
const LUT5: [BasisBlade; 0] = [];
#[cfg(feature = "exploratory")]
use n6::LUT6;
#[cfg(not(feature = "exploratory"))]
const LUT6: [BasisBlade; 0] = [];
#[cfg(feature = "exploratory")]
use n7::LUT7;
#[cfg(not(feature = "exploratory"))]
const LUT7: [BasisBlade; 0] = [];

#[cfg(feature = "rudimentary")]
use n0::TAB0;
#[cfg(not(feature = "rudimentary"))]
const TAB0: [BasisBlade; 0] = [];
#[cfg(feature = "rudimentary")]
use n1::TAB1;
#[cfg(not(feature = "rudimentary"))]
const TAB1: [BasisBlade; 0] = [];
use n2::TAB2;
use n3::TAB3;
use n4::TAB4;
#[cfg(feature = "exploratory")]
use n5::TAB5;
#[cfg(not(feature = "exploratory"))]
const TAB5: [BasisBlade; 0] = [];
#[cfg(feature = "exploratory")]
use n6::TAB6;
#[cfg(not(feature = "exploratory"))]
const TAB6: [BasisBlade; 0] = [];
#[cfg(feature = "exploratory")]
use n7::TAB7;
#[cfg(not(feature = "exploratory"))]
const TAB7: [BasisBlade; 0] = [];

const TAB: [&[BasisBlade]; 8] = [&TAB0, &TAB1, &TAB2, &TAB3, &TAB4, &TAB5, &TAB6, &TAB7];
const LUT: [&[BasisBlade]; 8] = [&LUT0, &LUT1, &LUT2, &LUT3, &LUT4, &LUT5, &LUT6, &LUT7];

#[test]
fn tabs() {
    let s = |l: &[BasisBlade], h: &[BasisBlade]| l.iter().all(|l| h.iter().any(|h| h.sym == l.sym));

    #[cfg(feature = "rudimentary")]
    assert!(s(TAB0.as_slice(), TAB1.as_slice()));
    #[cfg(feature = "rudimentary")]
    assert!(s(TAB1.as_slice(), TAB2.as_slice()));

    assert!(!s(TAB2.as_slice(), TAB3.as_slice()));
    assert!(!s(TAB3.as_slice(), TAB4.as_slice()));
    #[cfg(feature = "exploratory")]
    assert!(!s(TAB4.as_slice(), TAB5.as_slice()));
    #[cfg(feature = "exploratory")]
    assert!(!s(TAB5.as_slice(), TAB6.as_slice()));
    #[cfg(feature = "exploratory")]
    assert!(!s(TAB6.as_slice(), TAB7.as_slice()));
}

#[cfg(test)]
fn tab<const N: u32>() {
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

macro_rules! tab {
    ($($f:ident: $n:expr;)*) => {
        $(
            #[test]
            fn $f() {
                tab::<$n>();
            }
        )*
    };
}

#[cfg(feature = "rudimentary")]
tab! {
    tab0: 0;
    tab1: 1;
}

tab! {
    tab2: 2;
    tab3: 3;
    tab4: 4;
}

#[cfg(feature = "exploratory")]
tab! {
    tab5: 5;
    tab6: 6;
    tab7: 7;
}
