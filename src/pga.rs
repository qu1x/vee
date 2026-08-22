// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Plane-Based Pistachio Flavor -- Projective Geometric Algebra (PGA)

use super::{Algebra, Choose, Factor, Multivector, Symbol};
use core::{
    cmp::Ordering,
    fmt::{self, Debug, Display, Error, Write},
    ops::{Mul, Not},
};
use std::fmt::Alignment;

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
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Pga<const M: i8, const N: u32> {
    idx: u8,
}

/// Flavor-specific methods.
impl<const M: i8, const N: u32> Pga<M, N> {
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

impl<const M: i8, const N: u32> Default for Pga<M, N> {
    fn default() -> Self {
        const { Self::const_assert() };
        Self { idx: 0 }
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
        const { Self::const_assert() };
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
        let sig = if cnt.is_odd() { -1 } else { 1 };
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

mod n0;
mod n1;
mod n2;
mod n3;
mod n4;
mod n5;
mod n6;
mod n7;

use n0::LUT0;
use n1::LUT1;
use n2::LUT2;
use n3::LUT3;
use n4::LUT4;
use n5::LUT5;
use n6::LUT6;
use n7::LUT7;

use n0::TAB0;
use n1::TAB1;
use n2::TAB2;
use n3::TAB3;
use n4::TAB4;
use n5::TAB5;
use n6::TAB6;
use n7::TAB7;

const TAB: [&[BasisBlade]; 8] = [&TAB0, &TAB1, &TAB2, &TAB3, &TAB4, &TAB5, &TAB6, &TAB7];
const LUT: [&[BasisBlade]; 8] = [&LUT0, &LUT1, &LUT2, &LUT3, &LUT4, &LUT5, &LUT6, &LUT7];

#[test]
fn dim() {
    let s = |l: &[BasisBlade], h: &[BasisBlade]| l.iter().all(|l| h.iter().any(|h| h.sym == l.sym));

    assert!(s(TAB0.as_slice(), TAB1.as_slice()));
    assert!(s(TAB1.as_slice(), TAB2.as_slice()));

    assert!(!s(TAB2.as_slice(), TAB3.as_slice()));
    assert!(!s(TAB3.as_slice(), TAB4.as_slice()));
    assert!(!s(TAB4.as_slice(), TAB5.as_slice()));
    assert!(!s(TAB5.as_slice(), TAB6.as_slice()));
    assert!(!s(TAB6.as_slice(), TAB7.as_slice()));
}

#[test]
#[allow(clippy::too_many_lines, clippy::needless_raw_string_hashes)]
fn cnv() {
    use super::Tree;
    {
        use super::PgaP4 as Vee;

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
    {
        use super::PgaP3 as Vee;

        let norm = Vee::norm().eval([(('V', "e0123"), 2)]);
        assert_eq!(
            Tree::from(norm.clone()),
            Tree::with_factorization(norm, false)
        );

        format_eq!("{}", Vee::norm().eval([(('v', "e"), 1)]), ["+1", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(('v', "e"), -1)]), ["-1", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(('v', "e"), 2)]), ["+2", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(('V', "e0123"), 1)]), ["+v", "+I"]);
        format_eq!("{}", Vee::norm().eval([(('V', "e0123"), -1)]), ["+v", "-I"]);
        format_eq!("{}", Vee::norm().eval([(('V', "e0123"), 2)]), ["+v", "+2I"]);

        format_eq!("{:0}\n", Vee::norm().eval([(('v', "e"), 1)]), ["+1+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(('v', "e"), -1)]), ["-1+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(('v', "e"), 2)]), ["+2+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(('V', "e0123"), 1)]), ["+v+I"]);
        format_eq!("{:0}\n", Vee::norm().eval([(('V', "e0123"), -1)]), ["+v-I"]);
        format_eq!("{:0}\n", Vee::norm().eval([(('V', "e0123"), 2)]), ["+v+2I"]);

        format_eq!("{:<}", Vee::norm().eval([(('v', "e"), 1)]), ["1", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(('v', "e"), -1)]), ["-1", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(('v', "e"), 2)]), ["2", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(('V', "e0123"), 1)]), ["v", "+I"]);
        format_eq!(
            "{:<}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            ["v", "-I"]
        );
        format_eq!(
            "{:<}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            ["v", "+2I"]
        );

        format_eq!(
            "{:#}",
            Vee::norm().eval([(('v', "e"), 1)]),
            ["+1", "+v0123*I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(('v', "e"), -1)]),
            ["-1", "+v0123*I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(('v', "e"), 2)]),
            ["+2", "+v0123*I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(('V', "e0123"), 1)]),
            ["+v", "+I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            ["+v", "-I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            ["+v", "+2*I"]
        );

        format_eq!(
            "    {:$^4}",
            Vee::norm().eval([(('v', "e"), 1)]),
            [
                r"    \begin{aligned}[t]",
                r"      1 & \\",
                r"      + v_{0123} & \boldsymbol{I}",
                r"    \end{aligned}",
            ]
        );
        format_eq!(
            "    {:$<4}",
            Vee::norm().eval([(('v', "e"), 1)]),
            [r"      1 & \\", r"      + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<0}\n",
            Vee::norm().eval([(('v', "e"), 1)]),
            [r"1 + v_{0123} \boldsymbol{I}"]
        );

        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('v', "e"), 1)]),
            [r"  1 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('v', "e"), -1)]),
            [r"  -1 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('v', "e"), 2)]),
            [r"  2 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('V', "e0123"), 1)]),
            [r"  v & \\", r"  + & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            [r"  v & \\", r"  - & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            [r"  v & \\", r"  + 2 & \boldsymbol{I}"]
        );

        format_eq!(
            "{:>#4x}",
            Vee::norm().eval([(('v', "e"), 1)]),
            ["    let e = 1.0;", "    let e0123 = v0123;"]
        );

        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), 1)]),
            ["let e = 1.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), -1)]),
            ["let e = -1.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), 2)]),
            ["let e = 2.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('V', "e0123"), 1)]),
            ["let e = v;", "let e0123 = 1.0;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            ["let e = v;", "let e0123 = -1.0;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            ["let e = v;", "let e0123 = 2.0;"]
        );

        format_eq!(
            "{:#x}",
            Vee::plane()
                .norm_squared()
                .eval([(('x', "e1"), 2), (('y', "e2"), 3)]),
            ["let e = 13.0 + v3 * v3;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), 0)]),
            ["let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), (1, 2))]),
            ["let e = 1.0 / 2.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(('v', "e"), (-1, 2))]),
            ["let e = -1.0 / 2.0;", "let e0123 = v0123;"]
        );

        format_eq!(
            "{:^4x}",
            Vee::norm().eval([(('v', "e"), 1)]),
            ["    o.e = 1.0;", "    o.e0123 = v.e0123;"]
        );

        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('v', "e"), 1)]),
            ["o.e = 1.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('v', "e"), -1)]),
            ["o.e = -1.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('v', "e"), 2)]),
            ["o.e = 2.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('V', "e0123"), 1)]),
            ["o.e = v.e;", "o.e0123 = 1.0;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            ["o.e = v.e;", "o.e0123 = -1.0;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            ["o.e = v.e;", "o.e0123 = 2.0;"]
        );

        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('v', "e"), 1)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="1" shape=diamond]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="*" shape=box]"#,
                r#"  n0 -> n2"#,
                r#"  n3 [label="v0123" shape=ellipse]"#,
                r#"  n2 -> n3"#,
                r#"  n4 [label="I" shape=diamond]"#,
                r#"  n2 -> n4"#,
                r#"}"#,
            ]
        );
        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('v', "e"), -1)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="*" shape=box]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="-1" shape=circle]"#,
                r#"  n1 -> n2"#,
                r#"  n3 [label="1" shape=diamond]"#,
                r#"  n1 -> n3"#,
                r#"  n4 [label="*" shape=box]"#,
                r#"  n0 -> n4"#,
                r#"  n5 [label="v0123" shape=ellipse]"#,
                r#"  n4 -> n5"#,
                r#"  n6 [label="I" shape=diamond]"#,
                r#"  n4 -> n6"#,
                r#"}"#,
            ]
        );
        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('v', "e"), 2)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="*" shape=box]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="2" shape=circle]"#,
                r#"  n1 -> n2"#,
                r#"  n3 [label="1" shape=diamond]"#,
                r#"  n1 -> n3"#,
                r#"  n4 [label="*" shape=box]"#,
                r#"  n0 -> n4"#,
                r#"  n5 [label="v0123" shape=ellipse]"#,
                r#"  n4 -> n5"#,
                r#"  n6 [label="I" shape=diamond]"#,
                r#"  n4 -> n6"#,
                r#"}"#,
            ]
        );
        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('V', "e0123"), 1)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="*" shape=box]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="v" shape=ellipse]"#,
                r#"  n1 -> n2"#,
                r#"  n3 [label="1" shape=diamond]"#,
                r#"  n1 -> n3"#,
                r#"  n4 [label="I" shape=diamond]"#,
                r#"  n0 -> n4"#,
                r#"}"#,
            ]
        );
        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('V', "e0123"), -1)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="*" shape=box]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="v" shape=ellipse]"#,
                r#"  n1 -> n2"#,
                r#"  n3 [label="1" shape=diamond]"#,
                r#"  n1 -> n3"#,
                r#"  n4 [label="*" shape=box]"#,
                r#"  n0 -> n4"#,
                r#"  n5 [label="-1" shape=circle]"#,
                r#"  n4 -> n5"#,
                r#"  n6 [label="I" shape=diamond]"#,
                r#"  n4 -> n6"#,
                r#"}"#,
            ]
        );
        format_eq!(
            "{:#o}",
            Vee::norm().eval([(('V', "e0123"), 2)]),
            [
                r#"digraph vee {"#,
                r#"  n0 [label="+" shape=box]"#,
                r#"  n1 [label="*" shape=box]"#,
                r#"  n0 -> n1"#,
                r#"  n2 [label="v" shape=ellipse]"#,
                r#"  n1 -> n2"#,
                r#"  n3 [label="1" shape=diamond]"#,
                r#"  n1 -> n3"#,
                r#"  n4 [label="*" shape=box]"#,
                r#"  n0 -> n4"#,
                r#"  n5 [label="2" shape=circle]"#,
                r#"  n4 -> n5"#,
                r#"  n6 [label="I" shape=diamond]"#,
                r#"  n4 -> n6"#,
                r#"}"#,
            ]
        );
    }
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
                let cnt = if not.cnt(*b).is_odd() && sym.len() > 2 {
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
