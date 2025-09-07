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
    /// use vee::PgaP1 as Vee;
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// assert_eq!(format!("{translator:#}"), concat!(
    ///     "+1LwRw\n",
    ///     "+(+1LWRw-1LwRW)I\n",
    /// ));
    ///
    /// let point = Vee::point().pin() << Vee::translator();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+(+[+1vv]~W+[+2Vv]~w)e0\n",
    ///     "+[+1vv]~we1\n",
    /// ));
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
    /// use vee::PgaP2 as Vee;
    ///
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// assert_eq!(format!("{line:#}"), concat!(
    ///     "+(+1LXRY-1LYRX)e0\n",
    ///     "+(+1LYRw-1LwRY)e1\n",
    ///     "+(-1LXRw+1LwRX)e2\n",
    /// ));
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
    /// use vee::PgaP2 as Vee;
    ///
    /// let point = Vee::line().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LWRy+1LyRW)e20\n",
    ///     "+(+1LWRx-1LxRW)e01\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + P_0`$.
    ///
    /// ```
    /// use vee::PgaP2 as Vee;
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// assert_eq!(format!("{rotator:#}"), concat!(
    ///     "+1LxRx+1LyRy\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + P_\infty`$.
    ///
    /// ```
    /// use vee::PgaP2 as Vee;
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// assert_eq!(format!("{translator:#}"), concat!(
    ///     "-1LwRw\n",
    ///     "+(-1LYRw+1LwRY)e20\n",
    ///     "+(+1LXRw-1LwRX)e01\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::direction()
    }
    /// The multivector of motor $`m \equiv s + P`$.
    ///
    /// ```
    /// use vee::PgaP2 as Vee;
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// assert_eq!(format!("{motor:#}"), concat!(
    ///     "+1LxRx+1LyRy\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LWRy+1LyRW)e20\n",
    ///     "+(+1LWRx-1LxRW)e01\n",
    /// ));
    ///
    /// let point = Vee::point().pin() << Vee::motor();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+(+[+1vv+1ww]~w)e12\n",
    ///     "+(+[+2vw]~Y+[+2Xw-2Yv]~w+[+1vv-1ww]~X)e20\n",
    ///     "+(+[+2Xv+2Yw]~w+[-2vw]~X+[+1vv-1ww]~Y)e01\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::scalar() + Self::point()
    }
    /// The multivector of flector $`f \equiv \ell + S`$.
    ///
    /// ```
    /// use vee::PgaP2 as Vee;
    ///
    /// let flector = Vee::line().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// assert_eq!(format!("{flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRY+1LyRX)e0\n",
    ///     "+(+1LxRv-1LyRw)e1\n",
    ///     "+(+1LxRw+1LyRv)e2\n",
    ///     "+(+1LWRw+1LxRX+1LyRY)I\n",
    /// ));
    ///
    /// let point = Vee::point().pin() << Vee::flector();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+(+[-1xx-1yy]~w)e12\n",
    ///     "+(+[+2xy]~Y+[+2Vy+2Wx]~w+[+1xx-1yy]~X)e20\n",
    ///     "+(+[-2Vx+2Wy]~w+[+2xy]~X+[-1xx+1yy]~Y)e01\n",
    /// ));
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
    /// use vee::PgaP3 as Vee;
    ///
    /// let squared_norm = Vee::line().squared_norm();
    ///
    /// assert_eq!(squared_norm.basis_blades(), Vee::norm().basis_blades());
    /// assert_eq!(format!("{squared_norm:#}"), concat!(
    ///     "+1xx+1yy+1zz\n",
    ///     "+(-2Xx-2Yy-2Zz)I\n",
    /// ));
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
    /// use vee::PgaP3 as Vee;
    ///
    /// let plane = Vee::line().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    /// assert_eq!(format!("{plane:#}"), concat!(
    ///     "+(-1LXRX-1LYRY-1LZRZ)e0\n",
    ///     "+(+1LXRw+1LyRZ-1LzRY)e1\n",
    ///     "+(+1LYRw-1LxRZ+1LzRX)e2\n",
    ///     "+(+1LZRw+1LxRY-1LyRX)e3\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn plane() -> Self {
        Self::bias() + Self::normal()
    }
    /// The multivector of displacement $`\ell_0 \equiv x\e_{23} + y\e_{31} + z\e_{12}`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// // A line through the origin as the join of a point and the origin.
    /// let displacement = Vee::point().lhs() & Vee::weight().rhs();
    ///
    /// assert_eq!(displacement.basis_blades(), Vee::displacement().basis_blades());
    /// assert_eq!(format!("{displacement:#}"), concat!(
    ///     "-1LXRwe23\n",
    ///     "-1LYRwe31\n",
    ///     "-1LZRwe12\n",
    /// ));
    ///
    /// // A line through the origin as the meet of two planes through the origin.
    /// let displacement = Vee::normal().lhs() ^ Vee::normal().rhs();
    ///
    /// assert_eq!(displacement.basis_blades(), Vee::displacement().basis_blades());
    /// assert_eq!(format!("{displacement:#}"), concat!(
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    /// ));
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
    /// use vee::PgaP3 as Vee;
    ///
    /// // A line as the join of two points.
    /// let line = Vee::point().lhs() & Vee::point().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// assert_eq!(format!("{line:#}"), concat!(
    ///     "+(+1LYRZ-1LZRY)e01\n",
    ///     "+(-1LXRZ+1LZRX)e02\n",
    ///     "+(+1LXRY-1LYRX)e03\n",
    ///     "+(-1LXRw+1LwRX)e23\n",
    ///     "+(-1LYRw+1LwRY)e31\n",
    ///     "+(-1LZRw+1LwRZ)e12\n",
    /// ));
    ///
    /// // A line as the meet of two planes.
    /// let line = Vee::plane().lhs() ^ Vee::plane().rhs();
    ///
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// assert_eq!(format!("{line:#}"), concat!(
    ///     "+(+1LWRx-1LxRW)e01\n",
    ///     "+(+1LWRy-1LyRW)e02\n",
    ///     "+(+1LWRz-1LzRW)e03\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    /// ));
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
    /// use vee::PgaP3 as Vee;
    ///
    /// let point = Vee::plane().lhs() ^ Vee::line().rhs();
    ///
    /// assert_eq!(point.basis_blades(), Vee::point().basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+(+1LxRx+1LyRy+1LzRz)e123\n",
    ///     "+(-1LWRx+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LWRy-1LxRZ+1LzRX)e013\n",
    ///     "+(-1LWRz+1LxRY-1LyRX)e021\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn point() -> Self {
        Self::weight() + Self::direction()
    }
    /// The multivector of rotator $`r \equiv s + \ell_0`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// assert_eq!(format!("{rotator:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    /// ));
    ///
    /// let rotator = Vee::displacement().lhs() * Vee::displacement().rhs();
    ///
    /// assert_eq!(rotator.basis_blades(), Vee::rotator().basis_blades());
    /// assert_eq!(format!("{rotator:#}"), concat!(
    ///     "-1LxRx-1LyRy-1LzRz\n",
    ///     "+(-1LyRz+1LzRy)e23\n",
    ///     "+(+1LxRz-1LzRx)e31\n",
    ///     "+(-1LxRy+1LyRx)e12\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn rotator() -> Self {
        Self::scalar() + Self::displacement()
    }
    /// The multivector of translator $`t \equiv s + \ell_\infty`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// assert_eq!(format!("{translator:#}"), concat!(
    ///     "-1LwRw\n",
    ///     "+(+1LXRw-1LwRX)e01\n",
    ///     "+(+1LYRw-1LwRY)e02\n",
    ///     "+(+1LZRw-1LwRZ)e03\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + \ell`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let simple_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(simple_motor.basis_blades(), Vee::simple_motor().basis_blades());
    /// assert_eq!(format!("{simple_motor:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz\n",
    ///     "+(+1LWRx-1LxRW)e01\n",
    ///     "+(+1LWRy-1LyRW)e02\n",
    ///     "+(+1LWRz-1LzRW)e03\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::line()
    }
    /// The multivector of motor $`m \equiv s + \ell + S`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// assert_eq!(format!("{motor:#}"), concat!(
    ///     "-1LxRx-1LyRy-1LzRz\n",
    ///     "+(-1LYRz+1LZRy-1LyRZ+1LzRY)e01\n",
    ///     "+(+1LXRz-1LZRx+1LxRZ-1LzRX)e02\n",
    ///     "+(-1LXRy+1LYRx-1LxRY+1LyRX)e03\n",
    ///     "+(-1LyRz+1LzRy)e23\n",
    ///     "+(+1LxRz-1LzRx)e31\n",
    ///     "+(-1LxRy+1LyRx)e12\n",
    ///     "+(+1LXRx+1LYRy+1LZRz+1LxRX+1LyRY+1LzRZ)I\n",
    /// ));
    ///
    /// let motor = Vee::rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(motor.basis_blades(), Vee::motor().basis_blades());
    /// assert_eq!(format!("{motor:#}"), concat!(
    ///     "+1LvRv\n",
    ///     "+(+1LvRX-1LyRZ+1LzRY)e01\n",
    ///     "+(+1LvRY+1LxRZ-1LzRX)e02\n",
    ///     "+(+1LvRZ-1LxRY+1LyRX)e03\n",
    ///     "+1LxRve23\n",
    ///     "+1LyRve31\n",
    ///     "+1LzRve12\n",
    ///     "+(+1LxRX+1LyRY+1LzRZ)I\n"
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn motor() -> Self {
        Self::norm() + Self::line()
    }
    /// The multivector of rotoreflector $`f_r \equiv p_0 + P_0`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// assert_eq!(format!("{rotoreflector:#}"), concat!(
    ///     "+(+1LxRv-1LyRz+1LzRy)e1\n",
    ///     "+(+1LxRz+1LyRv-1LzRx)e2\n",
    ///     "+(-1LxRy+1LyRx+1LzRv)e3\n",
    ///     "+(+1LxRx+1LyRy+1LzRz)e123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv p + P_\infty`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// assert_eq!(format!("{transflector:#}"), concat!(
    ///     "+(-1LxRX-1LyRY-1LzRZ)e0\n",
    ///     "+1LxRve1\n+1LyRve2\n",
    ///     "+1LzRve3\n",
    ///     "+(+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LxRZ+1LzRX)e013\n",
    ///     "+(+1LxRY-1LyRX)e021\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::plane() + Self::direction()
    }
    /// The multivector of flector $`f \equiv p + P`$.
    ///
    /// ```
    /// use vee::PgaP3 as Vee;
    ///
    /// let flector = Vee::plane().lhs() * Vee::motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// assert_eq!(format!("{flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ)e0\n",
    ///     "+(+1LxRv-1LyRz+1LzRy)e1\n",
    ///     "+(+1LxRz+1LyRv-1LzRx)e2\n",
    ///     "+(-1LxRy+1LyRx+1LzRv)e3\n",
    ///     "+(+1LxRx+1LyRy+1LzRz)e123\n",
    ///     "+(-1LWRx+1LxRV+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LWRy-1LxRZ+1LyRV+1LzRX)e013\n",
    ///     "+(-1LWRz+1LxRY-1LyRX+1LzRV)e021\n",
    /// ));
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
    /// use vee::PgaP4 as Vee;
    ///
    /// let quadvector_squared_norm = Vee::point().squared_norm();
    ///
    /// assert_eq!(format!("{quadvector_squared_norm:#}"), concat!(
    ///     "+1ww\n",
    /// ));
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
    /// $`\ell_0 \equiv x\e_{234} + y\e_{314} + z\e_{124} + þ\e_{123}`$.
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
    /// $`P_\infty \equiv X\e_{0324} + Y\e_{0134} + Z\e_{0214} + Þ\e_{0123}`$.
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
    /// use vee::PgaP4 as Vee;
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// assert_eq!(format!("{single_rotator:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LþRþ\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    /// ));
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// assert_eq!(format!("{single_rotator:#}"), concat!(
    ///     "-1LxRx-1LyRy-1LzRz-1LþRþ\n",
    ///     "+(-1LyRz+1LzRy)e23\n",
    ///     "+(+1LxRz-1LzRx)e31\n",
    ///     "+(-1LxRy+1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    /// ));
    ///
    /// let squared_norm = Vee::single_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// assert_eq!(format!("{squared_norm:#}"), concat!(
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv\n",
    ///     "+(+2ad+2be+2cf)e1234\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::plane_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + \ell_0 + P_0`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// assert_eq!(format!("{double_rotator:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf+1LvRv\n",
    ///     "+(+1LaRv-1LbRc+1LcRb-1LeRf+1LfRe+1LvRa)e23\n",
    ///     "+(+1LaRc+1LbRv-1LcRa+1LdRf-1LfRd+1LvRb)e31\n",
    ///     "+(-1LaRb+1LbRa+1LcRv-1LdRe+1LeRd+1LvRc)e12\n",
    ///     "+(-1LbRf+1LcRe+1LdRv-1LeRc+1LfRb+1LvRd)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc+1LeRv-1LfRa+1LvRe)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa+1LfRv+1LvRf)e43\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    /// ));
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// assert_eq!(format!("{double_rotator:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf\n",
    ///     "+(-1LbRc+1LcRb-1LeRf+1LfRe)e23\n",
    ///     "+(+1LaRc-1LcRa+1LdRf-1LfRd)e31\n",
    ///     "+(-1LaRb+1LbRa-1LdRe+1LeRd)e12\n",
    ///     "+(-1LbRf+1LcRe-1LeRc+1LfRb)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc-1LfRa)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa)e43\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    /// ));
    ///
    /// let squared_norm = Vee::double_rotator().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// assert_eq!(format!("{squared_norm:#}"), concat!(
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv+1ww\n",
    ///     "+(+2ad+2be+2cf+2vw)e1234\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of translator $`t \equiv s + p_\infty`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// assert_eq!(format!("{translator:#}"), concat!(
    ///     "+1LwRw\n",
    ///     "+(-1LXRw+1LwRX)e01\n",
    ///     "+(-1LYRw+1LwRY)e02\n",
    ///     "+(-1LZRw+1LwRZ)e03\n",
    ///     "+(-1LwRÞ+1LÞRw)e40\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::plane_moment()
    }
    /// The multivector of simple motor $`m_s \equiv s + p`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let simple_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(simple_motor.basis_blades(), Vee::simple_motor().basis_blades());
    /// assert_eq!(format!("{simple_motor:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LþRþ\n",
    ///     "+(+1LWRx-1LxRW)e01\n",
    ///     "+(+1LWRy-1LyRW)e02\n",
    ///     "+(+1LWRz-1LzRW)e03\n",
    ///     "+(-1LWRþ+1LþRW)e40\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    /// ));
    ///
    /// let squared_norm = Vee::simple_motor().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), Vee::norm().basis_blades());
    /// assert_eq!(format!("{squared_norm:#}"), concat!(
    ///     // Scalar condition.
    ///     "+1aa+1bb+1cc+1dd+1ee+1ff+1vv\n",
    ///     // Point condition.
    ///     "+(+2ad+2be+2cf)e1234\n", // Weight condition.
    ///     "+(-2Yf+2Ze-2aÞ)e0324\n", // Direction condition.
    ///     "+(+2Xf-2Zd-2bÞ)e0134\n", // Direction condition.
    ///     "+(-2Xe+2Yd-2cÞ)e0214\n", // Direction condition.
    ///     "+(-2Xa-2Yb-2Zc)e0123\n", // Direction condition.
    /// ));
    ///
    /// let point = Vee::point().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+[+2ad+2be+2cf]~w\n", // Vanishes with weight condition.
    ///     "+(+[+1aa+1bb+1cc+1dd+1ee+1ff+1vv]~w)e1234\n",
    ///     "+(+[+2ac-2bv-2df]~Z+[-2Xv-2Yc+2Zb-2dÞ]~w+[+2bf-2ce-2dv]~Þ",
    ///       "+[+1aa-1bb-1cc-1dd+1ee+1ff+1vv]~X+[+2ab+2cv-2de]~Y)e0324\n",
    ///     "+(+[+2Xc-2Yv-2Za-2eÞ]~w+[-2af+2cd-2ev]~Þ+[+2ab-2cv-2de]~X",
    ///       "+[-1aa+1bb-1cc+1dd-1ee+1ff+1vv]~Y+[+2av+2bc-2ef]~Z)e0134\n",
    ///     "+(+[+2ae-2bd-2fv]~Þ+[+2ac+2bv-2df]~X+[-2av+2bc-2ef]~Y",
    ///       "+[-1aa-1bb+1cc+1dd+1ee-1ff+1vv]~Z+[-2Xb+2Ya-2Zv-2fÞ]~w)e0214\n",
    ///     "+(+[+2bf-2ce+2dv]~X+[-2af+2cd+2ev]~Y+[+2ae-2bd+2fv]~Z",
    ///       "+[-2Xd-2Ye-2Zf+2vÞ]~w+[+1aa+1bb+1cc-1dd-1ee-1ff+1vv]~Þ)e0123\n",
    /// ));
    ///
    /// let line = Vee::line().pin() << Vee::simple_motor();
    /// assert_eq!(line.basis_blades(), Vee::line().basis_blades());
    /// let plane = Vee::plane().pin() << Vee::simple_motor();
    /// assert_eq!(plane.basis_blades(), Vee::plane().basis_blades());
    ///
    /// let volume = Vee::volume().pin() << Vee::simple_motor();
    ///
    /// assert_eq!(volume.basis_blades(), (Vee::pseudoscalar() + Vee::volume()).basis_blades());
    /// assert_eq!(format!("{volume:#}"), concat!(
    ///     "+(+[+1aa+1bb+1cc+1dd+1ee+1ff+1vv]~W+[+2Xv-2Yc+2Zb-2dÞ]~x",
    ///       "+[+2Xc+2Yv-2Za-2eÞ]~y+[-2Xb+2Ya+2Zv-2fÞ]~z+[-2Xd-2Ye-2Zf-2vÞ]~þ)e0\n",
    ///     "+(+[+2ab+2cv-2de]~y+[+2ac-2bv-2df]~z+[+2bf-2ce-2dv]~þ",
    ///       "+[+1aa-1bb-1cc-1dd+1ee+1ff+1vv]~x)e1\n",
    ///     "+(+[+2av+2bc-2ef]~z+[-2af+2cd-2ev]~þ+[+2ab-2cv-2de]~x",
    ///       "+[-1aa+1bb-1cc+1dd-1ee+1ff+1vv]~y)e2\n",
    ///     "+(+[+2ae-2bd-2fv]~þ+[+2ac+2bv-2df]~x+[-2av+2bc-2ef]~y",
    ///       "+[-1aa-1bb+1cc+1dd+1ee-1ff+1vv]~z)e3\n",
    ///     "+(+[+2bf-2ce+2dv]~x+[-2af+2cd+2ev]~y+[+2ae-2bd+2fv]~z",
    ///       "+[+1aa+1bb+1cc-1dd-1ee-1ff+1vv]~þ)e4\n",
    ///     "+(+[+2ad+2be+2cf]~W+[-2Yf+2Ze-2aÞ]~x+[+2Xf-2Zd-2bÞ]~y",
    ///       "+[-2Xe+2Yd-2cÞ]~z+[-2Xa-2Yb-2Zc]~þ)I\n", // Vanishes with point condition.
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_motor() -> Self {
        Self::scalar() + Self::plane()
    }
    /// The multivector of single motor $`m_1 \equiv s + p + P_\infty`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// assert_eq!(format!("{single_motor:#}"), concat!(
    ///     "+1LvRv\n",
    ///     "+(-1LbRZ+1LcRY+1LdRÞ+1LvRX)e01\n",
    ///     "+(+1LaRZ-1LcRX+1LeRÞ+1LvRY)e02\n",
    ///     "+(-1LaRY+1LbRX+1LfRÞ+1LvRZ)e03\n",
    ///     "+(-1LdRX-1LeRY-1LfRZ+1LvRÞ)e40\n",
    ///     "+1LaRve23\n",
    ///     "+1LbRve31\n",
    ///     "+1LcRve12\n",
    ///     "+1LdRve41\n",
    ///     "+1LeRve42\n",
    ///     "+1LfRve43\n",
    ///     "+(+1LaRÞ-1LeRZ+1LfRY)e0324\n",
    ///     "+(+1LbRÞ+1LdRZ-1LfRX)e0134\n",
    ///     "+(+1LcRÞ-1LdRY+1LeRX)e0214\n",
    ///     "+(+1LaRX+1LbRY+1LcRZ)e0123\n",
    /// ));
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// assert_eq!(format!("{single_motor:#}"), concat!(
    ///     "-1LxRx-1LyRy-1LzRz-1LþRþ\n",
    ///     "+(-1LBRz+1LCRy+1LDRþ-1LyRC+1LzRB-1LþRD)e01\n",
    ///     "+(+1LARz-1LCRx+1LERþ+1LxRC-1LzRA-1LþRE)e02\n",
    ///     "+(-1LARy+1LBRx+1LFRþ-1LxRB+1LyRA-1LþRF)e03\n",
    ///     "+(-1LDRx-1LERy-1LFRz+1LxRD+1LyRE+1LzRF)e40\n",
    ///     "+(-1LyRz+1LzRy)e23\n",
    ///     "+(+1LxRz-1LzRx)e31\n",
    ///     "+(-1LxRy+1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    ///     "+(+1LARþ-1LERz+1LFRy+1LyRF-1LzRE+1LþRA)e0324\n",
    ///     "+(+1LBRþ+1LDRz-1LFRx-1LxRF+1LzRD+1LþRB)e0134\n",
    ///     "+(+1LCRþ-1LDRy+1LERx+1LxRE-1LyRD+1LþRC)e0214\n",
    ///     "+(+1LARx+1LBRy+1LCRz+1LxRA+1LyRB+1LzRC)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::plane() + Self::direction()
    }
    /// The multivector of double motor $`m_2 \equiv s + p + P`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// assert_eq!(format!("{double_motor:#}"), concat!(
    ///     "+1LvRv\n",
    ///     "+(-1LbRZ+1LcRY+1LdRÞ+1LvRX)e01\n",
    ///     "+(+1LaRZ-1LcRX+1LeRÞ+1LvRY)e02\n",
    ///     "+(-1LaRY+1LbRX+1LfRÞ+1LvRZ)e03\n",
    ///     "+(-1LdRX-1LeRY-1LfRZ+1LvRÞ)e40\n",
    ///     "+1LaRve23\n",
    ///     "+1LbRve31\n",
    ///     "+1LcRve12\n",
    ///     "+1LdRve41\n",
    ///     "+1LeRve42\n",
    ///     "+1LfRve43\n",
    ///     "+1LwRve1234\n",
    ///     "+(+1LaRÞ-1LeRZ+1LfRY+1LwRX)e0324\n",
    ///     "+(+1LbRÞ+1LdRZ-1LfRX+1LwRY)e0134\n",
    ///     "+(+1LcRÞ-1LdRY+1LeRX+1LwRZ)e0214\n",
    ///     "+(+1LaRX+1LbRY+1LcRZ-1LwRÞ)e0123\n",
    /// ));
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// assert_eq!(format!("{double_motor:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf\n",
    ///     "+(-1LYRc+1LZRb-1LbRZ+1LcRY+1LdRÞ-1LÞRd)e01\n",
    ///     "+(+1LXRc-1LZRa+1LaRZ-1LcRX+1LeRÞ-1LÞRe)e02\n",
    ///     "+(-1LXRb+1LYRa-1LaRY+1LbRX+1LfRÞ-1LÞRf)e03\n",
    ///     "+(+1LXRd+1LYRe+1LZRf-1LdRX-1LeRY-1LfRZ)e40\n",
    ///     "+(-1LbRc+1LcRb-1LeRf+1LfRe)e23\n",
    ///     "+(+1LaRc-1LcRa+1LdRf-1LfRd)e31\n",
    ///     "+(-1LaRb+1LbRa-1LdRe+1LeRd)e12\n",
    ///     "+(-1LbRf+1LcRe-1LeRc+1LfRb)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc-1LfRa)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa)e43\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    ///     "+(+1LYRf-1LZRe+1LaRÞ-1LeRZ+1LfRY+1LÞRa)e0324\n",
    ///     "+(-1LXRf+1LZRd+1LbRÞ+1LdRZ-1LfRX+1LÞRb)e0134\n",
    ///     "+(+1LXRe-1LYRd+1LcRÞ-1LdRY+1LeRX+1LÞRc)e0214\n",
    ///     "+(+1LXRa+1LYRb+1LZRc+1LaRX+1LbRY+1LcRZ)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::plane() + Self::point()
    }
    /// The multivector of rotoreflector $`f_r \equiv v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(rotoreflector.basis_blades(), Vee::rotoreflector().basis_blades());
    /// assert_eq!(format!("{rotoreflector:#}"), concat!(
    ///     "+(+1LxRv-1LyRc+1LzRb+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf+1LþRv)e4\n",
    ///     "+(-1LyRf+1LzRe+1LþRa)e234\n",
    ///     "+(+1LxRf-1LzRd+1LþRb)e314\n",
    ///     "+(-1LxRe+1LyRd+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc)e123\n",
    /// ));
    ///
    /// let squared_norm = Vee::rotoreflector().squared_norm();
    /// assert_eq!(squared_norm.basis_blades(), (Vee::scalar() + Vee::weight()).basis_blades());
    /// assert_eq!(format!("{squared_norm:#}"), concat!(
    ///     "+2xx+2yy+2zz+2þþ\n",
    ///     "+(-2xx-2yy-2zz+2þþ)e1234\n", // Weight condition.
    /// ));
    ///
    /// let point = Vee::point().pin() << Vee::rotoreflector();
    ///
    /// assert_eq!(point.basis_blades(), (Vee::scalar() + Vee::point()).basis_blades());
    /// assert_eq!(format!("{point:#}"), concat!(
    ///     "+[+2xx+2yy+2zz-2þþ]~w\n", // Vanishes with weight condition.
    ///     "+(+[-2xx-2yy-2zz-2þþ]~w)e1234\n",
    ///     "+(+[+4xþ]~Þ+[-4zþ]~Y+[+4yþ]~Z)e0324\n",
    ///     "+(+[+4zþ]~X+[-4xþ]~Z+[+4yþ]~Þ)e0134\n",
    ///     "+(+[+4xþ]~Y+[+4zþ]~Þ+[-4yþ]~X)e0214\n",
    ///     "+(+[+4zþ]~Z+[+4xþ]~X+[+4yþ]~Y)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn rotoreflector() -> Self {
        Self::normal() + Self::line_displacement()
    }
    /// The multivector of transflector $`f_t \equiv v + P_\infty`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// assert_eq!(format!("{transflector:#}"), concat!(
    ///     "+(-1LxRX-1LyRY-1LzRZ+1LþRÞ)e0\n",
    ///     "+1LxRve1\n",
    ///     "+1LyRve2\n",
    ///     "+1LzRve3\n",
    ///     "+1LþRve4\n",
    ///     "+(+1LxRÞ+1LþRX)e014\n",
    ///     "+(+1LyRÞ+1LþRY)e024\n",
    ///     "+(+1LzRÞ+1LþRZ)e034\n",
    ///     "+(+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LxRZ+1LzRX)e013\n",
    ///     "+(+1LxRY-1LyRX)e021\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::volume() + Self::line_moment()
    }
    /// The multivector of flector $`f \equiv v + \ell + S`$.
    ///
    /// ```
    /// use vee::PgaP4 as Vee;
    ///
    /// let flector = Vee::volume().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(flector.basis_blades(), Vee::flector().basis_blades());
    /// assert_eq!(format!("{flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ+1LþRÞ)e0\n",
    ///     "+(+1LxRv-1LyRc+1LzRb+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf+1LþRv)e4\n",
    ///     "+(-1LyRf+1LzRe+1LþRa)e234\n",
    ///     "+(+1LxRf-1LzRd+1LþRb)e314\n",
    ///     "+(-1LxRe+1LyRd+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc)e123\n",
    ///     "+(-1LWRd+1LxRÞ-1LyRZ+1LzRY+1LþRX)e014\n",
    ///     "+(-1LWRe+1LxRZ+1LyRÞ-1LzRX+1LþRY)e024\n",
    ///     "+(-1LWRf-1LxRY+1LyRX+1LzRÞ+1LþRZ)e034\n",
    ///     "+(-1LWRa+1LxRÞ+1LyRZ-1LzRY-1LþRX)e032\n",
    ///     "+(-1LWRb-1LxRZ+1LyRÞ+1LzRX-1LþRY)e013\n",
    ///     "+(-1LWRc+1LxRY-1LyRX+1LzRÞ-1LþRZ)e021\n",
    ///     "+(+1LxRX+1LyRY+1LzRZ+1LþRÞ)I\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn flector() -> Self {
        Self::volume() + Self::line() + Self::pseudoscalar()
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
    /// use vee::PgaP5 as Vee;
    ///
    /// let quadvector_squared_norm = Vee::line().squared_norm();
    ///
    /// assert_eq!(quadvector_squared_norm.basis_blades(),
    ///     (Vee::scalar() + Vee::line_moment()).basis_blades());
    /// assert_eq!(format!("{quadvector_squared_norm:#}"), concat!(
    ///     "+1xx+1yy+1zz+1ðð+1þþ\n",
    ///     "+(+2Dð-2Gþ-2Jx)e0145\n",
    ///     "+(+2Eð-2Hþ-2Jy)e0245\n",
    ///     "+(+2Fð-2Iþ-2Jz)e0345\n",
    ///     "+(+2Að-2Hz+2Iy)e0325\n",
    ///     "+(+2Bð+2Gz-2Ix)e0135\n",
    ///     "+(+2Cð-2Gy+2Hx)e0215\n",
    ///     "+(-2Aþ+2Ez-2Fy)e0324\n",
    ///     "+(-2Bþ-2Dz+2Fx)e0134\n",
    ///     "+(-2Cþ+2Dy-2Ex)e0214\n",
    ///     "+(-2Ax-2By-2Cz)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn norm() -> Self {
        Self::scalar() + Self::line()
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
    /// The multivector of single rotator $`r_1 \equiv s + v_0`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let single_rotator = Vee::normal().lhs() * Vee::normal().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// assert_eq!(format!("{single_rotator:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LðRð+1LþRþ\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    ///     "+(+1LxRð-1LðRx)e15\n",
    ///     "+(+1LyRð-1LðRy)e25\n",
    ///     "+(+1LzRð-1LðRz)e35\n",
    ///     "+(-1LðRþ+1LþRð)e45\n",
    /// ));
    ///
    /// let single_rotator = Vee::line_displacement().lhs() * Vee::line_displacement().rhs();
    ///
    /// assert_eq!(single_rotator.basis_blades(), Vee::single_rotator().basis_blades());
    /// assert_eq!(format!("{single_rotator:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LðRð+1LþRþ\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(+1LxRþ-1LþRx)e41\n",
    ///     "+(+1LyRþ-1LþRy)e42\n",
    ///     "+(+1LzRþ-1LþRz)e43\n",
    ///     "+(+1LxRð-1LðRx)e15\n",
    ///     "+(+1LyRð-1LðRy)e25\n",
    ///     "+(+1LzRð-1LðRz)e35\n",
    ///     "+(+1LðRþ-1LþRð)e45\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotator() -> Self {
        Self::scalar() + Self::volume_displacement()
    }
    /// The multivector of double rotator $`r_2 \equiv s + v_0 + \ell_0`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let double_rotator = Vee::single_rotator().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// assert_eq!(format!("{double_rotator:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf-1LgRg-1LhRh-1LiRi-1LjRj+1LvRv\n",
    ///     "+(+1LaRv-1LbRc+1LcRb-1LeRf+1LfRe-1LhRi+1LiRh+1LvRa)e23\n",
    ///     "+(+1LaRc+1LbRv-1LcRa+1LdRf-1LfRd+1LgRi-1LiRg+1LvRb)e31\n",
    ///     "+(-1LaRb+1LbRa+1LcRv-1LdRe+1LeRd-1LgRh+1LhRg+1LvRc)e12\n",
    ///     "+(-1LbRf+1LcRe+1LdRv-1LeRc+1LfRb+1LgRj-1LjRg+1LvRd)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc+1LeRv-1LfRa+1LhRj-1LjRh+1LvRe)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa+1LfRv+1LiRj-1LjRi+1LvRf)e43\n",
    ///     "+(-1LbRi+1LcRh-1LdRj+1LgRv-1LhRc+1LiRb+1LjRd+1LvRg)e15\n",
    ///     "+(+1LaRi-1LcRg-1LeRj+1LgRc+1LhRv-1LiRa+1LjRe+1LvRh)e25\n",
    ///     "+(-1LaRh+1LbRg-1LfRj-1LgRb+1LhRa+1LiRv+1LjRf+1LvRi)e35\n",
    ///     "+(+1LdRg+1LeRh+1LfRi-1LgRd-1LhRe-1LiRf+1LjRv+1LvRj)e45\n",
    ///     "+(+1LaRj+1LeRi-1LfRh-1LhRf+1LiRe+1LjRa)e2345\n",
    ///     "+(+1LbRj-1LdRi+1LfRg+1LgRf-1LiRd+1LjRb)e3145\n",
    ///     "+(+1LcRj+1LdRh-1LeRg-1LgRe+1LhRd+1LjRc)e1245\n",
    ///     "+(+1LaRg+1LbRh+1LcRi+1LgRa+1LhRb+1LiRc)e1235\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    /// ));
    ///
    /// let double_rotator = Vee::volume_displacement().lhs() * Vee::volume_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// assert_eq!(format!("{double_rotator:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf-1LgRg-1LhRh-1LiRi-1LjRj\n",
    ///     "+(-1LbRc+1LcRb-1LeRf+1LfRe-1LhRi+1LiRh)e23\n",
    ///     "+(+1LaRc-1LcRa+1LdRf-1LfRd+1LgRi-1LiRg)e31\n",
    ///     "+(-1LaRb+1LbRa-1LdRe+1LeRd-1LgRh+1LhRg)e12\n",
    ///     "+(-1LbRf+1LcRe-1LeRc+1LfRb+1LgRj-1LjRg)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc-1LfRa+1LhRj-1LjRh)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa+1LiRj-1LjRi)e43\n",
    ///     "+(-1LbRi+1LcRh-1LdRj-1LhRc+1LiRb+1LjRd)e15\n",
    ///     "+(+1LaRi-1LcRg-1LeRj+1LgRc-1LiRa+1LjRe)e25\n",
    ///     "+(-1LaRh+1LbRg-1LfRj-1LgRb+1LhRa+1LjRf)e35\n",
    ///     "+(+1LdRg+1LeRh+1LfRi-1LgRd-1LhRe-1LiRf)e45\n",
    ///     "+(+1LaRj+1LeRi-1LfRh-1LhRf+1LiRe+1LjRa)e2345\n",
    ///     "+(+1LbRj-1LdRi+1LfRg+1LgRf-1LiRd+1LjRb)e3145\n",
    ///     "+(+1LcRj+1LdRh-1LeRg-1LgRe+1LhRd+1LjRc)e1245\n",
    ///     "+(+1LaRg+1LbRh+1LcRi+1LgRa+1LhRb+1LiRc)e1235\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    /// ));
    ///
    /// let double_rotator = Vee::plane_displacement().lhs() * Vee::plane_displacement().rhs();
    ///
    /// assert_eq!(double_rotator.basis_blades(), Vee::double_rotator().basis_blades());
    /// assert_eq!(format!("{double_rotator:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf-1LgRg-1LhRh-1LiRi-1LjRj\n",
    ///     "+(+1LbRc-1LcRb+1LfRg-1LgRf-1LiRj+1LjRi)e23\n",
    ///     "+(+1LaRc-1LcRa+1LeRg-1LgRe+1LhRj-1LjRh)e31\n",
    ///     "+(+1LaRb-1LbRa+1LeRf-1LfRe-1LhRi+1LiRh)e12\n",
    ///     "+(-1LaRd+1LdRa+1LfRj+1LgRi-1LiRg-1LjRf)e41\n",
    ///     "+(+1LbRd-1LdRb+1LeRj-1LgRh+1LhRg-1LjRe)e42\n",
    ///     "+(-1LcRd+1LdRc-1LeRi-1LfRh+1LhRf+1LiRe)e43\n",
    ///     "+(-1LbRj-1LcRi+1LdRe-1LeRd+1LiRc+1LjRb)e15\n",
    ///     "+(-1LaRj+1LcRh-1LdRf+1LfRd-1LhRc+1LjRa)e25\n",
    ///     "+(+1LaRi+1LbRh+1LdRg-1LgRd-1LhRb-1LiRa)e35\n",
    ///     "+(+1LaRe+1LbRf+1LcRg-1LeRa-1LfRb-1LgRc)e45\n",
    ///     "+(-1LbRg+1LcRf+1LdRh+1LfRc-1LgRb+1LhRd)e2345\n",
    ///     "+(-1LaRg+1LcRe+1LdRi+1LeRc-1LgRa+1LiRd)e3145\n",
    ///     "+(-1LaRf+1LbRe+1LdRj+1LeRb-1LfRa+1LjRd)e1245\n",
    ///     "+(-1LaRh+1LbRi-1LcRj-1LhRa+1LiRb-1LjRc)e1235\n",
    ///     "+(-1LeRh+1LfRi-1LgRj-1LhRe+1LiRf-1LjRg)e1234\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotator() -> Self {
        Self::scalar() + Self::volume_displacement() + Self::line_displacement()
    }
    /// The multivector of translator $`t \equiv s + v_\infty`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let translator = Vee::point().lhs() * Vee::point().rhs();
    ///
    /// assert_eq!(translator.basis_blades(), Vee::translator().basis_blades());
    /// assert_eq!(format!("{translator:#}"), concat!(
    ///     "+1LwRw\n",
    ///     "+(-1LXRw+1LwRX)e01\n",
    ///     "+(-1LYRw+1LwRY)e02\n",
    ///     "+(-1LZRw+1LwRZ)e03\n",
    ///     "+(-1LwRÞ+1LÞRw)e40\n",
    ///     "+(+1LwRÐ-1LÐRw)e05\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn translator() -> Self {
        Self::scalar() + Self::volume_moment()
    }
    /// The multivector of simple single motor $`m_{s1} \equiv s + v`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let simple_single_motor = Vee::hypervolume().lhs() * Vee::hypervolume().rhs();
    ///
    /// assert_eq!(simple_single_motor.basis_blades(), Vee::simple_single_motor().basis_blades());
    /// assert_eq!(format!("{simple_single_motor:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LðRð+1LþRþ\n",
    ///     "+(+1LWRx-1LxRW)e01\n",
    ///     "+(+1LWRy-1LyRW)e02\n",
    ///     "+(+1LWRz-1LzRW)e03\n",
    ///     "+(-1LWRþ+1LþRW)e40\n",
    ///     "+(+1LWRð-1LðRW)e05\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(-1LxRþ+1LþRx)e41\n",
    ///     "+(-1LyRþ+1LþRy)e42\n",
    ///     "+(-1LzRþ+1LþRz)e43\n",
    ///     "+(+1LxRð-1LðRx)e15\n",
    ///     "+(+1LyRð-1LðRy)e25\n",
    ///     "+(+1LzRð-1LðRz)e35\n",
    ///     "+(-1LðRþ+1LþRð)e45\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_motor() -> Self {
        Self::scalar() + Self::volume()
    }
    /// The multivector of single motor $`m_1 \equiv s + v + \ell_\infty`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let single_motor = Vee::single_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// assert_eq!(format!("{single_motor:#}"), concat!(
    ///     "+1LvRv\n",
    ///     "+(-1LbRZ+1LcRY+1LdRÞ+1LgRÐ+1LvRX)e01\n",
    ///     "+(+1LaRZ-1LcRX+1LeRÞ+1LhRÐ+1LvRY)e02\n",
    ///     "+(-1LaRY+1LbRX+1LfRÞ+1LiRÐ+1LvRZ)e03\n",
    ///     "+(-1LdRX-1LeRY-1LfRZ-1LjRÐ+1LvRÞ)e40\n",
    ///     "+(-1LgRX-1LhRY-1LiRZ+1LjRÞ+1LvRÐ)e05\n",
    ///     "+1LaRve23\n",
    ///     "+1LbRve31\n",
    ///     "+1LcRve12\n",
    ///     "+1LdRve41\n",
    ///     "+1LeRve42\n",
    ///     "+1LfRve43\n",
    ///     "+1LgRve15\n",
    ///     "+1LhRve25\n",
    ///     "+1LiRve35\n",
    ///     "+1LjRve45\n",
    ///     "+(-1LdRÐ+1LgRÞ+1LjRX)e0145\n",
    ///     "+(-1LeRÐ+1LhRÞ+1LjRY)e0245\n",
    ///     "+(-1LfRÐ+1LiRÞ+1LjRZ)e0345\n",
    ///     "+(-1LaRÐ+1LhRZ-1LiRY)e0325\n",
    ///     "+(-1LbRÐ-1LgRZ+1LiRX)e0135\n",
    ///     "+(-1LcRÐ+1LgRY-1LhRX)e0215\n",
    ///     "+(+1LaRÞ-1LeRZ+1LfRY)e0324\n",
    ///     "+(+1LbRÞ+1LdRZ-1LfRX)e0134\n",
    ///     "+(+1LcRÞ-1LdRY+1LeRX)e0214\n",
    ///     "+(+1LaRX+1LbRY+1LcRZ)e0123\n",
    /// ));
    ///
    /// let single_motor = Vee::line().lhs() * Vee::line().rhs();
    ///
    /// assert_eq!(single_motor.basis_blades(), Vee::single_motor().basis_blades());
    /// assert_eq!(format!("{single_motor:#}"), concat!(
    ///     "+1LxRx+1LyRy+1LzRz+1LðRð+1LþRþ\n",
    ///     "+(+1LBRz-1LCRy-1LDRþ-1LGRð+1LyRC-1LzRB+1LðRG+1LþRD)e01\n",
    ///     "+(-1LARz+1LCRx-1LERþ-1LHRð-1LxRC+1LzRA+1LðRH+1LþRE)e02\n",
    ///     "+(+1LARy-1LBRx-1LFRþ-1LIRð+1LxRB-1LyRA+1LðRI+1LþRF)e03\n",
    ///     "+(+1LDRx+1LERy+1LFRz+1LJRð-1LxRD-1LyRE-1LzRF-1LðRJ)e40\n",
    ///     "+(+1LGRx+1LHRy+1LIRz-1LJRþ-1LxRG-1LyRH-1LzRI+1LþRJ)e05\n",
    ///     "+(+1LyRz-1LzRy)e23\n",
    ///     "+(-1LxRz+1LzRx)e31\n",
    ///     "+(+1LxRy-1LyRx)e12\n",
    ///     "+(+1LxRþ-1LþRx)e41\n",
    ///     "+(+1LyRþ-1LþRy)e42\n",
    ///     "+(+1LzRþ-1LþRz)e43\n",
    ///     "+(+1LxRð-1LðRx)e15\n",
    ///     "+(+1LyRð-1LðRy)e25\n",
    ///     "+(+1LzRð-1LðRz)e35\n",
    ///     "+(+1LðRþ-1LþRð)e45\n",
    ///     "+(+1LDRð-1LGRþ-1LJRx-1LxRJ+1LðRD-1LþRG)e0145\n",
    ///     "+(+1LERð-1LHRþ-1LJRy-1LyRJ+1LðRE-1LþRH)e0245\n",
    ///     "+(+1LFRð-1LIRþ-1LJRz-1LzRJ+1LðRF-1LþRI)e0345\n",
    ///     "+(+1LARð-1LHRz+1LIRy+1LyRI-1LzRH+1LðRA)e0325\n",
    ///     "+(+1LBRð+1LGRz-1LIRx-1LxRI+1LzRG+1LðRB)e0135\n",
    ///     "+(+1LCRð-1LGRy+1LHRx+1LxRH-1LyRG+1LðRC)e0215\n",
    ///     "+(-1LARþ+1LERz-1LFRy-1LyRF+1LzRE-1LþRA)e0324\n",
    ///     "+(-1LBRþ-1LDRz+1LFRx+1LxRF-1LzRD-1LþRB)e0134\n",
    ///     "+(-1LCRþ+1LDRy-1LERx-1LxRE+1LyRD-1LþRC)e0214\n",
    ///     "+(-1LARx-1LBRy-1LCRz-1LxRA-1LyRB-1LzRC)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line_moment()
    }
    /// The multivector of simple double motor $`m_{s2} \equiv s + v + \ell`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let simple_double_motor = Vee::volume().lhs() * Vee::volume().rhs();
    ///
    /// assert_eq!(simple_double_motor.basis_blades(), Vee::simple_double_motor().basis_blades());
    /// assert_eq!(format!("{simple_double_motor:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf-1LgRg-1LhRh-1LiRi-1LjRj\n",
    ///     "+(-1LYRc+1LZRb-1LbRZ+1LcRY+1LdRÞ+1LgRÐ-1LÐRg-1LÞRd)e01\n",
    ///     "+(+1LXRc-1LZRa+1LaRZ-1LcRX+1LeRÞ+1LhRÐ-1LÐRh-1LÞRe)e02\n",
    ///     "+(-1LXRb+1LYRa-1LaRY+1LbRX+1LfRÞ+1LiRÐ-1LÐRi-1LÞRf)e03\n",
    ///     "+(+1LXRd+1LYRe+1LZRf-1LdRX-1LeRY-1LfRZ-1LjRÐ+1LÐRj)e40\n",
    ///     "+(+1LXRg+1LYRh+1LZRi-1LgRX-1LhRY-1LiRZ+1LjRÞ-1LÞRj)e05\n",
    ///     "+(-1LbRc+1LcRb-1LeRf+1LfRe-1LhRi+1LiRh)e23\n",
    ///     "+(+1LaRc-1LcRa+1LdRf-1LfRd+1LgRi-1LiRg)e31\n",
    ///     "+(-1LaRb+1LbRa-1LdRe+1LeRd-1LgRh+1LhRg)e12\n",
    ///     "+(-1LbRf+1LcRe-1LeRc+1LfRb+1LgRj-1LjRg)e41\n",
    ///     "+(+1LaRf-1LcRd+1LdRc-1LfRa+1LhRj-1LjRh)e42\n",
    ///     "+(-1LaRe+1LbRd-1LdRb+1LeRa+1LiRj-1LjRi)e43\n",
    ///     "+(-1LbRi+1LcRh-1LdRj-1LhRc+1LiRb+1LjRd)e15\n",
    ///     "+(+1LaRi-1LcRg-1LeRj+1LgRc-1LiRa+1LjRe)e25\n",
    ///     "+(-1LaRh+1LbRg-1LfRj-1LgRb+1LhRa+1LjRf)e35\n",
    ///     "+(+1LdRg+1LeRh+1LfRi-1LgRd-1LhRe-1LiRf)e45\n",
    ///     "+(+1LaRj+1LeRi-1LfRh-1LhRf+1LiRe+1LjRa)e2345\n",
    ///     "+(+1LbRj-1LdRi+1LfRg+1LgRf-1LiRd+1LjRb)e3145\n",
    ///     "+(+1LcRj+1LdRh-1LeRg-1LgRe+1LhRd+1LjRc)e1245\n",
    ///     "+(+1LaRg+1LbRh+1LcRi+1LgRa+1LhRb+1LiRc)e1235\n",
    ///     "+(-1LaRd-1LbRe-1LcRf-1LdRa-1LeRb-1LfRc)e1234\n",
    ///     "+(+1LXRj-1LdRÐ+1LgRÞ+1LjRX-1LÐRd+1LÞRg)e0145\n",
    ///     "+(+1LYRj-1LeRÐ+1LhRÞ+1LjRY-1LÐRe+1LÞRh)e0245\n",
    ///     "+(+1LZRj-1LfRÐ+1LiRÞ+1LjRZ-1LÐRf+1LÞRi)e0345\n",
    ///     "+(-1LYRi+1LZRh-1LaRÐ+1LhRZ-1LiRY-1LÐRa)e0325\n",
    ///     "+(+1LXRi-1LZRg-1LbRÐ-1LgRZ+1LiRX-1LÐRb)e0135\n",
    ///     "+(-1LXRh+1LYRg-1LcRÐ+1LgRY-1LhRX-1LÐRc)e0215\n",
    ///     "+(+1LYRf-1LZRe+1LaRÞ-1LeRZ+1LfRY+1LÞRa)e0324\n",
    ///     "+(-1LXRf+1LZRd+1LbRÞ+1LdRZ-1LfRX+1LÞRb)e0134\n",
    ///     "+(+1LXRe-1LYRd+1LcRÞ-1LdRY+1LeRX+1LÞRc)e0214\n",
    ///     "+(+1LXRa+1LYRb+1LZRc+1LaRX+1LbRY+1LcRZ)e0123\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line()
    }
    /// The multivector of double motor $`m_2 \equiv s + v + \ell + S`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let double_motor = Vee::double_rotator().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// assert_eq!(format!("{double_motor:#}"), concat!(
    ///     "+1LvRv\n",
    ///     "+(-1LbRZ+1LcRY+1LdRÞ+1LgRÐ+1LvRX)e01\n",
    ///     "+(+1LaRZ-1LcRX+1LeRÞ+1LhRÐ+1LvRY)e02\n",
    ///     "+(-1LaRY+1LbRX+1LfRÞ+1LiRÐ+1LvRZ)e03\n",
    ///     "+(-1LdRX-1LeRY-1LfRZ-1LjRÐ+1LvRÞ)e40\n",
    ///     "+(-1LgRX-1LhRY-1LiRZ+1LjRÞ+1LvRÐ)e05\n",
    ///     "+1LaRve23\n",
    ///     "+1LbRve31\n",
    ///     "+1LcRve12\n",
    ///     "+1LdRve41\n",
    ///     "+1LeRve42\n",
    ///     "+1LfRve43\n",
    ///     "+1LgRve15\n",
    ///     "+1LhRve25\n",
    ///     "+1LiRve35\n",
    ///     "+1LjRve45\n",
    ///     "+1LxRve2345\n",
    ///     "+1LyRve3145\n",
    ///     "+1LzRve1245\n",
    ///     "+1LþRve1235\n",
    ///     "+1LðRve1234\n",
    ///     "+(-1LdRÐ+1LgRÞ+1LjRX-1LyRZ+1LzRY)e0145\n",
    ///     "+(-1LeRÐ+1LhRÞ+1LjRY+1LxRZ-1LzRX)e0245\n",
    ///     "+(-1LfRÐ+1LiRÞ+1LjRZ-1LxRY+1LyRX)e0345\n",
    ///     "+(-1LaRÐ+1LhRZ-1LiRY-1LxRÞ+1LþRX)e0325\n",
    ///     "+(-1LbRÐ-1LgRZ+1LiRX-1LyRÞ+1LþRY)e0135\n",
    ///     "+(-1LcRÐ+1LgRY-1LhRX-1LzRÞ+1LþRZ)e0215\n",
    ///     "+(+1LaRÞ-1LeRZ+1LfRY-1LxRÐ+1LðRX)e0324\n",
    ///     "+(+1LbRÞ+1LdRZ-1LfRX-1LyRÐ+1LðRY)e0134\n",
    ///     "+(+1LcRÞ-1LdRY+1LeRX-1LzRÐ+1LðRZ)e0214\n",
    ///     "+(+1LaRX+1LbRY+1LcRZ-1LðRÞ+1LþRÐ)e0123\n",
    ///     "+(+1LxRX+1LyRY+1LzRZ+1LðRÐ+1LþRÞ)I\n",
    /// ));
    ///
    /// let double_motor = Vee::plane().lhs() * Vee::plane().rhs();
    ///
    /// assert_eq!(double_motor.basis_blades(), Vee::double_motor().basis_blades());
    /// assert_eq!(format!("{double_motor:#}"), concat!(
    ///     "-1LaRa-1LbRb-1LcRc-1LdRd-1LeRe-1LfRf-1LgRg-1LhRh-1LiRi-1LjRj\n",
    ///     "+(-1LBRg+1LCRf+1LDRh+1LFRc-1LGRb+1LHRd+1LbRG-1LcRF-1LdRH-1LfRC+1LgRB-1LhRD)e01\n",
    ///     "+(-1LARg+1LCRe+1LDRi+1LERc-1LGRa+1LIRd+1LaRG-1LcRE-1LdRI-1LeRC+1LgRA-1LiRD)e02\n",
    ///     "+(-1LARf+1LBRe+1LDRj+1LERb-1LFRa+1LJRd+1LaRF-1LbRE-1LdRJ-1LeRB+1LfRA-1LjRD)e03\n",
    ///     "+(-1LARh+1LBRi-1LCRj-1LHRa+1LIRb-1LJRc+1LaRH-1LbRI+1LcRJ+1LhRA-1LiRB+1LjRC)e40\n",
    ///     "+(-1LERh+1LFRi-1LGRj-1LHRe+1LIRf-1LJRg+1LeRH-1LfRI+1LgRJ+1LhRE-1LiRF+1LjRG)e05\n",
    ///     "+(+1LbRc-1LcRb+1LfRg-1LgRf-1LiRj+1LjRi)e23\n",
    ///     "+(+1LaRc-1LcRa+1LeRg-1LgRe+1LhRj-1LjRh)e31\n",
    ///     "+(+1LaRb-1LbRa+1LeRf-1LfRe-1LhRi+1LiRh)e12\n",
    ///     "+(-1LaRd+1LdRa+1LfRj+1LgRi-1LiRg-1LjRf)e41\n",
    ///     "+(+1LbRd-1LdRb+1LeRj-1LgRh+1LhRg-1LjRe)e42\n",
    ///     "+(-1LcRd+1LdRc-1LeRi-1LfRh+1LhRf+1LiRe)e43\n",
    ///     "+(-1LbRj-1LcRi+1LdRe-1LeRd+1LiRc+1LjRb)e15\n",
    ///     "+(-1LaRj+1LcRh-1LdRf+1LfRd-1LhRc+1LjRa)e25\n",
    ///     "+(+1LaRi+1LbRh+1LdRg-1LgRd-1LhRb-1LiRa)e35\n",
    ///     "+(+1LaRe+1LbRf+1LcRg-1LeRa-1LfRb-1LgRc)e45\n",
    ///     "+(-1LbRg+1LcRf+1LdRh+1LfRc-1LgRb+1LhRd)e2345\n",
    ///     "+(-1LaRg+1LcRe+1LdRi+1LeRc-1LgRa+1LiRd)e3145\n",
    ///     "+(-1LaRf+1LbRe+1LdRj+1LeRb-1LfRa+1LjRd)e1245\n",
    ///     "+(-1LaRh+1LbRi-1LcRj-1LhRa+1LiRb-1LjRc)e1235\n",
    ///     "+(-1LeRh+1LfRi-1LgRj-1LhRe+1LiRf-1LjRg)e1234\n",
    ///     "+(-1LBRc+1LCRb-1LFRg+1LGRf+1LIRj-1LJRi+1LbRC-1LcRB+1LfRG-1LgRF-1LiRJ+1LjRI)e0145\n",
    ///     "+(-1LARc+1LCRa-1LERg+1LGRe-1LHRj+1LJRh+1LaRC-1LcRA+1LeRG-1LgRE+1LhRJ-1LjRH)e0245\n",
    ///     "+(-1LARb+1LBRa-1LERf+1LFRe+1LHRi-1LIRh+1LaRB-1LbRA+1LeRF-1LfRE-1LhRI+1LiRH)e0345\n",
    ///     "+(+1LARd-1LDRa-1LFRj-1LGRi+1LIRg+1LJRf-1LaRD+1LdRA+1LfRJ+1LgRI-1LiRG-1LjRF)e0325\n",
    ///     "+(-1LBRd+1LDRb-1LERj+1LGRh-1LHRg+1LJRe+1LbRD-1LdRB+1LeRJ-1LgRH+1LhRG-1LjRE)e0135\n",
    ///     "+(+1LCRd-1LDRc+1LERi+1LFRh-1LHRf-1LIRe-1LcRD+1LdRC-1LeRI-1LfRH+1LhRF+1LiRE)e0215\n",
    ///     "+(+1LBRj+1LCRi-1LDRe+1LERd-1LIRc-1LJRb-1LbRJ-1LcRI+1LdRE-1LeRD+1LiRC+1LjRB)e0324\n",
    ///     "+(+1LARj-1LCRh+1LDRf-1LFRd+1LHRc-1LJRa-1LaRJ+1LcRH-1LdRF+1LfRD-1LhRC+1LjRA)e0134\n",
    ///     "+(-1LARi-1LBRh-1LDRg+1LGRd+1LHRb+1LIRa+1LaRI+1LbRH+1LdRG-1LgRD-1LhRB-1LiRA)e0214\n",
    ///     "+(-1LARe-1LBRf-1LCRg+1LERa+1LFRb+1LGRc+1LaRE+1LbRF+1LcRG-1LeRA-1LfRB-1LgRC)e0123\n",
    ///     "+(-1LARa-1LBRb-1LCRc-1LDRd-1LERe-1LFRf-1LGRg-1LHRh-1LIRi-1LJRj",
    ///       "+1LaRA+1LbRB+1LcRC+1LdRD+1LeRE+1LfRF+1LgRG+1LhRH+1LiRI+1LjRJ)I\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_motor() -> Self {
        Self::scalar() + Self::volume() + Self::line() + Self::pseudoscalar()
    }
    /// The multivector of single rotoreflector $`f_{r1} \equiv h_0 + p_0`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let single_rotoreflector = Vee::normal().lhs() * Vee::single_rotator().rhs();
    ///
    /// assert_eq!(single_rotoreflector.basis_blades(), Vee::single_rotoreflector().basis_blades());
    /// assert_eq!(format!("{single_rotoreflector:#}"), concat!(
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(-1LyRf+1LzRe+1LþRa)e234\n",
    ///     "+(-1LxRf+1LzRd-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc)e123\n",
    ///     "+(-1LyRi+1LzRh-1LðRa)e253\n",
    ///     "+(-1LxRi+1LzRg+1LðRb)e315\n",
    ///     "+(-1LxRh+1LyRg-1LðRc)e152\n",
    ///     "+(+1LxRj-1LðRd-1LþRg)e145\n",
    ///     "+(+1LyRj-1LðRe-1LþRh)e245\n",
    ///     "+(+1LzRj-1LðRf-1LþRi)e345\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement()
    }
    /// The multivector of double rotoreflector $`f_{r2} \equiv h_0 + p_0 + P_0`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let double_rotoreflector = Vee::normal().lhs() * Vee::double_rotator().rhs();
    ///
    /// assert_eq!(double_rotoreflector.basis_blades(), Vee::double_rotoreflector().basis_blades());
    /// assert_eq!(format!("{double_rotoreflector:#}"), concat!(
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(+1LxRð-1LyRf+1LzRe-1LðRx+1LþRa)e234\n",
    ///     "+(-1LxRf-1LyRð+1LzRd+1LðRy-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LzRð-1LðRz+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc-1LðRþ-1LþRð)e123\n",
    ///     "+(-1LxRþ-1LyRi+1LzRh-1LðRa-1LþRx)e253\n",
    ///     "+(-1LxRi+1LyRþ+1LzRg+1LðRb+1LþRy)e315\n",
    ///     "+(-1LxRh+1LyRg-1LzRþ-1LðRc-1LþRz)e152\n",
    ///     "+(+1LxRj-1LyRz+1LzRy-1LðRd-1LþRg)e145\n",
    ///     "+(+1LxRz+1LyRj-1LzRx-1LðRe-1LþRh)e245\n",
    ///     "+(-1LxRy+1LyRx+1LzRj-1LðRf-1LþRi)e345\n",
    ///     "+(+1LxRx+1LyRy+1LzRz+1LðRð-1LþRþ)e12345\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_rotoreflector() -> Self {
        Self::normal() + Self::plane_displacement() + Self::weight()
    }
    /// The multivector of transflector $`f_t \equiv h + p_\infty`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let transflector = Vee::normal().lhs() * Vee::translator().rhs();
    ///
    /// assert_eq!(transflector.basis_blades(), Vee::transflector().basis_blades());
    /// assert_eq!(format!("{transflector:#}"), concat!(
    ///     "+(-1LxRX-1LyRY-1LzRZ-1LðRÐ+1LþRÞ)e0\n",
    ///     "+1LxRve1\n",
    ///     "+1LyRve2\n",
    ///     "+1LzRve3\n",
    ///     "+1LþRve4\n",
    ///     "+1LðRve5\n",
    ///     "+(-1LxRÐ+1LðRX)e015\n",
    ///     "+(+1LyRÐ-1LðRY)e052\n",
    ///     "+(-1LzRÐ+1LðRZ)e035\n",
    ///     "+(+1LðRÞ+1LþRÐ)e054\n",
    ///     "+(+1LxRÞ+1LþRX)e014\n",
    ///     "+(-1LyRÞ-1LþRY)e042\n",
    ///     "+(+1LzRÞ+1LþRZ)e034\n",
    ///     "+(+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LxRZ+1LzRX)e013\n",
    ///     "+(+1LxRY-1LyRX)e021\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn transflector() -> Self {
        Self::hypervolume() + Self::plane_moment()
    }
    /// The multivector of simple single flector $`f_{s1} \equiv h + p`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let simple_single_flector = Vee::hypervolume().lhs() * Vee::simple_single_motor().rhs();
    ///
    /// assert_eq!(simple_single_flector.basis_blades(),
    ///     Vee::simple_single_flector().basis_blades());
    /// assert_eq!(format!("{simple_single_flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ-1LðRÐ+1LþRÞ)e0\n",
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(+1LWRg-1LxRÐ+1LðRX)e015\n",
    ///     "+(-1LWRh+1LyRÐ-1LðRY)e052\n",
    ///     "+(+1LWRi-1LzRÐ+1LðRZ)e035\n",
    ///     "+(-1LWRj+1LðRÞ+1LþRÐ)e054\n",
    ///     "+(-1LWRd+1LxRÞ+1LþRX)e014\n",
    ///     "+(+1LWRe-1LyRÞ-1LþRY)e042\n",
    ///     "+(-1LWRf+1LzRÞ+1LþRZ)e034\n",
    ///     "+(-1LWRa+1LyRZ-1LzRY)e032\n",
    ///     "+(-1LWRb-1LxRZ+1LzRX)e013\n",
    ///     "+(-1LWRc+1LxRY-1LyRX)e021\n",
    ///     "+(-1LyRf+1LzRe+1LþRa)e234\n",
    ///     "+(-1LxRf+1LzRd-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc)e123\n",
    ///     "+(-1LyRi+1LzRh-1LðRa)e253\n",
    ///     "+(-1LxRi+1LzRg+1LðRb)e315\n",
    ///     "+(-1LxRh+1LyRg-1LðRc)e152\n",
    ///     "+(+1LxRj-1LðRd-1LþRg)e145\n",
    ///     "+(+1LyRj-1LðRe-1LþRh)e245\n",
    ///     "+(+1LzRj-1LðRf-1LþRi)e345\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn simple_single_flector() -> Self {
        Self::hypervolume() + Self::plane()
    }
    /// The multivector of single flector $`f_1 \equiv h + p + P_\infty`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let single_flector = Vee::hypervolume().lhs() * Vee::single_motor().rhs();
    ///
    /// assert_eq!(single_flector.basis_blades(), Vee::single_flector().basis_blades());
    /// assert_eq!(format!("{single_flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ-1LðRÐ+1LþRÞ)e0\n",
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(+1LWRg-1LxRÐ-1LyRF+1LzRE+1LðRX+1LþRA)e015\n",
    ///     "+(-1LWRh-1LxRF+1LyRÐ+1LzRD-1LðRY-1LþRB)e052\n",
    ///     "+(+1LWRi-1LxRE+1LyRD-1LzRÐ+1LðRZ+1LþRC)e035\n",
    ///     "+(-1LWRj+1LxRA+1LyRB+1LzRC+1LðRÞ+1LþRÐ)e054\n",
    ///     "+(-1LWRd+1LxRÞ-1LyRI+1LzRH-1LðRA+1LþRX)e014\n",
    ///     "+(+1LWRe-1LxRI-1LyRÞ+1LzRG+1LðRB-1LþRY)e042\n",
    ///     "+(-1LWRf-1LxRH+1LyRG+1LzRÞ-1LðRC+1LþRZ)e034\n",
    ///     "+(-1LWRa+1LxRJ+1LyRZ-1LzRY-1LðRD-1LþRG)e032\n",
    ///     "+(-1LWRb-1LxRZ+1LyRJ+1LzRX-1LðRE-1LþRH)e013\n",
    ///     "+(-1LWRc+1LxRY-1LyRX+1LzRJ-1LðRF-1LþRI)e021\n",
    ///     "+(-1LyRf+1LzRe+1LþRa)e234\n",
    ///     "+(-1LxRf+1LzRd-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc)e123\n",
    ///     "+(-1LyRi+1LzRh-1LðRa)e253\n",
    ///     "+(-1LxRi+1LzRg+1LðRb)e315\n",
    ///     "+(-1LxRh+1LyRg-1LðRc)e152\n",
    ///     "+(+1LxRj-1LðRd-1LþRg)e145\n",
    ///     "+(+1LyRj-1LðRe-1LþRh)e245\n",
    ///     "+(+1LzRj-1LðRf-1LþRi)e345\n",
    ///     "+(+1LyRC-1LzRB+1LðRG-1LþRD)e03245\n",
    ///     "+(-1LxRC+1LzRA+1LðRH-1LþRE)e01345\n",
    ///     "+(+1LxRB-1LyRA+1LðRI-1LþRF)e02145\n",
    ///     "+(+1LxRD+1LyRE+1LzRF+1LðRJ)e01235\n",
    ///     "+(-1LxRG-1LyRH-1LzRI-1LþRJ)e01243\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn single_flector() -> Self {
        Self::hypervolume() + Self::plane() + Self::direction()
    }
    /// The multivector of double flector $`f_2 \equiv h + p + P`$.
    ///
    /// ```
    /// use vee::PgaP5 as Vee;
    ///
    /// let double_flector = Vee::hypervolume().lhs() * Vee::simple_double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// assert_eq!(format!("{double_flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ-1LðRÐ+1LþRÞ)e0\n",
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(+1LWRg-1LxRÐ-1LyRF+1LzRE+1LðRX+1LþRA)e015\n",
    ///     "+(-1LWRh-1LxRF+1LyRÐ+1LzRD-1LðRY-1LþRB)e052\n",
    ///     "+(+1LWRi-1LxRE+1LyRD-1LzRÐ+1LðRZ+1LþRC)e035\n",
    ///     "+(-1LWRj+1LxRA+1LyRB+1LzRC+1LðRÞ+1LþRÐ)e054\n",
    ///     "+(-1LWRd+1LxRÞ-1LyRI+1LzRH-1LðRA+1LþRX)e014\n",
    ///     "+(+1LWRe-1LxRI-1LyRÞ+1LzRG+1LðRB-1LþRY)e042\n",
    ///     "+(-1LWRf-1LxRH+1LyRG+1LzRÞ-1LðRC+1LþRZ)e034\n",
    ///     "+(-1LWRa+1LxRJ+1LyRZ-1LzRY-1LðRD-1LþRG)e032\n",
    ///     "+(-1LWRb-1LxRZ+1LyRJ+1LzRX-1LðRE-1LþRH)e013\n",
    ///     "+(-1LWRc+1LxRY-1LyRX+1LzRJ-1LðRF-1LþRI)e021\n",
    ///     "+(+1LxRð-1LyRf+1LzRe-1LðRx+1LþRa)e234\n",
    ///     "+(-1LxRf-1LyRð+1LzRd+1LðRy-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LzRð-1LðRz+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc-1LðRþ-1LþRð)e123\n",
    ///     "+(-1LxRþ-1LyRi+1LzRh-1LðRa-1LþRx)e253\n",
    ///     "+(-1LxRi+1LyRþ+1LzRg+1LðRb+1LþRy)e315\n",
    ///     "+(-1LxRh+1LyRg-1LzRþ-1LðRc-1LþRz)e152\n",
    ///     "+(+1LxRj-1LyRz+1LzRy-1LðRd-1LþRg)e145\n",
    ///     "+(+1LxRz+1LyRj-1LzRx-1LðRe-1LþRh)e245\n",
    ///     "+(-1LxRy+1LyRx+1LzRj-1LðRf-1LþRi)e345\n",
    ///     "+(+1LxRx+1LyRy+1LzRz+1LðRð-1LþRþ)e12345\n",
    ///     "+(-1LWRx+1LyRC-1LzRB+1LðRG-1LþRD)e03245\n",
    ///     "+(-1LWRy-1LxRC+1LzRA+1LðRH-1LþRE)e01345\n",
    ///     "+(-1LWRz+1LxRB-1LyRA+1LðRI-1LþRF)e02145\n",
    ///     "+(+1LWRþ+1LxRD+1LyRE+1LzRF+1LðRJ)e01235\n",
    ///     "+(-1LWRð-1LxRG-1LyRH-1LzRI-1LþRJ)e01243\n",
    /// ));
    ///
    /// let double_flector = Vee::hypervolume().lhs() * Vee::double_motor().rhs();
    ///
    /// assert_eq!(double_flector.basis_blades(), Vee::double_flector().basis_blades());
    /// assert_eq!(format!("{double_flector:#}"), concat!(
    ///     "+(+1LWRv-1LxRX-1LyRY-1LzRZ-1LðRÐ+1LþRÞ)e0\n",
    ///     "+(+1LxRv-1LyRc+1LzRb-1LðRg+1LþRd)e1\n",
    ///     "+(+1LxRc+1LyRv-1LzRa-1LðRh+1LþRe)e2\n",
    ///     "+(-1LxRb+1LyRa+1LzRv-1LðRi+1LþRf)e3\n",
    ///     "+(-1LxRd-1LyRe-1LzRf-1LðRj+1LþRv)e4\n",
    ///     "+(+1LxRg+1LyRh+1LzRi+1LðRv+1LþRj)e5\n",
    ///     "+(+1LWRg-1LxRÐ-1LyRF+1LzRE+1LðRX+1LþRA)e015\n",
    ///     "+(-1LWRh-1LxRF+1LyRÐ+1LzRD-1LðRY-1LþRB)e052\n",
    ///     "+(+1LWRi-1LxRE+1LyRD-1LzRÐ+1LðRZ+1LþRC)e035\n",
    ///     "+(-1LWRj+1LxRA+1LyRB+1LzRC+1LðRÞ+1LþRÐ)e054\n",
    ///     "+(-1LWRd+1LxRÞ-1LyRI+1LzRH-1LðRA+1LþRX)e014\n",
    ///     "+(+1LWRe-1LxRI-1LyRÞ+1LzRG+1LðRB-1LþRY)e042\n",
    ///     "+(-1LWRf-1LxRH+1LyRG+1LzRÞ-1LðRC+1LþRZ)e034\n",
    ///     "+(-1LWRa+1LxRJ+1LyRZ-1LzRY-1LðRD-1LþRG)e032\n",
    ///     "+(-1LWRb-1LxRZ+1LyRJ+1LzRX-1LðRE-1LþRH)e013\n",
    ///     "+(-1LWRc+1LxRY-1LyRX+1LzRJ-1LðRF-1LþRI)e021\n",
    ///     "+(+1LxRð-1LyRf+1LzRe-1LðRx+1LþRa)e234\n",
    ///     "+(-1LxRf-1LyRð+1LzRd+1LðRy-1LþRb)e134\n",
    ///     "+(-1LxRe+1LyRd+1LzRð-1LðRz+1LþRc)e124\n",
    ///     "+(+1LxRa+1LyRb+1LzRc-1LðRþ-1LþRð)e123\n",
    ///     "+(-1LxRþ-1LyRi+1LzRh-1LðRa-1LþRx)e253\n",
    ///     "+(-1LxRi+1LyRþ+1LzRg+1LðRb+1LþRy)e315\n",
    ///     "+(-1LxRh+1LyRg-1LzRþ-1LðRc-1LþRz)e152\n",
    ///     "+(+1LxRj-1LyRz+1LzRy-1LðRd-1LþRg)e145\n",
    ///     "+(+1LxRz+1LyRj-1LzRx-1LðRe-1LþRh)e245\n",
    ///     "+(-1LxRy+1LyRx+1LzRj-1LðRf-1LþRi)e345\n",
    ///     "+(+1LxRx+1LyRy+1LzRz+1LðRð-1LþRþ)e12345\n",
    ///     "+(-1LWRx+1LxRV+1LyRC-1LzRB+1LðRG-1LþRD)e03245\n",
    ///     "+(-1LWRy-1LxRC+1LyRV+1LzRA+1LðRH-1LþRE)e01345\n",
    ///     "+(-1LWRz+1LxRB-1LyRA+1LzRV+1LðRI-1LþRF)e02145\n",
    ///     "+(+1LWRþ+1LxRD+1LyRE+1LzRF+1LðRJ+1LþRV)e01235\n",
    ///     "+(-1LWRð-1LxRG-1LyRH-1LzRI+1LðRV-1LþRJ)e01243\n",
    /// ));
    /// ```
    #[must_use]
    #[inline]
    pub fn double_flector() -> Self {
        Self::hypervolume() + Self::plane() + Self::point()
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
