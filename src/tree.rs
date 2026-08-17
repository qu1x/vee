// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use super::{Algebra, Factorization, Monomial, Multivector, Polynomial, Rational, Symbol};
use core::{iter::repeat_n, ops::Not};

/// Non-binary algebraic expression tree up to symbolic [`Multivector`] expressions.
///
/// Unifies and simplifies following implementations due to being an intermediate and recursive data
/// structure:
///
///   * <code>impl [Display] for [Multivector]</code>
///   * <code>impl [Octal] for [Multivector]</code>
///
/// Convert from:
///
///   * <code>Tree::from([Multivector])</code> or <code>[Tree::with_factorization](v, true)</code>
///   * <code>Tree::from([Factorization])</code>
///   * <code>Tree::from([Polynomial])</code>
///   * <code>Tree::from([Monomial])</code>
///   * <code>Tree::from([Rational])</code>
///   * <code>Tree::from([Symbol])</code>
///
/// Convert to:
///
///   * <code>[Multivector]::try_from(tree)</code> where <code>Error = [Symbol]</code> indicates
///     that [`Symbol`] is not part of the [`Algebra`]. This is due [`Tree`] and [`Symbol`] being
///     non-generic and hence agnostic of the [`Algebra`].
///
/// [Display]: `core::fmt::Display`
/// [Octal]: `core::fmt::Octal`
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tree {
    /// Sum of subtrees.
    Add(Vec<Self>),
    /// Product of subtrees.
    Mul(Vec<Self>),
    /// Rational leaf node.
    Num(Rational),
    /// Symbolic leaf node.
    Sym(Symbol),
}

impl Tree {
    /// The zero constant.
    pub const ZERO: Self = Self::Num(Rational::ZERO);
    /// The one constant.
    pub const ONE: Self = Self::Num(Rational::ONE);

    /// Performs factorization on `v` and creates tree.
    ///
    /// Optionally, the GCDs are `signed` comprising the factored predominant sign.
    #[must_use]
    #[allow(clippy::missing_panics_doc)]
    pub fn with_factorization<B: Algebra>(v: Multivector<B>, signed: bool) -> Self {
        let mut add = Vec::with_capacity(v.map.len());
        for (b, p) in v.map {
            let f = Factorization::new(p, signed);
            if f.is_one() {
                add.push(b.into().into());
            } else {
                add.push(Self::Mul(vec![f.into(), b.into().into()]));
            }
        }
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.pop().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
    /// Returns [`Rational`] if [`Self::Num`].
    #[must_use]
    #[inline]
    pub const fn as_num(&self) -> Option<Rational> {
        if let Self::Num(num) = self {
            Some(*num)
        } else {
            None
        }
    }
    /// Returns [`Symbol`] if [`Self::Sym`].
    #[must_use]
    #[inline]
    pub const fn as_sym(&self) -> Option<Symbol> {
        if let Self::Sym(sym) = self {
            Some(*sym)
        } else {
            None
        }
    }
}

impl Default for Tree {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl<B: Algebra> From<Multivector<B>> for Tree {
    fn from(v: Multivector<B>) -> Self {
        let add = v
            .map
            .into_iter()
            .filter(|(_b, p)| !p.is_zero())
            .map(|(b, p)| {
                if p.is_one() {
                    b.into().into()
                } else {
                    Self::Mul(vec![p.into(), b.into().into()])
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.into_iter().next().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
}

impl From<Factorization> for Tree {
    fn from(f: Factorization) -> Self {
        let add = f
            .map
            .into_iter()
            .filter(|(_m, (p, c))| !(p.is_zero() || c.is_zero()))
            .map(|(m, (p, c))| {
                let is_mul = [c.is_one(), p.is_one(), m.is_one()].map(bool::not);
                let len = is_mul.into_iter().map(usize::from).sum::<usize>();
                let mut mul = Vec::with_capacity(len);
                if is_mul[0] {
                    mul.push(c.into());
                }
                if is_mul[1] {
                    mul.push(p.into());
                }
                if is_mul[2] {
                    mul.push(m.into());
                }
                if len == 0 {
                    Self::ONE
                } else if len == 1 {
                    mul.pop().expect("unreachable")
                } else {
                    Self::Mul(mul)
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            let any = add.into_iter().next().expect("unreachable");
            if f.gcd.is_one() {
                any
            } else if let Self::Num(num) = any {
                Self::Num(f.gcd * num)
            } else {
                Self::Mul(vec![f.gcd.into(), any])
            }
        } else if f.gcd.is_one() {
            Self::Add(add)
        } else {
            Self::Mul(vec![f.gcd.into(), Self::Add(add)])
        }
    }
}

impl From<Polynomial> for Tree {
    fn from(p: Polynomial) -> Self {
        let add = p
            .map
            .into_iter()
            .filter(|(_m, c)| !c.is_zero())
            .map(|(m, c)| {
                if c.is_one() {
                    m.into()
                } else {
                    let m = Self::from(m);
                    if let Self::Mul(mul) = m {
                        let mut vec = vec![Self::from(c)];
                        vec.extend(mul);
                        Self::Mul(vec)
                    } else if let Self::Num(num) = m {
                        Self::Num(c * num)
                    } else {
                        Self::Mul(vec![c.into(), m])
                    }
                }
            })
            .collect::<Vec<Self>>();
        if add.is_empty() {
            Self::ZERO
        } else if add.len() == 1 {
            add.into_iter().next().expect("unreachable")
        } else {
            Self::Add(add)
        }
    }
}

impl From<Monomial> for Tree {
    fn from(m: Monomial) -> Self {
        let mul = m
            .map
            .into_iter()
            .flat_map(|(s, e)| {
                repeat_n(s, e.get().try_into().expect("negative exponent")).map(From::from)
            })
            .collect::<Vec<Self>>();
        if mul.is_empty() {
            Self::ONE
        } else if mul.len() == 1 {
            mul.into_iter().next().expect("unreachable")
        } else {
            Self::Mul(mul)
        }
    }
}

impl From<Rational> for Tree {
    #[inline]
    fn from(r: Rational) -> Self {
        Self::Num(r)
    }
}

impl From<Symbol> for Tree {
    #[inline]
    fn from(s: Symbol) -> Self {
        Self::Sym(s)
    }
}
