// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

use core::{
    cmp::Ordering,
    fmt::{self, Alignment, Display},
    ops::Not,
};

/// Symbol as Unicode character with optional *combining diacritical mark*.
#[derive(Debug, Clone, Copy)]
pub struct Symbol {
    var: char,
    alt: char,
    pub(crate) cdm: char,
    pub(crate) lab: &'static str,
}

impl PartialEq for Symbol {
    fn eq(&self, other: &Self) -> bool {
        if self.is_vec() && other.is_vec() {
            self.lab == other.lab
        } else {
            self.var == other.var && self.alt == other.alt && self.cdm == other.cdm
        }
    }
}

impl Eq for Symbol {}

impl Ord for Symbol {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.is_vec() && other.is_vec() {
            self.lab.cmp(other.lab)
        } else {
            self.var
                .cmp(&other.var)
                .then(self.alt.cmp(&other.alt))
                .then(self.cdm.cmp(&other.cdm))
        }
    }
}

impl PartialOrd for Symbol {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Symbol {
    /// Unicode null (i.e., `'\0'`).
    pub const NIL: char = '\0';
    /// Unicode zero-width space.
    pub const VEC: char = '\u{200b}';

    /// Unicode *combining dot above* (i.e., `"◌̇"`).
    pub const ALT: char = '\u{0307}';
    /// Unicode *combining x below* (i.e., `"◌͓"`).
    pub const PIN: char = '\u{0353}';
    /// Unicode *combining left arrowhead below* (i.e., `"◌͔"`).
    pub const LHS: char = '\u{0354}';
    /// Unicode *combining right arrowhead below* (i.e., `"◌͕"`).
    pub const RHS: char = '\u{0355}';

    /// Whether this symbol is a basis blade.
    #[must_use]
    #[inline]
    pub const fn is_vec(&self) -> bool {
        self.var == Self::VEC
    }
    /// Whether this symbol is the scalar.
    #[must_use]
    #[inline]
    pub const fn is_scalar(&self) -> bool {
        self.is_vec() && self.lab.len() == 1
    }
    /// Whether this symbol is the pseudoscalar.
    #[must_use]
    #[inline]
    pub const fn is_pseudoscalar(&self) -> bool {
        self.is_vec() && self.alt == Self::ALT
    }
    /// Whether this symbol is pinned.
    #[must_use]
    #[inline]
    pub const fn is_pin(&self) -> bool {
        self.cdm == Self::PIN
    }
    /// Whether this symbol is alternative.
    #[must_use]
    #[inline]
    pub const fn is_alt(&self) -> bool {
        self.alt == Self::ALT
    }
    /// Creates symbol for `variable` and non-empty `label` in $`\e`$-notation.
    ///
    /// # Panics
    ///
    /// Panics on empty symbol `label`.
    #[must_use]
    #[inline]
    pub const fn new(variable: char, label: &'static str) -> Self {
        assert!(!label.is_empty(), "empty symbol label");
        Self {
            var: variable,
            alt: Self::NIL,
            cdm: Self::NIL,
            lab: label,
        }
    }
    /// Marks this symbol with [`Self::ALT`].
    #[must_use]
    #[inline]
    pub const fn alt(mut self) -> Self {
        self.alt = Self::ALT;
        self
    }
    /// Marks this symbol with Unicode *combining diacritical mark*.
    #[must_use]
    #[inline]
    pub const fn cdm(mut self, cdm: char) -> Self {
        self.cdm = cdm;
        self
    }
    /// Marks this symbol with [`Self::PIN`].
    #[must_use]
    #[inline]
    pub const fn pin(self) -> Self {
        self.cdm(Self::PIN)
    }
    /// Marks this symbol with [`Self::LHS`] as left-hand side.
    #[must_use]
    #[inline]
    pub const fn lhs(self) -> Self {
        self.cdm(Self::LHS)
    }
    /// Marks this symbol with [`Self::RHS`] as right-hand side.
    #[must_use]
    #[inline]
    pub const fn rhs(self) -> Self {
        self.cdm(Self::RHS)
    }
}

impl From<(&str, &'static str)> for Symbol {
    #[inline]
    fn from((vars, lab): (&str, &'static str)) -> Self {
        let mut vars = vars.chars();
        let var = vars.next().unwrap_or_default();
        assert_eq!(vars.next(), None, "multi-character symbol");
        Self::new(var, lab)
    }
}

impl From<(char, &'static str)> for Symbol {
    #[inline]
    fn from((var, lab): (char, &'static str)) -> Self {
        Self::new(var, lab)
    }
}

impl Not for Symbol {
    type Output = Self;

    /// Swaps lowercase and uppercase character.
    fn not(self) -> Self {
        let var = if self.var.is_lowercase() {
            let mut iter = self.var.to_uppercase();
            assert_eq!(iter.len(), 1, "no uppercase for {}", self.var);
            iter.next().unwrap()
        } else {
            let mut iter = self.var.to_lowercase();
            assert_eq!(iter.len(), 1, "no lowercase for {}", self.var);
            iter.next().unwrap()
        };
        Self {
            var,
            alt: self.alt,
            cdm: self.cdm,
            lab: self.lab,
        }
    }
}

impl Not for &Symbol {
    type Output = Symbol;

    /// Swaps lowercase and uppercase character.
    #[inline]
    fn not(self) -> Self::Output {
        !*self
    }
}

impl Display for Symbol {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        if self.is_vec() {
            if fmt.fill() == '$' {
                if self.is_pseudoscalar() {
                    write!(fmt, "\\boldsymbol{{I}}")?;
                } else {
                    let lab = &self.lab[1..];
                    if lab.is_empty() {
                        write!(fmt, "\\boldsymbol{{e}}")?;
                    } else if lab.len() == 1 {
                        write!(fmt, "\\boldsymbol{{e}}_{lab}")?;
                    } else {
                        write!(fmt, "\\boldsymbol{{e}}_{{{lab}}}")?;
                    }
                }
            } else if self.is_scalar() {
                write!(
                    fmt,
                    "{}",
                    if fmt.precision().is_some() {
                        "1.0"
                    } else {
                        "1"
                    }
                )?;
            } else if fmt.align() == Some(Alignment::Center) {
                write!(fmt, "o.{}", self.lab)?;
            } else if self.is_pseudoscalar() {
                write!(fmt, "I")?;
            } else {
                write!(fmt, "{}", self.lab)?;
            }
        } else {
            if fmt.alternate() || fmt.align() == Some(Alignment::Center) || fmt.fill() == '$' {
                let var = match self.cdm {
                    Self::PIN => 'p',
                    Self::LHS => 'l',
                    Self::RHS => 'r',
                    _ => 'v',
                };
                let lab = &self.lab[1..];
                if fmt.fill() == '$' {
                    if lab.is_empty() {
                        write!(fmt, "{var}")?;
                    } else if lab.len() == 1 {
                        write!(fmt, "{var}_{lab}")?;
                    } else {
                        write!(fmt, "{var}_{{{lab}}}")?;
                    }
                } else if fmt.align() == Some(Alignment::Center) {
                    write!(fmt, "{var}.e{lab}")?;
                } else {
                    write!(fmt, "{var}{lab}")?;
                }
            } else {
                write!(fmt, "{}", self.var)?;
                if self.alt != Self::NIL {
                    write!(fmt, "{}", self.alt)?;
                }
                if self.cdm != Self::NIL {
                    write!(fmt, "{}", self.cdm)?;
                }
            }
        }
        Ok(())
    }
}
