# Version 0.1.6 (2025-11-05)

  * Fix KaTeX, maybe?
  * Fix text form for `"{:-}"`.

# Version 0.1.5 (2025-11-04)

  * Add DOT form (i.e., `text/vnd.graphviz`) via `Octal` trait.
  * Add alternative symbols labelled after basis blades (`"{:#}"`).
  * Implement factorization of pinned symbols, GCD coefficients, and predominant sign.
      * Expand (`"{:+}"`) or reduce (i.e., factorized) (`"{}"`).
      * Factor predominant sign (`"{:-}"`).
  * Eliminate orthonormalization conditions from expressions using reflection/projection operator.
  * Evaluate symbols as rationals.
  * Count operations (i.e., multiplications and additions).
  * Fix KaTeX rendering on large pages.

# Version 0.1.4 (2025-10-16)

  * Update KaTeX.
  * Explore 6D/7D PGA.
  * Assert entities have unique symbols and exactly one per basis blade fixing 4D PGA.
  * Introduce `Symbol::alt()` to extend symbol space.
  * Rechoose basis blades following recipe. After a refactoring, they can be genrated this way.
  * Refactor rational numbers without dependency and find their GCD/LCM.
  * Refactor symbol without dependency.
  * Swap case without dependency.
  * Always export `format_eq!` by using non-pretty fallback.
  * Let `format_eq!` omit unit coefficient for more compact text form.
  * Introduce `{:+}` for not omitting unit coefficient.
  * Fix invalid coefficient omission.
  * Fix typo and link.
  * Construct Cayley table with pre-allocations only.

# Version 0.1.3 (2025-09-18)

  * Add projection operator. Supports rejection, i.e., projecting lower- onto higher-grade blade.
  * Use Unicode *combining diacritical marks* instead of ~/L/R.
  * Define 6D/7D PGA blades.

# Version 0.1.2 (2025-09-07)

  * Generate constructor methods of basis blades with compile-time blocks.
  * Distinguish bias/weight (W/w) from scalar/pseudoscalar (v/V) and relabel norm.
  * Distinguish motors.
  * Add more examples.
  * Add polarity operator.
  * Fix bugs, typos, and improve documentation.
  * Explore 4D/5D PGA.

# Version 0.1.1 (2025-09-02)

  * Fix KaTeX documentation.

# Version 0.1.0 (2025-09-02)

  * Implement symbolic reduction of multivector expressions.
  * Implement generation of expressions in text form.
  * Define the metric-agnostic basis (i.e., elliptic, parabolic, hyperbolic) and the multivector
    entities of the projective geometric algebra (PGA) for dimensions N < 6. Assert the Cayley
    tables in unit tests.

# Version 0.0.0 (2021-06-27)

  * Reserve name.
