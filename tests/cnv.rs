// Copyright © 2025-2026 Rouven Spreckels <rs@qu1x.dev>
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of
// the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(missing_docs)]

use vee::format_eq;

#[test]
#[allow(clippy::too_many_lines, clippy::needless_raw_string_hashes)]
fn cnv() {
    use vee::Tree;
    {
        use vee::PgaP4 as Vee;

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
        use vee::{PgaP3 as Vee, pga::PgaP3 as Bee};

        let norm = Vee::norm().eval([(Bee::e0123(), 2)]);
        assert_eq!(
            Tree::from(norm.clone()),
            Tree::with_factorization(norm, false)
        );

        format_eq!("{}", Vee::norm().eval([(Bee::e(), 1)]), ["+1", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(Bee::e(), -1)]), ["-1", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(Bee::e(), 2)]), ["+2", "+VI"]);
        format_eq!("{}", Vee::norm().eval([(Bee::e0123(), 1)]), ["+v", "+I"]);
        format_eq!("{}", Vee::norm().eval([(Bee::e0123(), -1)]), ["+v", "-I"]);
        format_eq!("{}", Vee::norm().eval([(Bee::e0123(), 2)]), ["+v", "+2I"]);

        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e(), 1)]), ["+1+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e(), -1)]), ["-1+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e(), 2)]), ["+2+VI"]);
        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e0123(), 1)]), ["+v+I"]);
        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e0123(), -1)]), ["+v-I"]);
        format_eq!("{:0}\n", Vee::norm().eval([(Bee::e0123(), 2)]), ["+v+2I"]);

        format_eq!("{:<}", Vee::norm().eval([(Bee::e(), 1)]), ["1", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(Bee::e(), -1)]), ["-1", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(Bee::e(), 2)]), ["2", "+VI"]);
        format_eq!("{:<}", Vee::norm().eval([(Bee::e0123(), 1)]), ["v", "+I"]);
        format_eq!("{:<}", Vee::norm().eval([(Bee::e0123(), -1)]), ["v", "-I"]);
        format_eq!("{:<}", Vee::norm().eval([(Bee::e0123(), 2)]), ["v", "+2I"]);

        format_eq!(
            "{:#}",
            Vee::norm().eval([(Bee::e(), 1)]),
            ["+1", "+v0123*I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(Bee::e(), -1)]),
            ["-1", "+v0123*I"]
        );
        format_eq!(
            "{:#}",
            Vee::norm().eval([(Bee::e(), 2)]),
            ["+2", "+v0123*I"]
        );
        format_eq!("{:#}", Vee::norm().eval([(Bee::e0123(), 1)]), ["+v", "+I"]);
        format_eq!("{:#}", Vee::norm().eval([(Bee::e0123(), -1)]), ["+v", "-I"]);
        format_eq!(
            "{:#}",
            Vee::norm().eval([(Bee::e0123(), 2)]),
            ["+v", "+2*I"]
        );

        format_eq!(
            "    {:$^4}",
            Vee::norm().eval([(Bee::e(), 1)]),
            [
                r"    \begin{aligned}[t]",
                r"      1 & \\",
                r"      + v_{0123} & \boldsymbol{I}",
                r"    \end{aligned}",
            ]
        );
        format_eq!(
            "    {:$<4}",
            Vee::norm().eval([(Bee::e(), 1)]),
            [r"      1 & \\", r"      + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<0}\n",
            Vee::norm().eval([(Bee::e(), 1)]),
            [r"1 + v_{0123} \boldsymbol{I}"]
        );

        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e(), 1)]),
            [r"  1 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e(), -1)]),
            [r"  -1 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e(), 2)]),
            [r"  2 & \\", r"  + v_{0123} & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e0123(), 1)]),
            [r"  v & \\", r"  + & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e0123(), -1)]),
            [r"  v & \\", r"  - & \boldsymbol{I}"]
        );
        format_eq!(
            "{:$<}",
            Vee::norm().eval([(Bee::e0123(), 2)]),
            [r"  v & \\", r"  + 2 & \boldsymbol{I}"]
        );

        format_eq!(
            "{:>#4x}",
            Vee::norm().eval([(Bee::e(), 1)]),
            ["    let e = 1.0;", "    let e0123 = v0123;"]
        );

        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), 1)]),
            ["let e = 1.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), -1)]),
            ["let e = -1.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), 2)]),
            ["let e = 2.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e0123(), 1)]),
            ["let e = v;", "let e0123 = 1.0;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e0123(), -1)]),
            ["let e = v;", "let e0123 = -1.0;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e0123(), 2)]),
            ["let e = v;", "let e0123 = 2.0;"]
        );

        format_eq!(
            "{:#x}",
            Vee::plane()
                .norm_squared()
                .eval([(Bee::e1(), 2), (Bee::e2(), 3)]),
            ["let e = 13.0 + v3 * v3;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), 0)]),
            ["let e = 0.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e0123(), 0)]),
            ["let e = v;", "let e0123 = 0.0;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), (1, 2))]),
            ["let e = 1.0 / 2.0;", "let e0123 = v0123;"]
        );
        format_eq!(
            "{:#x}",
            Vee::norm().eval([(Bee::e(), (-1, 2))]),
            ["let e = -1.0 / 2.0;", "let e0123 = v0123;"]
        );

        format_eq!(
            "{:^4x}",
            Vee::norm().eval([(Bee::e(), 1)]),
            ["    o.e = 1.0;", "    o.e0123 = v.e0123;"]
        );

        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e(), 1)]),
            ["o.e = 1.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e(), -1)]),
            ["o.e = -1.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e(), 2)]),
            ["o.e = 2.0;", "o.e0123 = v.e0123;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e0123(), 1)]),
            ["o.e = v.e;", "o.e0123 = 1.0;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e0123(), -1)]),
            ["o.e = v.e;", "o.e0123 = -1.0;"]
        );
        format_eq!(
            "{:^x}",
            Vee::norm().eval([(Bee::e0123(), 2)]),
            ["o.e = v.e;", "o.e0123 = 2.0;"]
        );

        format_eq!(
            "{:#o}",
            Vee::norm().eval([(Bee::e(), 1)]),
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
            Vee::norm().eval([(Bee::e(), -1)]),
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
            Vee::norm().eval([(Bee::e(), 2)]),
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
            Vee::norm().eval([(Bee::e0123(), 1)]),
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
            Vee::norm().eval([(Bee::e0123(), -1)]),
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
            Vee::norm().eval([(Bee::e0123(), 2)]),
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
