use std::{iter::Peekable, str::Chars};

#[derive(Debug)]
enum Token {
    LBrack,
    RBrack,
    Colon,
    Atom(usize),
    HCount(usize),
    Digit(usize),
    SingleBond,
    TripleBond,
    End,
}

fn get_digits(chars: &mut Peekable<Chars<'_>>) -> String {
    let mut digits = String::new();
    while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
        digits.push(chars.next().unwrap());
    }
    digits
}

fn scan(s: String) -> Vec<Token> {
    use Token as T;
    let mut chars = s.chars().peekable();
    let mut ret = Vec::new();
    while let Some(c) = chars.next() {
        let got = match c {
            '[' => T::LBrack,
            ']' => T::RBrack,
            '#' => {
                // # can either be a number inside of an atom, eg [#6], or a
                // triple bond. at some point we might have to improve this
                // check
                let digits = get_digits(&mut chars);
                if digits.is_empty() {
                    T::TripleBond
                } else {
                    T::Atom(digits.parse().unwrap())
                }
            }
            'H' => T::HCount(get_digits(&mut chars).parse().unwrap_or(1)),
            ':' => T::Colon,
            '0'..='9' => T::Digit(
                // combine the digit in c with any following digits
                format!("{c}{}", get_digits(&mut chars)).parse().unwrap(),
            ),
            '-' => T::SingleBond,
            _ => panic!("unrecognized token {c}"),
        };
        ret.push(got);
    }
    ret
}

#[test]
fn test_scan() {
    let s = r#"[#6H3:1]-[#6H2:2]-[#7H:3]-[#7H:4]-[#6H3:5]"#;
    dbg!(scan(s.to_owned()));
}
