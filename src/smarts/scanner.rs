use std::{iter::Peekable, str::Chars};

#[derive(Clone, Debug, PartialEq)]
pub(super) enum Token {
    // punctuation
    LBrack,
    RBrack,
    LParen,
    RParen,
    Colon,
    Dash, // could be bond or charge at this point
    At,
    // counts
    Atom(usize),
    HCount(usize),
    Digit(usize),
    Plus(usize),
    // bonds
    DoubleBond,
    TripleBond,
    UpBond,
    DownBond,
    // end
    End,
}

impl Token {
    /// Returns `true` if the token is [`End`].
    ///
    /// [`End`]: Token::End
    #[must_use]
    pub(super) fn is_end(&self) -> bool {
        matches!(self, Self::End)
    }
}

fn get_digits(chars: &mut Peekable<Chars<'_>>) -> String {
    let mut digits = String::new();
    while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
        digits.push(chars.next().unwrap());
    }
    digits
}

pub(super) fn scan(s: String) -> Vec<Token> {
    use Token as T;
    let mut chars = s.chars().peekable();
    let mut ret = Vec::new();
    while let Some(c) = chars.next() {
        let got = match c {
            '[' => T::LBrack,
            ']' => T::RBrack,
            '(' => T::LParen,
            ')' => T::RParen,
            ':' => T::Colon,
            '-' => T::Dash,
            '@' => T::At,
            '=' => T::DoubleBond,
            '\\' => T::DownBond,
            '/' => T::UpBond,
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
            '+' => T::Plus(get_digits(&mut chars).parse().unwrap_or(1)),
            '0'..='9' => T::Digit(
                // combine the digit in c with any following digits
                format!("{c}{}", get_digits(&mut chars)).parse().unwrap(),
            ),
            _ => panic!("unrecognized token {c} in \n{s}"),
        };
        ret.push(got);
    }
    ret.push(T::End);
    ret
}

#[cfg(test)]
mod tests {
    use crate::{rdkit::to_smarts, Dataset};

    use super::*;

    #[test]
    fn simple_scan() {
        let s = r#"[#6H3:1]-[#6H2:2]-[#7H:3]-[#7H:4]-[#6H3:5]"#;
        scan(s.to_owned());
    }

    #[test]
    fn big_scan() {
        let mut smiles =
            Dataset::load("testfiles/opt.json").unwrap().to_smiles();
        smiles.dedup();
        for smile in smiles {
            scan(to_smarts(smile));
        }
    }
}
