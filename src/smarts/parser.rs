//! Parser for SMARTS. Grammar:
//!
//! smarts -> atom | atom bond smarts
//! atom -> "[" "#" DIGIT+ ("H" DIGIT*)* ":" DIGIT+ "]"
//! bond -> DIGIT? grouping* ( "-" | "/" | "\" | "=" | "#" | ":" | "@" )
//! grouping -> "(" bond smarts ")"
//!
//! a smarts is either just an atom or an atom followed by a bond and further
//! smarts
//!
//! an atom is something inside of square brackets. again, this can be a lot
//! more complicated but maybe not for the concrete patterns I intend to parse
//!
//! a bond is just one of the bond symbols optionally prefixed by a digit for
//! connecting back later in the string. it can also be a grouping. digit and
//! grouping might be mutually exclusive, but I'm not sure
//!
//! I'm pretty sure there are actually further complications to this. for
//! example, I think you can apply modifiers to bonds, and there are logical
//! operators in here somewhere too. I'm not running into these yet because I'm
//! labeling concrete molecules. So this might only be of concern for generated
//! Smarts structs, not actually parsing
//!
//! Usually you would turn the sequence of tokens into an AST, which you can
//! then evaluate, but I think I can turn my tokens directly into my desired
//! Smarts struct

use super::{scanner::Token, Atom, Bond, BondOrder};

pub(super) struct Parser {
    /// `tokens` represents a single input SMARTS string decomposed into a
    /// sequence of tokens. We turn this sequence back into a [Smarts] struct
    tokens: Vec<Token>,
    cur: usize,
}

impl Parser {
    pub(super) fn new(tokens: Vec<Token>) -> Self {
        Self { tokens, cur: 0 }
    }

    #[inline]
    fn at_end(&self) -> bool {
        self.peek().is_end()
    }

    fn peek(&self) -> &Token {
        &self.tokens[self.cur]
    }

    fn advance(&mut self) -> Token {
        let ret = self.peek().clone();
        if !self.at_end() {
            self.cur += 1;
        }
        ret
    }

    pub(super) fn parse(mut self) -> (Vec<Atom>, Vec<Bond>) {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        while !self.at_end() {
            match self.peek() {
                Token::LBrack => atoms.push(self.atom()),
                Token::LParen => {
                    let (a, b) = self.grouping();
                    atoms.extend(a);
                    bonds.extend(b);
                }
                _ => bonds.push(self.bond()),
            }
        }
        (atoms, bonds)
    }

    fn atom(&mut self) -> Atom {
        self.advance(); // discard LBrack signaling we're in here
        let mut atomic_number = 0;
        let mut n_hydrogens = 0;
        let mut mol_index = 0;
        loop {
            match self.advance() {
                Token::Atom(n) => atomic_number = n,
                Token::HCount(n) => n_hydrogens = n,
                Token::Colon => {
                    let Token::Digit(i) = self.advance() else {
                        unreachable!();
                    };
                    mol_index = i;
                }
                Token::RBrack => break,
                Token::End => panic!("EOF while parsing atom"),
                x => panic!("unknown atom component: {x:?}"),
            };
        }
        Atom {
            atomic_number,
            n_hydrogens,
            mol_index,
        }
    }

    fn grouping(&mut self) -> (Vec<Atom>, Vec<Bond>) {
        self.advance(); // discard LParen
        let mut tokens = Vec::new();
        while !matches!(self.peek(), Token::RParen) {
            tokens.push(self.advance());
        }
        self.advance(); // discard closing RParen
        Parser::new(tokens).parse()
    }

    fn bond(&mut self) -> Bond {
        // we can't actually know atom2 when we call this. we'll have to tidy
        // that up after the call. in fact, we can only determine the bond order
        // from what we're parsing here
        let mut order = BondOrder::Single;
        loop {
            if matches!(self.peek(), Token::LBrack) {
                break; // signals start of next atom, don't consume
            }

            match self.advance() {
                Token::Dash => order = BondOrder::Single,
                x => panic!("unknown bond component: {x:?}"),
            }
        }
        Bond {
            atom1: 0,
            atom2: 0,
            order,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::smarts::scanner::scan;

    use super::*;

    #[test]
    fn parse_single() {
        let s = r#"[#6H3:1]-[#6H2:2]-[#7H:3]-[#7H:4]-[#6H3:5]"#;
        let tokens = scan(s.to_owned());
        let ret = Parser::new(tokens).parse();
        dbg!(ret);
    }
}
