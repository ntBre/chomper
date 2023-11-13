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

    fn context(&self, n: usize) -> &[Token] {
        let beg = self.cur.saturating_sub(n);
        let end = self.cur + n;
        &self.tokens[beg..end.min(self.tokens.len())]
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

    pub(super) fn parse(&mut self) -> (Vec<Atom>, Vec<Bond>) {
        let mut atoms = Vec::new();
        let mut bonds: Vec<Bond> = Vec::new();
        // for updating bonds. this is not going to work at all for groupings
        let mut atom1 = 0;
        while !self.at_end() {
            match self.peek() {
                Token::LBrack => {
                    let atom = self.atom();
                    if !atoms.is_empty() {
                        // a bond must have come before
                        bonds.last_mut().unwrap().atom2 = atom.mol_index;
                    }
                    atom1 = atom.mol_index;
                    atoms.push(atom);
                }
                Token::LParen => {
                    let (a, b) = self.grouping();
                    atoms.extend(a);
                    bonds.extend(b);
                }
                Token::RParen => break, // for recursive calls from grouping
                _ => {
                    let mut bond = self.bond();
                    bond.atom1 = atom1;
                    bonds.push(bond)
                }
            }
        }
        (atoms, bonds)
    }

    fn atom(&mut self) -> Atom {
        self.advance(); // discard LBrack signaling we're in here
        let mut atomic_number = 0;
        let mut n_hydrogens = 0;
        let mut mol_index = 0;
        let mut charge = 0;
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
                Token::Plus(n) => charge = n as isize,
                Token::Dash => {
                    let t = self.peek().clone();
                    let n = if let Token::Digit(n) = t {
                        self.advance();
                        n
                    } else {
                        1
                    };
                    charge = -(n as isize);
                }
                Token::RBrack => break,
                Token::End => panic!("EOF while parsing atom"),
                x => self.error("atom", x),
            };
        }
        Atom {
            atomic_number,
            n_hydrogens,
            charge,
            mol_index,
        }
    }

    fn error(&self, label: &str, t: Token) -> ! {
        panic!(
            "unknown {label} component {t:?} at {}: {:?}",
            self.cur,
            self.context(5)
        )
    }

    fn grouping(&mut self) -> (Vec<Atom>, Vec<Bond>) {
        self.advance(); // discard LParen
        let ret = self.parse();
        self.advance(); // discard closing RParen
        ret
    }

    fn bond(&mut self) -> Bond {
        // we can't actually know atom2 when we call this. we'll have to tidy
        // that up after the call. in fact, we can only determine the bond order
        // from what we're parsing here
        let mut order = BondOrder::Single;
        let mut ring_marker = None;
        loop {
            if matches!(self.peek(), Token::LBrack | Token::RParen) {
                // signals start of next atom or end of current group, don't
                // consume in either case
                break;
            }

            match self.advance() {
                Token::Dash => order = BondOrder::Single,
                Token::Colon => {
                    // ring_marker can come on either side of the order symbols.
                    // TODO handle dash case too
                    if let Token::Digit(n) = self.peek() {
                        ring_marker = Some(*n);
                        self.advance();
                    }
                    order = BondOrder::Aromatic;
                }
                Token::DoubleBond => order = BondOrder::Double,
                Token::Digit(n) => {
                    ring_marker = Some(n);
                    match self.advance() {
                        Token::Dash => order = BondOrder::Single,
                        Token::Colon => order = BondOrder::Aromatic,
                        Token::DoubleBond => order = BondOrder::Double,
                        Token::End => break,
                        x => self.error("inner bond", x),
                    }
                }
                Token::End => break,
                x => self.error("bond", x),
            }
        }
        Bond {
            atom1: 0,
            atom2: 0,
            order,
            ring_marker,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{rdkit::to_smarts, smarts::scanner::scan, Dataset};

    use super::*;

    #[test]
    fn parse_single() {
        let s = r#"[#6H3:1]-[#6H2:2]-[#7H:3]-[#7H:4]-[#6H3:5]"#;
        let tokens = scan(s.to_owned());
        let got = Parser::new(tokens).parse();
        let want_atoms = vec![
            Atom::new(6, 3, 0, 1),
            Atom::new(6, 2, 0, 2),
            Atom::new(7, 1, 0, 3),
            Atom::new(7, 1, 0, 4),
            Atom::new(6, 3, 0, 5),
        ];
        use BondOrder as B;
        let want_bonds = vec![
            Bond::new(1, 2, B::Single, None),
            Bond::new(2, 3, B::Single, None),
            Bond::new(3, 4, B::Single, None),
            Bond::new(4, 5, B::Single, None),
        ];
        assert_eq!(got, (want_atoms, want_bonds));
    }

    #[test]
    fn parse_all() {
        let mut smiles =
            Dataset::load("testfiles/opt.json").unwrap().to_smiles();
        smiles.dedup();
        for smile in smiles {
            let smarts = to_smarts(smile);
            let tokens = scan(smarts);
            Parser::new(tokens).parse();
        }
    }
}
