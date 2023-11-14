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

use super::{scanner::Token, Atom, BondOrder, Chiral};

#[derive(Clone, PartialEq)]
pub enum Expr {
    Atom(Atom),
    Bond(BondOrder),
    Grouping(Vec<Expr>),
    Connect(usize),
}

impl std::fmt::Debug for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Expr::Atom(a) => write!(
                f,
                "[#{}H{}{:+}:{}]",
                a.atomic_number, a.n_hydrogens, a.charge, a.mol_index
            ),
            Expr::Bond(order) => write!(f, "{order:?}"),
            Expr::Grouping(g) => write!(f, "Grouping({g:?})"),
            Expr::Connect(n) => write!(f, "Connect({n})"),
        }
    }
}

impl Expr {
    /// Returns `true` if the expr is [`Connect`].
    ///
    /// [`Connect`]: Expr::Connect
    #[must_use]
    pub fn is_connect(&self) -> bool {
        matches!(self, Self::Connect(..))
    }

    /// Returns `true` if the expr is [`Grouping`].
    ///
    /// [`Grouping`]: Expr::Grouping
    #[must_use]
    pub fn is_grouping(&self) -> bool {
        matches!(self, Self::Grouping(..))
    }

    /// Returns `true` if the expr is [`Atom`].
    ///
    /// [`Atom`]: Expr::Atom
    #[must_use]
    pub fn is_atom(&self) -> bool {
        matches!(self, Self::Atom(..))
    }

    /// Returns `true` if the expr is [`Bond`].
    ///
    /// [`Bond`]: Expr::Bond
    #[must_use]
    pub fn is_bond(&self) -> bool {
        matches!(self, Self::Bond(..))
    }
}

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

    pub(super) fn parse(&mut self) -> Vec<Expr> {
        let mut ret = Vec::new();
        while !self.at_end() {
            match self.peek() {
                Token::LBrack => ret.push(self.atom()),
                Token::LParen => ret.push(self.grouping()),
                Token::RParen => break, // for recursive calls from grouping
                Token::Digit(n) => {
                    ret.push(Expr::Connect(*n));
                    self.advance();
                }
                _ => ret.push(self.bond()),
            }
        }
        ret
    }

    fn atom(&mut self) -> Expr {
        self.advance(); // discard LBrack signaling we're in here
        let mut chirality = Chiral::None;
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
                Token::At => chirality = Chiral::Acw,
                Token::AtAt => chirality = Chiral::Cw,
                Token::RBrack => break,
                Token::End => panic!("EOF while parsing atom"),
                x => self.error("atom", x),
            };
        }
        Expr::Atom(Atom {
            atomic_number,
            n_hydrogens,
            charge,
            chirality,
            mol_index,
        })
    }

    fn error(&self, label: &str, t: Token) -> ! {
        panic!(
            "unknown {label} component {t:?} at {}: {:?}",
            self.cur,
            self.context(5)
        )
    }

    fn grouping(&mut self) -> Expr {
        self.advance(); // discard LParen
        let ret = self.parse();
        self.advance(); // discard closing RParen
        Expr::Grouping(ret)
    }

    fn bond(&mut self) -> Expr {
        match self.advance() {
            Token::Dash => Expr::Bond(BondOrder::Single),
            Token::DoubleBond => Expr::Bond(BondOrder::Double),
            Token::Colon => Expr::Bond(BondOrder::Aromatic),
            Token::TripleBond => Expr::Bond(BondOrder::Triple),
            Token::At => Expr::Bond(BondOrder::Ring),
            Token::DownBond => Expr::Bond(BondOrder::Down),
            Token::UpBond => Expr::Bond(BondOrder::Up),
            x => self.error("bond", x),
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
        let want = vec![
            Expr::Atom(Atom::new(6, 3, 0, Chiral::None, 1)),
            Expr::Bond(BondOrder::Single),
            Expr::Atom(Atom::new(6, 2, 0, Chiral::None, 2)),
            Expr::Bond(BondOrder::Single),
            Expr::Atom(Atom::new(7, 1, 0, Chiral::None, 3)),
            Expr::Bond(BondOrder::Single),
            Expr::Atom(Atom::new(7, 1, 0, Chiral::None, 4)),
            Expr::Bond(BondOrder::Single),
            Expr::Atom(Atom::new(6, 3, 0, Chiral::None, 5)),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn parse_problems() {
        let smiles = [
            "[H:12][C:1]([H:13])([H:14])[C:2]([H:15])([C:3](=[O:4])[C:5]1=[N:6]\
             [S:7](=[O:8])[O:9][N:10]1[H:16])[O:11][H:17]"
        ];
        use BondOrder as B;
        let wants = [vec![
            Expr::Atom(Atom::new(6, 3, 0, Chiral::None, 1)),
            Expr::Bond(B::Single),
            Expr::Atom(Atom::new(6, 1, 0, Chiral::None, 2)),
            Expr::Grouping(vec![
                Expr::Bond(B::Single),
                Expr::Atom(Atom::new(6, 0, 0, Chiral::None, 3)),
                Expr::Grouping(vec![
                    Expr::Bond(B::Double),
                    Expr::Atom(Atom::new(8, 0, 0, Chiral::None, 4)),
                ]),
                Expr::Bond(B::Single),
                Expr::Atom(Atom::new(6, 0, 0, Chiral::None, 5)),
                Expr::Connect(1),
                Expr::Bond(B::Double),
                Expr::Atom(Atom::new(7, 0, 0, Chiral::None, 6)),
                Expr::Bond(B::Single),
                Expr::Atom(Atom::new(16, 0, 0, Chiral::None, 7)),
                Expr::Grouping(vec![
                    Expr::Bond(B::Double),
                    Expr::Atom(Atom::new(8, 0, 0, Chiral::None, 8)),
                ]),
                Expr::Bond(B::Single),
                Expr::Atom(Atom::new(8, 0, 0, Chiral::None, 9)),
                Expr::Bond(B::Single),
                Expr::Atom(Atom::new(7, 1, 0, Chiral::None, 10)),
                Expr::Bond(B::Single),
                Expr::Connect(1),
            ]),
            Expr::Bond(B::Single),
            Expr::Atom(Atom::new(8, 1, 0, Chiral::None, 11)),
        ]];
        for (i, smile) in smiles.into_iter().enumerate() {
            let smarts = to_smarts(smile.to_owned());
            let got = Parser::new(scan(smarts)).parse();
            let want = wants[i].clone();
            assert_eq!(got, want);
        }
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
