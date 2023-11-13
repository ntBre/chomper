//! SMARTS pattern parser

use crate::smarts::parser::Parser;

use self::scanner::scan;

mod parser;
mod scanner;

#[derive(Debug)]
pub struct Atom {
    pub atomic_number: usize,
    pub n_hydrogens: usize,
    pub mol_index: usize,
}

#[derive(Debug)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Up,
    Down,
}

#[derive(Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
}

pub struct Smarts {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl Smarts {
    pub fn parse(s: String) -> Self {
        let tokens = scan(s);
        let parser = Parser::new(tokens);
        let (atoms, bonds) = parser.parse();
        Self { atoms, bonds }
    }
}
