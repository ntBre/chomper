//! SMARTS pattern parser

use crate::smarts::parser::Parser;

use self::scanner::scan;

mod parser;
mod scanner;

#[derive(Debug, PartialEq)]
pub struct Atom {
    pub atomic_number: usize,
    pub n_hydrogens: usize,
    pub mol_index: usize,
}

impl Atom {
    pub fn new(
        atomic_number: usize,
        n_hydrogens: usize,
        mol_index: usize,
    ) -> Self {
        Self {
            atomic_number,
            n_hydrogens,
            mol_index,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Up,
    Down,
}

#[derive(Debug, PartialEq)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
}

impl Bond {
    pub fn new(atom1: usize, atom2: usize, order: BondOrder) -> Self {
        Self {
            atom1,
            atom2,
            order,
        }
    }
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
