//! SMARTS pattern parser

use crate::smarts::parser::Parser;

use self::scanner::scan;

mod parser;
mod scanner;

#[derive(Default, Debug, PartialEq)]
pub enum Chiral {
    Cw,
    Acw,
    #[default]
    None,
}

#[derive(Debug, PartialEq)]
pub struct Atom {
    pub atomic_number: usize,
    pub n_hydrogens: usize,
    pub charge: isize,
    pub chirality: Chiral,
    pub mol_index: usize,
}

impl Atom {
    pub fn new(
        atomic_number: usize,
        n_hydrogens: usize,
        charge: isize,
        chirality: Chiral,
        mol_index: usize,
    ) -> Self {
        Self {
            atomic_number,
            n_hydrogens,
            charge,
            chirality,
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
    Ring,
    Up,
    Down,
}

#[derive(Debug, PartialEq)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
    /// optional digit for specifying a backreference to a connection
    pub ring_marker: Option<usize>,
}

impl Bond {
    pub fn new(
        atom1: usize,
        atom2: usize,
        order: BondOrder,
        ring_marker: Option<usize>,
    ) -> Self {
        Self {
            atom1,
            atom2,
            order,
            ring_marker,
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
        let mut parser = Parser::new(tokens);
        let (atoms, bonds) = parser.parse();
        Self { atoms, bonds }
    }
}
