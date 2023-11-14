//! SMARTS pattern parser

use crate::smarts::parser::Parser;

use self::{evaluator::Evaluator, scanner::scan};

mod evaluator;
mod parser;
mod scanner;

#[derive(Clone, Default, Debug, PartialEq)]
pub enum Chiral {
    Cw,
    Acw,
    #[default]
    None,
}

#[derive(Clone, Debug, PartialEq)]
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

#[derive(Clone, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Ring,
    Up,
    Down,
}

impl std::fmt::Debug for BondOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                BondOrder::Single => "-",
                BondOrder::Double => "=",
                BondOrder::Triple => "#",
                BondOrder::Aromatic => ":",
                BondOrder::Ring => "@",
                BondOrder::Up => "/",
                BondOrder::Down => "\\",
            }
        )
    }
}
impl std::fmt::Debug for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{:?}{}", self.atom1, self.order, self.atom2)
    }
}

#[derive(PartialEq)]
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
        let mut parser = Parser::new(tokens);
        let exprs = parser.parse();
        let eval = Evaluator::new(exprs);
        let (atoms, bonds) = eval.eval();
        Self { atoms, bonds }
    }
}
