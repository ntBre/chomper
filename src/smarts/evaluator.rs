use std::collections::HashMap;

use super::{parser::Expr, Atom, Bond, BondOrder};

/// search `exprs` for the first atom occuring before `cur`
fn prev_atom(exprs: &[Expr], cur: usize) -> Option<&Expr> {
    if cur == 0 {
        return None;
    }
    for p in (0..cur - 1).rev() {
        let r = &exprs[p];
        if r.is_atom() {
            return Some(r);
        }
    }
    None
}

fn next_atom(exprs: &[Expr], cur: usize) -> Option<&Expr> {
    for p in cur..exprs.len() {
        let r = &exprs[p];
        if !r.is_connect() && !r.is_grouping() {
            return Some(r);
        }
    }
    None
}

pub(super) struct Evaluator {
    exprs: Vec<Expr>,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    cur: usize,
    /// connection table for ring bongs
    ctab: HashMap<usize, Vec<usize>>,
}

impl Evaluator {
    pub(super) fn new(exprs: Vec<Expr>) -> Self {
        Self {
            exprs,
            atoms: Vec::new(),
            bonds: Vec::new(),
            cur: 0,
            ctab: HashMap::new(),
        }
    }

    fn at_end(&self) -> bool {
        self.cur == self.exprs.len()
    }

    fn prev(&self) -> &Expr {
        &self.exprs[self.cur - 2]
    }

    fn peek(&self) -> Option<&Expr> {
        self.exprs.get(self.cur)
    }

    fn prev_atom(&self) -> &Expr {
        prev_atom(&self.exprs, self.cur).unwrap()
    }

    fn next(&mut self) -> Expr {
        let ret = self.exprs[self.cur].clone();
        if self.cur < self.exprs.len() {
            self.cur += 1;
        };
        ret
    }

    // I think we're actually going to need to do the next/prev stuff from the
    // parser so we can look ahead and behind as neede
    pub(crate) fn eval(mut self) -> (Vec<Atom>, Vec<Bond>) {
        while !self.at_end() {
            let expr = self.next();
            self.inner(expr);
        }
        let Evaluator { atoms, bonds, .. } = self;
        (atoms, bonds)
    }

    fn inner(&mut self, expr: Expr) {
        match expr {
            Expr::Atom(a) => {
                self.atom(a);
            }
            Expr::Bond(order) => {
                self.bond(order);
            }
            Expr::Grouping(g) => {
                self.grouping(g);
            }
            Expr::Connect(_) => {
                todo!()
            }
        }
    }

    fn bond(&mut self, order: BondOrder) {
        let atom1 = match self.prev_atom() {
            Expr::Atom(a) => a.mol_index,
            Expr::Bond(_) => todo!("{}", self.cur),
            Expr::Grouping(_) => todo!(),
            Expr::Connect(_) => todo!(),
        };
        // cloning so we can remove below
        let atom2 = match self.peek().unwrap().clone() {
            Expr::Atom(a) => a.mol_index,
            Expr::Bond(_) => todo!("{}", self.cur),
            Expr::Grouping(_) => todo!(),
            Expr::Connect(n) => {
                self.next(); // advance over connection
                self.get_connection(n)
            }
        };
        self.bonds.push(Bond {
            atom1,
            atom2,
            order,
        });
    }

    fn get_connection(&mut self, n: usize) -> usize {
        self.ctab.get_mut(&n).unwrap().pop().unwrap()
    }

    fn atom(&mut self, a: Atom) {
        if let Some(&Expr::Connect(n)) = self.peek() {
            self.next();
            // treat ctab as a stack so we can reuse labels
            self.ctab.entry(n).or_insert(Vec::new()).push(a.mol_index);
        }
        self.atoms.push(a);
    }

    fn grouping(&mut self, g: Vec<Expr>) {
        let mut giter = g.iter().enumerate().peekable();
        while let Some((i, expr)) = giter.next() {
            match expr {
                Expr::Atom(a) => {
                    if let Some((_, &Expr::Connect(n))) = giter.peek() {
                        giter.next();
                        // treat ctab as a stack so we can reuse labels
                        self.ctab
                            .entry(n)
                            .or_insert(Vec::new())
                            .push(a.mol_index);
                    }
                    self.atoms.push(a.clone());
                }
                Expr::Bond(order) => {
                    let atom1 = match prev_atom(&g, i + 1) {
                        Some(Expr::Atom(a)) => a.mol_index,
                        None => {
                            // bond is to the atom before the group. can't use
                            // self.prev_atom because we could be in a recursive
                            // group
                            self.atoms.last().unwrap().mol_index
                        }
                        _ => unreachable!(),
                    };
                    let atom2 = match giter.peek().unwrap().1 {
                        Expr::Atom(a) => a.mol_index,
                        Expr::Connect(n) => {
                            giter.next(); // discard Connect expr
                            self.get_connection(*n)
                        }
                        _ => unreachable!(),
                    };
                    let (atom1, atom2) = (atom1.min(atom2), atom1.max(atom2));
                    let bond = Bond {
                        atom1,
                        atom2,
                        order: order.clone(),
                    };
                    self.bonds.push(bond);
                }
                Expr::Grouping(h) => self.grouping(h.clone()),
                Expr::Connect(_) => panic!("{i}"),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        rdkit::to_smarts,
        smarts::{parser::Parser, scanner::scan, Chiral, Smarts},
        Dataset,
    };

    use super::*;

    #[test]
    fn problems() {
        let smiles = [
            "[H:12][C:1]([H:13])([H:14])[C:2]([H:15])([C:3](=[O:4])[C:5]1=[N:6]\
             [S:7](=[O:8])[O:9][N:10]1[H:16])[O:11][H:17]"
        ];
        use BondOrder as B;
        let wants = [Smarts {
            atoms: vec![
                Atom::new(6, 3, 0, Chiral::None, 1),
                Atom::new(6, 1, 0, Chiral::None, 2),
                Atom::new(6, 0, 0, Chiral::None, 3),
                Atom::new(8, 0, 0, Chiral::None, 4),
                Atom::new(6, 0, 0, Chiral::None, 5),
                Atom::new(7, 0, 0, Chiral::None, 6),
                Atom::new(16, 0, 0, Chiral::None, 7),
                Atom::new(8, 0, 0, Chiral::None, 8),
                Atom::new(8, 0, 0, Chiral::None, 9),
                Atom::new(7, 1, 0, Chiral::None, 10),
                Atom::new(8, 1, 0, Chiral::None, 11),
            ],
            bonds: vec![
                Bond::new(1, 2, B::Single),
                Bond::new(2, 3, B::Single),
                Bond::new(3, 4, B::Double),
                Bond::new(3, 5, B::Single),
                Bond::new(5, 6, B::Double),
                Bond::new(6, 7, B::Single),
                Bond::new(7, 8, B::Double),
                Bond::new(7, 9, B::Single),
                Bond::new(9, 10, B::Single),
                Bond::new(5, 10, B::Single),
                Bond::new(2, 11, B::Single),
            ],
        }];
        for (smile, want) in smiles.into_iter().zip(wants.into_iter()) {
            let smarts = to_smarts(dbg!(smile).to_owned());
            let tokens = scan(dbg!(smarts));
            let p = Parser::new(tokens).parse();
            let (atoms, bonds) = Evaluator::new(dbg!(p)).eval();
            assert_eq!(atoms, want.atoms);
            assert_eq!(bonds, want.bonds);
        }
    }

    #[test]
    fn all() {
        let mut smiles =
            Dataset::load("testfiles/opt.json").unwrap().to_smiles();
        smiles.dedup();
        for smile in smiles {
            let smarts = to_smarts(dbg!(smile));
            let tokens = scan(dbg!(smarts));
            let p = Parser::new(tokens).parse();
            Evaluator::new(dbg!(p)).eval();
        }
    }
}
