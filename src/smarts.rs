//! SMARTS pattern parser

use self::scanner::scan;

mod scanner;
mod parser;

pub struct Smarts;

impl Smarts {
    pub fn parse(s: String) -> Self {
        let tokens = scan(s);
        todo!()
    }
}
