use chomper::{rdkit::to_smarts, Dataset};

fn main() {
    let mut smiles = Dataset::load("testfiles/opt.json").unwrap().to_smiles();
    smiles.dedup();
    for smile in smiles {
        println!("{:.80}", to_smarts(smile));
    }
}

// main idea is to read a dataset to get SMILES, convert the smiles to smarts
// with rdkit (MolToSmarts), then start processing the smarts
