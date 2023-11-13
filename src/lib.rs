use std::{collections::HashMap, error::Error, fs::File, path::Path};

use serde::Deserialize;

#[derive(Deserialize)]
struct Record {
    cmiles: String,
}

#[derive(Deserialize)]
pub struct Dataset {
    entries: HashMap<String, Vec<Record>>,
}

impl Dataset {
    pub fn load(path: impl AsRef<Path>) -> Result<Dataset, Box<dyn Error>> {
        let f = File::open(path)?;
        let r: Self = serde_json::from_reader(f)?;
        Ok(r)
    }

    /// consume `self` and return the contained vector of canonical SMILES
    /// strings
    pub fn to_smiles(self) -> Vec<String> {
        self.entries
            .into_values()
            .flatten()
            .map(|v| v.cmiles)
            .collect()
    }
}
