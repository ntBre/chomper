use pyo3::{types::PyModule, Python};

pub fn to_smarts(smiles: String) -> String {
    Python::with_gil(|py| {
        let chem = PyModule::import(py, "rdkit.Chem").unwrap();
        let mol = chem.call_method1("MolFromSmiles", (smiles,)).unwrap();
        chem.call_method1("MolToSmarts", (mol, true))
            .unwrap()
            .extract()
            .unwrap()
    })
}
