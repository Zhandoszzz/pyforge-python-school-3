from rdkit import Chem


def substructure_search(mols, mol):
    mol = Chem.MolFromSmiles(mol)
    result = [m for m in mols if Chem.MolFromSmiles(m).HasSubstructMatch(mol)]
    return result
