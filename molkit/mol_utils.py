from rdkit.Chem import AllChem, Draw


def smiles_to_2d_mol(smiles):
    try:
        mol = AllChem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        coords = mol.GetConformer().GetPositions()
        max_x = max(coords, key=lambda x: x[0])[0]
        min_x = min(coords, key=lambda x: x[0])[0]
        max_y = max(coords, key=lambda x: x[1])[1]
        min_y = min(coords, key=lambda x: x[1])[1]
        return Draw.MolToImage(mol), max_x - min_x, max_y - min_y
    except:
        return None, None, None


def smiles_to_3d_mol(smiles):
    try:
        mol = AllChem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return Draw.MolToImage(mol)
    except:
        return None
