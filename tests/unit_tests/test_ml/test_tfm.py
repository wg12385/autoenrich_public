
import autoenrich.ml.features.BCAI_calc.mol_graph_setup as BCAI
from test_generators.dummy_mol import get_ethane_mol, get_random_mol, get_random_ethane

from autoenrich.reference.periodic_table import Get_periodic_table

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

def get_mols():
    p_table = Get_periodic_table()

    mols = []
    for i in range(4):
        mol = get_random_ethane()
        mol.molid = 'random' + str(i)
        mol.get_distance_matrix(heavy_only=False)
        mols.append(mol)
    mol = get_ethane_mol()
    mol.get_distance_matrix(heavy_only=False)
    mols.append(mol)

    return mols

def test_make_atom_df():

    mols = get_mols()
    ats = 0
    for mol in mols:
        ats += len(mol.types)

    atoms = BCAI.make_atom_df(mols)
    assert len(atoms["molecule_name"].unique()) == len(mols)

    counted = 0
    for i, idx in enumerate(atoms['atom_index'].values):
        for mol in mols:
            if atoms['molecule_name'][i] == mol.molid:
                counted += 1
                assert np.array_equal(mol.conn[idx], atoms['conn'][i])
                assert mol.xyz[idx][0] == atoms['x'][i]
                assert mol.xyz[idx][1] == atoms['y'][i]
                assert mol.xyz[idx][2] == atoms['z'][i]

    assert counted == ats


def test_make_structure_dict():

    p_table = Get_periodic_table()

    mols = get_mols()
    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    assert len(structure_dict.keys()) == len(mols)

    for mol in mols:
        assert structure_dict[mol.molid]['symbols'] == [p_table[type] for type in mol.types]
        assert np.array_equal(structure_dict[mol.molid]['positions'], mol.xyz)
        assert np.array_equal(structure_dict[mol.molid]['conn'], mol.conn)


def test_enhance_structure_dict():

    p_table = Get_periodic_table()
    #####
    mols = get_mols()
    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    #####

    BCAI.enhance_structure_dict(structure_dict)

    for mol in mols:
        assert structure_dict[mol.molid]['symbols'] == [p_table[type] for type in mol.types]
        assert np.array_equal(structure_dict[mol.molid]['positions'], mol.xyz)
        assert np.array_equal(structure_dict[mol.molid]['conn'], mol.conn)
        assert np.array_equal(structure_dict[mol.molid]['distances'], mol.dist)


def test_enhance_atoms():

    p_table = Get_periodic_table()

    #####
    mols = get_mols()
    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    BCAI.enhance_structure_dict(structure_dict)
    ###########

    BCAI.enhance_atoms(atoms, structure_dict)

    for i, idx in enumerate(atoms['atom_index'].values):

        molid = atoms['molecule_name'][i]
        mol = 0
        for ml_fnd in mols:
            if ml_fnd.molid == molid:
                mol = ml_fnd

        atid = atoms['atom_index'][i]
        assert p_table.index(atoms['atom'][i]) == mol.types[atid]
        assert np.array_equal(atoms['conn'][i], mol.conn[atid])
        assert np.array_equal(atoms['distance'][i], mol.dist[atid])


def test_make_bonds_df():

    p_table = Get_periodic_table()

    #####
    mols = get_mols()
    #####

    bonds = BCAI.make_bonds_df(mols)
    assert len(bonds["molecule_name"].unique()) == len(mols)

    for idx, bond in enumerate(bonds):

        molid = bonds['molecule_name'][idx]

        at1 = bonds['atom_index_0'][idx]
        at2 = bonds['atom_index_1'][idx]

        mol = 0
        for ml_fnd in mols:
            if ml_fnd.molid == molid:
                mol = ml_fnd

        assert mol.coupling_len[at1][at2] == int(bonds['type'][idx][0])
        assert mol.coupling[at1][at2] == bonds['scalar_coupling_constant'][idx]


def test_enhance_bonds():

    p_table = Get_periodic_table()

    #####
    mols = get_mols()

    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    BCAI.enhance_structure_dict(structure_dict)
    BCAI.enhance_atoms(atoms, structure_dict)

    bonds = BCAI.make_bonds_df(mols)
    ############

    BCAI.enhance_bonds(bonds, structure_dict, flag='3JHH')

    for idx, bond in enumerate(bonds):
        molid = bonds['molecule_name'][idx]
        at1 = bonds['atom_index_0'][idx]
        at2 = bonds['atom_index_1'][idx]

        mol = 0
        for ml_fnd in mols:
            if ml_fnd.molid == molid:
                mol = ml_fnd

        assert mol.coupling_len[at1][at2] == int(bonds['type'][idx][0])
        assert mol.coupling[at1][at2] == bonds['scalar_coupling_constant'][idx]
        if bonds['labeled_type'][idx] == '3JHH':
            assert bonds['predict'][idx] == 1
        else:
            assert bonds['predict'][idx] == 0

    for mol in mols:
        for atom1, type1 in enumerate(mol.types):
            for atom2, type2 in enumerate(mol.types):

                if atom1 == atom2:
                    continue

                row = bonds.loc[(bonds['molecule_name'] == mol.molid)
                        & (bonds['atom_index_0'] == atom1)
                        & (bonds['atom_index_1'] == atom2)]

                cpl = row['scalar_coupling_constant'].values

                if type1 == 1 and type2 == 1 and mol.coupling_len[atom1][atom2] == 3:
                    assert row.predict.values == 1
                    assert mol.coupling[atom1][atom2] == cpl[0]



def test_make_triplets():

    p_table = Get_periodic_table()

    #####
    mols = get_mols()

    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    BCAI.enhance_structure_dict(structure_dict)
    BCAI.enhance_atoms(atoms, structure_dict)

    bonds = BCAI.make_bonds_df(mols)
    BCAI.enhance_bonds(bonds, structure_dict, flag='3JHH')
    ############

    triplets = BCAI.make_triplets(bonds["molecule_name"].unique(), structure_dict)
    assert len(triplets["molecule_name"].unique()) == len(mols)

    count = 0
    for mol in mols:
        for atom1, type1 in enumerate(mol.types):
            for atom2, type2 in enumerate(mol.types):
                if atom1 == atom2:
                    continue

                for atom3, type3 in enumerate(mol.types):

                    if atom3 in [atom1, atom2] or atom3 < atom2:
                        continue


                    if mol.conn[atom1][atom2] != 1 or mol.conn[atom1][atom3] != 1:
                        continue

                    row = triplets.loc[(triplets.molecule_name == mol.molid)
                                & (triplets.atom_index_0 == atom1)
                                & (triplets.atom_index_1 == atom2)
                                & (triplets.atom_index_2 == atom3)]

                    assert len(row.index) == 1

                    ba = mol.xyz[atom2] - mol.xyz[atom1]
                    bc = mol.xyz[atom3] - mol.xyz[atom1]

                    angle = np.sum(ba * bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                    angle = np.arccos(np.clip(angle, -1.0, 1.0))

                    assert angle == row.angle.values
                    count += 1

    assert count == len(triplets.index)


def test_add_embedding():

    p_table = Get_periodic_table()

    #####
    mols = get_mols()

    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    BCAI.enhance_structure_dict(structure_dict)
    BCAI.enhance_atoms(atoms, structure_dict)

    bonds = BCAI.make_bonds_df(mols)
    BCAI.enhance_bonds(bonds, structure_dict, flag='3JHH')

    triplets = BCAI.make_triplets(bonds["molecule_name"].unique(), structure_dict)
    #####

    embeddings, atoms, bonds, triplets = BCAI.add_embedding(atoms, bonds, triplets)


    # Not sure what to check . . .


def test_get_scaling():

    #####
    mols = get_mols()

    atoms = BCAI.make_atom_df(mols)
    structure_dict = BCAI.make_structure_dict(atoms)
    BCAI.enhance_structure_dict(structure_dict)
    BCAI.enhance_atoms(atoms, structure_dict)

    bonds = BCAI.make_bonds_df(mols)
    BCAI.enhance_bonds(bonds, structure_dict, flag='3JHH')

    #####

    means, stds = BCAI.get_scaling(bonds)


















###
