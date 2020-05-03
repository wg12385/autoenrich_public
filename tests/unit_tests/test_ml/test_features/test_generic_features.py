# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.
import sys
import os

import numpy as np

from autoenrich.reference.periodic_table import Get_periodic_table
import autoenrich.ml.features.GNR_features as GNR

import test_generators.dummy_mol as dmy


def test_get_dummy_features():

    mol = dmy.get_ethane_mol()

    mols = [mol]

    x, y, r = GNR.get_dummy_features(mols, targetflag='CCS')
    assert x == []
    assert np.array_equal(y, mol.shift[np.where(mol.types==6)])
    assert np.array_equal(r, [['ethane', 0], ['ethane', 1]])

    x, y, r = GNR.get_dummy_features(mols, targetflag='HCS')
    assert x == []
    assert np.array_equal(y, mol.shift[np.where(mol.types==1)])
    assert np.array_equal(r, [['ethane', 2], ['ethane', 3], ['ethane', 4],
                            ['ethane', 5], ['ethane', 6], ['ethane', 7]])

    x, y, r = GNR.get_dummy_features(mols, targetflag='1JCH')
    for i, index in enumerate(r):
        assert y[i] == mol.coupling[index[1]][index[2]]
    assert np.array_equal(r, [['ethane', 0, 2], ['ethane', 0, 3], ['ethane', 0, 4],
                            ['ethane', 1, 5], ['ethane', 1, 6], ['ethane', 1, 7]])

    x, y, r = GNR.get_dummy_features(mols, targetflag='3JHH')
    for i, index in enumerate(r):
        assert y[i] == mol.coupling[index[1]][index[2]]
    assert np.array_equal(r, [['ethane', 2, 5], ['ethane', 2, 6], ['ethane', 2, 7],
                            ['ethane', 3, 5], ['ethane', 3, 6], ['ethane', 3, 7],
                            ['ethane', 4, 5], ['ethane', 4, 6], ['ethane', 4, 7]])



def test_make_atom_df():

    mols = dmy.get_rndethane_mols()
    ats = 0
    for mol in mols:
        ats += len(mol.types)

    atoms = GNR.make_atom_df(mols)
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


def test_make_struc_dict():

    p_table = Get_periodic_table()

    mols = dmy.get_rndethane_mols()
    atoms = GNR.make_atom_df(mols)
    structure_dict = GNR.make_struc_dict(atoms)
    assert len(structure_dict.keys()) == len(mols)

    for mol in mols:
        assert structure_dict[mol.molid]['typesstr'] == [p_table[type] for type in mol.types]
        assert np.array_equal(structure_dict[mol.molid]['positions'], mol.xyz)
        assert np.array_equal(structure_dict[mol.molid]['conn'], mol.conn)

def test_make_bonds_df():

    p_table = Get_periodic_table()

    #####
    mols = dmy.get_rndethane_mols()
    #####

    bonds = GNR.make_bonds_df(mols)
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
