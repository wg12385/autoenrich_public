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
import glob

from autoenrich.molecule.dataset import dataset
from autoenrich.file_creation.structure_formats.nmredata import nmrmol_to_nmredata

from test_generators.dummy_mol import get_random_ethane, get_random_mol, mol_is_ethane
from test_generators.dummy_dataset import get_test_dataset



def test_dataset():

    dset = dataset()

    assert dset.mols == []
    assert dset.atoms.empty
    assert dset.bonds.empty
    assert dset.struc.empty
    assert dset.mol_order == []
    assert dset.big_data == False
    assert dset.files == []
    assert dset.params == {}

def test_get_mols():

    files = glob.glob('tests/test_store/dataset/eth_rnd_mol_*.nmredata.sdf')
    dset = dataset()
    dset.get_mols(files)

    assert dset.files == files
    assert dset.big_data == False
    assert len(dset.mols) == 5
    for mol in dset.mols:
        assert mol_is_ethane(mol)


def test_get_featuresCMAT():

    dset = get_test_dataset(size = 1)
    args = {'featureflag': 'CMAT',
            'targetflag': 'HCS'}
    params = {'cutoff': 5.0}

    dset.get_features_frommols(args, params=params)

    assert len(dset.atoms['atomic_rep'].values[0]) == 50*(50+1)/2


def test_get_featuresaSLATM():

    dset = get_test_dataset(size = 1)
    args = {'featureflag': 'aSLATM',
            'targetflag': 'HCS'}
    params = {'cutoff': 5.0}

    dset.get_features_frommols(args, params=params)

    assert len(dset.atoms['atomic_rep'].values[0]) == 20105


def test_get_featuresFCHL():

    dset = get_test_dataset(size = 1)
    args = {'featureflag': 'FCHL',
            'targetflag': 'HCS'}
    params = {'cutoff': 5.0}

    dset.get_features_frommols(args, params=params)

    assert dset.atoms['atomic_rep'].values[0].shape == (5, 50)






#
