


from autoenrich.ml.features import QML_features as QML
from autoenrich.ml.features import GNR_features as GNR

from test_generators.dummy_dataset import get_test_dataset
import test_generators.dummy_mol as dmy

def test_get_atomic_cmat():

    mols = dmy.get_rndethane_mols()
    atoms = GNR.make_atom_df(mols)
    struc = GNR.make_struc_dict(atoms)
    bonds = GNR.make_bonds_df(mols)

    atoms = QML.get_atomic_qml_features(atoms, bonds, struc, featureflag='CMAT')

    assert len(atoms['atomic_rep'].values[0]) == 50*(50+1)/2
    assert len(atoms['atomic_rep'].values[0].nonzero()[0]) == len(mols[0].types)*(len(mols[0].types)+1)/2


def test_get_atomic_aSLATM():

    mols = dmy.get_rndethane_mols()
    atoms = GNR.make_atom_df(mols)
    struc = GNR.make_struc_dict(atoms)
    bonds = GNR.make_bonds_df(mols)

    atoms = QML.get_atomic_qml_features(atoms, bonds, struc, featureflag='aSLATM')

    assert len(atoms['atomic_rep'].values[0]) == 20105

def test_get_atomic_FCHL():

    mols = dmy.get_rndethane_mols()
    atoms = GNR.make_atom_df(mols)
    struc = GNR.make_struc_dict(atoms)
    bonds = GNR.make_bonds_df(mols)

    atoms = QML.get_atomic_qml_features(atoms, bonds, struc, featureflag='FCHL')

    assert atoms['atomic_rep'].values[0].shape == (5, 50)

def test_get_atomic_ACSF():

    mols = dmy.get_rndethane_mols()
    atoms = GNR.make_atom_df(mols)
    struc = GNR.make_struc_dict(atoms)
    bonds = GNR.make_bonds_df(mols)

    atoms = QML.get_atomic_qml_features(atoms, bonds, struc, featureflag='ACSF')
    






#
