




from autoenrich.ml.models.model import genericmodel
from autoenrich.ml.models.KRRmodel import KRRmodel
from autoenrich.ml.models.FCHLmodel import FCHLmodel

from autoenrich.molecule.dataset import dataset

from test_generators.dummy_mol import get_random_ethane
from test_generators.dummy_dataset import get_test_dataset

from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target


import numpy as np

import os



def test_trainFCHLmodel():

    dset = get_test_dataset()

    for target in ['HCS', '1JCH']:
        args = {'featureflag': 'FCHL',
                'targetflag': target}
        params = {'cutoff': 5.0}

        dset.get_features_frommols(args, params)
        model = FCHLmodel()
        model.get_x(dset, args['targetflag'], assign_train=True)

        model.train()

        assert model.trained == True


def test_predictFCHLmodel():

    dset = get_test_dataset()

    for target in ['HCS', '1JCH']:
        args = {'featureflag': 'FCHL',
                'targetflag': target}
        params = {'cutoff': 5.0}

        dset.get_features_frommols(args, params)
        model = FCHLmodel()

        target = flag_to_target(args['targetflag'])

        model.get_x(dset, args['targetflag'], assign_train=True)
        model.train()

        assert model.trained == True
        train_x, train_y = model.get_x(dset, args['targetflag'], assign_train=False)
        test_x = np.asarray(train_x[:, :4])
        y_pred = model.predict(test_x)

        assert np.allclose(y_pred, train_y[:4], atol=0.1, rtol=0.0)
        assert sum(y_pred) != 0.0


def test_cvpredict_FCHLmodel_shift():

    dset = get_test_dataset(size=6)

    for target in ['HCS', '1JCH']:
        args = {'featureflag': 'FCHL',
                'targetflag': target}
        params = {'cutoff': 5.0}

        dset.get_features_frommols(args, params)
        model = KRRmodel()
        model.get_x(dset, args['targetflag'], assign_train=True)

        model.train()

        assert model.trained == True

        y_pred = model.cv_predict(fold=3)

        assert y_pred.shape == model.train_y.shape

def test_save_model():
    dset = get_test_dataset()
    args = {'featureflag': 'FCHL',
                'targetflag': 'HCS'}
    params = {'cutoff': 5.0}

    dset.get_features_frommols(args, params)
    model = FCHLmodel()
    model.get_x(dset, args['targetflag'], assign_train=True)
    model.train()

    filename='tests/test_tmp/tmp_model.pkl'
    model.save_model(filename=filename)

    assert os.path.isfile(filename)

    ld_model = FCHLmodel()
    ld_model.load_model(filename)

    train_x, train_y = model.get_x(dset, args['targetflag'], assign_train=False)
    test_x = np.asarray(train_x[:, :4])
    y_pred1 = ld_model.predict(test_x)
    y_pred2 = model.predict(test_x)

    assert np.array_equal(y_pred1, y_pred2)

    os.remove(filename)














#
