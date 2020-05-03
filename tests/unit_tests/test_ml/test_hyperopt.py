# somehow write testing for the hyper-parameter optimisation


import os

from autoenrich.ml.hyperparameter_tuning import HPS_grid, HPS_random, HPS_gaussian
from test_generators.dummy_dataset import get_test_dataset
from autoenrich.ml.models.KRRmodel import KRRmodel



def test_grid_hyperopt():

    dataset = get_test_dataset()
    args = {'featureflag': 'CMAT',
            'targetflag': 'HCS',
            'feature_optimisation': 'True',
            'logfile': 'tests/test_tmp/tmp_logfile',
            'searchflag': 'test_search',
            'modelflag': 'KRR',
            'param_list': ['cutoff', 'sigma', 'lamda'],
            'param_logs': {'cutoff': 'lin', 'sigma': 'log', 'lamda': 'log'},
            'param_ranges': {'cutoff': [0, 10], 'sigma': [-5, -2], 'lamda': [-10, -4]},
            'output_dir': './'}

    args['grid_density'] = 2


    model = KRRmodel()
    dataset, BEST_SCORE = HPS_grid.grid_search(model, dataset, args)


def test_random_hyperopt():

    dataset = get_test_dataset()
    args = {'featureflag': 'CMAT',
            'targetflag': 'HCS',
            'feature_optimisation': 'True',
            'logfile': 'tests/test_tmp/tmp_logfile',
            'searchflag': 'test_search',
            'modelflag': 'KRR',
            'param_list': ['cutoff', 'sigma', 'lamda'],
            'param_logs': {'cutoff': 'lin', 'sigma': 'log', 'lamda': 'log'},
            'param_ranges': {'cutoff': [0, 10], 'sigma': [-5, -2], 'lamda': [-10, -4]},
            'output_dir': './'}

    args['epochs'] = 3

    model = KRRmodel()
    dataset, BEST_SCORE = HPS_random.random_search(model, dataset, args)


def test_gaussian_hyperopt():

    dataset = get_test_dataset()
    args = {'featureflag': 'CMAT',
            'targetflag': 'HCS',
            'feature_optimisation': 'True',
            'logfile': 'tests/test_tmp/tmp_logfile',
            'searchflag': 'test_search',
            'modelflag': 'KRR',
            'param_list': ['cutoff', 'sigma', 'lamda'],
            'param_logs': {'cutoff': 'lin', 'sigma': 'log', 'lamda': 'log'},
            'param_ranges': {'cutoff': [0, 10], 'sigma': [-5, -2], 'lamda': [-10, -4]},
            'output_dir': './'}

    args['epochs'] = 3
    args['kappa'] = 5
    args['xi'] = 0.1
    args['load'] = False
    args['random'] = 0

    model = KRRmodel()
    dataset, BEST_SCORE = HPS_gaussian.gaussian_search(model, dataset, args)











##
