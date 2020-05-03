# somehow write testing for the hyper-parameter optimisation


import os

from autoenrich.ml.hyperparameter_tuning import HPS_generic
from test_generators.dummy_dataset import get_test_dataset
from autoenrich.ml.models.KRRmodel import KRRmodel


def test_setup_logfile():

    args = {}
    args['searchflag'] = 'gaussian'
    args['modelflag'] = 'FCHL'
    args['featureflag'] = 'FCHL'
    args['targetflag'] = 'HCS'
    args['param_list'] = ['cutoff', 'lambda', 'sigma']
    args['param_ranges'] = {'cutoff': [0, 1],
                            'lambda': [-1, -5],
                            'sigma': [-4, -10]}
    args['param_logs'] = {'cutoff': 'lin',
                            'lambda': 'log',
                            'sigma': 'log'}

    args['logfile'] = 'tests/test_tmp/tmp.log'
    HPS_generic.setup_logfile(args)
    assert os.path.isfile(args['logfile'])
    os.remove(args['logfile'])


def test_update_logfile():

    args = {}
    args['searchflag'] = 'gaussian'
    args['modelflag'] = 'FCHL'
    args['featureflag'] = 'FCHL'
    args['targetflag'] = 'HCS'
    args['param_list'] = ['cutoff', 'lamda', 'sigma']
    args['param_ranges'] = {'cutoff': [0, 1],
                            'lamda': [-1, -5],
                            'sigma': [-4, -10]}
    args['param_logs'] = {'cutoff': 'lin',
                            'lamda': 'log',
                            'sigma': 'log'}

    args['logfile'] = 'tests/test_tmp/tmp.log'
    HPS_generic.setup_logfile(args)
    assert os.path.isfile(args['logfile'])

    iter = 1
    score = 40.04
    next_point_to_probe = {'sigma': 0.1, 'lamda': 20, 'cutoff': 5}
    time_taken = 0.42
    HPS_generic.update_logfile(args, iter, score, next_point_to_probe, time_taken)
    os.remove(args['logfile'])


def test_iteration():

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

    params = {'cutoff': 5.0, 'sigma': 0.1, 'lamda': 0.01}

    model = KRRmodel()

    BEST_SCORE = 99999.9999
    BEST_PARAMS ={}

    HPS_generic.setup_logfile(args)
    for iter in range(2):
        next_point_to_probe = params
        score, BEST_SCORE, BEST_PARAMS = HPS_generic.HPS_iteration(model, iter, dataset, args,
                        next_point_to_probe=next_point_to_probe, BEST_SCORE=BEST_SCORE, BEST_PARAMS=BEST_PARAMS)

    print(score, BEST_SCORE, BEST_PARAMS)

    os.remove(args['logfile'])

















##
