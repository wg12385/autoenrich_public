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

from sklearn.model_selection import KFold
import pickle
import numpy as np
import gzip
import copy
import time

def setup_logfile(args):
	strings = []
	strings.append('HPS ' + args['searchflag'] + ' SEARCH')
	strings.append(args['modelflag'] + '   ' + args['featureflag'] + '   ' + args['targetflag'])
	for param in args['param_list']:
		if args['param_logs'][param] == 'no':
			continue
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}'.format(param=param,
																	low=args['param_ranges'][param][0],
																	high=args['param_ranges'][param][1]))
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in args['param_list']:
		if args['param_logs'][param] == 'no':
			continue
		string = string + '\t{param:<20s}'.format(param=param)
	string = string + '\t{time:<10s}'.format(time='Mins')
	strings.append(string)

	if args['logfile'] == "":
		args['logfile'] = args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + args['searchflag'] + '.log'

	with open(args['logfile'], 'w') as f:
		for string in strings:
			print(string, file=f)

def update_logfile(args, iter, score, next_point_to_probe, time_taken):

	with open(args['logfile'], 'a') as f:
		string = '{i:<10d}\t{score:<10.5f}'.format(i=iter, score=score)
		for param in args['param_list']:
			if args['param_logs'][param] == 'no':
				continue
			string = string + '\t{param:<20.4g}'.format(param=next_point_to_probe[param])

		string = string + '\t{time:<10.4f}'.format(time=time_taken)
		print(string, file=f)

def HPS_iteration(model, iter, dataset, args, next_point_to_probe={}, BEST_SCORE=100000.00, BEST_PARAMS={}):

	#print('HPS_iteration. . .')
	time0 = time.time()

	args['max_size'] = 200
	args['load_dataset'] = 'false'

	if args['featureflag'] != 'BCAI' and args['feature_optimisation'] == 'True':
		dataset.get_features_frommols(args, params=next_point_to_probe)

		#print('Number of molecules in training set: ', len(dataset.mols))



	model.reset()
	model.get_x(dataset, args['targetflag'], assign_train=True)
	model.params = next_point_to_probe
	#print('Number of envs in training set: ', len(model.train_x[0]))
	y_pred = model.cv_predict()
	print(y_pred, file=open('preds.txt', 'w'))
	print(model.train_y, file=open('train_y.txt', 'w'))
	score = np.mean(np.absolute(y_pred - model.train_y))

	if score > 99999.99 or np.isnan(score):
		score = 99999.99

	time_taken = (time0 - time.time())/60

	update_logfile(args, iter, score, next_point_to_probe, time_taken)

	if score < BEST_SCORE:
		BEST_SCORE = score
		BEST_PARAMS = next_point_to_probe
		outname = args['output_dir'] +  args['modelflag'] + '_' + args['targetflag'] + '_' + args['featureflag'] + '_' + args['searchflag'] + '_model.pkl'
		model.save_model(outname)
		#print('BEST_SCORE = ', score)

	if BEST_PARAMS == {}:
		BEST_PARAMS = next_point_to_probe

	return score, BEST_SCORE, BEST_PARAMS

'''
def save_models(dataset, BEST_PARAMS, args):

	# create model
	if args['modelflag'] == 'KRR':
		from autoenrich.ml.models import KRRmodel
		model = KRRmodel.KRRmodel(id, np.asarray(dataset.x), np.asarray(dataset.y), params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'FCHL':
		from autoenrich.ml.models import FCHLmodel
		model = FCHLmodel.FCHLmodel(id, np.asarray(dataset.x), np.asarray(dataset.y), params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'NN':
		from autoenrich.ml.models import NNmodel
		model = NNmodel.NNmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'TFM':
		from autoenrich.ml.models import TFMmodel
		model = TFMmodel.TFMmodel(id, dataset.x, dataset.y, dataset.r, params=BEST_PARAMS, model_args=args)

	model.train()

	outname = args['output_dir'] +  args['modelflag'] + '_' + args['targetflag'] + '_' + args['featureflag'] + '_' + args['searchflag'] + '_model.pkl'

	pickle.dump(model, open(outname, "wb"))

	if args['modelflag'] != 'TFM':
		kf = KFold(n_splits=args['cv_steps'])
		kf.get_n_splits(dataset.x)

		i = 0
		for train_index, _ in kf.split(dataset.x):
			i += 1

			tmp_model = copy.deepcopy(model)

			tmp_dataset.x = np.asarray(dataset.x)[train_index]
			tmp_model.train_y = np.asarray(model.train_y)[train_index]
			tmp_model.train()

			outfile = args['output_dir'] + outname.split('/')[-1].split('.')[0] + '_' + str(i) + '.pkl'
			pickle.dump(tmp_model, open(outfile, "wb"))

	else:
		import torch
		molnames = list(set([dataset.r[i][0] for i in range(len(dataset.r))]))

		kf = KFold(n_splits=args['cv_steps'])
		kf.get_n_splits(dataset.x)

		i = 0
		for train_index, _ in kf.split(dataset.x[0]):
			i += 1
			train_x_list = []
			train_y_list = []
			for ii in range(len(dataset.x)):
				train_x_list.append(torch.index_select(dataset.x[ii], 0, torch.tensor(train_index)))
			for r, ref in enumerate(dataset.r):
				if ref[0] in [molnames[idx] for idx in train_index]:
					train_y_list.append(dataset.y[r])

			model.train(train_x=train_x_list, train_y=train_y_list)

			outfile = args['output_dir'] + outname.split('/')[-1].split('.')[0] + '_' + str(i) + '.pkl'
			pickle.dump(model, open(outfile, "wb"))

	return outname

'''










#
