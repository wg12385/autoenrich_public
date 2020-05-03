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

import numpy as np
from autoenrich.ml.models.model import genericmodel
import qml
from sklearn.model_selection import KFold
import copy
import sys
import pickle

class KRRmodel(genericmodel):

	def __init__(self, id='KRRmodel', x=[], y=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)
		self.alpha = []

		if 'sigma' not in self.params.keys():
			self.params['sigma'] = 0.1

		if 'lamda' not in self.params.keys():
			self.params['lamda'] = 0.0001

	def reset(self):
		genericmodel.__init__(self, self.id, [], [], {}, {})
		self.alpha = []

		if 'sigma' not in self.params.keys():
			self.params['sigma'] = 0.1

		if 'lamda' not in self.params.keys():
			self.params['lamda'] = 0.0001



	def train(self, train_x=[], train_y=[]):

		# train_x must be of shape (i, j, k)
		# where i is number of atoms in each 'enrivonment'
		# 		j is the number of atoms (training cases)
		#		k is the size of the atomic representation

		k = []

		if len(train_x) == 0:
			train_x = self.train_x
		if len(train_y) == 0:
			train_y = self.train_y

		dimensions = train_x.shape[0]

		Xr = train_x


		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			k.append(qml.kernels.laplacian_kernel(Xr[d], Xr[d], self.params['sigma']) + self.params['lamda'] * np.identity(Xr[d].shape[0]))

		# loop through kernels
		K = k[0]
		if len(k) > 1:
			for i in range(1, len(k)):
				# multiply kernels so result is k1 * k2 * k3 ...
				K = K * k[i]

		# get the KRR prefactors, i.e. this is the training of the network
		self.alpha = qml.math.cho_solve(K, train_y)

		# report training state
		self.trained = True


	def predict(self, test_x, train_x=[]):

		if len(train_x) == 0:
			train_x = self.train_x

		ks = []
		assert test_x.shape[0] == train_x.shape[0], print(test_x.shape, train_x.shape)
		assert test_x.shape[2] == train_x.shape[2], print(test_x.shape, train_x.shape)

		dimensions = train_x.shape[0]
		Xe = test_x
		Xr = train_x
		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			ks.append(qml.kernels.laplacian_kernel(Xe[d], Xr[d], self.params['sigma']))

		# loop through kernels
		Ks = ks[0]
		if len(ks) > 1:
			for i in range(1, len(ks)):
				# multiply kernels so result is k1 * k2 * k3 ...
				Ks = Ks * ks[i]

		# predict values of y
		y_pred = np.dot(Ks, self.alpha)

		return y_pred


	def get_x(self, dataset, targetflag, assign_train=False):

		if len(targetflag) == 3:
			x = dataset.atoms.loc[(dataset.atoms['typestr'] == targetflag[0])]['atomic_rep'].to_numpy()
			train_x = []
			for rep in x:
				train_x.append(rep)

			train_x = np.asarray([train_x])
			train_y = dataset.atoms.loc[(dataset.atoms['typestr'] == targetflag[0])]['shift'].to_numpy()


		elif len(targetflag) == 4:
			cpl_idx = dataset.bonds.loc[(dataset.bonds['type'] == targetflag)][['molecule_name', 'atom_index_0', 'atom_index_1']].to_numpy()
			x1 = []
			x2 = []
			for cpl in cpl_idx:
				x1.append(dataset.atoms.loc[(dataset.atoms['molecule_name'] == cpl[0])
									& (dataset.atoms['atom_index'] == cpl[1])]['atomic_rep'].to_numpy()[0])
				x2.append(dataset.atoms.loc[(dataset.atoms['molecule_name'] == cpl[0])
									& (dataset.atoms['atom_index'] == cpl[2])]['atomic_rep'].to_numpy()[0])
			train_x = [x1, x2]

			train_x = np.asarray(train_x)
			train_y = np.squeeze(dataset.bonds.loc[(dataset.bonds['type'] == targetflag)][['scalar_coupling_constant']].to_numpy())


		assert train_x.shape[1] == train_y.shape[0], print(train_x.shape, train_y.shape)

		if assign_train:
			self.train_x = train_x
			self.train_y = train_y
		else:
			return train_x, train_y

	def save_model(self, filename='model.pkl'):

		model_store = {
				"alpha": self.alpha,
				"train_x": self.train_x,
				"train_y": self.train_y,
				"params": self.params,
				"id": self.id
		}

		pickle.dump(model_store, open(filename, "wb"))

	def load_model(self, filename):

		model_store = pickle.load(open(filename, "rb"))
		self.alpha = model_store["alpha"]
		self.train_x = model_store["train_x"]
		self.train_y = model_store["train_y"]
		self.params = model_store["params"]
		self.id = model_store["id"]


	def cv_predict(self, fold=5):

		kf = KFold(n_splits=fold)
		kf.get_n_splits(self.train_x)
		pred_y = []

		for train_index, test_index in kf.split(self.train_x[0]):

			self.train(train_x=self.train_x[:, train_index], train_y=self.train_y[train_index])

			pred_y.extend(self.predict(self.train_x[:, test_index], train_x=self.train_x[:, train_index]))

		pred_y = np.asarray(pred_y)

		return pred_y



##
