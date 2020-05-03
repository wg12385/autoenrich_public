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

import numpy as np
import copy
from qml.math import cho_solve
from qml.kernels import laplacian_kernel
from sklearn.model_selection import KFold
from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target

import tracemalloc
# top level model class

class genericmodel(object):

	def __init__(self, id='genericmodel', train_x=[], train_y=[], params={},
									model_args={}):
		self.id = id

		self.train_x = train_x
		self.train_y = train_y

		self.params = params
		self.args = model_args

		self.trained = False

	def train(self):
		print('Stub function, why are you running this ??')

	def predict(self):
		print('Stub function, why are you running this ??')


	def cv_predict(self, fold):
		print('Stub function, why are you running this ??')


























##
