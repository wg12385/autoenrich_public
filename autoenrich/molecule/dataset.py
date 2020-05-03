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


from .nmrmol import nmrmol
from autoenrich.util.file_gettype import get_type
from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target
from autoenrich.util.filename_utils import get_unique_part
from autoenrich.ml.features import GNR_features as GNR

from tqdm import tqdm
import gzip
import pickle

import pandas as pd

import numpy as np
np.set_printoptions(threshold=99999999999)

class dataset(object):

	def __init__(self):

		self.mols = []
		self.atoms = pd.DataFrame()
		self.bonds = pd.DataFrame()
		self.struc = pd.DataFrame()
		self.BCAI = []
		self.mol_order = []
		self.big_data = False
		self.files = []

		self.params = {}


	def get_mols(self, files, type='', label_part=-1, fallback=False):
		self.mols = []
		self.files = files

		if len(files) > 2000:
			self.big_data = True

		if label_part == -1:
			label_part = get_unique_part(files)
			if label_part == -1:
				fallback = True

		for f, file in enumerate(files):
			if fallback:
				id = str(f)
			else:
				id = file.split('/')[-1].split('.')[0].split('_')[label_part]

			try:
				int(id)
				id = 'nmrmol' + str(id)
			except:
				pass

			if self.big_data:
				self.mols.append([file, id, type])
			else:
				mol = nmrmol(molid=id)

				if type == '':
					ftype = get_type(file)
				else:
					ftype = type
				mol.read_nmr(file, ftype)
				self.mols.append(mol)


	def get_features_frommols(self, args, params={}, molcheck_run=False, training=True, max=200):

		self.params = params

		target = flag_to_target(args['targetflag'])
		self.remove_mols(target)
		if molcheck_run:
			return


		for mol in self.mols:
			if len(mol.types) > max:
				max = len(mol.types) + 1
				print('WARNING, setting max atoms to ', max)

		self.params['max'] = max

		self.atoms = GNR.make_atom_df(self.mols)
		self.struc = GNR.make_struc_dict(self.atoms)
		if len(args['targetflag']) == 4:
			self.bonds = GNR.make_bonds_df(self.mols)

		if args['featureflag'] in ['aSLATM', 'CMAT', 'FCHL', 'ACSF']:
			from autoenrich.ml.features import QML_features
			self.atoms = QML_features.get_atomic_qml_features(self.atoms, self.bonds,
															self.struc, featureflag=args['featureflag'],
															cutoff=params['cutoff'], max=max)

		elif args['featureflag'] == 'BCAI':
			from autoenrich.ml.features import TFM_features
			self.BCAI, self.atoms, self.bonds, self.struc, xfiles, rfiles, yfiles = TFM_features.get_BCAI_features(self.atoms,
															self.bonds, self.struc, targetflag=args['targetflag'], training=training)

		elif args['featureflag'] != 'dummy':
			return

		else:
			print('Feature flag not recognised, no feature flag: ', args['featureflag'])
			return 0


	def assign_from_ml(self, pred_y, var, zero=True):
		assert len(self.r) > 0, print('Something went wrong, nothing to assign')

		try:
			with gzip.open(self.r[0], "rb") as f:
				self.r = pickle.load(f)
		except Exception as e:
			print(e)
			pass

		for molrf in tqdm(self.mols, desc='Assigning predictions'):

			if type(molrf) == list:
					mol = nmrmol(molid=molrf[1])

					if molrf[2] == '':
						ftype = get_type(molrf[2])
					else:
						ftype = molrf[2]

					mol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			if zero:
				mol.coupling = np.zeros((len(mol.types), len(mol.types)), dtype=np.float64)
				mol.shift = np.zeros((len(mol.types)), dtype=np.float64)

			for m, refs in enumerate(self.r):
				for r, ref in enumerate(refs):
					if mol.molid != ref[0]:
						continue

					if len(ref) == 2:
						for t in range(len(mol.types)):
							if ref[1] == t:
								mol.shift[t] = pred_y[r]
								mol.shift_var[t] = var[r]

					elif len(ref) == 3:
						for t1 in range(len(mol.types)):
							for t2 in range(len(mol.types)):
								if ref[1] == t1 and ref[2] == t2:
									mol.coupling[t1][t2] = pred_y[r]
									mol.coupling_var[t1][t2] = var[r]

	def remove_mols(self, target, progress=False):
		## Discard useless molecules:
		to_remove = []

		if progress:
			pbar = tqdm(self.mols)
		else:
			pbar = self.mols

		for molrf in pbar:
			if self.big_data:
				mol = nmrmol(molid=molrf[1])

				if molrf[2] == '':
					ftype = get_type(molrf[0])
				else:
					ftype = molrf[2]

				mol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			if len(target) == 1:
				if not target in mol.types:
					to_remove.append(mol.molid)
			elif len(target) == 3:
				found = False
				for i, it in enumerate(mol.types):
					for j, jt in enumerate(mol.types):
						if i >= j:
							continue

						if mol.coupling_len[i][j] == target[0]:
							if it == target[1] and jt == target[2]:
								found = True
							elif it == target[2] and jt == target[1]:
								found = True

				if not found:
					to_remove.append(mol.molid)

		if len(to_remove) > 1:
			keep = []
			for i in range(len(self.mols)):
				if self.big_data:
					if self.mols[i][1] in to_remove:
						continue
				else:
					if self.mols[i].molid in to_remove:
						continue
				keep.append(self.mols[i])

			self.mols = keep

			#print('REMOVED ', len(to_remove), '/', len(to_remove)+len(keep), ' molecules due to lack of features')
			#print(to_remove[:10])
			assert len(keep) > 1
