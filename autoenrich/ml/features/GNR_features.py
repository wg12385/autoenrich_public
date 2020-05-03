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

from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target
from autoenrich.reference.periodic_table import Get_periodic_table
import numpy as np
import pandas as pd
from tqdm import tqdm

import collections

def make_atom_df(mols, progress=False):

	p_table = Get_periodic_table()

	# construct dataframe as BCAI requires from mols
	# atoms has: molecule_name, atom, labeled atom,
	molecule_name = [] 	# molecule name
	atom_index = []		# atom index
	typestr = []		# atom type (string)
	typeint = []		# atom type (integer)
	x = []				# x coordinate
	y = []				# y coordinate
	z = []				# z coordinate
	conns = []
	shifts = []

	mol_order = []
	m = -1

	if progress:
		pbar = tqdm(mols, desc='Constructing atom dictionary')
	else:
		pbar = mols

	for molrf in pbar:
		m += 1
		if len(mols) > 2000:
			mol = nmrmol(molid=molrf[1])

			if molrf[2] == '':
				ftype = get_type(molrf[2])
			else:
				ftype = molrf[2]
			mol.read_nmr(molrf[0], ftype)
		else:
			mol = molrf
		mol_order.append(mol.molid)
		for t, type in enumerate(mol.types):
			molecule_name.append(mol.molid)
			atom_index.append(t)
			typestr.append(p_table[type])
			typeint.append(type)
			x.append(mol.xyz[t][0])
			y.append(mol.xyz[t][1])
			z.append(mol.xyz[t][2])
			conns.append(mol.conn[t])
			shifts.append(mol.shift[t])

	atoms = {	'molecule_name': molecule_name,
				'atom_index': atom_index,
				'typestr': typestr,
				'typeint': typeint,
				'x': x,
				'y': y,
				'z': z,
				'conn': conns,
				'shift': shifts
			}
	atoms = pd.DataFrame(atoms)

	return atoms


def make_bonds_df(mols, progress=False):

		p_table = Get_periodic_table()

		# construct dataframe as BCAI requires from mols
		# atoms has: molecule_name, atom, labeled atom,
		id = []				# number
		molecule_name = [] 	# molecule name
		atom_index_0 = []	# atom index for atom 1
		atom_index_1 = []	# atom index for atom 2
		cpltype = []			# coupling type
		coupling = []	# coupling value
		r = []
		y = []

		i = -1
		m = -1

		if progress:
			pbar = tqdm(mols, desc='Constructing bonds dictionary')
		else:
			pbar = mols

		for molrf in pbar:
			m += 1
			if len(mols) > 2000:
				mol = nmrmol(molid=molrf[1])

				if molrf[2] == '':
					ftype = get_type(molrf[2])
				else:
					ftype = molrf[2]
				mol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			moly = []
			molr = []

			for t, type in enumerate(mol.types):
				for t2, type2 in enumerate(mol.types):
					if t == t2:
						continue

					TFM_flag = str(mol.coupling_len[t][t2]) + 'J' + p_table[type] + p_table[type2]

					#if TFM_flag != flag and flag != 'all':
					#	continue

					i += 1
					id.append(i)
					molecule_name.append(mol.molid)
					atom_index_0.append(t)
					atom_index_1.append(t2)


					cpltype.append(TFM_flag)

					coupling.append(mol.coupling[t][t2])

					moly.append(mol.coupling[t][t2])
					molr.append([mol.molid, t, t2])

			y.append(moly)
			r.append(molr)

		bonds = {	'id': id,
					'molecule_name': molecule_name,
					'atom_index_0': atom_index_0,
					'atom_index_1': atom_index_1,
					'type': cpltype,
					'scalar_coupling_constant': coupling
				}

		bonds = pd.DataFrame(bonds)

		return bonds

def make_struc_dict(atoms_dataframe):

	atoms = atoms_dataframe.sort_values(["molecule_name", "atom_index"])  # ensure ordering is consistent
	# Make a molecule-based dictionary of the information
	structure_dict = collections.defaultdict(lambda: {"typesstr":[], "typesint": [], "positions":[],"conn":[]})
	for index,row in atoms.iterrows():
		structure_dict[row["molecule_name"]]["typesstr"].append(row["typestr"])
		structure_dict[row["molecule_name"]]["typesint"].append(row["typeint"])
		structure_dict[row["molecule_name"]]["positions"].append([row["x"],row["y"],row["z"]])
		structure_dict[row["molecule_name"]]["conn"].append(row["conn"])

	return structure_dict

# Make 'dummy features' for a given target, used to get list of references
def get_dummy_features(mols, targetflag='CCS'):
	# Input:
	#	mols: list of autoenrich nmrmol objects
	#	targetflag: flag corresponding to nmr parameter (string)

	# Returns:
	#	x: empty list (where features would normally go)
	#	y: 1D list of NMR parameters
	#	r: 1D list of NMR parameter references (molid, atom(1, atom2))

	target = flag_to_target(targetflag)

	x = []
	y = []
	r = []

	for mol in mols:
		if len(target) == 1:
			for i in range(len(mol.types)):
				if mol.types[i] == target[0]:
					y.append(mol.shift[i])
					r.append([mol.molid, i])

		if len(target) == 3:
			for i in range(len(mol.types)):
				for j in range(i+1, len(mol.types)):

					if not ( mol.types[i] == target[1] and mol.types[j] == target[2] ):
						continue

					if mol.coupling_len[i][j] != target[0]:
						continue

					y.append(mol.coupling[i][j])
					r.append([mol.molid, i, j])

	return x, y, r
