# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

from autoenrich.file_creation import HPC_submission as HPCsub
import autoenrich.file_read.orca_read as orcaread
import autoenrich.file_read.g09_read as g09read
import autoenrich.file_read.g16_read as g16read
import autoenrich.file_creation.orca_submission as orcasub
import autoenrich.file_creation.g09_submission as g09sub
import autoenrich.file_creation.g16_submission as g16sub
import autoenrich.file_creation.structure_formats.nmredata as nmredata
import sys
import os.path

def setup_nmr(molecule, prefs, path='', ids=[], max=50):

	nmr_files = []
	for conformer in molecule.conformers:
		if (conformer.nmr_status == 'None' or conformer.nmr_status == 'pre-submission') and conformer.opt_status == 'successful':
			if prefs['NMR']['software'] == 'orca':
				conformer.nmr_in = orcasub.make_nmrin(prefs, conformer.molid, conformer.xyz, conformer.types,
															path + 'NMR/')
			elif prefs['NMR']['software'] == 'g09':
				conformer.nmr_in = g09sub.make_nmrcom(prefs, conformer.molid, conformer.xyz, conformer.types,
															path + 'NMR/')
			elif prefs['NMR']['software'] == 'g16':
				conformer.nmr_in = g16sub.make_nmrcom(prefs, conformer.molid, conformer.xyz, conformer.types,
															path + 'NMR/')

			conformer.nmr_log = conformer.nmr_in.split('.')[0] + '.log'
			conformer.nmr_status = 'pre-submission'
			nmr_files.append('NMR/' + conformer.nmr_in.split('/')[-1])

	if len(nmr_files) == 0:
		print('No files to submit. . .')
		sys.exit(0)

	IN_ARRAY = 'NMR/NMR_IN_ARRAY.txt'
	with open(path + IN_ARRAY, 'w') as f:
		for file in nmr_files:
			print(file, file=f)

	system = prefs['comp']['system']
	memory = prefs['NMR']['memory']
	processors = prefs['NMR']['processors']
	walltime = prefs['NMR']['walltime']
	software = prefs['NMR']['software']

	files = len(nmr_files)
	chunks = HPCsub.get_chunks(files)
	qsub_names = []
	for ck in range(chunks):
		start = (ck * max) + 1
		end = ((ck + 1) * max)
		if end > files:
			end = files


		#header = HPCsub.make_HPC_header(jobname=jobname, system=system, nodes=1, ppn=processors, walltime=walltime, mem=memory)
		jobname = 'aE_' + molecule.molid + '_' + str(ck) + '_NMR'
		strings = HPCsub.make_HPC_batch_submission(prefs, molecule.molid, IN_ARRAY, start, end, software=software,
									jobname=jobname, nodes=1, ppn=processors, walltime=walltime, mem=memory)

		if prefs['comp']['system'] == 'PBS':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.qsub'
		elif prefs['comp']['system'] == 'slurm':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.slurm'
		elif prefs['comp']['system'] == 'local':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.sh'
		with open(filename, 'w') as f:
			for string in strings:
				print(string, file=f)
		qsub_names.append(filename)

	print('Created ', chunks, ' submission files. . .')
	if prefs['comp']['system'] == 'PBS':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('qsub ', file)
	elif prefs['comp']['system'] == 'PBS':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('bash ', file)
	elif prefs['comp']['system'] == 'slurm':
		print('Havent finished this yet, good luck pal. . . .')


def process_nmr(molecule, prefs, path=''):

	good = 0
	bad = 0
	process = False

	for conformer in molecule.conformers:
		if os.path.isfile(conformer.nmr_log):
			if prefs['NMR']['software'] == 'orca':
				status = orcaread.get_nmr_status(conformer.nmr_log)
			elif prefs['NMR']['software'] == 'g09':
				status = g09read.get_nmr_status(conformer.nmr_log)
			elif prefs['NMR']['software'] == 'g16':
				status = g16read.get_nmr_status(conformer.nmr_log)

			if status == 'successful':
				good +=1
				conformer.read_nmr(conformer.nmr_log, type='orca')
			else:
				bad += 1
		else:
			bad += 1
			status = 'pre-submission'

		if conformer.nmr_status != status and status == 'successful':
			process = True
			outname = path + 'OUTPUT/' + molecule.molid + '_' + conformer.molid + '.nmredata.sdf'
			nmredata.nmrmol_to_nmredata(conformer, outname)

		conformer.nmr_status = status
		string = 'Conformer {molid:^3s} status: {status:^10s}'.format(molid=str(conformer.molid),
																								status=conformer.nmr_status)
		print(string)

	if process:
		molecule.boltzmann_average()
		outname = path + 'OUTPUT/' + molecule.molid + '.nmredata.sdf'
		nmredata.nmrmol_to_nmredata(molecule, outname)

	#print(good, ' successful NMR calculations, ', bad, ' failed, out of ', len(status))

	print('NMR Processing completed.')














###
