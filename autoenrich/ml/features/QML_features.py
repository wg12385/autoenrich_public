import qml
import numpy as np
import sys


def get_atomic_qml_features(atoms, bonds, struc, featureflag='CMAT', cutoff=5.0, max=50):

    atoms['atomic_rep'] = np.empty((len(atoms), 0)).tolist()

    for molid in atoms['molecule_name'].unique():

        xyz = struc[molid]['positions']
        types = struc[molid]['typesint']

        if featureflag == 'aSLATM':
            mbtypes = [[1],[1,1], [1,1,1], [1,1,6], [1,1,7], [1,1,8], [1,1,9], [1,6], [1,6,1], [1,6,6], [1,6,7], [1,6,8], [1,6,9], [1,7], [1,7,1], [1,7,6], [1,7,7], [1,7,8], [1,7,9], [1,8], [1,8,1], [1,8,6], [1,8,7], [1,8,8], [1,8,9], [1,9], [1,9,1], [1,9,6], [1,9,7], [1,9,8], [1,9,9], [6], [6,1], [6,1,1], [6,1,6], [6,1,7], [6,1,8], [6,1,9], [6,6], [6,6,1], [6,6,6], [6,6,7], [6,6,8], [6,6,9], [6,7], [6,7,1], [6,7,6], [6,7,7], [6,7,8], [6,7,9], [6,8], [6,8,1], [6,8,6], [6,8,7], [6,8,8], [6,8,9], [6,9], [6,9,1], [6,9,6], [6,9,7], [6,9,8], [6,9,9], [7],[7,1], [7,1,1], [7,1,6], [7,1,7], [7,1,8], [7,1,9], [7,6], [7,6,1], [7,6,6], [7,6,7], [7,6,8], [7,6,9], [7,7], [7,7,1], [7,7,6], [7,7,7], [7,7,8], [7,7,9], [7,8], [7,8,1], [7,8,6], [7,8,7], [7,8,8], [7,8,9], [7,9], [7,9,1], [7,9,6], [7,9,7], [7,9,8], [7,9,9], [8], [8,1], [8,1,1], [8,1,6], [8,1,7], [8,1,8], [8,1,9], [8,6], [8,6,1], [8,6,6], [8,6,7], [8,6,8], [8,6,9], [8,7], [8,7,1], [8,7,6], [8,7,7], [8,7,8], [8,7,9], [8,8], [8,8,1], [8,8,6], [8,8,7], [8,8,8], [8,8,9], [8,9], [8,9,1], [8,9,6], [8,9,7], [8,9,8], [8,9,9], [9], [9,1], [9,1,1], [9,1,6], [9,1,7], [9,1,8], [9,1,9], [9,6], [9,6,1], [9,6,6], [9,6,7], [9,6,8], [9,6,9], [9,7], [9,7,1], [9,7,6], [9,7,7], [9,7,8], [9,7,9], [9,8], [9,8,1], [9,8,6], [9,8,7], [9,8,8], [9,8,9], [9,9], [9,9,1], [9,9,6], [9,9,7], [9,9,8], [9,9,9]]
            '''
            nuclear_charges = []
            for molid in atoms['molecule_name'].unique():
                nuclear_charges.append(struc[molid]['typesint'])
            mbtypes = qml.representations.get_slatm_mbtypes(nuclear_charges)
            '''
            reps = qml.representations.generate_slatm(np.asarray(xyz), np.asarray(types), mbtypes, rcut=cutoff, local=True)

        elif featureflag == 'CMAT':
            reps = qml.representations.generate_atomic_coulomb_matrix(np.asarray(types), xyz, size=max, central_cutoff = cutoff)

        elif featureflag == 'FCHL':
            reps = qml.fchl.generate_representation(np.asarray(xyz), np.asarray(types), max_size=max, cut_distance=cutoff)

        elif featureflag == 'ACSF':
            try:
                reps = qml.representations.generate_acsf(types, xyz, elements=[1, 6, 7, 8, 9, 14, 15, 16, 17, 35],
                                                        nRs2=3, nRs3=3,
                                                        nTs=2, eta2=1, eta3=1, zeta=1, rcut=cutoff, acut=cutoff,
                                                        bin_min=0.1, gradients=False)
            except AttributeError as e:
                print(e)
                print('No ACSF function available')
                return 0

        else:
            print('feature not found !!! not such feature: ', featureflag)
            return

        for i in range(len(types)):
            idx = atoms.loc[(atoms['molecule_name'] == molid)
                            & (atoms['atom_index'] == i)].index[0]
            atoms.at[idx, 'atomic_rep'] = np.asarray(reps[i])

    return atoms












###
