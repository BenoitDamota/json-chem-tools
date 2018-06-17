#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Functions for calculating molecular orbitals and electron density differences.

from physics import A_to_a0
from orbkit import grid, core
## Patched version of orbkit.read
import read
import numpy as np



## Get grid parameters and initialize grid
def _init_ORB_grid(data, grid_par=-6, over_s=7):
	u"""Initializes ORBKIT's grid and returns necessary data for visualization.
	**Parameters:**

	  data : dict
	An instance of ORBKIT's QC class.
	  grid_par : int
	Parameter controlling grid. If positive or zero, specifies number of points; else specifies the resolution in atomic units.
	  over_s : int|float
	Oversizing, to be tuned

	**Returns**
	  X, Y, Z 
	  out
	List of numpy.ndarrays of the same shape as X, Y, Z
	"""

	## Spacing/Number of points
	if grid_par > 0:
		grid.N_ = [grid_par]*3
	elif grid_par == 0:
		grid.N_ = [80]*3
	else:
		grid.delta_ = [1.0/(-grid_par)]*3

	grid.max_ = np.amax(data.geo_spec, axis=0) + over_s
	grid.min_ = np.amin(data.geo_spec, axis=0) - over_s

	grid.init()

	return np.mgrid[grid.min_[0]:grid.max_[0]:1j*len(grid.x),   # X
                        grid.min_[1]:grid.max_[1]:1j*len(grid.y),   # Y
                        grid.min_[2]:grid.max_[2]:1j*len(grid.z)]   # Z



## Calculations
def MO(j_data, MO_list, grid_par=-6):
	u"""Calculates the voxels representing the requested molecular orbitals of a molecule.
	** Parameters **
	  j_data : dict
	Data on the molecule, deserialized from the scanlog format.
	  MO_list : list(str) or list(int)
	A list containing either:
	  - integers designating the molecular orbitals (starting from 1)
	  - strings designating the molecular orbitals either through:
	    - keywords (homo, lumo)
	    - ranges (e.g. "homo-1:lumo+2")
	    - symmetries (e.g. "1.A2")
	  grid_par : int, optional
	Governs the grid to be used for the voxels.
	  - If positive, it indicates the number of points in all 3 dimensions.
	  - If zero, it indicates the default number of points (80).
	  - If negative, it indicates the resolution, in Á^-1.
	** Returns **
	  out: list
	List of sequences of scalar values at each voxel of the grid for each molecular orbital, of the same shape as the grid.
	  X, Y, Z
	Meshgrids (as generated by numpy.mgrid), required for placing the voxel values contained in out.
	"""

	## Prepare data for ORBKIT
	qc = read.convert_json(j_data, all_mo=True)
	X, Y, Z = _init_ORB_grid(qc, grid_par=grid_par)

	## Get list of orbitals
	qc.mo_spec = read.mo_select(qc.mo_spec, MO_list)["mo_spec"]

	## Calculate
	out = core.rho_compute(qc, calc_mo=True, numproc=4)

	return out, X, Y, Z



def EDD(j_data, transitions, grid_par=-6):
	u"""Calculates the voxels representing the density differences of the requested transitions of a molecule.
	** Parameters **
	  j_data : dict
	Data on the molecule, as deserialized from the scanlog format.
	  transitions : list
	A list of integers designating the transitions to consider (starting from 0)
	  grid_par : int, optional
	Governs the grid to be used for the voxels.
	  - If positive, it indicates the number of points in all 3 dimensions.
	  - If zero, it indicates the default number of points (80).
	  - If negative, it indicates the resolution, in Á^-1.
	** Returns ** 
	  out: list
	List of sequences of scalar values at each voxel of the grid for each molecular orbital, of the same shape as the grid.
	  X, Y, Z
	Meshgrids (as generated by numpy.mgrid), required for placing the voxel values contained in out.
	"""

	## Prepare data for ORBKIT
	qc = read.convert_json(j_data, all_mo=True)

	## To save time, the calculation is done in two phases:

	## 1. Get all MOs involved in transitions and calculate them once
	MO_list = [STT[0] for T in transitions for ST in T for STT in ST[:2]]
	MO_set = set(MO_list)
	## The dictionary is needed because the MO numbers do not correspond
	## to their indices in the output of rho_compute, and we can only access
	## them through indices (assuming ORBKIT conserves the order of the MO numbers
	## in qc.mo_spec)
	tab = dict([(MO, i) for i, MO in enumerate(MO_set)])

	qc.mo_spec = read.mo_select(qc.mo_spec, [str(MO + 1) for MO in MO_set])["mo_spec"]
	X, Y, Z = _init_ORB_grid(qc, grid_par=grid_par)
	MOs = core.rho_compute(qc, calc_mo=True, numproc=4)

	## 2. Combine MOs according to info in `et_transitions`
	out = []
	for i, T in enumerate(transitions):
		series = np.zeros(MOs[0].shape)
		for j, ST in enumerate(T):
			## Dp_i = S_j(C_ij**2*(MO2_ij**2 - MO1_ij**2))
			print "Calculating transition {}.{}".format(i, j)
			series += ST[2]**2*(np.square(MOs[tab[ST[1][0]]]) - np.square(MOs[tab[ST[0][0]]]))
		out += [series]

	return out, X, Y, Z
