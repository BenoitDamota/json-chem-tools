#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-

from sys import argv
from physics import A_to_a0
from calc_orb import *
from viz_orb_mayavi import *
import numpy as np
from sys import stderr, exit
from json import load, dump
import os

## Usage
usage = ( "Usage: {0} TOPO                         $npts ${{input}}.json                      ${{output}}\n"
	+ "     | {0} MO[=$value]                  $npts ${{input}}.json                      ${{output}}\n"
        + "     | {0} POTENTIAL                    $npts ${{input}}.json                      ${{output}}\n"
        + "     | {0} (TD|EDD|BARY|Tozer)[=$value] $npts ${{OPT-input}}.json ${{TD-input}}.json ${{output}}\n"
        + "     | {0} FUKUI                        $npts ${{OPT-input}}.json ${{SP-input}}.json ${{output}}\n"
        + "     | {0} DUAL                         $npts ${{OPT-input}}.json ${{SP-input}}.json ${{SP-input}}.json ${{output}}").format(argv[0])

## Get job type and parameters
val = argv[1].split("=")
job = val[0].upper()   ## .upper() added to not take the letter case into account

## Verifying usage

if len(argv) < 4:
    print('Not enough arguments. Check usage below :')
    stderr.write(usage)
    exit(0)

elif job not in {'TOPO', 'MO', 'POTENTIAL', 'TD', 'EDD', 'BARY', 'TOZER', 'FUKUI', 'DUAL'}:
    print(job,' is not an expected keyword for calculations.')
    stderr.write(usage)
    exit(1)

elif isinstance(int(argv[2]), int) != True:
    print('Numbers of points should be an integer (positive or negative)')
    stderr.write(usage)
    exit(1)

## Output file name
if ".json" in argv[-1]:
    (filepath, file_name) = os.path.split(os.getcwd())
    #file_name = None
else:
    file_name = argv[-1]

## Input
if ".json" not in argv[3]:
    print('A .json file should be given as an input')
    stderr.write(usage)
    exit(1)

with open(argv[3], "r") as f:
        data = load(f)

	## Code of uncertain status due to instability in specs
	## Kept for later
	#from cclib.parser import ccopen
	#cc_data = ccopen(argv[4]).parse()
	##cc_data.mocoeffs[0] = [[0.0 if np.log(np.abs(x)) < -9 else x for x in S] for S in cc_data.mocoeffs[0]]
	#qc = read.convert_cclib(cc_data, all_mo=True)

## Get grid parameters and initialize grid

## Oversizing in Bohr radii
## ORBKIT uses 5 by default, tune this as required
over_s = 7

## Spacing/Number of points
par = int(argv[2])

## Calculate
if job == "TOPO":
    ## Creating the figure in viz_orb_mayavi.py
        topo(data, file_name)
        print('Job done')

elif job == "MO":

    ## Get list of orbitals
	MO_list = val[1].split(",")
    ## Calculations of MO
	out, X, Y, Z = MO(data, MO_list, grid_par=par)
	for series in out:
		## The length product works because all voxels of the ORBKIT grid have the same dimensions
		print np.sum(np.square(series))*(X[1,0,0] - X[0,0,0])*(Y[0,1,0] - Y[0,0,0])*(Z[0,0,1] - Z[0,0,0])
    ## Visulation of the MO
	viz_MO(out, X, Y, Z, data, file_name=file_name, labels=MO_list)

elif job in {"TD", "EDD", "BARY"}:
    ## Check if there is a .json associated with the OPT
	if ".json" not in argv[4]:
		stderr.write(usage)
		exit(1)
    ## Load the data
	with open(argv[4], "r") as f:
		TD_data = load(f)
    ## List existing transitions
	try:
		T_list = map(int, val[1].split(","))
		transitions = [TD_data[0]["results"]["excited_states"]["et_transitions"][i] for i in T_list]
	except IndexError:
		transitions = TD_data[0]["results"]["excited_states"]["et_transitions"]
    ## Calculations
	out, X, Y, Z = TD(data, transitions, grid_par=par)

	if job == "TD":
        ## Returns the calculated values of the tozer_lambda, d_CT, Q_CT, Mu_CT and e- barycenter and hole barycenter
            TD_data[0]["results"]["excited_states"]["Tozer_lambda"] = []
            TD_data[0]["results"]["excited_states"]["d_ct"] = []
            TD_data[0]["results"]["excited_states"]["q_ct"] = []
            TD_data[0]["results"]["excited_states"]["mu_ct"] = []
            TD_data[0]["results"]["excited_states"]["e-_barycenter"] = []
            TD_data[0]["results"]["excited_states"]["hole_barycenter"] = []
            for e in out:
                TD_data[0]["results"]["excited_states"]["Tozer_lambda"].append(e[1])
                TD_data[0]["results"]["excited_states"]["d_ct"].append(e[2][0])
                TD_data[0]["results"]["excited_states"]["q_ct"].append(e[2][1])
                TD_data[0]["results"]["excited_states"]["mu_ct"].append(e[2][2])
                TD_data[0]["results"]["excited_states"]["e-_barycenter"].append(e[2][3].tolist())
                TD_data[0]["results"]["excited_states"]["hole_barycenter"].append(e[2][4].tolist())
            # for l in [(e[1], e[2]) for e in out]:
	    #     print l
            #with open("./test_bdm.json", "w+") as f:
            with open(argv[4], "w+") as f:
                dump(TD_data, f)

	elif job == "EDD":
		viz_EDD([e[0] for e in out], X, Y, Z, data, file_name=file_name, labels=T_list)

	elif job == "BARY":
		B = [e[2][3:] for e in out]
		for P in B:
			print P
		viz_BARY(B, data, file_name=file_name)


elif job == "POTENTIAL":
    ## Calculations
	out_r, out_V, X, Y, Z = Potential(data, grid_par=par)

	P = np.array([X, Y, Z])

	if file_name is not None:
		np.save("{}-P".format(file_name), P)
		np.save("{}-rho".format(file_name), out_r)
		np.save("{}-V".format(file_name), out_V)

	viz_Potential(out_r, out_V, X, Y, Z, data, file_name=file_name)

elif job == "FUKUI":

    """Usage : {0} FUKUI $npts ${{OPT-input}}.json ${{SP-input}}.json ${{output}}"""
            
## Check if there is two .json files given to discrete.py
    if ".json" not in argv[4]:
        stderr.write(usage)
        exit(1)
            
## Open and load the datas from .json files
    with open(argv[4], "r") as f_sp:
        data_sp = load(f_sp)

## Detection if SP is SP_minus of SP_plus, label is given to make the calculations in correct order to calc_orb.py and viz_orb_mayavi.py
    charge_opt = data[0]["molecule"]["charge"]
    charge_sp = data_sp[0]["molecule"]["charge"]

    if charge_opt - charge_sp == +1:
        print ("f_plus")
        label = "plus"
    elif charge_opt - charge_sp == -1:
        print("f_minus")
        label = "minus"
    else:
        print("Error, no correct json file")
        exit(1)

## Calculation of delta_rho plus or minus for visualisation
    delta_rho, X, Y, Z = Fukui(data, data_sp, label=label, grid_par=par)
## Visualisation of fukui plus or minus
    viz_Fukui(delta_rho, X, Y, Z, data, file_name=file_name, labels=label)

## Calculation of the Reactivity_indices depending on the SP (plus or minus)
    if label == 'plus':
        A, fplus_lambda_mulliken, fplus_lambda_hirshfeld = CDFT_plus_Indices(data, data_sp)
    ## Initialisation of the new entry for .json dict
        data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_affinity"] = []
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_plus"] = []
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_plus"] = []

    ## Data entry
        data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_affinity"] = A
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_plus"] = fplus_lambda_mulliken.tolist()
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_plus"] = fplus_lambda_hirshfeld.tolist()

    ## Writing the .json
        #with open("./test_bdm.json", "w+") as f:
        with open(argv[4], "w+") as f:
            dump(data_sp, f)

    elif label == 'minus':
        I, fminus_lambda_mulliken, fminus_lambda_hirshfeld = CDFT_minus_Indices(data, data_sp)
     ## Initialisation of the new entry for .json dict
        data[0]["results"]["wavefunction"]["reactivity_indices"]["ionization_potential"] = []
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_minus"] = []
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_minus"] = []

    ## Data entry
        data[0]["results"]["wavefunction"]["reactivity_indices"]["ionization_potential"] = I
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_minus"] = fminus_lambda_mulliken.tolist()
        data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_minus"] = fminus_lambda_hirshfeld.tolist()

    ## Writing the .json
        #with open("./test_bdm.json", "w+") as f:
        with open(argv[4], "w+") as f:
            dump(data_sp, f)

elif job == 'DUAL':

    """{0} DUAL $npts ${{OPT-input}}.json ${{SP-input}}.json ${{SP-input}}.json ${{output}}"""

## Check if there is three .json files given to discrete.py
    if ".json" not in argv[4] and argv[5]:
        stderr.write(usage)
        exit(1)

## Loading datas from SP json
    with open(argv[4], "r") as f_sp1:
        data_sp1 = load(f_sp1)
    with open(argv[5], "r") as f_sp2:
        data_sp2 = load(f_sp2)

    charge_opt = data[0]["molecule"]["charge"]
    charge_sp1 = data_sp1[0]["molecule"]["charge"]
    charge_sp2 = data_sp2[0]["molecule"]["charge"]

## Assignation of json file to SP_minus or SP_plus
    if charge_opt - charge_sp1 == +1:
        data_spminus = data_sp1
    if charge_opt - charge_sp1 == -1:
        data_spplus = data_sp1
    if charge_opt - charge_sp2 == +1:
        data_spminus = data_sp2
    if charge_opt - charge_sp2 == -1:
        data_spplus = data_sp2
    if data_spminus is None or data_spplus is None:
        print("Error, .json files are not correct for Fukui indices calculations")
        exit(1)

## Calculation of delta_rho for the visualisation of the dual reactivity descriptor
    delta_rho_2, X, Y, Z = Fdual(data, data_spplus, data_spminus)
## Visualisation of fukui plus or minus
    viz_Fdual(delta_rho_2, X, Y, Z, data, file_name=file_name)

## Calculation of the reactivity indices
    A, I, Khi, Eta, Omega, DeltaN, fplus_lambda_mulliken, fminus_lambda_mulliken, fdual_lambda_mulliken, fplus_lambda_hirshfeld, fminus_lambda_hirshfeld, fdual_lambda_hirshfeld, c_Omega_plus_mulliken, c_Omega_minus_mulliken, c_Omega_dual_mulliken, c_Omega_plus_hirshfeld, c_Omega_minus_hirshfeld, c_Omega_dual_hirshfeld = CDFT_Indices(data, data_spplus, data_spminus)

## Initialization of the new data in .json file
    data[0]["results"]["wavefunction"]["reactivity_indices"] = {}
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_affinity"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["ionization_potential"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electronegativity"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["hardness"] = [] 
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electrophilicity"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_flow"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_plus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_plus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_minus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_minus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["dual_mulliken"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["dual_hirshfeld"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_plus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_minus"] = []
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_dual"] = []

## Data writing    
    
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_affinity"] = A
    data[0]["results"]["wavefunction"]["reactivity_indices"]["ionization_potential"] = I
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electronegativity"] = Khi
    data[0]["results"]["wavefunction"]["reactivity_indices"]["hardness"] = Eta
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electrophilicity"] = Omega
    data[0]["results"]["wavefunction"]["reactivity_indices"]["electron_flow"] = DeltaN
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_plus"] = fplus_lambda_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_plus"] = fplus_lambda_hirshfeld.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_mulliken_minus"] = fminus_lambda_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_fukui_hirshfeld_minus"] = fminus_lambda_hirshfeld.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["fukui_dual_mulliken"] = fdual_lambda_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["fukui_dual_hirshfeld"] = fdual_lambda_hirshfeld.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_plus_mulliken"] = c_Omega_plus_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_minus_mulliken"] = c_Omega_minus_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_dual_mulliken"] = c_Omega_dual_mulliken.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_plus_hirshfeld"] = c_Omega_plus_hirshfeld.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_minus_hirshfeld"] = c_Omega_minus_hirshfeld.tolist()
    data[0]["results"]["wavefunction"]["reactivity_indices"]["condensed_electrophilicity_dual_hirshfeld"] = c_Omega_dual_hirshfeld.tolist()

## Writing the .json
   #with open("./test_bdm.json", "w+") as f:
    with open(argv[3], "w+") as f:
        dump(data, f)


    pass
