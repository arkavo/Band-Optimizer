# NOTES: Since all the structures we use in this is tetragonal, we're not changing the
# atomic lattice positions from (0, 0, 0) to (0.25, 0.25, 0.25).

# Start with decent upper and lower limits for k_points and ecut.
# if error "Cannot read output" appears, it means that the scf.in file has failed.
# k_cut must be in integer format(9 instead of 9.0) or else the scf.in file will fail.


import subprocess as sp
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#===============================================================================
# Change whatever you want here

raw_input_directory = "./Raw_Inputs/"
raw_data_directory = "./Raw_data/"
raw_output_directory = "./Raw_Outputs/"
# UPF files

# Used limits

# Sn = 4 15 10 120 5. 15. 3
# C  = 4 32 10 120 4. 20. 3

# Change this to your UPF file

# upf = "Sn.pbe-dn-rrkjus_psl.1.0.0.UPF"
# upf = "C.pbe-n-rrkjus_psl.1.0.0.UPF"
upf = "Ge.pbe-dn-kjpaw_psl.1.0.0.UPF"
# Set K points limits (integers ONLY)
K_POINT_LOWER_LIMIT = 4
K_POINT_UPPER_LIMIT = 32
# Set ecutwfc limits
ECUT_LOWER_LIMIT = 10
ECUT_UPPER_LIMIT = 120
# Set lattice_parameter limits
LATT_K_LOWER_LIMIT = 4.
LATT_K_UPPER_LIMIT = 20.
# Set accuracy level
ACCURACY_LEVEL = 3
# Device specific
# ESPRESSO_path = "~/espresso-5.0.2/build_ompi/bin/"
# default
ESPRESSO_path = "pw.x"
# Don't change anything below this line
#===============================================================================
cwd = os.getcwd()
cl_ele = upf.split(".")[0]
element = "\'"+upf.split(".")[0]+"\'"
mass = 0
table = json.loads(open('Table.json').read())
for elem in table['elements']:
    if str(elem['symbol']).lower() == cl_ele.lower():
        name = elem['name']
        mass = float(elem['atomic_mass'])
        ec = elem['shells']
        print(f"Element: {name}")
        print(f"Atomic mass: {mass}")
        print(f"Shells: {ec}")
        print("\n________________________________________________________________________________\n\n")
        break

def clean():
    print("Cleaning... in/out files")
    os.system("rm *.in *.out")
    print("Done")
    print("Cleaning... raw data")
    os.chdir(path=raw_data_directory)
    os.system("rm -rf *")
    print("Clean complete\n")
    os.chdir(path=cwd)
    
def SCF_INPUT(upf, latt_k, ecut, k_pts):
    
    element = upf.split(".")[0]
    fname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_" + "scf.in"
    oname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_" + "scf.out"
    infile = open(fname, "w")
    automatic = "{automatic}"
    
    blocktext = f"""\
    &control
        calculation = 'scf',
        restart_mode = 'from_scratch',
        pseudo_dir = './Raw_Inputs/',
        outdir = './Raw_data/',
        prefix = '{element}'
    /
    &system    
        ibrav=2, celldm(1)={latt_k}, nat=2, ntyp=1,
        ecutwfc ={ecut}
    /
    &electrons
        conv_thr =  1.0d-8
        mixing_beta = 0.7
        diagonalization = 'david'
    /
    
    ATOMIC_SPECIES
    {element}   {mass}   {upf} 
    ATOMIC_POSITIONS
    {element}  0.00 0.00 0.00
    {element}  0.25 0.25 0.25
    K_POINTS  {automatic}
    {k_pts} {k_pts} {k_pts} 0 0 0"""

    infile.write(blocktext)
    infile.close()
    print("{fname} created".format(fname=fname))
    cmd = ESPRESSO_path + " < " + fname + " > " + oname
    os.system(cmd)
    
    cmd2 = "grep ! " + oname
    a = os.popen(cmd2).read()
    a = a.split(" ")
    
    return [float(latt_k), float(ecut), float(k_pts), float(a[-2])]



def k_pt_opt(lower_lt=K_POINT_LOWER_LIMIT, upper_lt=K_POINT_UPPER_LIMIT):
    data = np.array([])
    
    for k_points in range(lower_lt, upper_lt):
        data = np.append(data, SCF_INPUT(upf, 10.0, 20, k_points))
        print(f"K_POINTS: {k_points} = DONE")
    
    data = data.reshape(-1, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title("Energy vs K points")
    ax.set_xlabel("K point")
    ax.set_ylabel("Energy (eV)")
    ax = plt.plot(data[:,2], data[:,3])
    ind = np.argmin(data[:,3])
    k_point_min = int(data[:,2][ind])
    plt.savefig(raw_output_directory + name + "_k-points.png")
    print(data)
    print(f"\n\nK-point optimized = {k_point_min}\n")
    KPT = pd.DataFrame(data, columns=["Lattice Parameter", "Ecutwfc", "K points", "Energy"])
    KPT = KPT[["K points", "Energy"]]
    KPT.to_csv(raw_output_directory + name + "_k points.csv",index=False,sep=" ",header=False)
    print("\nData written to csv\n")
    return k_point_min



def ecut_opt(lower_ecut=ECUT_LOWER_LIMIT, upper_ecut=ECUT_UPPER_LIMIT):
    data = np.array([0,lower_ecut,res_k_pt,0])
    ec_prev = 0
    for ecut in range(lower_ecut, upper_ecut, 10):
        data = np.append(data, SCF_INPUT(upf, 10.0, ecut, int(res_k_pt)))
        data = data.reshape(-1, 4)
        ec_next = float(data[-1,3])
        if (np.abs(ec_prev - ec_next) / np.abs(ec_next)) < 0.0001:
            print("Next value within 0.01 of previous value")
            print(f"terminate ECUT = {ecut}")
            break
        print(f"Ecut: {ecut} = DONE")
        ec_prev = ec_next
    
    data = data.reshape(-1, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title("Energy vs Ecut")
    ax.set_xlabel("Ecut")
    ax.set_ylabel("Energy (eV)")
    ax = plt.plot(data[1:,1], data[1:,3])
    ind = np.argmin(data[1:,3])
    ecut_min = data[1:,1][ind]
    plt.savefig(raw_output_directory + name + "_ecut.png")
    print(data)
    print(f"\n\nEcut optimized = {ecut_min}\n")
    ECUT = pd.DataFrame(data, columns=[
                        "Lattice Parameter", "Ecutwfc", "K points", "Energy"])
    ECUT = ECUT[["Ecutwfc", "Energy"]]
    ECUT.to_csv(raw_output_directory + name + "_ecut.csv",index=False,sep=" ",header=False)
    print("\nData written to csv\n")
    return ecut_min



def latt_opt(lower_lt=LATT_K_LOWER_LIMIT, upper_lt=LATT_K_UPPER_LIMIT,LEVEL=ACCURACY_LEVEL):
    #data = np.array([])
    fig = plt.figure(figsize=((LEVEL+1)*10,10))
    for level in range(LEVEL):
        data = np.array([])
        print(f"\nBeginning lattice optimization level {level}\n")
        for latt in np.arange(lower_lt, upper_lt,(2*0.1**level)):
            latt = round(latt, level)
            data = np.append(data, SCF_INPUT(upf, latt, int(res_ecut), int(res_k_pt)))
            print(f"Latt: {latt} = DONE")
        
        data = data.reshape(-1, 4)
        ind = np.argmin(data[1:,3])
        latt_min = data[1:,0][ind]
        print(f"\n\nLattice minimum at level {level} = {latt_min}\n")
        lower_lt = latt_min - (0.1**level)
        upper_lt = latt_min + (0.1**level)
        ax = fig.add_subplot(1 ,LEVEL ,level+1)
        plt.title(label=f"Lattice optimization level {level}")
        ax.set_xlabel("Lattice parameter")
        ax.set_ylabel("Energy (eV)")
        ax = plt.plot(data[1:,0], data[1:,3])
    plt.savefig(raw_output_directory + name + "_lattice.png")
    
    LATT = pd.DataFrame(data, columns=["Lattice Parameter", "Ecutwfc", "K points", "Energy"])
    LATT = LATT[["Lattice Parameter", "Energy"]]
    LATT.to_csv(raw_output_directory + name + "_lattice_level_{level}.csv",index=False,sep=" ",header=False)
    print("\nData written to csv\n")
    return latt_min

res_k_pt = k_pt_opt()
res_ecut = ecut_opt()
res_latt = latt_opt()

print("\nCreate final .in file with optimized parameters\n")
SCF_INPUT(upf, res_latt, res_ecut, res_k_pt)
fname = element + "_" + str(res_latt) + "_" + str(res_ecut) + "_" + str(res_k_pt) + "_" + "scf.in"
cmd2 = "mv " + fname + " " + raw_output_directory
os.system(cmd2)
print("Final .in file created with optimized parameters\n")

clean()

