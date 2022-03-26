import subprocess as sp
import os
import json
import numpy as np
import matplotlib.pyplot as plt

prefix = ""
raw_input_directory = "./Raw_Inputs/"
raw_data_directory = "./Raw_data/"
cwd = os.getcwd()
#upf = "C.pbe-n-rrkjus_psl.1.0.0.UPF"
upf = "Sn.pbe-dn-rrkjus_psl.1.0.0.UPF"
cl_ele = upf.split(".")[0]
element = "\'"+upf.split(".")[0]+"\'"
mass = 0
table = json.loads(open('Table.json').read())
for elem in table['elements']:
    if str(elem['symbol']).lower() == cl_ele.lower():
        name = elem['name']
        mass = float(elem['atomic_mass'])
        ec = elem['shells']
        print(name, mass)
        print(ec)
        break

# device specific
# ESPRESSO_path = "~/espresso-5.0.2/build_ompi/bin/"
# default
ESPRESSO_path = "pw.x"

def clean():
    print("Cleaning... in/out files")
    os.system("rm *.in *.out")
    print("Done")
    print("Cleaning... raw data")
    os.chdir(path=raw_data_directory)
    os.system("rm -rf *")
    print("Clean complete")
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



def k_pt_opt(lower_lt=8, upper_lt=9):
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
    k_point_min = data[:,2][ind]
    plt.savefig(name + " k points.png")
    print(data)
    print(f"K-point minimum = {k_point_min}")

    return k_point_min

res_k_pt = k_pt_opt()

def ecut_opt(lower_ecut=30, upper_ecut=120):
    data = np.array([0,lower_ecut,res_k_pt,0])
    ec_prev = 0
    for ecut in range(lower_ecut, upper_ecut, 10):
        data = np.append(data, SCF_INPUT(upf, 10.0, ecut, int(res_k_pt)))
        data = data.reshape(-1, 4)
        #print(data)
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
    plt.savefig(name + " ecut.png")
    print(data)
    print(f"Ecut minimum = {ecut_min}")
    
    return ecut_min

res_ecut = ecut_opt()

def latt_opt(lower_lt=1, upper_lt=10,LEVEL=3):
    #data = np.array([])
    fig = plt.figure(figsize=((LEVEL+1)*10,10))
    for level in range(LEVEL):
        data = np.array([])
        print(f"Beginning lattice optimization level {level}")
        for latt in np.arange(lower_lt, upper_lt,(2*0.1**level)):
            data = np.append(data, SCF_INPUT(upf, latt, int(res_ecut), int(res_k_pt)))
            print(f"Latt: {latt} = DONE")
        
        data = data.reshape(-1, 4)
        ind = np.argmin(data[1:,3])
        latt_min = data[1:,0][ind]
        print(f"Lattice minimum at level {level} = {latt_min}")
        lower_lt = latt_min - (0.1**level)
        upper_lt = latt_min + (0.1**level)
        ax = fig.add_subplot(1 ,LEVEL ,level+1)
        plt.title(label=f"Lattice optimization level {level}")
        ax.set_xlabel("Lattice parameter")
        ax.set_ylabel("Energy (eV)")
        ax = plt.plot(data[1:,0], data[1:,3])
    plt.savefig(name + " lattice.png")
    return latt_min

res_latt = latt_opt(lower_lt=5, upper_lt=19, LEVEL=3)
        
clean()


