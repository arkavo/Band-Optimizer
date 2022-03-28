import os
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

upf = "Ge.pbe-dn-kjpaw_psl.1.0.0.UPF"
file = "Ge_11.12_50.0_24_scf.in"
ESPRESSO_path = "pw.x"
BANDS_path = "bands.x"
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

def READER(file):
    opt_pts = file.split("_")
    latt_k = float(opt_pts[1])
    ecut = float(opt_pts[2])
    k_pts = int(opt_pts[3])
    return latt_k, ecut, k_pts

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
    fname = element + "_" + str(latt_k) + "_" + \
        str(ecut) + "_" + str(k_pts) + "_" + "scf.in"
    oname = element + "_" + str(latt_k) + "_" + \
        str(ecut) + "_" + str(k_pts) + "_" + "scf.out"
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




def BAND_INPUT(upf, latt_k, ecut, k_pts):
    element = upf.split(".")[0]
    fname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_band.in"
    oname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_band.out"
    infile = open(fname, "w")
    crystal_b = "{crystal_b}"
    
    
    blocktext = f"""\
    &control
        calculation='bands'
        pseudo_dir = './Raw_Inputs/',
        outdir='./Raw_data/',
        prefix='{element}'
    /
    &system
        ibrav=2, celldm(1) ={latt_k}, nat=2, ntyp=1,
        ecutwfc={ecut}, nbnd=4
    /
    &electrons
        conv_thr =  1.0d-8
        mixing_beta = 0.7
        diagonalization='david'
    /
    ATOMIC_SPECIES
    {element}  {mass}  {upf}
    
    ATOMIC_POSITIONS
    {element} 0.00 0.00 0.00
    {element} 0.25 0.25 0.25
    K_POINTS {crystal_b}
    4
    0.000  0.500  0.000  40 !L
    0.000  0.000  0.000  40 !G 
    0.500  0.500  0.000  40 !K
    0.625  0.650  0.250  40 !U
    """
    infile.write(blocktext)
    infile.close()
    print("Input file created")
    cmd = ESPRESSO_path + " < " + fname + " > " + oname
    print("Output file created")
    os.system(cmd)
    
def BAND_OUTPUT(oname):
    element = upf.split(".")[0]
    fname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_bandX.in"
    oname = element + "_" + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_bandX.out"
    blocktext = f"""\
    &BANDS
        prefix='{element}'
        outdir="./Raw_data/"
        filband="BandX.dat"
    /
    """
    infile = open(fname, "w")
    infile.write(blocktext)
    infile.close()
    cmd = BANDS_path + " < " + fname + " > " + oname
    os.system(cmd)
    print("Band output file created")
    return oname
dataline = READER(file)
print(dataline)
latt_k = dataline[0]
ecut = dataline[1]
k_pts = dataline[2]
# SCF_INPUT(upf, latt_k, ecut, k_pts)
ONAME = BAND_INPUT(upf, latt_k, ecut, k_pts)
BAND_OUTPUT(ONAME)

#clean()