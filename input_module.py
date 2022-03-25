import subprocess as sp
import os
import json

prefix = ""
raw_input_directory = "./Raw_Inputs/"
directory = "./Raw_data/"
#upf = "C.pbe-n-rrkjus_psl.1.0.0.UPF"
upf = "Sn.pbe-dn-rrkjus_psl.1.0.0.UPF"
element = "\'"+upf.split(".")[0]+"\'"
mass = 0
table = json.loads(open('Table.json').read())
for elem in table['elements']:
    if elem['symbol'] == element:
        name = elem['name']
        mass = float(elem['atomic_mass'])
        print(name, mass)
        break
        
# device specific
# ESPRESSO_path = "~/espresso-5.0.2/build_ompi/bin/"
# default
ESPRESSO_path = "pw.x"


def SCF_INPUT(upf, latt_k, ecut, k_pts):
    
    element = upf.split(".")[0]
    fname = element + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_" + "scf.in"
    oname = element + str(latt_k) + "_" + str(ecut) + "_" + str(k_pts) + "_" + "scf.out"
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
    
    cmd = ESPRESSO_path + " < " + fname + " > " + oname
    os.system(cmd)
    
    cmd2 = "grep ! " + oname
    a = os.popen(cmd2).read()
    a = a.split(" ")
    
    #print(a)
    print(element,latt_k,ecut,k_pts, a[-2])


SCF_INPUT(upf, 11, 12, 4)
#os.system("rm *.in *.out")
#os.system(["cd ./Raw_data", "rm -rf *"])


