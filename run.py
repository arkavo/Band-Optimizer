import subprocess as sp
import os
import json

a = os.popen("grep ! Sn11_12_4_scf.out").read()
a=a.split(" ")
print(a[-2]) 