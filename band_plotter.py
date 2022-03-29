import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

IN_FILE ="Ge_bandX.dat"
OUT_DIR = "./Band_plots/"

with open(IN_FILE) as f:
    data = f.readlines()
    dataline = data[0].split()
    nband = int(dataline[2][:-1])
    nkpt = int(dataline[4])
    length = len(data)
    unit = int((length-1)/nkpt)
    print(f"nbands = {nband}")
    print(f"nkpts = {nkpt}")
    print(f"unit data length = {unit}")
    
    x0 = 0.
    y0 = 0.
    z0 = 0.
    
    x1 = 0.
    y1 = 0.
    z1 = 0.
    band_data = np.array([])
    k_pts = np.array([])
    for i in range(nkpt):
        k_line = data[i*unit+1].split()
        x0 = float(k_line[0])
        y0 = float(k_line[1])
        z0 = float(k_line[2])
        k_pts = np.append(k_pts, [np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)])
        bandline = np.array([])
        #band_data = np.array([])
        for j in range(1,unit):
            bandline = data[1+i*unit+j].split()
            for item in bandline:
                item = float(item)
                #print(item)
                band_data = np.append(band_data, item)
            #print(bandline)
            #np.append(band_data, bandline)
    band_data = band_data.reshape(nkpt, nband)
        #print(f"bandline = {bandline}")
    #print(band_data)
    #print(k_pts)
    symm_lines = k_pts[0], k_pts[40], k_pts[80], k_pts[120]
    print(f"symmertry lines = {symm_lines}")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(nband):
        ax.plot(k_pts, band_data[:,i])
    for i in range(len(symm_lines)):
        ax.axvline(x=symm_lines[i], color='r', linestyle='--')
    element = IN_FILE.split('_')[0]
    print("DONE")
    outfile = OUT_DIR + element+"band_plot.png"
    plt.savefig(outfile)
    print(f"Saved to {outfile}")
    cmd = "mv *.dat *.gnu *.rap "+OUT_DIR
    os.system(cmd)
    print(f"Raw files Moved to {OUT_DIR}")