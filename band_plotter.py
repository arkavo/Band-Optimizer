import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


OUT_DIR = "./Band_plots/"
IN_FILE = "Ge_bandx.dat"
IN_DIR = OUT_DIR+IN_FILE

element = IN_FILE.split('_')[0]
with open(IN_DIR) as f:
    data = f.readlines()
    dataline = data[0].split()
    nband = int(dataline[2][:-1])
    nkpt = int(dataline[4])
    length = len(data)
    unit = int((length-1)/nkpt)
    print(f"nbands = {nband}")
    print(f"nkpts = {nkpt}")
    print(f"unit data length = {unit}\n\n")
    
    x0 = 0.
    y0 = 0.
    z0 = 0.
    
    x1 = 0.
    y1 = 0.
    z1 = 0.
    band_data = np.array([])
    k_pts = np.array([])
    r = 0
    for i in range(nkpt):
        k_line = data[i*unit+1].split()
        x0 = float(k_line[0])
        y0 = float(k_line[1])
        z0 = float(k_line[2])
        r += np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        k_pts = np.append(k_pts, r)
        x1 = x0
        y1 = y0
        z1 = z0
        bandline = np.array([])
        #band_data = np.array([])
        for j in range(1,unit):
            bandline = data[1+i*unit+j].split()
            for item in bandline:
                item = float(item)
                band_data = np.append(band_data, item)
    band_data = band_data.reshape(nkpt, nband)
    symm_lines = k_pts[0], k_pts[40], k_pts[80], k_pts[120]
    symm_points= ["L","G","K","U"]
    print(f"symmertry lines = {symm_lines}\n")
    ymax = np.max(band_data)
    ymin = np.min(band_data)
    absmax = np.max([np.abs(ymax), np.abs(ymin)]) + 0.5
    xmin = np.min(k_pts)
    xmax = np.max(k_pts)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-absmax, absmax)
    ax.set_xlabel("k-points")
    ax.set_ylabel("Energy")
    ax.set_title(element + " Band structure")
    ax.set_xticks(symm_lines)
    ax.set_xticklabels(symm_points)
    for i in range(nband):
        ax.plot(k_pts, band_data[:,i])
    for i in range(len(symm_lines)):
        ax.axvline(x=symm_lines[i], color='r', linestyle='--')
    ax.axhline(y=0, color='k', linestyle=':')
    
    print("DONE")
    outfile = OUT_DIR + element+"band_plot.png"
    plt.savefig(outfile)
    print(f"Saved to {outfile}\n")
    cmd = "mv *.dat *.gnu *.rap " + OUT_DIR
    os.system(cmd)
    print(f"Raw files Moved to {OUT_DIR}")
