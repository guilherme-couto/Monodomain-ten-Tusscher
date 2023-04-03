import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import imageio.v2

if len(sys.argv) != 6:
    print('Usage: python3 plot.py <method> <dt_ode> <dt_pde> <total_frames> <frame_rate>')
    sys.exit(1)

method = sys.argv[1]
dt_ode = sys.argv[2]
dt_pde = sys.argv[3]
totalframes = int(sys.argv[4])
framerate = int(sys.argv[5])

U = np.zeros((totalframes, 200, 200))
t = np.zeros(totalframes)

filename = f'tnnp-{method}-{dt_ode}-{dt_pde}.txt'
timesfile = f'sim-times-{method}-{dt_ode}-{dt_pde}.txt'

# Read data from files
with open(filename, 'r') as f:
    for n in range(totalframes):
        for i in range(len(U[0])):
            for j in range(len(U[0])):
                U[n][i][j] = float(f.readline())

with open(timesfile, 'r') as f:
    for n in range(totalframes):
        t[n] = float(f.readline())

# Make plots
plots = []
for n in range(len(U)):
    if n % framerate == 0:
        plotname = 'plot-' + str(n) + '.png'
        plots.append(plotname)
        
        plt.imshow(U[n], cmap='plasma', vmin=-100, vmax=50)
        plt.colorbar(label='V (mV)')
        plt.title(f'Monodomain TNNP {method} dtode = {dt_ode} dtpde = {dt_pde} t = {t[n]:.2f}')
        plt.xticks([])
        plt.yticks([])
        
        plt.savefig(plotname)
        plt.close()

# Build gif
with imageio.v2.get_writer(f'{method}-{dt_ode}-{dt_pde}.gif', mode='I') as writer:
    for plot in plots:
        image = imageio.v2.imread(plot)
        writer.append_data(image)
        
# Remove files
for png in set(plots):
    os.remove(png)

    
    