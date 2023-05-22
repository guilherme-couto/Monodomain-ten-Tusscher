import os
import numpy as np
import matplotlib.pyplot as plt

# Execute the pacing for different BCLs

# Initialize the BCLs
bcls = np.arange(1000, 200, -10)

# Clean APD file
os.system('rm -f ./simulation-files/apd90.txt')

# Iterate over the BCLs
for bcl in bcls:
    # Execute the pacing protocol
    os.system(f'./pacing {bcl}')
    
    # Build the graph
    os.system(f'cd plots && python plot_pacing.py {bcl:.1f}')

# Plot the APD90 vs BCL graph
bcl = []
apd90 = []
with open('./simulation-files/apd90.txt', 'r') as f:
    for line in f:
        data = line.split(' | ')
        bcl.append(float(data[0]))
        apd90.append(float(data[1]))
        
plt.plot(bcl, apd90, '-o')
plt.xlabel('BCL (ms)')
plt.ylabel('APD90 (ms)')
plt.title('APD90 vs BCL')
plt.savefig('./png/apd90_vs_bcl.png')