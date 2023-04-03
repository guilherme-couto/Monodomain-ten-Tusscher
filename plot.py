import numpy as np
import matplotlib.pyplot as plt

U = np.zeros((50, 200, 200))


# Read data from file
with open('tnnp-ADI-0.020-0.040.txt', 'r') as f:
    for k in range(len(U)):
        for i in range(len(U[0])):
            for j in range(len(U[0])):
                U[k][i][j] = float(f.readline())

# Fill variable
plt.imshow(U[5], cmap='hot', vmin=-100, vmax=50)
plt.colorbar()
plt.savefig('tnnp-ADI-0.020-0.040.png')
    
    