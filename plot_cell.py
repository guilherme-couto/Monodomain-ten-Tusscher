import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print('Usage: python3 plot_cell.py <dt>')
    sys.exit(1)

dt = sys.argv[1]

# Read data from files
t = []
timesfile = f'./simulation-files/sim-times-cell-{dt}.txt'
with open(timesfile, 'r') as f:
    for line in f:
        t.append(float(line))

filename = f'./simulation-files/tnnp-cell-{dt}.txt'
U = []
with open(filename, 'r') as f:
    for line in f:
        U.append(float(line))

# Make plot
plt.plot(t, U)
plt.title(f'Cell TNNP Epi')
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')

plt.savefig('cell-epi.png')
plt.close()


    
    