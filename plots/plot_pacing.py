import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print('Usage: python3 plot_cell.py <BCL (ms)>')
    sys.exit(1)

dt = 0.05
BCL = sys.argv[1]

# Read data from files
t = []
timesfile = f'../simulation-files/times-pacing-{dt}-{BCL}.txt'
with open(timesfile, 'r') as f:
    for line in f:
        t.append(float(line))

filename = f'../simulation-files/tnnp-pacing-{dt}-{BCL}.txt'
U = []
with open(filename, 'r') as f:
    for line in f:
        U.append(float(line))

# Make plot
plt.plot(t, U)
plt.title(f'TNNP EPI Pacing (BCL = {BCL} ms)')
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')

plt.savefig(f'../png/pacing-{BCL}.png')
plt.close()


    
    