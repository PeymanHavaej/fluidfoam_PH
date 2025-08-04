from fluidfoam import readfield, readmesh
import numpy as np
import matplotlib.pyplot as plt

# Load mesh and field
sol = './'
timename = '1.09998478e-05'
x, y, z = readmesh(sol, structured=False)
T = readfield(sol, timename, 'T')

# Stack coordinates and T values
coords = np.vstack((x, y, z)).T
T = np.array(T)

# Define x bins (adjust number of bins as needed)
num_bins = 100
x_bins = np.linspace(x.min(), x.max(), num_bins + 1)
x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

# Bin-averaged T along x
T_avg_x = np.zeros(num_bins)
for i in range(num_bins):
    in_bin = (x >= x_bins[i]) & (x < x_bins[i + 1])
    if np.any(in_bin):
        T_avg_x[i] = np.mean(T[in_bin])
    else:
        T_avg_x[i] = np.nan  # Or leave as 0, depending on preference

# Plot
plt.figure()
plt.plot(x_centers, T_avg_x, '-o')
plt.xlabel('x')
plt.ylabel('Average Temperature')
plt.title('1D Profile of Averaged Temperature along x')
plt.grid(True)
plt.tight_layout()
plt.show()
