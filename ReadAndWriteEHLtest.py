from fluidfoam import readfield, readmesh, OpenFoamFileWriter
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Load mesh and field
sol = './EHL_output'
timename = '0'
x, y, z = readmesh(sol, structured=False)
p = readfield(sol, timename, 'p')
p_new = p - 101325


# Load x and p arrays from Reynolds_output directory
x_reynolds = np.load('Reynolds_output/x.npy')
p_reynolds = np.load('Reynolds_output/p.npy')

# Interpolate 1D Reynolds pressure onto 2D EHL mesh

# Create 1D interpolator for Reynolds pressure
interp_func = interp1d(x_reynolds, p_reynolds, bounds_error=False, fill_value="extrapolate")

# Interpolate onto x-coordinates of EHL mesh
p_reynolds_interp = interp_func(x)

# Optionally, reshape if needed to match p_new's shape
if p_new.shape == p_reynolds_interp.shape:
    p_new = p_reynolds_interp
else:
    p_new = p_reynolds_interp.reshape(p_new.shape)







writer = OpenFoamFileWriter(sol, timename, is_ascii=True, verbose=True)

writer.writefield(
    field_data=p_new,
    field_name="p_new",
    field_type="scalar",
    source_file="EHL_output/0/p"
)
