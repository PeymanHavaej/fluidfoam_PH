from fluidfoam import readfield, readmesh
import numpy as np
import matplotlib.pyplot as plt
from fluidfoam_write import OpenFoamFileWriter

# Load mesh and field
sol = './EHL_output'
timename = '0'
x, y, z = readmesh(sol, structured=False)
T = readfield(sol, timename, 'T')

T_new = T - 313.15

writer = OpenFoamFileWriter(sol, timename, is_ascii=True, verbose=True)

writer.writefield(
    field_data=T_new,
    field_name="T_new",
    field_type="scalar",
    source_file="output_samples/ascii/0/T"
)
