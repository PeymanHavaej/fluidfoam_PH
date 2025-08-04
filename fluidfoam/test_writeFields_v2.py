import fluidfoam
from fluidfoam_write import OpenFoamFileWriter
import numpy as np

# Test scalar field with boundary conditions from source file
writer = OpenFoamFileWriter("output_samples/ascii", time_name="0", is_ascii=True, verbose=True)
read_data = fluidfoam.readscalar("output_samples/ascii", "0", "alpha")
alpha_new = read_data 
writer.writefield(
    field_data=alpha_new,
    field_name="alpha_write_2",
    field_type="scalar",
    source_file="output_samples/ascii/0/alpha"
)

# Verify by reading back the written field
written_data = fluidfoam.readscalar("output_samples/ascii", "0", "alpha_write_2")
assert np.allclose(written_data, alpha_new), "Written and read-back data do not match"
print("Test passed: Written scalar field matches input data scaled by 1.5")

# Test scalar field without boundary conditions
writer.writefield(
    field_data=alpha_new,
    field_name="alpha_write_3",
    field_type="scalar"
)
print("Test passed: Written scalar field without boundaryField")
