import fluidfoam
from writeOF import OpenFoamFileWriter
import numpy as np

# Test scalar field
writer = OpenFoamFileWriter("output_samples/ascii", time_name="0", is_ascii=True)
read_data = fluidfoam.readscalar("output_samples/ascii", "0", "alpha")
alpha_new = read_data * 1.5  # Scale the field by 1.5
writer.writefield(alpha_new, "alpha_write_2", "scalar")
