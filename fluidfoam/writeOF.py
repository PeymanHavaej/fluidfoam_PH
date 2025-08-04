"""
Writing OpenFOAM Files with Python
=================================
This module provides a class to write OpenFOAM files, complementing the reading
capabilities of the fluidfoam module.

.. autoclass:: OpenFoamFileWriter
   :members:
"""

import os
import gzip
import struct
import numpy as np

# Define colors for verbose output (consistent with fluidfoam.py)
W = '\033[0m'  # white (normal)
R = '\033[31m'  # red
G = '\033[32m'  # green
O = '\033[33m'  # orange
B = '\033[34m'  # blue
P = '\033[35m'  # purple

def _make_path(path, time_name=None, name=None, region=None):
    """Construct file path for OpenFOAM case, supporting time and region directories."""
    if region is not None:
        path = os.path.join(path, "constant", region, name)
    elif time_name is not None and name is not None:
        path = os.path.join(path, time_name, name)
    elif name is not None:
        path = os.path.join(path, name)
    return path

class OpenFoamFileWriter:
    """Class to write OpenFOAM files, such as field data, in ASCII or binary format."""

    def __init__(
        self,
        path,
        time_name="0",
        region=None,
        is_ascii=True,
        is_compressed=False,
        precision=15,
        verbose=True
    ):
        """
        Initialize the OpenFOAM file writer.

        Args:
            path: str
                Path to the OpenFOAM case directory.
            time_name: str
                Time directory name (default: '0').
            region: str, optional
                Region name for multi-region cases (default: None).
            is_ascii: bool
                Write files in ASCII format if True, binary if False (default: True).
            is_compressed: bool
                Compress files using gzip if True (default: False).
            precision: int
                Number of decimal places for floating-point values in ASCII mode (default: 15).
            verbose: bool
                Print progress messages if True (default: True).
        """
        self.pathcase = path
        self.time_name = time_name
        self.region = region
        self.is_ascii = is_ascii
        self.is_compressed = is_compressed
        self.precision = precision
        self.verbose = verbose
        self._open = gzip.open if is_compressed else open
        self.file_ext = ".gz" if is_compressed else ""
        self.file_mode = "wt" if is_ascii else "wb"

    def _write_header(self, f, class_name, object_name, location=None):
        """Write OpenFOAM file header."""
        if location is None:
            location = self.time_name if self.region is None else f"constant/{self.region}"
        header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Web:      www.OpenFOAM.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    {class_name} file {object_name} written by fluidfoam

\\*---------------------------------------------------------------------------*/

FoamFile
{{
    version     2.0;
    format      {"ascii" if self.is_ascii else "binary"};
    class       {class_name};
    location    "{location}";
    object      {object_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
        f.write(header.encode() if not self.is_ascii else header)

    def writefield(
        self,
        field_data,
        field_name,
        field_type,
        structured=False,
        boundary=None,
        order="F"
    ):
        """
        Write OpenFOAM field data (scalar, vector, symmtensor, or tensor).

        Args:
            field_data: ndarray
                Array containing field data. Shape depends on field_type:
                - scalar: (n_cells,) or (nx, ny, nz) if structured
                - vector: (3, n_cells) or (3, nx, ny, nz) if structured
                - symmtensor: (6, n_cells) or (6, nx, ny, nz) if structured
                - tensor: (9, n_cells) or (9, nx, ny, nz) if structured
            field_name: str
                Name of the field (e.g., 'U', 'alpha', 'sigma').
            field_type: str
                Type of field ('scalar', 'vector', 'symmtensor', 'tensor').
            structured: bool
                If True, assumes field_data is shaped for a structured mesh (default: False).
            boundary: str, optional
                Boundary patch name if writing boundary field data (default: None).
            order: str
                Array ordering for structured meshes, 'F' (Fortran) or 'C' (C) (default: 'F').

        Returns:
            None

        A way you might use me is:
            writer = OpenFoamFileWriter('path_of_OpenFoam_case')
            field_data = np.ones(1000)
            writer.writefield(field_data, 'alpha', 'scalar')
        """
        # Validate field_type
        valid_types = ["scalar", "vector", "symmtensor", "tensor"]
        if field_type not in valid_types:
            raise ValueError(f"field_type must be one of {valid_types}")

        # Construct file path
        field_path = _make_path(self.pathcase, self.time_name, field_name, self.region)
        os.makedirs(os.path.dirname(field_path), exist_ok=True)
        field_path += self.file_ext

        # Determine class name and number of components
        class_map = {
            "scalar": ("scalarField", 1, "[0 0 0 0 0 0 0]"),
            "vector": ("vectorField", 3, "[0 1 -1 0 0 0 0]"),
            "symmtensor": ("symmTensorField", 6, "[0 2 -2 0 0 0 0]"),
            "tensor": ("tensorField", 9, "[0 2 -2 0 0 0 0]")
        }
        class_name, n_components, dimensions = class_map[field_type]

        # Prepare field data
        field_data = np.asarray(field_data)
        if structured:
            if field_type == "scalar":
                if field_data.ndim != 3:
                    raise ValueError("Structured scalar field must have 3D shape (nx, ny, nz)")
                field_data = field_data.ravel(order=order)
            else:
                if field_data.ndim != 4 or field_data.shape[0] != n_components:
                    raise ValueError(f"Structured {field_type} field must have shape ({n_components}, nx, ny, nz)")
                field_data = field_data.reshape(n_components, -1, order=order).T.flatten()
        else:
            if field_type == "scalar":
                if field_data.ndim != 1:
                    raise ValueError("Unstructured scalar field must be 1D")
            else:
                if field_data.ndim != 2 or field_data.shape[0] != n_components:
                    raise ValueError(f"Unstructured {field_type} field must have shape ({n_components}, n_cells)")

        n_points = field_data.size // n_components if field_type != "scalar" else field_data.size

        if self.verbose:
            print(f"Writing {field_type} field to {field_path}")

        with self._open(field_path, self.file_mode) as f:
            # Write header and dimensions
            self._write_header(f, class_name, field_name)
            f.write(f"dimensions      {dimensions};\n\n".encode() if not self.is_ascii else f"dimensions      {dimensions};\n\n")
            f.write(f"internalField   nonuniform {class_name}({n_points})\n(\n".encode() if not self.is_ascii else f"internalField   nonuniform {class_name}({n_points})\n(\n")

            # Write field data
            if self.is_ascii:
                if field_type == "scalar":
                    for val in field_data:
                        f.write(f"{val:.{self.precision}f}\n")
                else:
                    for i in range(n_points):
                        vals = field_data[i * n_components:(i + 1) * n_components]
                        f.write(f"({' '.join(f'{v:.{self.precision}f}' for v in vals)})\n")
            else:
                if field_type == "scalar":
                    f.write(struct.pack(f"{n_points}d", *field_data))
                else:
                    for i in range(n_points):
                        vals = field_data[i * n_components:(i + 1) * n_components]
                        f.write(struct.pack(f"{n_components}d", *vals))

            f.write(b")\n" if not self.is_ascii else ")\n")

            # Write boundary field if specified
            if boundary is not None:
                boundary_field = f"""
boundaryField
{{
    {boundary}
    {{
        type            calculated;
        value           uniform (0{' 0' * (n_components - 1)});
    }}
}}
"""
                f.write(boundary_field.encode() if not self.is_ascii else boundary_field)

        if self.verbose:
            print(f"Field {field_name} written successfully to {field_path}")
