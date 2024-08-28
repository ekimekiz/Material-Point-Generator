"""
STL to Tetrahedral Mesh Data Extraction Script

This script processes an input file containing volumes, centroids, and surface triangle nodes extracted from a tetrahedral mesh. It generates two output files:

1. A file containing the volume and centroid information for each tetrahedron, which become material points.
2. A file containing the centroid and normal vector information for each surface triangle associated with the tetrahedrons, which become surface points.

### Usage:
To run the script, use the following command:

    python sp_and_normals.py <input_filename>

Replace `<input_filename>` with the name of your input file, which should be a text file containing the data. For example:

    python sp_and_normals.py some_mesh_tetra_surface_data.txt

### Output:
The script will generate two output files:
1. `<input_filename>_centroids_volumes.txt`: Contains the volume and centroid coordinates (x, y, z) for each tetrahedron.
2. `<input_filename>_centroids_normals.txt`: Contains the centroid coordinates (x, y, z) and the normal vector components (x, y, z) for each surface triangle.

### Example:
If your input file is named `some_mesh_tetra_surface_data.txt`, the script will produce:
- `some_mesh_tetra_surface_data_centroids_volumes.txt`
- `some_mesh_tetra_surface_data_centroids_normals.txt`

These files will contain the processed data that can be used for further analysis or visualization.

### Dependencies:
- Python 3.x
- NumPy library (install using `pip install numpy`)

Make sure to run the script from the directory containing the input file, or provide the full path to the input file.
"""

import os
import sys
import numpy as np

# Check if a filename was provided as a command-line argument
if len(sys.argv) < 2:
    raise ValueError("Please provide the filename as a command-line argument.")

# Get the filename from the command-line arguments
filename = sys.argv[1]

# Change the working directory to the location of this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Attempt to open the file for reading
try:
    file = open(filename, 'r')
except IOError:
    raise IOError(f'Failed to open file: {filename}')

# Initialize arrays and lists to store the data
volumes = []       # To store the volumes of the tetrahedrons
centroids = []     # To store the centroids of the tetrahedrons
triangle_nodes = []  # To store the surface triangle nodes

# Initialize an index for tracking the current centroid
current_index = -1

# Read through the file line by line
line_number = 0  # To keep track of line numbers
for line in file:
    line_number += 1
    if line.strip():  # Ensure the line is not empty
        # Split the line into numbers
        values = np.fromstring(line, dtype=float, sep=' ')
        num_values = len(values)

        if num_values == 4:
            # This line contains volume and centroid data
            current_index += 1
            volumes.append(values[0])  # The first value is the volume
            centroids.append(values[1:4])  # The next three values are the centroid coordinates (x, y, z)
            triangle_nodes.append([])  # Initialize as an empty list since no triangles are on this line
        elif num_values == 13:
            # This line contains volume, a centroid, and one surface triangle
            current_index += 1
            volumes.append(values[0])  # The first value is the volume
            centroids.append(values[1:4])  # The next three values are the centroid coordinates (x, y, z)
            triangle_nodes.append([values[4:13].reshape(3, 3)])  # The last nine values are the triangle vertices
        elif num_values == 22:
            # This line contains volume, a centroid, and two surface triangles
            current_index += 1
            volumes.append(values[0])  # The first value is the volume
            centroids.append(values[1:4])  # The next three values are the centroid coordinates (x, y, z)
            # The last 18 values are the two triangles' vertices
            triangle_nodes.append([values[4:13].reshape(3, 3), values[13:22].reshape(3, 3)])
        elif num_values == 31:
            # This line contains volume, a centroid, and three surface triangles
            current_index += 1
            volumes.append(values[0])  # The first value is the volume
            centroids.append(values[1:4])  # The next three values are the centroid coordinates (x, y, z)
            # The last 27 values are the three triangles' vertices
            triangle_nodes.append([values[4:13].reshape(3, 3), values[13:22].reshape(3, 3), values[22:31].reshape(3, 3)])

# Close the file after reading all lines
file.close()

# Convert lists to numpy arrays for easier manipulation and numerical operations
volumes = np.array(volumes)
centroids = np.array(centroids)

# Generate output filenames based on the input filename
tetra_output_filename = f'{os.path.splitext(filename)[0]}_centroids_volumes.txt'
triangle_output_filename = f'{os.path.splitext(filename)[0]}_centroids_normals.txt'

# Open files to save the tetrahedron centroid and volume information and triangle centroid and normal information
with open(tetra_output_filename, 'w') as tetra_file, open(triangle_output_filename, 'w') as triangle_file:

    # Loop through each centroid and corresponding surface triangles
    for i in range(centroids.shape[0]):
        # Save the tetrahedron centroid and volume information
        tetra_file.write(f'{volumes[i]:.4f} {centroids[i, 0]:.4f} {centroids[i, 1]:.4f} {centroids[i, 2]:.4f}\n')

        # If there are surface triangles associated with this tetrahedron, process them
        if len(triangle_nodes[i]) > 0:
            for triangle in triangle_nodes[i]:
                v1 = triangle[0, :]  # First vertex of the triangle
                v2 = triangle[1, :]  # Second vertex of the triangle
                v3 = triangle[2, :]  # Third vertex of the triangle

                # Calculate the vectors for two edges of the triangle
                edge1 = v2 - v1
                edge2 = v3 - v1

                # Compute the normal vector using the cross product
                normal = np.cross(edge1, edge2)

                # Calculate the centroid of the triangle (average of the three vertices)
                centroid_of_triangle = (v1 + v2 + v3) / 3

                # Normalize the normal vector to unit length
                normal = normal / np.linalg.norm(normal)

                # Save the triangle centroid and normal information
                triangle_file.write(f'{centroid_of_triangle[0]:.4f} {centroid_of_triangle[1]:.4f} {centroid_of_triangle[2]:.4f} '
                                    f'{normal[0]:.4f} {normal[1]:.4f} {normal[2]:.4f}\n')

# Print confirmation messages indicating the completion of file writing
print(f'Tetrahedron centroids and volumes saved to {tetra_output_filename}')
print(f'Triangle centroids and normals saved to {triangle_output_filename}')
