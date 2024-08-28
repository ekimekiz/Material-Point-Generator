"""
STL to Tetrahedral Mesh Volume and Surface Triangle Data Extraction using Gmsh

This script processes an STL file and a corresponding tetrahedral mesh (.msh file) to:
1. Calculate the total volume of the object from the STL file.
2. Generate a tetrahedral mesh using Gmsh and extract the volume and centroids of each tetrahedron.
3. Cross-check surface triangles with their corresponding tetrahedrons to identify which tetrahedrons have surface elements.
4. Output the volume, centroid, and surface triangle coordinates for each tetrahedron into a text file.

Output:
    The script will generate a text file named 'file_name_mesh_tetra_surface_data.txt' containing the following data for each tetrahedron in the mesh:
    1. **Volume:** The volume of each tetrahedron.
    2. **Centroid:** The location of the centroid of the tetrahedron, saved as x, y, z coordinates.
    3. **Surface Triangles:** The faces of the tetrahedral elements that lie on the surface of the object. The data for these triangles is saved next to the centroid data as the coordinates of the three corners of each triangle in the format:
       - Corner 1: x, y, z
       - Corner 2: x, y, z
       - Corner 3: x, y, z

    A tetrahedron can at most have three surface faces for a closed object made of multiple elements.
    Therefore each line will have entries ranging from 4 (volume, centroid_x, centroid_y, centroid_z) to 31 (3 triangles, 3 nodes for each triangle, x, y, z for each node).

Note:
    Ensure that Gmsh is installed and accessible in your system's PATH before running this script.
    The script expects that the STL file and the corresponding .msh file have been generated beforehand.

Usage:
    Run the script using the following command:

        python script_name.py <stl_filename> <msh_filename>

    Replace `<stl_filename>` and `<msh_filename>` with the appropriate file names, including their extensions.

Example:
    If you have an STL file named "file_name.stl" and a corresponding mesh file "file_name_mesh.msh", use the command:

        python script_name.py file_name.stl file_name_mesh.msh

Dependencies:
    - Python 3.x
    - numpy
    - gmsh
    - numpy-stl (install using pip: pip install numpy-stl)

The script will output a text file named 'file_name_mesh_tetra_surface_data.txt' containing the volume, centroid, and surface triangle data for each tetrahedron in the mesh.
Volume is the volume of each tetrahedron, next to it the location of the centroid of the tetrahedron is saved as x, y, z coordinates.
Surface triangles are the faces of the tetrahedral elements that lie on the surface of the object. Their data is saved next to the centroid data as the coordinates of the three corners of each triangle: Corner 1 x, y, z; Corner 2 x, y, z.
"""

import os
import sys
import gmsh
import numpy as np
from stl import mesh

# Ensure the script runs in the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script_name.py <stl_filename> <msh_filename>")
    sys.exit(1)

# Get the STL and mesh file names from the command-line arguments
stl_filename = sys.argv[1]
msh_filename = sys.argv[2]

# Load the STL file
your_mesh = mesh.Mesh.from_file(stl_filename)

# Function to calculate the signed volume of a triangle
def signed_volume_of_triangle(v1, v2, v3):
    return np.dot(v1, np.cross(v2, v3)) / 6.0

# Initialize the volume of the object from the STL file
stl_volume = 0.0

# Iterate through each triangle in the mesh and accumulate volume
for i in range(len(your_mesh)):
    v1, v2, v3 = your_mesh.vectors[i]
    stl_volume += signed_volume_of_triangle(v1, v2, v3)

# Take the absolute value since the volume should be positive
stl_volume = abs(stl_volume)

# Initialize the Gmsh API
gmsh.initialize()

# Create a new model for visualization
gmsh.model.add("Mesh and Surface Triangle Data")

# Load the mesh from the saved .msh file
gmsh.merge(msh_filename)

# Generate mesh data
gmsh.model.mesh.generate(3)

# Get the nodes and elements
element_types, element_tags, node_tags = gmsh.model.mesh.getElements()

# Function to calculate the centroid of a tetrahedron
def calculate_centroid(coords):
    return np.mean(coords, axis=0)

# Function to calculate the signed volume of a tetrahedron given its vertices
def calculate_signed_volume(coords):
    v1, v2, v3, v4 = coords
    return np.dot(v1 - v4, np.cross(v2 - v4, v3 - v4)) / 6.0

# Initialize a list to store the mapping between 2D elements and tetrahedra
surface_to_tetra_map = {}

# Initialize a variable to accumulate the total signed volume
total_signed_volume = 0.0

# Extract 2D and 3D elements separately
tetra_elements = []
surface_elements = []

for etype, etags, n_tags in zip(element_types, element_tags, node_tags):
    if etype == 4:  # Gmsh type 4 corresponds to tetrahedral elements
        for i in range(len(etags)):
            nodes = n_tags[i * 4:(i + 1) * 4]
            if len(nodes) == 4:
                coords = [np.array(gmsh.model.mesh.getNode(node)[0]) for node in nodes]
                volume = calculate_signed_volume(coords)
                total_signed_volume += volume
                tetra_elements.append({
                    'tetra_tag': etags[i],
                    'nodes': tuple(sorted(nodes)),
                    'coords': coords,
                    'volume': volume
                })
    elif etype == 2:  # Gmsh type 2 corresponds to triangular elements (surface)
        for i in range(len(etags)):
            nodes = n_tags[i * 3:(i + 1) * 3]
            surface_elements.append({
                'surface_tag': etags[i],
                'nodes': tuple(sorted(nodes)),
                'coords': [np.array(gmsh.model.mesh.getNode(node)[0]) for node in nodes]
            })

# Cross-check each surface triangle to determine which tetrahedron it belongs to
for surface in surface_elements:
    surface_nodes = surface['nodes']
    for tetra in tetra_elements:
        tetra_nodes = tetra['nodes']
        faces = [
            tuple(sorted([tetra_nodes[0], tetra_nodes[1], tetra_nodes[2]])),
            tuple(sorted([tetra_nodes[0], tetra_nodes[1], tetra_nodes[3]])),
            tuple(sorted([tetra_nodes[0], tetra_nodes[2], tetra_nodes[3]])),
            tuple(sorted([tetra_nodes[1], tetra_nodes[2], tetra_nodes[3]])),
        ]
        if surface_nodes in faces:
            if tetra['tetra_tag'] not in surface_to_tetra_map:
                surface_to_tetra_map[tetra['tetra_tag']] = {
                    'volume': tetra['volume'],
                    'centroid': calculate_centroid(tetra['coords']),
                    'surface_coords': []
                }
            surface_to_tetra_map[tetra['tetra_tag']]['surface_coords'].append(surface['coords'])

# Check and print if any tetrahedra have more than one surface triangle
for tetra_tag, data in surface_to_tetra_map.items():
    if len(data['surface_coords']) > 1:
        print(f"Tetrahedron {tetra_tag} has {len(data['surface_coords'])} surface triangles.")

# Save the combined tetrahedron volumes, centroids, and surface triangle node coordinates
output_filename = f"{os.path.splitext(msh_filename)[0]}_tetra_surface_data.txt"
with open(output_filename, 'w') as f:
    for tetra in tetra_elements:
        volume = tetra['volume']
        centroid = calculate_centroid(tetra['coords'])
        tetra_tag = tetra['tetra_tag']

        f.write(f"{volume} {centroid[0]} {centroid[1]} {centroid[2]}")

        if tetra_tag in surface_to_tetra_map:
            surface_triangles = surface_to_tetra_map[tetra_tag]['surface_coords']
            for triangle in surface_triangles:
                for node in triangle:
                    f.write(f" {node[0]} {node[1]} {node[2]}")
        f.write("\n")

# Print the number of tetrahedral elements generated
print(f"Total number of tetrahedral elements generated: {len(tetra_elements)}")


# Print the volume calculated from the STL mesh as a check
print(f"Volume calculated from STL: {stl_volume}")

# Finalize the Gmsh API
gmsh.finalize()

# Sum the absolute value of the signed volume to avoid orientation issues
total_absolute_volume = sum(abs(tetra['volume']) for tetra in tetra_elements)

# Print the total absolute volume
print(f"Total absolute volume of all tetrahedra: {total_absolute_volume}")

# Calculate the difference with the STL volume
print(f"Difference between STL volume and absolute tetrahedral volume: {stl_volume - total_absolute_volume}")

print(f"Tetrahedron centroids and surface triangles saved to ‘{output_filename}’")
