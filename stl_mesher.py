"""
STL to 3D Mesh Generator using Gmsh

This script processes an STL file to generate a 3D tetrahedral mesh using Gmsh.
It also calculates and compares the volume from the original STL file with the volume of the generated mesh to check the accuracy of the mesh.

Output:
    The script generates a .msh file with the same name as the input STL file, containing the tetrahedral mesh.
    It also prints the total volume calculated from the STL and the generated mesh, along with their difference.

Note:
    Ensure that Gmsh is installed and accessible in your system's PATH before running this script.
    
Usage:
    python stl_mesher.py <STL_FILE> [options]

Arguments:
    <STL_FILE>             The path to the STL file to be processed.

Options:
    -char_length=<value>   Set the characteristic length for mesh generation. This determines the size of the mesh elements.
                           Default is 1.0. The smaller the value, the finer the mesh.
                           The effect of this value depends on the size of the object in the STL file. A value of 1.0 may be too small
                           for an object with a size of 1000 units.
                           STL files do not store size units, so it depends on how the object was created.
                           (e.g., a 1m x 1m x 1m cube created in CAD software with default units in mm would have dimensions
                           of 1000 x 1000 x 1000 units in the STL file).
                           If mesh generation takes too long, try increasing the characteristic length.
                           If the characteristic length is too large, Gmsh will produce the largest elements it can without error.

    -refinement=<level>    Set the mesh refinement level. Higher levels generate finer meshes.
                           Default is 0 (no refinement).

Examples:
    1. Generate a mesh with default settings:
       python stl_mesher.py model.stl

    2. Generate a mesh with a specific characteristic length of 0.5:
       python stl_mesher.py model.stl -char_length=0.5

    3. Generate a mesh with a specific characteristic length of 0.5 and a refinement level of 2:
       python stl_mesher.py model.stl -char_length=0.5 -refinement=2
"""

import os
import gmsh
import sys
import numpy as np
from stl import mesh

# Ensure the script runs in the directory where the script is located
# This ensures that any file paths provided will be relative to the script's location.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Check for the STL filename passed as a command-line argument
# The script expects the user to provide the STL file as the first command-line argument.
if len(sys.argv) < 2 or sys.argv[1].startswith('-'):
    print("Error: Please provide the STL file as a command-line argument.")
    sys.exit(1)

# Assign the provided STL filename to a variable for use in the script
stl_filename = sys.argv[1]

# Set default characteristic length and refinement level
# The characteristic length controls the mesh resolution, and the refinement level determines how many times the mesh will be refined.
uniform_size = 1.0  # Default characteristic length
refinement_level = 0  # Default refinement level

# Parse command-line arguments for optional settings
# Users can override the default characteristic length and refinement level using additional command-line arguments.
for arg in sys.argv[2:]:
    if arg.startswith('-char_length='):
        uniform_size = float(arg.split('=')[1])
    if arg.startswith('-refinement='):
        refinement_level = int(arg.split('=')[1])

# Load the STL file using the numpy-stl library
your_mesh = mesh.Mesh.from_file(stl_filename)

# Extract all the unique vertices from the mesh
all_vertices = your_mesh.vectors.reshape(-1, 3)
unique_vertices = np.unique(all_vertices, axis=0)

# Function to calculate the signed volume of a triangle
# This function is used to compute the volume of the mesh by summing the volumes of all triangles in the STL file.
def signed_volume_of_triangle(v1, v2, v3):
    return np.dot(v1, np.cross(v2, v3)) / 6.0

# Initialize the volume of the object from the STL file
# This will accumulate the volume of each triangle in the STL to get the total volume.
volume = 0.0

# Iterate through each triangle in the mesh and accumulate volume
# The script calculates the total volume of the object by summing the volumes of all the triangles in the STL file.
for i in range(len(your_mesh)):
    v1, v2, v3 = your_mesh.vectors[i]
    volume += signed_volume_of_triangle(v1, v2, v3)

# Take the absolute value since the volume should be positive
# The calculated volume could be negative due to orientation, so we take the absolute value.
volume = abs(volume)

# Initialize the Gmsh API for meshing and further processing
gmsh.initialize()

# Increase the tolerance to allow for more flexibility with overlapping facets
# This setting helps in handling STL files with minor overlaps in facets.
gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap", 180.0)  # Adjust the tolerance as needed

# Enable the automatic removal of duplicate entities
# gmsh.option.setNumber("Mesh.RemoveAllDuplicates", 1)

# Load the STL file into Gmsh for meshing
gmsh.model.add("STL_Model")
gmsh.merge(stl_filename)

# Set a uniform characteristic length across the entire geometry
# This determines the size of the elements in the generated mesh.
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", uniform_size)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", uniform_size)


# Create a surface loop and a volume from the STL surface
# Gmsh requires surfaces to be defined as loops to generate volumes.
surfaces = gmsh.model.getEntities(2)  # Get all surface entities
if len(surfaces) == 0:
    raise ValueError("No surfaces found in the STL file.")

# Create a volume by defining a surface loop
surface_loop = gmsh.model.geo.addSurfaceLoop([surf[1] for surf in surfaces])
volume_entity = gmsh.model.geo.addVolume([surface_loop])

# Synchronize the CAD kernel
# This ensures all geometrical operations are updated before mesh generation.
gmsh.model.geo.synchronize()

# Generate the 3D mesh (tetrahedrons)
# This generates a mesh from the geometry, which is essential for further numerical analysis.
gmsh.model.mesh.generate(3)

# Get the nodes and translate them to ensure proper alignment
# This retrieves the coordinates of all nodes in the mesh.
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = np.array(node_coords).reshape(-1, 3)

# Find the minimum coordinates along each axis
# This helps in translating the mesh so that it starts at the origin.
min_coords = np.min(node_coords, axis=0)

# Translate the coordinates so that the minimum corner is at the origin
# Shifting the mesh to the origin can be useful for certain simulations or analyses.
shifted_coords = node_coords - min_coords

# Update the node coordinates in the mesh with the shifted coordinates
# This updates the mesh in Gmsh with the newly shifted node coordinates.
for i, node_tag in enumerate(node_tags):
    gmsh.model.mesh.setNode(node_tag, shifted_coords[i], [])

# Apply mesh refinement according to the specified level (if any)
# Users can request multiple levels of refinement to make the mesh finer.
if refinement_level > 0:
    for _ in range(refinement_level):
        gmsh.model.mesh.refine()

# Save the generated mesh to a .msh file named after the STL file
# The mesh is saved with a filename that matches the input STL file but with a .msh extension.
mesh_filename = f"{os.path.splitext(stl_filename)[0]}_mesh.msh"
gmsh.write(mesh_filename)

# Get the elements and nodes from the mesh
# This retrieves the mesh elements, particularly tetrahedrons, for further processing.
element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim=3)

# Function to calculate the signed volume of a tetrahedron
# This is used to calculate the volume of each tetrahedron in the mesh.
def signed_volume_of_tetrahedron(v1, v2, v3, v4):
    return np.abs(np.dot(v1 - v4, np.cross(v2 - v4, v3 - v4))) / 6.0

# Initialize a counter for tetrahedral elements and total volume
# These variables will keep track of the total number of tetrahedrons and their combined volume.
tetrahedral_count = 0
total_tetrahedral_volume = 0.0

# Loop through the tetrahedral elements to calculate their volumes
# The script calculates the volume of each tetrahedron and sums them to get the total volume.
for i in range(len(node_tags[0]) // 4):  # Each tetrahedron has 4 nodes
    # Extract the nodes corresponding to the current tetrahedron
    tetra_nodes = node_tags[0][i * 4:(i + 1) * 4]
    
    # Ensure you're correctly using the element tag for this tetrahedron
    tetra_tag = element_tags[0][i]  # Correctly aligns element tag with the current tetrahedron

    # Get the coordinates of the 4 nodes
    v1, v2, v3, v4 = [np.array(gmsh.model.mesh.getNode(node)[0]) for node in tetra_nodes]
    
    # Calculate the volume of the tetrahedron
    tetra_volume = signed_volume_of_tetrahedron(v1, v2, v3, v4)
    
    # Accumulate the total tetrahedral volume
    total_tetrahedral_volume += tetra_volume
    
    # Increment the tetrahedral element count
    tetrahedral_count += 1

# Print the number of tetrahedral elements processed
print(f"Total number of tetrahedral elements processed: {tetrahedral_count}")

# Print the total volume of tetrahedral elements
print(f"Total volume of tetrahedral elements: {total_tetrahedral_volume}")

# Compare the calculated volume with the total tetrahedral volume
# This provides a comparison between the volume from the STL file and the volume calculated from the mesh.
print(f"Volume calculated from STL: {volume}")
print(f"Difference between calculated volume and tetrahedral volume: {abs(volume - total_tetrahedral_volume)}")

# Finalize the Gmsh API
# This ensures that Gmsh is properly shut down, freeing any resources it used.
gmsh.finalize()

# Output the mesh filename to the user
# The script informs the user where the output mesh file is saved.
output_filename = f"{os.path.splitext(stl_filename)[0]}_mesh.msh"
print(f"Mesh file saved to '{output_filename}'")
