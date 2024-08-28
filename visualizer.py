import sys
import numpy as np
import os
import matplotlib.pyplot as plt

# Change the working directory to the location of this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python visualizer.py <tetra_filename> <triangle_filename>")
    sys.exit(1)

# Get the filenames from the command-line arguments
tetra_filename = sys.argv[1]
triangle_filename = sys.argv[2]

# Load the tetrahedron centroids and volumes
tetra_data = np.loadtxt(tetra_filename)

# Extract tetrahedron centroids (columns 2-4)
tetra_centroids = tetra_data[:, 1:4]

# Load the triangle centroids and normals
triangle_data = np.loadtxt(triangle_filename)

# Extract triangle centroids (columns 1-3)
triangle_centroids = triangle_data[:, 0:3]

# Extract triangle normals (columns 4-6)
triangle_normals = triangle_data[:, 3:6]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot tetrahedron centroids in blue
ax.plot3D(tetra_centroids[:, 0], tetra_centroids[:, 1], tetra_centroids[:, 2], 'bo', markersize=5, label='Tetrahedra Centroids')

# Plot triangle centroids in red
ax.plot3D(triangle_centroids[:, 0], triangle_centroids[:, 1], triangle_centroids[:, 2], 'ro', markersize=5, label='Triangle Centroids')

# Plot normals at triangle centroids
ax.quiver(triangle_centroids[:, 0], triangle_centroids[:, 1], triangle_centroids[:, 2],
          triangle_normals[:, 0], triangle_normals[:, 1], triangle_normals[:, 2],
          length=0.5, color='k', linewidth=1, label='Normals')

# Add labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Tetrahedron and Triangle Centroids with Normals')
ax.legend(loc='best')

# Set axis properties for better visualization
ax.set_box_aspect([1, 1, 1])  # aspect ratio is 1:1:1
ax.grid(True)

# Display the plot
plt.show()
