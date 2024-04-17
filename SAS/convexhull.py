import numpy as np
import plotly.graph_objects as go
from scipy.spatial import ConvexHull

# Let's read the content of the XYZ file to understand its format
file_path = './../working-code/test-calc/clusters/c1.xyz'

# Read the first few lines of the file to get an idea of its structure
with open(file_path, 'r') as file:
    lines = [next(file) for _ in range(5)]

lines

# Function to read XYZ file and extract atomic coordinates
def read_xyz_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        atom_count = int(lines[0].strip())
        atoms = []
        coordinates = []
        for line in lines[2:2+atom_count]:
            parts = line.split()
            atoms.append(parts[0])  # Atom type
            coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
        coordinates = np.array(coordinates)
    return atoms, coordinates

# Read the coordinates from the XYZ file
atoms, coordinates = read_xyz_file(file_path)

# Calculate centroid of the cluster
centroid = np.mean(coordinates, axis=0)

# Calculate the furthest distance from the centroid to define the radius of the sphere
radius = np.max(np.linalg.norm(coordinates - centroid, axis=1))

# open the sphere using Plotly
fig = go.Figure()

# Add atoms to the figure
fig.add_trace(go.Scatter3d(
    x=coordinates[:, 0],
    y=coordinates[:, 1],
    z=coordinates[:, 2],
    mode='markers',
    marker=dict(size=5, color='blue'),
    name='Atoms'
))

# Add the centroid as a red sphere
fig.add_trace(go.Scatter3d(
    x=[centroid[0]],
    y=[centroid[1]],
    z=[centroid[2]],
    mode='markers',
    marker=dict(size=3, color='red'),
    name='Centroid'
))

# Add a transparent sphere to represent the cluster
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)
x = centroid[0] + radius * np.outer(np.cos(theta), np.sin(phi))
y = centroid[1] + radius * np.outer(np.sin(theta), np.sin(phi))
z = centroid[2] + radius * np.outer(np.ones(np.size(theta)), np.cos(phi))

# Add the convex hull as a mesh3d trace
fig.add_trace(go.Surface(
    x=x,
    y=y,
    z=z,
    colorscale=[[0, 'red'], [1, 'red']],
    opacity=0.4,
    showscale=False,
    name='Sphere'
))

# Set plot layout details
fig.update_layout(
    title="Interactive Molecular Cluster with Sphere",
    scene=dict(
        xaxis_title='X Coordinates',
        yaxis_title='Y Coordinates',
        zaxis_title='Z Coordinates',
        aspectmode='auto'
    )
)

# Show figure
fig.show()

# Calculate the surface area of the sphere
surface_area = 4 * np.pi * radius**2
# Print the surface area
print(f"Surface Area of Sphere: {surface_area:.2f} square units")

# volume of the sphere
volume = (4/3) * np.pi * radius**3
# Print the volume
print(f"Volume of Sphere: {volume:.2f} cubic units")

centroid, radius, surface_area

# Your coordinates array from previous calculations
#coordinates = np.array([...])  # Replace with actual coordinates array

# Calculate the convex hull of the molecular coordinates
hull = ConvexHull(coordinates)

# Save the surface to a xyz file that is compatible with VMD
with open('convex_hull.xyz', 'w') as file:
    file.write(f"{len(hull.points)}\n")
    file.write("Convex Hull\n")
    for point in hull.points:
        file.write(f"C {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")

# Create a Plotly figure for the convex hull
fig = go.Figure()

# Add atoms to the figure
fig.add_trace(go.Scatter3d(
    x=coordinates[:, 0],
    y=coordinates[:, 1],
    z=coordinates[:, 2],
    mode='markers',
    marker=dict(size=5, color='blue'),
    name='Atoms'
))

# Add the convex hull as a mesh3d trace
fig.add_trace(go.Mesh3d(
    x=coordinates[:, 0],
    y=coordinates[:, 1],
    z=coordinates[:, 2],
    i=hull.simplices[:, 0],
    j=hull.simplices[:, 1],
    k=hull.simplices[:, 2],
    color='red',
    opacity=0.4,
    name='Convex Hull'
))

# Set plot layout details
fig.update_layout(
    title="Interactive Molecular Cluster with Convex Hull",
    scene=dict(
        xaxis_title='X Coordinates',
        yaxis_title='Y Coordinates',
        zaxis_title='Z Coordinates',
        aspectmode='auto'
    )
)

# Calculate the surface area of the convex hull
surface_area = hull.area
# Print the surface area
print(f"Surface Area of Convex Hull: {surface_area:.2f} square units")
# Calculate the volume of the convex hull
volume = hull.volume
# Print the volume
print(f"Volume of Convex Hull: {volume:.2f} cubic units")

# Show figure
fig.show()
