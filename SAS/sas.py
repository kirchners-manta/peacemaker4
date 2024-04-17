# Let's first read the contents of the van der Waals radii file and the XYZ file

vdw_radii_path = 'vdw_radii.txt'
xyz_file_path = 'c1m1w4-8.xyz'

# Load van der Waals radii
with open(vdw_radii_path, 'r') as file:
    vdw_radii = eval(file.read())

# Load the XYZ file
with open(xyz_file_path, 'r') as file:
    xyz_data = file.readlines()

import numpy as np

# Parse the XYZ data to extract atom types and coordinates
atom_types = []
coordinates = []

for line in xyz_data[2:]:  # Skip the first two lines (header and energy)
    parts = line.split()
    if len(parts) == 4:
        atom_type, x, y, z = parts
        atom_types.append(atom_type)
        coordinates.append([float(x), float(y), float(z)])

# Convert list of coordinates to numpy array
coordinates = np.array(coordinates)

# Associate each atom type with its van der Waals radius
radii = [vdw_radii[atom] if atom in vdw_radii else 0.0 for atom in atom_types]

# Convert radii list to a numpy array for element-wise operations
radii_np = np.array(radii)

# Generate points for the solvent-accessible surface
def generate_sphere_points(n_points=100):
    points = np.zeros((n_points, 3))
    inc = np.pi * (3 - np.sqrt(5))
    off = 2 / n_points
    for k in range(n_points):
        y = k * off - 1 + (off / 2)
        r = np.sqrt(1 - y**2)
        phi = k * inc
        points[k] = [np.cos(phi) * r, y, np.sin(phi) * r]
    return points

# Function that calculated the distance from the points on the sphere to the atoms
def cdist(sphere_points, coordinates):
    return np.sqrt(((sphere_points[:, None] - coordinates)**2).sum(axis=2))



# Redefine the function to use the numpy array for radii
def shrake_rupley_algorithm(atom_types, coordinates, radii_np, n_points=100000):
    #Calculate the SASA using the Shrake-Rupley algorithm
    
    total_sasa = 0.0
    points = generate_sphere_points(n_points)
    
    # Calculate SASA for each atom
    for i, (center, radius) in enumerate(zip(coordinates, radii_np)):
        # Scale points to the surface of the sphere including probe radius
        sphere_points = center + points * radius 
        
        # Calculate distances from these points to the centers of all other atoms
        distances = cdist(sphere_points, coordinates)
        accessible = np.all((distances > radii_np) | (np.isclose(distances, 0)), axis=1)
        
        # Calculate the fraction of points that are accessible
        accessible_fraction = np.sum(accessible) / n_points
        
        # Surface area of a sphere: 4 * pi * r^2
        atom_sasa = 4 * np.pi * radius**2 * accessible_fraction
        total_sasa += atom_sasa        
    
    return total_sasa

# Calculate the SASA again using the corrected function
sasa_corrected = shrake_rupley_algorithm(atom_types, coordinates, radii_np)
print(sasa_corrected)



# Create a simple Plotly visualization of atoms
import plotly.graph_objects as go
from scipy.spatial import Delaunay

# Assign a color to each atom type
atom_colors = {'Ac': 'grey',
               'Al': 'grey',
               'Am': 'grey',
               'Sb': 'grey',
               'Ar': 'grey',
               'As': 'grey',
               'At': 'grey',
               'Ba': 'grey',
               'Bk': 'grey',
               'Be': 'grey',
               'Bi': 'grey',
               'Bh': 'grey',
               'B' : 'grey',
               'Br': 'grey',
               'Cd': 'grey',
               'Ca': 'grey',
               'Cf': 'grey',
               'C' : 'darkslategray',
               'Ce': 'grey',
               'Cs': 'grey',
               'Cl': 'green',
               'Cr': 'grey',
               'Co': 'grey',
               'Cu': 'grey',
               'Cm': 'grey',
               'Db': 'grey',
               'Dy': 'grey',
               'Es': 'grey',
               'Er': 'grey',
               'Eu': 'grey',
               'Fm': 'grey',
               'F' : 'grey',
               'Fr': 'grey',
               'Gd': 'grey',
               'Ga': 'grey',
               'Ge': 'grey',
               'Au': 'grey',
               'Hf': 'grey',
               'Hs': 'grey',
               'He': 'grey',
               'Ho': 'grey',
               'H' : 'snow',
               'In': 'grey',
               'I' : 'grey',
               'Ir': 'grey',
               'Fe': 'grey',
               'Kr': 'grey',
               'La': 'grey',
               'Lr': 'grey',
               'Pb': 'grey',
               'Li': 'grey',
               'Lu': 'grey',
               'Mg': 'grey',
               'Mn': 'grey',
               'Mt': 'grey',
               'Md': 'grey',
               'Hg': 'grey',
               'Mo': 'grey',
               'Nd': 'grey',
               'Ne': 'grey',
               'Np': 'grey',
               'Ni': 'grey',
               'Nb': 'grey',
               'N' : 'grey',
               'No': 'grey',
               'Os': 'grey',
               'O' : 'red',
               'Pd': 'grey',
               'P' : 'grey',
               'Pt': 'grey',
               'Pu': 'grey',
               'Po': 'grey',
               'K' : 'grey',
               'Pr': 'grey',
               'Pm': 'grey',
               'Pa': 'grey',
               'Ra': 'grey',
               'Rn': 'grey',
               'Re': 'grey',
               'Rh': 'grey',
               'Rb': 'grey',
               'Ru': 'grey',
               'Rf': 'grey',
               'Sm': 'grey',
               'Sc': 'grey',
               'Sg': 'grey',
               'Se': 'grey',
               'Si': 'grey',
               'Ag': 'grey',
               'Na': 'grey',
               'Sr': 'grey',
               'S' : 'grey',
               'Ta': 'grey',
               'Tc': 'grey',
               'Te': 'grey',
               'Tb': 'grey',
               'Tl': 'grey',
               'Th': 'grey',
               'Tm': 'grey',
               'Sn': 'grey',
               'Ti': 'grey',
               'W' : 'grey',
               'U' : 'grey',
               'V' : 'grey',
               'Xe': 'grey',
               'Yb': 'grey',
               'Y' : 'grey',
               'Zn': 'grey',
               'Zr': 'grey',
               'Xx': 'grey',}
            
# Map each atom type to its corresponding color
colors_mapped = [atom_colors[atom] for atom in atom_types]

# Function to visualize the molecular structure with SASA points
def visualize_sasa(atom_types, coordinates, radii_np, colors_mapped, n_points=5000):
    points = generate_sphere_points(n_points)
    
    # Store all SASA points for visualization
    sasa_points = []
    sasa_colors = []
    
    # Calculate SASA for each atom
    for i, (center, radius) in enumerate(zip(coordinates, radii_np)):
        # Scale points to the surface of the sphere including probe radius
        sphere_points = center + points * radius 
        
        # Calculate distances from these points to the centers of all other atoms
        distances = cdist(sphere_points, coordinates)
        accessible = np.all((distances > radii_np) | (np.isclose(distances, 0)), axis=1)
        
        # Collect accessible points for visualization
        for point, is_accessible in zip(sphere_points, accessible):
            if is_accessible:
                sasa_points.append(point)
                sasa_colors.append('black') 
    
    return np.array(sasa_points), sasa_colors

# Calculate SASA points for visualization
sasa_points, sasa_colors = visualize_sasa(atom_types, coordinates, radii_np, colors_mapped)

# Prepare traces for Plotly visualization
atom_trace = go.Scatter3d(
    x=coordinates[:, 0],
    y=coordinates[:, 1],
    z=coordinates[:, 2],
    mode='markers',
    marker=dict(
        size=radii_np * 20,  # Scale radii for visibility
        color=colors_mapped,  # Use mapped colors
        opacity=0.9
    ),
    name='Atoms'
)

sasa_trace = go.Scatter3d(
    x=sasa_points[:, 0],
    y=sasa_points[:, 1],
    z=sasa_points[:, 2],
    mode='markers',
    marker=dict(
        size=2,  # Small size for surface points
        color=sasa_colors,  # Colors corresponding to atom types
        opacity=0.3
    ),
    name='SASA Points'
)

# Define layout settings
layout = go.Layout(
    title="Molecular Structure with SASA Visualization",
    scene=dict(
        xaxis=dict(title="X-axis"),
        yaxis=dict(title="Y-axis"),
        zaxis=dict(title="Z-axis")
    )
)

# Create the figure and display it
fig = go.Figure(data=[atom_trace, sasa_trace], layout=layout)
fig.show()
