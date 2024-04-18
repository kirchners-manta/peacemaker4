#!/usr/bin/env python3

#=========================================================================================
# Peacemaker -- A Quantum Cluster Equilibrium Code.
#
# Copyright 2004-2006 Barbara Kirchner, University of Bonn
# Copyright 2007-2012 Barbara Kirchner, University of Leipzig
# Copyright 2013-2024 Barbara Kirchner, University of Bonn
#
# This file is part of Peacemaker.
#
# Peacemaker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Peacemaker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
#=========================================================================================

# This python sript generates The solvent-accessible surface area (SASA) of a cluster using 
# the Shrake-Rupley algorithm.
# For the generation of the SASA the VdW radii of the atoms are used.
# For each atom a sphere with the corresponding VdW radius is generated. Random points on the
# surface of the sphere are generated and the distance to all other atoms is calculated.
# If the distance is larger than the VdW radius of the other atom, the point is considered, 
# otherwise it is discarded.
# The SASA is then calculated by the fraction of accessible points on the sphere.
# The SASA is visualized using plotly if wished.

#### User Input #################################################################################

# Please enter the path to the directory containing the .xyz files
# path = "/path/to/xyz_files"
#path = "../working-code/test-calc/clusters"
path = "."
# Do you want to visualize the SASA using plotly? (yes/no)
visualize = True

#################################################################################################


import numpy as np
import glob
import plotly.graph_objects as go


# The VdW radii are taken from Bondi's compilation.
vdw_radii = {'Ac': 2.47, 'Al': 1.84, 'Am': 2.44, 'Sb': 2.06, 'Ar': 1.88, 'As': 1.85, 'At': 2.02, 'Ba': 2.68,
             'Bk': 2.44, 'Be': 1.53, 'Bi': 2.07, 'Bh': 0.0, 'B': 1.92, 'Br': 1.85, 'Cd': 2.18, 'Ca': 2.31,
             'Cf': 2.45, 'C': 1.7, 'Ce': 2.42, 'Cs': 3.43, 'Cl': 1.75, 'Cr': 2.06, 'Co': 2.0, 'Cu': 1.96,
             'Cm': 2.45, 'Db': 0.0, 'Dy': 2.31, 'Es': 2.45, 'Er': 2.29, 'Eu': 2.35, 'Fm': 2.45, 'F': 1.47,
             'Fr': 3.48, 'Gd': 2.34, 'Ga': 1.87, 'Ge': 2.11, 'Au': 2.14, 'Hf': 2.23, 'Hs': 0.0, 'He': 1.4,
             'Ho': 2.3, 'H': 1.1, 'In': 1.93, 'I': 1.98, 'Ir': 2.13, 'Fe': 2.04, 'Kr': 2.02, 'La': 2.43,
             'Lr': 2.46, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2.24, 'Mg': 1.73, 'Mn': 2.05, 'Mt': 0.0, 'Md': 2.46,
             'Hg': 2.23, 'Mo': 2.17, 'Nd': 2.39, 'Ne': 1.54, 'Np': 2.39, 'Ni': 1.97, 'Nb': 2.18, 'N': 1.55,
             'No': 2.46, 'Os': 2.16, 'O': 1.52, 'Pd': 2.1, 'P': 1.8, 'Pt': 2.13, 'Pu': 2.43, 'Po': 1.97,
             'K': 2.75, 'Pr': 2.4, 'Pm': 2.38, 'Pa': 2.43, 'Ra': 2.83, 'Rn': 2.2, 'Re': 2.16, 'Rh': 2.1,
             'Rb': 3.03, 'Ru': 2.13, 'Rf': 0.0, 'Sm': 2.36, 'Sc': 2.15, 'Sg': 0.0, 'Se': 1.9, 'Si': 2.1,
             'Ag': 2.11, 'Na': 2.27, 'Sr': 2.49, 'S': 1.8, 'Ta': 2.22, 'Tc': 2.16, 'Te': 2.06, 'Tb': 2.33,
             'Tl': 1.96, 'Th': 2.45, 'Tm': 2.27, 'Sn': 2.17, 'Ti': 2.11, 'W': 2.18, 'U': 2.41, 'V': 2.07,
             'Xe': 2.16, 'Yb': 2.26, 'Y': 2.32, 'Zn': 2.01, 'Zr': 2.23, 'Xx': 0.0}

# Import the clusters (.xyz)
xyz_files = glob.glob( path +"/*.xyz")

# Parse the xyz data to extract the atom types and coordinates
atom_types = []
coordinates = []

# Loop over all .xyz files
for file in xyz_files:
    with open(file, 'r') as f:
        xyz_data = f.readlines()
        # Loop over all lines in the .xyz file
        for line in xyz_data[2:]: 
            parts = line.split()
            if len(parts) == 4:
                atom_type, x, y, z = parts
                atom_types.append(atom_type)
                coordinates.append([float(x), float(y), float(z)])
                
# Convert the atom types and coordinates to numpy arrays
atom_types = np.array(atom_types)
coordinates = np.array(coordinates)

# Associate each atom with its vdw radius
radii = [vdw_radii[atom] if atom in vdw_radii else 0.0 for atom in atom_types]
radii = np.array(radii)



## CALCULATION ##############################################################################

# Generate random points on the surface of a unit sphere
def random_sphere_points(n_points):
    points = np.zeros((n_points, 3))
    # x² + y² + z² = 1
    theta = np.random.uniform(0, 2*np.pi, n_points) # azimuthal angle
    phi = np.random.uniform(0, np.pi, n_points)     # polar angle
    # r = 1
    points[:,0] = np.sin(phi) * np.cos(theta) # x
    points[:,1] = np.sin(phi) * np.sin(theta) # y
    points[:,2] = np.cos(phi)                 # z
    
    return points

# Calculation of the distance between two points
def dist(p1, p2):
    return np.sqrt(((p1[:,None] - p2)**2).sum(axis=2))
    

# Calculate the SASA using the Shrake-Rupley algorithm 
def shrake_rupley(coordinates, radii, n_points=5000):
    
    # Initialize the total SASA
    total_sasa = 0.0
    
    # Generate random points on the surface of a sphere
    points = random_sphere_points(n_points)
    
    # Initialize accessible points and add counter for the accessible points
    accessible_points = []
    
    # Calculation of the SASA for each atom
    for (center, radius) in zip(coordinates, radii):
        
        # Scale the points according to the corresponding VdW radius
        sphere_points = center + points * radius
        
        # Calculate the distance from the points on the sphere to the atoms
        distance = dist(sphere_points, coordinates)
        
        # Check if the sphere point lies outside the VdW radius of the other atom (boolean array)
        accessible = np.all((distance > radii) | (np.isclose(distance, 0)), axis=1)
        
        # Save the accessible points
        accessible_points.append(sphere_points[accessible])
        
        # Calculate the ratio of accessible points
        accessible_ratio = np.sum(accessible) / n_points
        
        # Calculate the surface area of the sphere
        atom_sasa = 4 * np.pi * radius**2 * accessible_ratio
        total_sasa += atom_sasa
        
    return total_sasa, accessible_points

############################################################################################
    
# Calculate the SASA using the Shrake-Rupley algorithm
total_sasa, accessible_points = shrake_rupley(coordinates, radii)

# Print the calculated SASA
print(f'The calculated SASA is: {total_sasa:.2f} Å²')

# Visualize the SASA using plotly (optional) ###############################################

if (visualize == True):

    atom_colors = {'Ac': 'grey', 'Al': 'grey', 'Am': 'grey', 'Sb': 'grey', 'Ar': 'grey', 'As': 'grey', 'At': 'grey', 'Ba': 'grey',
                   'Bk': 'grey', 'Be': 'grey', 'Bi': 'grey', 'Bh': 'grey', 'B': 'grey', 'Br': 'grey', 'Cd': 'grey', 'Ca': 'grey',
                   'Cf': 'grey', 'C': 'black', 'Ce': 'grey', 'Cs': 'grey', 'Cl': 'green', 'Cr': 'grey', 'Co': 'blue', 'Cu': 'copper',
                   'Cm': 'grey', 'Db': 'grey', 'Dy': 'grey', 'Es': 'grey', 'Er': 'grey', 'Eu': 'grey', 'Fm': 'grey', 'F': 'grey',
                   'Fr': 'grey', 'Gd': 'grey', 'Ga': 'grey', 'Ge': 'grey', 'Au': 'grey', 'Hf': 'grey', 'Hs': 'grey', 'He': 'grey',
                   'Ho': 'grey', 'H': 'white', 'In': 'grey', 'I': 'grey', 'Ir': 'grey', 'Fe': 'brown', 'Kr': 'grey', 'La': 'grey',
                   'Lr': 'grey', 'Pb': 'grey', 'Li': 'grey', 'Lu': 'grey', 'Mg': 'grey', 'Mn': 'grey', 'Mt': 'grey', 'Md': 'grey',
                   'Hg': 'grey', 'Mo': 'grey', 'Nd': 'grey', 'Ne': 'grey', 'Np': 'grey', 'Ni': 'grey', 'Nb': 'grey', 'N': 'lightblue',
                   'No': 'grey', 'Os': 'grey', 'O': 'red', 'Pd': 'grey', 'P': 'grey', 'Pt': 'grey', 'Pu': 'grey', 'Po': 'grey',
                   'K': 'grey', 'Pr': 'grey', 'Pm': 'grey', 'Pa': 'grey', 'Ra': 'grey', 'Rn': 'grey', 'Re': 'grey', 'Rh': 'grey',
                   'Rb': 'grey', 'Ru': 'grey', 'Rf': 'grey', 'Sm': 'grey', 'Sc': 'grey', 'Sg': 'grey', 'Se': 'grey', 'Si': 'grey',
                   'Ag': 'grey', 'Na': 'grey', 'Sr': 'grey', 'S': 'yellow', 'Ta': 'grey', 'Tc': 'grey', 'Te': 'grey', 'Tb': 'grey',
                   'Tl': 'grey', 'Th': 'grey', 'Tm': 'grey', 'Sn': 'grey', 'Ti': 'grey', 'W': 'grey', 'U': 'grey', 'V': 'grey',
                   'Xe': 'grey', 'Yb': 'grey', 'Y': 'grey', 'Zn': 'grey', 'Zr': 'grey', 'Xx': 'grey'}

    # Map each atom type to its corresponding color
    color_map = [atom_colors[atom] for atom in atom_types]


    # Visualization of the cluster
    atom_trace = go.Scatter3d(
        x=coordinates[:, 0],
        y=coordinates[:, 1],
        z=coordinates[:, 2],
        mode='markers',
        marker=dict(
            size=radii * 20,  # Scale radii for visibility
            color=color_map,  # Use mapped colors
            opacity=0.9
        ),
        name='Atoms'
    )

    # Visualization of the SASA points
    sasa_trace = go.Scatter3d(
        x=np.concatenate([point[:, 0] for point in accessible_points]),
        y=np.concatenate([point[:, 1] for point in accessible_points]),
        z=np.concatenate([point[:, 2] for point in accessible_points]),
        mode='markers',
        marker=dict(
            size=2,
            color='lightblue',
            opacity=0.3
        ),
        name='SASA points'
    )

    # Define layout settings
    layout = go.Layout(
        title="Cluster with solvent-accessible surface area (SASA)",
        scene=dict(
            xaxis=dict(title="X-axis"),
            yaxis=dict(title="Y-axis"),
            zaxis=dict(title="Z-axis")
        )
    )

    # Create the figure and display it
    fig = go.Figure(data=[atom_trace, sasa_trace], layout=layout)
    fig.show()