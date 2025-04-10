import sys, os
for d in os.listdir( "build" ):
    sys.path.append( os.path.join( "build", d ) )
sys.path.append( "." )

from pysdot.domain_types import ConvexPolyhedraAssembly
from pysdot.radial_funcs import CompressibleFunc
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix
import pyvista as pv
from pysdot import OptimalTransport
from pysdot import PowerDiagram
import numpy as np
import matplotlib.pyplot as plt

def seed_transform(Z):
    """
    Given Z, an (N,2) array where each row is [z1, z2],
    compute and return an (N,2) array where the i-th row is
    y^i = (1/(2*z2)) * [z1, -1].
    """
    # Extract z1 and z2 from each row
    z1 = Z[:, 0]
    z2 = Z[:, 1]
    # Compute the two components for each row:
    y1 = z1 / (2 * z2)
    y2 = -1 / (2 * z2)
    return np.column_stack((y1, y2))

def weight_transform(w, Z, f=1.0, c_p=1.0, Pi_0=1.0):
    """
    Given w (an (N,) array) and Z (an (N,2) array, each row [z1, z2]),
    compute for each i:
    
        psi^i = w^i + (z1^i/(2*z2^i))^2 + (1/(2*z2^i))^2 - (f^2/(2*z2^i))*(z1^i)^2 + c_p * Pi_0
    
    Returns an array of the same shape as w.
    The parameters f, c_p, and Pi_0 can be adjusted or passed in.
    """
    # Extract z1 and z2 from Z
    z1 = Z[:, 0]
    z2 = Z[:, 1]
    
    term1 = (z1 / (2 * z2))**2
    term2 = (1 / (2 * z2))**2
    term3 = (f**2 / (2 * z2)) * (z1**2)
    
    psi = w + term1 + term2 - term3 + c_p * Pi_0
    return psi

def weight_transform_inv(psi, Z, f=1.0, c_p=1.0, Pi_0=1.0):
    """
    Given w (an (N,) array) and Z (an (N,2) array, each row [z1, z2]),
    compute for each i:
    
        w^i = psi^i - (z1^i/(2*z2^i))^2 - (1/(2*z2^i))^2 + (f^2/(2*z2^i))*(z1^i)^2 - c_p * Pi_0
    
    Returns an array of the same shape as w.
    The parameters f, c_p, and Pi_0 can be adjusted or passed in.
    """
    # Extract z1 and z2 from Z
    z1 = Z[:, 0]
    z2 = Z[:, 1]
    
    term1 = (z1 / (2 * z2))**2
    term2 = (1 / (2 * z2))**2
    term3 = (f**2 / (2 * z2)) * (z1**2)

    w = psi - term1 - term2 + term3 - c_p * Pi_0
    return w

def expmap_inverse_vec(points, y, f, g):
    """
    points: numpy array of shape (N, 2), where N is the number of boundary points.
    y: a 2-element array representing the mapping point [y1, y3].
    f, g: scalar parameters.
    Returns an array of mapped points of shape (N, 2).
    """
    # Extract columns from points array.
    p1 = points[:, 0]
    p3 = points[:, 1]
    
    # Unpack the mapping point.
    y1, y3 = y
    
    # Compute the mapped coordinates.
    x1 = (f ** 2 / y3) * (p1 - y1)
    x3 = (g / y3 ** 2) * p3 + (f ** 2 / (2 * y3 ** 2)) * (p1 - y1) ** 2 
    return np.stack((x1, x3), axis = 1)

# First we create a triangulation of a square
# Nx = 500
# Ny = 500
# Lx = 1 # width of the fundamental domain
# Ly = 1
# x = np.linspace(0, 1, Nx)
# y = np.linspace(0, 1, Ny)
# X, Y = np.meshgrid(x, y)
# points = np.array([X.flatten(), Y.flatten()]).T
# tri = Delaunay(points)
# distorted_points = expmap_inverse_vec(points, np.array([0, 1]), 1, 1)

# --- Plotting the triangulation ---

# fig, ax = plt.subplots()
# # Draw the triangulation edges
# plt.triplot(distorted_points[:, 0], distorted_points[:, 1], tri.simplices, 'k-', lw=0.5)
# # Draw the vertices as small circles
# plt.plot(distorted_points[:, 0], distorted_points[:, 1], 'o', markersize=2)
# plt.title("Distorted Triangulation")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.axis('equal')
# plt.show()

# domain
domain = ConvexPolyhedraAssembly()
domain.add_box( [0, -1], [2, 1] )
# numTri = np.shape(tri.simplices)[0] # number of triangles in the triangulation

# Add each triangle to the domain one by one
# for k in range(numTri):

#     T = tri.simplices[k] # indices of vertices in the triangle (ordered clockwise)
#     p = distorted_points[T,:] # coordinates of vertices in the triangle
#     domain.add_simplex(p) # add the triangle to the domain

Z = np.array([[1, 1], [1, 1.5]]);#, [0.5, 0.5], [1.5, 0.5]])
N = len(Z)
w0 = np.zeros(N) + 1

Y = seed_transform(Z)
psi0 = weight_transform(w0, Z)

err_tol = (1e-3 / 100) * (1 / N)

ot = OptimalTransport( positions = Y, weights = psi0, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = 1, gamma = 1.41, g = 1, f_cor = 1, pi_0 = 1, c_p = 1 ))
ot.set_stopping_criterion(err_tol, 'max delta masses')

print( "Pre-masses: ", ot.pd.integrals() )
mvs = ot.pd.der_integrals_wrt_weights()
m = csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets))
print( m.todense() )

ot.adjust_weights()
psi = ot.get_weights()
w = weight_transform_inv(psi, Z)

print( "Post-masses: ", ot.pd.integrals() )
print( "Internal Energy: ", ot.pd.internal_energy() )
print( "Centroids:", ot.pd.centroids() )
print( "Optimized weights: ", w )
