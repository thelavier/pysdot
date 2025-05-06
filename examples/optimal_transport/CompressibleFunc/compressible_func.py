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
from itertools import combinations
from collections import defaultdict

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
    x1 = points[:, 0]
    x3 = points[:, 1]
    
    # Unpack the mapping point.
    y1, y3 = y
    
    # Compute the mapped coordinates.
    p1 = (f ** 2 / y3) * (x1 - y1)
    p3 = (g / y3 ** 2) * x3 + (f ** 2 / (2 * y3 ** 2)) * (x1 - y1) ** 2 
    return np.stack((p1, p3), axis = 1)

def finite_difference_hessian(Y, w, domain, kappa, gamma, g, f, Pi_0, c_p, epsilon=1e-4):
    N = len(w)
    hessian_fd = np.zeros((N, N))

    for i in range(N):
        for j in range(i, N):
            if i == j:
                # Perturb forward and backward for w[i]
                w_perturbed_forward = np.copy(w)
                w_perturbed_backward = np.copy(w)

                w_perturbed_forward[i] += epsilon
                w_perturbed_backward[i] -= epsilon

                otf = OptimalTransport( positions = Y, weights = w_perturbed_forward, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p ))
                otb = OptimalTransport( positions = Y, weights = w_perturbed_backward, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p ))

                g_perturbed_forward = 1 / N - otf.pd.integrals()
                g_perturbed_backward = 1 / N - otb.pd.integrals()
                
                hessian_fd[i, i] = (g_perturbed_forward[i] - g_perturbed_backward[i]) / (2 * epsilon)
            else:
                # Off-diagonal elements: mixed partial derivatives w.r.t. w[i] and w[j]
                w_perturbed_forward_j = np.copy(w)
                w_perturbed_backward_j = np.copy(w)

                # Perturb in the j-th direction
                w_perturbed_forward_j[j] += epsilon
                w_perturbed_backward_j[j] -= epsilon

                otfj = OptimalTransport( positions = Y, weights = w_perturbed_forward_j, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p ))
                otbj = OptimalTransport( positions = Y, weights = w_perturbed_backward_j, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p ))

                # Compute g for both perturbed states in the j-th direction
                g_forward = 1 / N - otfj.pd.integrals()
                g_backward = 1 / N - otbj.pd.integrals()

                # First derivative of g[i] with respect to w[j]
                hessian_fd[i, j] = (g_forward[i] - g_backward[i]) / (2 * epsilon)

                # Since Hessian is symmetric
                hessian_fd[j, i] = hessian_fd[i, j]
    return hessian_fd

# System Parameters
f = 1 
g = 1 #10
c_p = 1 #1003.5
Pi_0 = 1 #0.864
gamma = 1.41
kappa = 1 #65.00526358589575

# First we create a triangulation of a square
Nx = 50
Ny = 50
Lx = 2 # width of the fundamental domain
Ly = 2
x = np.linspace(0, 2, Nx)
y = np.linspace(-1, 1, Ny)
X, Y = np.meshgrid(x, y)
points = np.array([X.flatten(), Y.flatten()]).T
tri = Delaunay(points)
distorted_points = expmap_inverse_vec(points, np.array([0, 1]), f, g)

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
# domain.add_box( [0, -1], [2, 1] )
numTri = np.shape(tri.simplices)[0] # number of triangles in the triangulation

# Add each triangle to the domain one by one
for T in tri.simplices:
    p = distorted_points[T,:] # coordinates of vertices in the triangle
    domain.add_simplex(p) # add the triangle to the domain

Z = np.array([[1, 1]])#, [1, 1.5], [0.5, 0.5], [1.5, 0.5]]) 
# Z = np.array([[1.94666077, 0.96859015], [0.0616024,  0.17076978]])
N = len(Z)

w0 = np.zeros(N)
Y = seed_transform(Z)
psi0 = weight_transform(w0, Z, f = f, c_p = c_p, Pi_0 = Pi_0)

err_tol = (1e-3 / 100) * (1 / N)

ot = OptimalTransport( positions = Y, weights = psi0, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True ), verbosity=0)
ot.set_stopping_criterion(err_tol, 'max delta masses')

# Adding replications based on periodicity
for x in range(-int(False), int(False) + 1):
    if x != 0 :
        ot.pd.add_replication([2 * x, 0])

print( "Pre-masses: ", ot.pd.integrals() )
# mvs = ot.pd.der_integrals_wrt_weights()
# m = csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets))
# Hess = finite_difference_hessian(Y, psi0, domain, kappa, gamma, g, f, Pi_0, c_p)
# print( "Finite Differences Hessian : \n", Hess)
# print( "Numerical Hessian : \n", -m.todense() )

ot.adjust_weights()
psi = ot.get_weights()
w = weight_transform_inv(psi, Z)

print( "Post-masses: ", ot.pd.integrals() )
print( "Internal Energy: ", ot.pd.internal_energy() )
print( "Centroids:", ot.pd.centroids() )
print( "Optimized weights: ", w )

filename = 'pb.vtk'
ot.pd.display_vtk( filename )

# Identify the vertices of the cell
grid = pv.read(filename)
grid = grid.clean(tolerance=1e-8)

# Grab triangle→cell IDs
cell_ids = grid.cell_data['num'].astype(int)

# Extract the two triangle-submeshes for cell 0 and cell 1
cells0 = np.where(cell_ids == 0)[0]
cells1 = np.where(cell_ids == 1)[0]
region0 = grid.extract_cells(cells0)
region1 = grid.extract_cells(cells1)

# Build the set of **domain** boundary edges (so we can ignore them later)
dom_edges = set()
dom_bnd = grid.extract_feature_edges(
    boundary_edges=True,
    feature_edges=False,
    manifold_edges=False,
    non_manifold_edges=False
)
# `dom_bnd.lines` is an Nx3 array [2, ptId0, ptId1] for each segment
for _, i0, i1 in dom_bnd.lines.reshape(-1, 3):
    p0 = tuple(dom_bnd.points[i0][:2])
    p1 = tuple(dom_bnd.points[i1][:2])
    dom_edges.add(tuple(sorted((p0, p1))))

# Helper to pull out the **interior** boundary edges of a region
def interior_edges(region):
    bnd = region.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False
    )
    segs = []
    for _, i0, i1 in bnd.lines.reshape(-1, 3):
        p0 = tuple(bnd.points[i0][:2])
        p1 = tuple(bnd.points[i1][:2])
        e = tuple(sorted((p0, p1)))
        # drop any edges that lie on the *domain* boundary
        if e not in dom_edges:
            segs.append(e)
    return set(segs)

edges0 = interior_edges(region0)
edges1 = interior_edges(region1)

# The **shared** edges between region0 & region1 are exactly their mutual face
shared = edges0 & edges1
if not shared:
    raise RuntimeError("No shared face found between cells 0 & 1!")
    
# Stitch those shared segments into a tiny graph to pick out the two degree-1 nodes
graph = defaultdict(list)
for a, b in shared:
    graph[a].append(b)
    graph[b].append(a)

endpoints = [pt for pt, nbrs in graph.items() if len(nbrs) == 1]
if len(endpoints) != 2:
    raise RuntimeError(f"Expected exactly two end-points, got {endpoints!r}")

print("Boundary between cell 0 ↔ 1 runs from", endpoints[0], "to", endpoints[1])

# Mass of cells
colours=ot.pd.integrals()

# Read the data
grid=pv.read(filename)

# create cell data that gives the cell volumes, this allows us to colour by cell volumes
cell_colours = colours[grid.cell_data['num'].astype(int)]
grid.cell_data['colours']=cell_colours

# plot the data with an automatically created plotter, for a static picture use backend='static'
plotter = pv.Plotter(window_size=[800,800])
plotter.add_mesh(grid) #, clim=[minvel, maxvel])

# Set the camera for 2D view
plotter.camera_position = 'xy'

# convert 2D→3D by dropping z=0
pts2d = np.array(endpoints)                # shape (2,2)
pts3d = np.hstack([pts2d, np.zeros((2,1))])  # shape (2,3)

# make a PyVista point cloud
cloud = pv.PolyData(pts3d)

# add to the existing plotter
plotter.add_mesh(
    cloud,
    color='red',
    point_size=10                     # bump up the radius
)

# Render the frame
plotter.show()