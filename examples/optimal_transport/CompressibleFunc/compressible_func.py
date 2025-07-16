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

## Various supporting functions to ensure that we transform to the c-exponential chart where the cost is quadratic

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

# This function now creates and returns the full periodic dataset.
def create_periodic_dataset(Z_initial, psi_initial, f, c_p, Pi_0, Lx, PeriodicX):
    """
    Calculates and returns the full set of diracs, including periodic replicas,
    to be treated as a single, larger non-periodic problem.
    """
    # If not periodic, just transform Z to Y and return the original data.
    if not PeriodicX:
        Y_initial = seed_transform(Z_initial)
        return Y_initial, psi_initial

    # --- Start with the original data ---
    # The OT object needs positions in Y-space, so we transform the initial Z points.
    Y_initial = seed_transform(Z_initial)
    
    # Create lists that we can append replica data to.
    Y_extended = list(Y_initial)
    psi_extended = list(psi_initial)

    # --- Calculate and append replica data ---
    w_base = weight_transform_inv(psi_initial, Z_initial, f, c_p, Pi_0)
    num_original_seeds = Z_initial.shape[0]

    for i in range(num_original_seeds):
        # Loop through periodic images (-1 for left, +1 for right)
        for n in [-1, 1]: 
            T1 = Lx * n # Assuming Lx is the half-width of your domain
            
            # 1. Calculate the true Z-replica position
            z_rep = np.array([Z_initial[i, 0] + T1, Z_initial[i, 1]])
            
            # 2. Calculate the corresponding Y-replica position
            y_rep = seed_transform(z_rep.reshape(1, -1))[0]
            
            # 3. Calculate the corresponding Y-replica weight from the base weight
            psi_rep = weight_transform(w_base[i], z_rep.reshape(1, -1), f, c_p, Pi_0)[0]
            
            # Append the new replica data to our lists
            Y_extended.append(y_rep)
            psi_extended.append(psi_rep)

    # Convert lists back to NumPy arrays and return them.
    return np.array(Y_extended), np.array(psi_extended)

def finite_difference_hessian(Y, w, domain, kappa, gamma, g, f, Pi_0, c_p, mea, res, PeriodicX, epsilon = 1e-4):
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

                if mea == True:
                    otf = OptimalTransport( positions = Y, weights = w_perturbed_forward, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True, Int_res = res ))
                    otb = OptimalTransport( positions = Y, weights = w_perturbed_backward, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True, Int_res = res ))

                    for x in range(-int(PeriodicX), int(PeriodicX) + 1):
                        if x != 0 :
                            otf.pd.add_replication([2 * x, 0])
                            otb.pd.add_replication([2 * x, 0])
                else:
                    otf = OptimalTransport( positions = Y, weights = w_perturbed_forward, masses = np.ones(N) / N, domain = domain)
                    otb = OptimalTransport( positions = Y, weights = w_perturbed_backward, masses = np.ones(N) / N, domain = domain)

                    for x in range(-int(PeriodicX), int(PeriodicX) + 1):
                        if x != 0 :
                            otf.pd.add_replication([2 * x, 0])
                            otb.pd.add_replication([2 * x, 0])

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

                if mea == True:
                    otfj = OptimalTransport( positions = Y, weights = w_perturbed_forward_j, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True, Int_res = res ))
                    otbj = OptimalTransport( positions = Y, weights = w_perturbed_backward_j, masses = np.ones(N) / N, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True, Int_res = res ))
                
                    for x in range(-int(PeriodicX), int(PeriodicX) + 1):
                        if x != 0 :
                            otfj.pd.add_replication([2 * x, 0])
                            otbj.pd.add_replication([2 * x, 0])
                else:
                    otfj = OptimalTransport( positions = Y, weights = w_perturbed_forward_j, masses = np.ones(N) / N, domain = domain)
                    otbj = OptimalTransport( positions = Y, weights = w_perturbed_backward_j, masses = np.ones(N) / N, domain = domain)

                    for x in range(-int(PeriodicX), int(PeriodicX) + 1):
                        if x != 0 :
                            otfj.pd.add_replication([2 * x, 0])
                            otbj.pd.add_replication([2 * x, 0])

                # Compute g for both perturbed states in the j-th direction
                g_forward = 1 / N - otfj.pd.integrals()
                g_backward = 1 / N - otbj.pd.integrals()

                # First derivative of g[i] with respect to w[j]
                hessian_fd[i, j] = (g_forward[i] - g_backward[i]) / (2 * epsilon)

                # Since Hessian is symmetric
                hessian_fd[j, i] = hessian_fd[i, j]
    return hessian_fd

#Define the rescaling function to improve the inital guess for the Damped Newton Solver
def rescale_weights(box, Z, PeriodicX):
    """
    Calculates initial weights for a 2D Laguerre tessellation by rescaling and
    centering the generator points within the domain. This helps ensure a
    good initial guess for the optimal transport solver.

    This function handles edge cases, including single-particle systems or
    collinear particles, without errors.

    Parameters:
        box (list or tuple): Domain bounds in the format [xmin, ymin, xmax, ymax].
        Z (numpy.ndarray): Seed points (shape: n x 2).
        PeriodicX (bool): Flag indicating if the x-dimension is periodic.

    Returns:
        numpy.ndarray: A 1D array of calculated weights, one for each seed point.
    """
    # For this 2D function, we assume the Y-dimension is not periodic.
    PeriodicY = False
    
    # Check for the trivial case of no points
    if Z.shape[0] == 0:
        return np.array([])
        
    min_Z, max_Z = np.min(Z, axis=0), np.max(Z, axis=0)
    lambda_vals = []
    
    # An array indicating which dimensions (x, y) are non-periodic
    non_periodic_dims = ~np.array([PeriodicX, PeriodicY])

    # Calculate the scaling factor lambda, considering only non-periodic dimensions
    for i, periodic in enumerate([PeriodicX, PeriodicY]):
        if not periodic:
            box_width = box[i + 2] - box[i]
            point_spread = max_Z[i] - min_Z[i]

            # --- FIX: Check for zero spread to prevent division by zero ---
            if point_spread > 1e-9:  # Use a small tolerance for floating point safety
                lambda_dim = box_width / point_spread
                lambda_vals.append(lambda_dim)

    # --- FIX: Handle case where points have no spread in any non-periodic dimension ---
    if lambda_vals:
        # Choose the most constrained scaling factor and add a small buffer
        lambda_ = min(lambda_vals) * (1 - 1e-2)
    else:
        # This occurs with a single particle or collinear points.
        # No scaling is necessary or meaningful.
        lambda_ = 1.0

    # Calculate the translation vector t to center the rescaled points
    translation = np.zeros(2)
    for i, periodic in enumerate([PeriodicX, PeriodicY]):
        if not periodic:
            center_of_box = (box[i + 2] + box[i]) / 2
            center_of_points = (min_Z[i] + max_Z[i]) / 2
            center_of_rescaled_points = lambda_ * center_of_points
            translation[i] = center_of_box - center_of_rescaled_points

    # Filter the translation vector for only non-periodic dimensions
    t = translation[non_periodic_dims]

    # Select only the non-periodic coordinate data from Z
    Z_mod = Z[:, non_periodic_dims]
    
    # If Z_mod is empty (all dimensions are periodic), return zero weights
    if Z_mod.shape[1] == 0:
        return np.zeros(Z.shape[0])

    # Calculate the final weights
    w = (1 - lambda_**2) * np.sum(np.square(Z_mod), axis=1) - 2 * lambda_ * np.dot(Z_mod, t)

    return w

## Set the System Parameters

f = 1 
g = 1 #10
c_p = 1 #1003.5
Pi_0 = 1 #0.864
gamma = 2
kappa = 1/2 #65.00526358589575
box = [-1, 0, 1, 1]
PeriodicX = True
PeriodicY = False

## Construct the Domain

domain = ConvexPolyhedraAssembly()

# If working on a box do the following 

Lx, Ly = [box[i+2] - box[i] for i in range(2)]

# Calculate the offset and size for each dimension based on periodicity
# size = [2 * Lx if PeriodicX else box[2], 2 * Ly if PeriodicY else box[3]]
# offset = [-Lx if PeriodicX else box[0], -Ly if PeriodicY else box[1]]

# domain.add_box(offset, size)

# # If working on a distorted domain do the following

Nx = 100
Ny = 100
x = np.linspace(-3, 3, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)
points = np.array([X.flatten(), Y.flatten()]).T
tri = Delaunay(points)
distorted_points = expmap_inverse_vec(points, np.array([0, 1]), f, g)

numTri = np.shape(tri.simplices)[0] # number of triangles in the triangulation

# Add each triangle to the domain one by one
for T in tri.simplices:
    p = distorted_points[T,:] # coordinates of vertices in the triangle
    domain.add_simplex(p) # add the triangle to the domain

    # # Add a periodic copy on the left
    # translation = [[Lx,0.],[Lx,0.],[Lx,0.]]
    # domain.add_simplex(p-translation)

    # # Add a periodic copy on the right
    # domain.add_simplex(p+translation)

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

## Set the seed positions

# Z = np.array([[1, 1], [1, 1.5], [0.5, 0.5], [1.5, 0.5]]) 
Z = np.array([[0.5, 5]])
# N = len(Z)

# NList = (np.rint(np.linspace(2, 250, 10)).astype(int)).tolist()  
NList = [1] 
Errors = []

petsc_opts_1 = {
    'ksp_type': 'gmres',      # Use Conjugate Gradient solver
    'pc_type': 'gamg',   # Use Jacobi (diagonal) preconditioner
    'ksp_monitor': None,   # Print solver progress
    'ksp_rtol': 1e-12      # Set a stricter relative tolerance
}

petsc_opts_2 = {
    'ksp_type': 'gmres',         # More robust than CG
    'pc_type': 'sor',            # Successive Over-Relaxation (a classic, robust choice)
    'ksp_monitor': None,         # Prints the residual at each iteration
    'ksp_converged_reason': None # Prints why the solver stopped (converged or failed)
}

petsc_opts_3 = {
    'ksp_type': 'gmres',
    'pc_type': 'ilu',            # Incomplete LU factorization
    'pc_factor_levels': 2,       # "Fill level": higher is more accurate but more expensive. Try 1, 2, or 3.
    'ksp_monitor': None,
    'ksp_converged_reason': None
}

petsc_opts_4 = {
    'ksp_type': 'cg',                # Go back to cg, as it's good for symmetric systems
    'pc_type': 'gamg',
    'pc_gamg_type': 'classical',     # Try a different aggregation strategy
    'mg_levels_pc_type': 'sor',      # Use SOR as the smoother on each multigrid level
    'ksp_monitor': None,
    'ksp_converged_reason': None
}

petsc_opts_5 = {
    'pc_type': 'lu',
    'pc_factor_mat_solver_type': 'mumps', # A robust parallel direct solver package
    'ksp_type': 'preonly'                # This tells KSP to just apply the PC (which is the direct solve)
}

solver_options = petsc_opts_4

for N in NList:
    print("This test with ", N, " seeds is periodic ?", PeriodicX)
    # print(solver_options)

    # # Generate x values between -1 and 1
    # x = np.random.uniform(-1, 1, N)

    # # Generate y values between 98900 and 102000
    # y = np.random.uniform(0.98900, 1.02000, N)

    # # Combine x and y into a single array of shape (N, 2)
    # Z = np.column_stack((x, y))

    # Transform the seeds (diracs) eto the exponential chart

    print("Initial seeds:", Z)

    Y_initial = seed_transform(Z)

    ## Set an inital guess for the weights

    psi_initial = rescale_weights(box, Y_initial, True) + 1 
    # psi_initial = weight_transform(np.zeros(N) - 11/30, Z, f, c_p, Pi_0)
    # psi0 = np.zeros(N)

    Y, psi0 = create_periodic_dataset(Z, psi_initial, f, c_p, Pi_0, Lx, True)

    print("Initial transformed seeds:", Y)
    print("Initial transformed weights:", psi0)

    target_masses = np.zeros(len(psi0))
    target_masses[:N] = 1.0 / N 

    ## Initialise the optimal transport problem setting the resolution of the integration with Int_res which corresponds to the number of Gaussian Quadrature points 

    ot = OptimalTransport( positions = Y, weights = psi0, masses = target_masses, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, Int = True, Int_res = 10 ), petsc_options=solver_options, verbosity=0)

    ## Set the error tolerance and the stopping criterion 

    err_tol = (1e-3 / 100) * (1 / N)
    ot.set_stopping_criterion(err_tol, 'max delta masses')

    ## Adding replications based on periodicity

    # for x in range(-int(PeriodicX), int(PeriodicX) + 1):
    #     if x != 0 :
    #         ot.pd.add_replication([2 * x, 0])

    ## Check that the optimal transport problem is set up correctly, all cells should initially have positive mass

    print( "Pre-masses: ", ot.pd.integrals() )
    # mvs = ot.pd.der_integrals_wrt_weights()
    # m = csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets))
    # NHess = -m.todense()
    # FDHess = finite_difference_hessian(Y, psi0, domain, kappa, gamma, g, f, Pi_0, c_p, mea = True, res = 100, PeriodicX=False)
    # rel_error = np.linalg.norm(NHess - FDHess) / np.linalg.norm(FDHess) 
    # print("Condition # of FD Hessian: ", np.linalg.cond(FDHess))
    # print("Condition # of Numerical Hessian: ", np.linalg.cond(NHess))
    # Errors.append(rel_error)
    # print("Done", N)
    # print( "Finite Differences Hessian : \n", FDHess)
    # print( "Numerical Hessian : \n", NHess)

    # 6. SOLVE THE OPTIMAL TRANSPORT PROBLEM
    ot.adjust_weights_periodic(Z, f, c_p, Pi_0, Lx) # Pass in Z_initial and other params
    psi_final_full = ot.get_weights()

    # 7. TRANSFORM FINAL WEIGHTS (Crucial)
    # Use only the weights from the REAL particles to transform back to 'w'.
    psi_final_real = psi_final_full[:N]
    w = weight_transform_inv(psi_final_real, Z) # Use original Z, not Y

    ## Check that the masses after solving the problem are all equal

    # print("Relative Error in the First Step Hessian", rel_error)
    print("Final Mass: ", ot.pd.integrals()[:N])
    print( "Error In Final Mass: ", np.linalg.norm(ot.get_masses()[:N] - ot.pd.integrals()[:N]) / np.linalg.norm(ot.get_masses()[:N]))
    print( "Internal Energy: ", ot.pd.internal_energy() )
    print( "Centroids:", ot.pd.centroids() )
    print( "Optimized weights: ", w )

# plt.figure()                        # optional: starts a fresh figure
# plt.plot(NList, Errors, marker='o') # line with dots at each point
# plt.xlabel('N')                     # x-axis label
# plt.ylabel('Error')                 # y-axis label
# plt.tight_layout()                  # tidy up margins
# plt.show()  

## Visualise the solution to the Optimal transport probelm

filename = 'pb.vtk'
ot.pd.display_vtk( filename )

# # Identify the vertices of the cell
# grid = pv.read(filename)
# grid = grid.clean(tolerance=1e-8)

# # Grab triangle→cell IDs
# cell_ids = grid.cell_data['num'].astype(int)

# # Extract the two triangle-submeshes for cell 0 and cell 1
# cells0 = np.where(cell_ids == 0)[0]
# cells1 = np.where(cell_ids == 1)[0]
# region0 = grid.extract_cells(cells0)
# region1 = grid.extract_cells(cells1)

# # Build the set of **domain** boundary edges (so we can ignore them later)
# dom_edges = set()
# dom_bnd = grid.extract_feature_edges(
#     boundary_edges=True,
#     feature_edges=False,
#     manifold_edges=False,
#     non_manifold_edges=False
# )
# # `dom_bnd.lines` is an Nx3 array [2, ptId0, ptId1] for each segment
# for _, i0, i1 in dom_bnd.lines.reshape(-1, 3):
#     p0 = tuple(dom_bnd.points[i0][:2])
#     p1 = tuple(dom_bnd.points[i1][:2])
#     dom_edges.add(tuple(sorted((p0, p1))))

# # Helper to pull out the **interior** boundary edges of a region
# def interior_edges(region):
#     bnd = region.extract_feature_edges(
#         boundary_edges=True,
#         feature_edges=False,
#         manifold_edges=False,
#         non_manifold_edges=False
#     )
#     segs = []
#     for _, i0, i1 in bnd.lines.reshape(-1, 3):
#         p0 = tuple(bnd.points[i0][:2])
#         p1 = tuple(bnd.points[i1][:2])
#         e = tuple(sorted((p0, p1)))
#         # drop any edges that lie on the *domain* boundary
#         if e not in dom_edges:
#             segs.append(e)
#     return set(segs)

# edges0 = interior_edges(region0)
# edges1 = interior_edges(region1)

# # The **shared** edges between region0 & region1 are exactly their mutual face
# shared = edges0 & edges1
# if not shared:
#     raise RuntimeError("No shared face found between cells 0 & 1!")
    
# # Stitch those shared segments into a tiny graph to pick out the two degree-1 nodes
# graph = defaultdict(list)
# for a, b in shared:
#     graph[a].append(b)
#     graph[b].append(a)

# endpoints = [pt for pt, nbrs in graph.items() if len(nbrs) == 1]
# if len(endpoints) != 2:
#     raise RuntimeError(f"Expected exactly two end-points, got {endpoints!r}")

# print("Boundary between cell 0 ↔ 1 runs from", endpoints[0], "to", endpoints[1])

# Mass of cells
colours=ot.pd.integrals()

# Read the data
grid=pv.read(filename)

# create cell data that gives the cell volumes, this allows us to colour by cell volumes
cell_colours = colours[grid.cell_data['num'].astype(int)]
grid.cell_data['Mass']=cell_colours

# plot the data with an automatically created plotter, for a static picture use backend='static'
plotter = pv.Plotter(window_size=[800,800])
plotter.add_mesh(grid) #, clim=[minvel, maxvel])
plotter.set_scale(xscale = 1, yscale = 1)

# Set the camera for 2D view
plotter.camera_position = 'xy'

# # convert 2D→3D by dropping z=0
# pts2d = np.array(endpoints)                # shape (2,2)
# pts3d = np.hstack([pts2d, np.zeros((2,1))])  # shape (2,3)

# # make a PyVista point cloud
# cloud = pv.PolyData(pts3d)

# # add to the existing plotter
# plotter.add_mesh(
#     cloud,
#     color='red',
#     point_size=10                     # bump up the radius
# )

# Render the frame
plotter.show()