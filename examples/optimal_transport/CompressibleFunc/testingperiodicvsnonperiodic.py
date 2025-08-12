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
import pandas as pd

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
    
    psi = w + term1 + term2 - term3
    return psi

def weight_transform_inv(psi, Z, f=1.0):
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

    w = psi - term1 - term2 + term3 
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
def create_periodic_dataset(Z, Lx):
    Z_ext = list(Z)
    for i in range(N):
        # Loop through periodic images (-1 for left, +1 for right)
        for n in [-1, 1]: 
            T1 = Lx * n 
                
            # 1. Calculate the true Z-replica position
            z_rep = np.array([Z[i, 0] + T1, Z[i, 1]])
                
            # Append the new replica data to our lists
            Z_ext.append(z_rep)

    return np.array(Z_ext)

#Define the rescaling function to improve the inital guess for the Damped Newton Solver
def rescale_weights(bx, Z, PeriodicX):
    """
    Rescales weights for 3D Laguerre tessellation, ensuring positive cell areas.

    Parameters:
        bx (list/tuple): Domain [xmin, ymin, zmin, xmax, ymax, zmax].
        Z (numpy.ndarray): Seed points (shape: n x 3).
        psi (float): Initial guess for weights.
        PeriodicX, PeriodicY (bool): Periodicity flags for each axis.

    Returns:
        tuple: Weights (numpy.ndarray), scaling factor (float), and translation vector (numpy.ndarray).
    """
    PeriodicY = False
    if PeriodicX and PeriodicY:
        return np.zeros(len(Z)), 0, 0

    min_Z, max_Z = np.min(Z, axis=0), np.max(Z, axis=0)
    lambda_vals = []

    # Calculate lambda_ for each non-periodic dimension
    for i, periodic in enumerate([PeriodicX, PeriodicY]):
        if not periodic:
            box_width = bx[i + 2] - bx[i]
            point_spread = max_Z[i] - min_Z[i]

            if point_spread > 1e-9:  # Use a small tolerance for floating point safety
                lambda_dim = box_width / point_spread
                lambda_vals.append(lambda_dim)

    if lambda_vals:
        # Choose the most constrained scaling factor and add a small buffer
        lambda_ = min(lambda_vals) * (1 - 1e-2)
    else:
        # This occurs with a single particle or collinear points.
        # No scaling is necessary or meaningful.
        lambda_ = 1.0

    # Calculate translation vector for non-periodic dimensions
    translation = []
    for i, periodic in enumerate([PeriodicX, PeriodicY]):
        if not periodic:
            center_dom = (bx[i + 2] + bx[i]) / 2
            center_rescaled = lambda_ * (min_Z[i] + max_Z[i]) / 2
            translation.append(center_dom - center_rescaled)
        else:
            translation.append(0)
    t = np.array(translation)[~np.array([PeriodicX, PeriodicY])]

    # Calculate weights
    Z_mod = Z[:, ~np.array([PeriodicX, PeriodicY])]
    w = (1 - lambda_) * np.square(np.linalg.norm(Z_mod, axis=1)) - 2 * np.dot(Z_mod, t)

    return w

## Set the System Parameters

f = 1 
g = 1 
c_p = 1003.5
Pi_0 = 0.864
gamma = 1.41
kappa = 65.00526358589575
box = [-0.1, 0.0, 0.1, 0.1]
PeriodicX = True
PeriodicY = False

## Construct the Domain

domain = ConvexPolyhedraAssembly()

# If working on a box do the following 

Lx, Ly = [box[i+2] - box[i] for i in range(2)]

# # If working on a distorted domain do the following

Nx = 25
Ny = 25
x = np.linspace(-0.3, 0.3, Nx)
y = np.linspace(0.0, 0.1, Ny)
X, Y = np.meshgrid(x, y)
points = np.array([X.flatten(), Y.flatten()]).T
tri = Delaunay(points)
distorted_points = expmap_inverse_vec(points, np.array([0, 1]), f, g)

numTri = np.shape(tri.simplices)[0] # number of triangles in the triangulation

# Add each triangle to the domain one by one
for T in tri.simplices:
    p = distorted_points[T,:] # coordinates of vertices in the triangle
    domain.add_simplex(p) # add the triangle to the domain

NList = np.linspace(1, 100, 20, True, dtype=int)
NumberofTrials = 10
TrialErrors = []

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

for M in range(0, NumberofTrials):
    Errors = []
    for N in NList:
        print("This test with ", N, " seeds is periodic ?", PeriodicX)

        # Generate x values between -1 and 1
        x = np.random.uniform(-0.1, 0.1, N)

        # Generate y values between 98900 and 102000
        y = np.random.uniform(900, 1100, N)

        # Combine x and y into a single array of shape (N, 2)
        Z = np.column_stack((x, y))

        # Transform the seeds (diracs) eto the exponential chart

        Z_ext = create_periodic_dataset(Z, Lx)
        Y = seed_transform(Z_ext)

        ## Set an inital guess for the weights

        psi0 = rescale_weights([-0.3, 0.045, 0.3, 0.1], Y, False)
        K_offset = kappa * gamma * (1 / 0.02) ** (gamma - 1)

        target_masses = np.ones(3 * N) / N 

        ## Initialise the optimal transport problem setting the resolution of the integration with Int_res which corresponds to the number of Gaussian Quadrature points 

        ot = OptimalTransport( positions = Y, weights = psi0, masses = target_masses, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, w_offset = K_offset, Int = True, Int_res = 9 ), petsc_options=solver_options, verbosity=0)
        pot = OptimalTransport( positions = Y, weights = psi0, masses = target_masses, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, w_offset = K_offset, Int = True, Int_res = 9 ), petsc_options=solver_options, verbosity=0)

        ## Set the error tolerance and the stopping criterion 

        err_tol = (1e-3 / 100) * (1 / N)
        ot.set_stopping_criterion(err_tol, 'max delta masses')
        pot.set_stopping_criterion(err_tol, 'max delta masses')

        # --- 3. SOLVE with both methods ---
        print("Solving with non-periodic (transform_aware) solver...")
        ot.adjust_weights_transform_aware(Z_ext, f)
        
        print("Solving with true periodic solver...")
        pot.adjust_weights_periodic(Z, f, Lx)

         # --- 3. CONVERGENCE CHECK ---
        target_mass_N = pot.get_masses()[:N]
        integrals_N = pot.pd.integrals()[:N]
        mass_error = np.linalg.norm(target_mass_N - integrals_N) / (np.linalg.norm(target_mass_N) + 1e-20)
        
        print(f"Periodic solver mass error: {mass_error:.4e}")
        
        # --- 4. STORE RESULT or NAN ---
        if mass_error > 1e-5:
            print(f"!! Convergence failed for N={N}. Storing nan. !!\n")
            Errors.append(np.nan)
        else:
            psi_final_non_periodic = ot.get_weights()
            psi_final_periodic = pot.get_weights()
            w_non_periodic = weight_transform_inv(psi_final_non_periodic + K_offset - c_p * Pi_0, Z_ext, f) 
            w_periodic = weight_transform_inv(psi_final_periodic + K_offset - c_p * Pi_0, Z_ext, f)

            w_real_non_periodic = w_non_periodic[:N]
            w_real_periodic = w_periodic[:N]
            epsilon = 1e-20
            relative_difference = np.linalg.norm(w_real_non_periodic - w_real_periodic) / (np.linalg.norm(w_real_periodic) + epsilon)
            
            print(f"-> Relative difference for N={N}: {relative_difference:.4e}\n")
            Errors.append(relative_difference)

    TrialErrors.append(Errors)

# --- Data Processing using nan-aware functions ---
trial_errors_array = np.array(TrialErrors, dtype=float)
mean_errors = np.nanmean(trial_errors_array, axis=0)
std_dev_errors = np.nanstd(trial_errors_array, axis=0)

# Filter out N-values where all trials failed (resulting in a nan mean)
valid_mask = ~np.isnan(mean_errors)
NList_final = np.array(NList)[valid_mask]
mean_errors_final = mean_errors[valid_mask]
std_dev_errors_final = std_dev_errors[valid_mask]

# --- Data Saving ---
save_df = pd.DataFrame({
    'N': NList_final,
    'Mean_Error': mean_errors_final,
    'Std_Dev_Error': std_dev_errors_final
})
csv_filename = 'solver_comparison_data_filtered.csv'
save_df.to_csv(csv_filename, index=False)
print(f"\nSaved processed data to {csv_filename}")

# --- Plotting ---
plt.figure(figsize=(12, 7))
plt.errorbar(NList_final, mean_errors_final, yerr=std_dev_errors_final,
             marker='o', linestyle='--', color='royalblue',
             label='Mean Relative Difference w/ Std Dev',
             capsize=5, ecolor='crimson', elinewidth=1.5)

# --- Analysis ---
positive_mask = mean_errors_final > 0
if np.any(positive_mask):
    coeffs = np.polyfit(np.log(NList_final[positive_mask]), np.log(mean_errors_final[positive_mask]), 1)
    slope = coeffs[0]
    fit_line = np.exp(coeffs[1]) * (NList_final**slope)
    plt.loglog(NList_final, fit_line, 'k-', label=f'Fit (Order ≈ {-slope:.2f})', alpha=0.7)

# --- Enhancements ---
plt.xscale('log')
plt.yscale('log')
plt.title('Convergence of Non-Periodic Solver to True Periodic Solution')
plt.xlabel('Number of Seeds (N)')
plt.ylabel('Mean Relative Difference Between Solver Weights')
plt.grid(True, which="both", linestyle=':', linewidth=0.5)
plt.legend()
plt.tight_layout()

# --- Save and Show ---
plot_filename = 'solver_comparison_plot_filtered.png'
plt.savefig(plot_filename, dpi=300)
print(f"Saved plot to {plot_filename}")
plt.show()