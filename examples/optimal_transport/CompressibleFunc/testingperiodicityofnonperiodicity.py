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

    w = psi - term1 - term2 + term3 #- c_p * Pi_0
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
g = 1 #10
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

# domain.add_box(offset, size)

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

NList = list(np.unique(np.logspace(0, 4, num=20, dtype=int)))
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
        # print(solver_options)
        # np.random.seed(42)

        # Generate x values between -1 and 1
        x = np.random.uniform(-0.1, 0.1, N)

        # Generate y values between 98900 and 102000
        y = np.random.uniform(900, 1100, N)

        # Combine x and y into a single array of shape (N, 2)
        Z = np.column_stack((x, y))

        # Transform the seeds (diracs) eto the exponential chart

        Z_ext = create_periodic_dataset(Z, Lx)
        Y = seed_transform(Z_ext)

        # print("Initial seeds:", Y)

        ## Set an inital guess for the weights

        psi0 = rescale_weights([-0.3, 0.045, 0.3, 0.1], Y, False)
        K_offset = kappa * gamma * (1 / 0.02) ** (gamma - 1)

        target_masses = np.ones(3 * N) / N 

        ## Initialise the optimal transport problem setting the resolution of the integration with Int_res which corresponds to the number of Gaussian Quadrature points 

        ot = OptimalTransport( positions = Y, weights = psi0, masses = target_masses, domain = domain, radial_func = CompressibleFunc( kappa = kappa, gamma = gamma, g = g, f_cor = f, pi_0 = Pi_0, c_p = c_p, w_offset = K_offset, Int = True, Int_res = 9 ), petsc_options=solver_options, verbosity=1)

        ## Set the error tolerance and the stopping criterion 

        err_tol = (1e-3 / 100) * (1 / N)
        ot.set_stopping_criterion(err_tol, 'max delta masses')

        # 6. SOLVE THE OPTIMAL TRANSPORT PROBLEM
        ot.adjust_weights_transform_aware(Z_ext, f)
        psi_final_full = ot.get_weights()

        # 7. TRANSFORM FINAL WEIGHTS 
        w = weight_transform_inv(psi_final_full + K_offset - c_p * Pi_0, Z_ext) 
        w_real = w[:N]
        w_replicas = w[N:]
        w_replicas_reshaped = w_replicas.reshape((N, 2))

        # 3. Calculate the difference between each real weight and its corresponding replicas
        # We use np.newaxis to broadcast the real weights across all replicas
        abs_difference = np.abs(w_replicas_reshaped - w_real[:, np.newaxis])

        # 4. Calculate the relative error (add a small epsilon for stability if a weight is zero)
        epsilon = 1e-20
        relative_error = abs_difference / (np.abs(w_real[:, np.newaxis]) + epsilon)

        # 5. Print summary statistics
        print("--- Periodicity Constraint Check ---")
        print(f"Max absolute difference:   {np.max(abs_difference):.4e}")
        print(f"Max relative error:        {np.max(relative_error):.4e}")
        print(f"Mean relative error:       {np.mean(relative_error):.4e}")
        Errors.append(np.mean(relative_error))

        ## Check that the masses after solving the problem are all equal

        print( "Error In Final Mass of the Periodic Solver: ", np.linalg.norm(ot.get_masses()[:N] - ot.pd.integrals()[:N]) / np.linalg.norm(ot.get_masses()[:N]))

    TrialErrors.append(Errors)

# --- Data Processing ---
# This section assumes you have just finished your experiment loops and have:
# NList: The list of resolutions, e.g., [1, 2, 4, ...]
# TrialErrors: A list of lists, e.g., [[err_n1_t1, err_n2_t1, ...], [err_n1_t2, err_n2_t2, ...], ...]

# Convert the list of lists into a NumPy array for efficient calculations.
# The shape will be (NumberofTrials, len(NList)).
trial_errors_array = np.array(TrialErrors)

# Calculate the mean error for each resolution N across all trials.
# axis=0 calculates the mean down each column.
mean_errors = np.mean(trial_errors_array, axis=0)

# Calculate the standard deviation for each resolution N.
std_dev_errors = np.std(trial_errors_array, axis=0)


# --- Data Saving ---
# Create a pandas DataFrame to store the processed data neatly.
results_df = pd.DataFrame({
    'N': NList,
    'Mean_Error': mean_errors,
    'Std_Dev_Error': std_dev_errors
})

# Save the DataFrame to a CSV file. This is great for records and sharing.
csv_filename = 'convergence_analysis_data.csv'
results_df.to_csv(csv_filename, index=False)
print(f"\nSaved processed data to {csv_filename}")


# --- Plotting ---
# Create a figure with a specific size for better readability.
plt.figure(figsize=(12, 7))

# Plot the mean error with error bars representing the standard deviation.
# This gives a clear picture of the result's variability.
plt.errorbar(NList, mean_errors, yerr=std_dev_errors,
             marker='o',          # Add markers at each point
             linestyle='--',       # Use a dashed line to connect points
             color='royalblue',    # Line color
             label='Mean Error w/ Std Dev',
             capsize=5,            # Adds caps to the error bars
             ecolor='crimson',     # Color of the error bars
             elinewidth=1.5)       # Thickness of the error bars

# --- Analysis: Fit a line to estimate convergence rate ---
# On a log-log plot, an error relationship like Error = C * N^(-p)
# becomes a straight line. The slope gives the order of convergence.
# We will use the mean errors for the fit.
coeffs = np.polyfit(np.log(NList), np.log(mean_errors), 1)
slope = coeffs[0]

# Create data for the fitted line to overlay on the plot.
fit_line = np.exp(coeffs[1]) * (np.array(NList)**slope)
plt.loglog(NList, fit_line, 'k-', label=f'Fit (Order ≈ {-slope:.2f})', alpha=0.7)


# --- Enhancements for Clarity ---
# Set both axes to a logarithmic scale.
plt.xscale('log')
plt.yscale('log')

# Add a title and labels.
plt.title('Convergence of Periodic Solver: Error vs. Resolution')
plt.xlabel('Number of Seeds (N)')
plt.ylabel('Mean Relative Error (Periodicity)')

# Add a grid for both major and minor ticks.
plt.grid(True, which="both", linestyle=':', linewidth=0.5)

# Add a legend.
plt.legend()

# Ensure all plot elements fit nicely.
plt.tight_layout()


# --- Save the Final Plot ---
# Save the figure with a high DPI for publication quality.
plot_filename = 'convergence_plot_with_stats.png'
plt.savefig(plot_filename, dpi=300)
print(f"Saved plot to {plot_filename}")

# Display the plot.
plt.show()