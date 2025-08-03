from .domain_types import ConvexPolyhedraAssembly
from .radial_funcs import RadialFuncEntropy
from .radial_funcs import RadialFuncInBall
from .radial_funcs import RadialFuncUnit
from .PowerDiagram import PowerDiagram
import numpy as np
import warnings
from scipy.optimize import line_search
from scipy.optimize.linesearch import LineSearchWarning
from scipy.sparse import csr_matrix
import importlib

# This is a new helper function to set PETSc options
def set_petsc_options(options_dict):
    """
    Sets PETSc options from a Python dictionary.
    Example: {'ksp_type': 'cg', 'pc_type': 'jacobi'}
    """
    import sys
    for key, value in options_dict.items():
        sys.argv.append(f'-{key}')
        if value is not None:
            sys.argv.append(str(value))

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

def recalculate_replica_weights(Z_initial, psi_extended, f, c_p, Pi_0, Lx):
    """
    Recalculates the psi weights of replicas based on the current psi weights of the real diracs.
    """
    num_real = Z_initial.shape[0]
    psi_real_current = psi_extended[:num_real]
    w_real_current = weight_transform_inv(psi_real_current, Z_initial, f, c_p, Pi_0)

    replica_index_offset = num_real
    for i in range(num_real):
        parent_w = w_real_current[i]
        for n in [-1, 1]:
            T1 = Lx * n
            z_rep = np.array([Z_initial[i, 0] + T1, Z_initial[i, 1]])
            psi_rep_new = weight_transform(parent_w, z_rep.reshape(1, -1), f, c_p, Pi_0)[0]
            psi_extended[replica_index_offset] = psi_rep_new
            replica_index_offset += 1

    return psi_extended

def get_periodic_jacobian(num_real, num_replicas_per_side=1):
    """
    Builds the Jacobian matrix J for the dependency dw_v = J * dw_r for
    a 1D periodic system.
    """
    num_virtual_per_real = 2 * num_replicas_per_side
    num_virtual = num_real * num_virtual_per_real
    
    J = np.zeros((num_virtual, num_real))
    
    # Each virtual particle's weight depends only on its corresponding
    # real particle's weight, with a derivative of 1.
    for i in range(num_real):
        # Identify the rows corresponding to the i-th real particle's replicas
        start_row = i * num_virtual_per_real
        end_row = start_row + num_virtual_per_real
        J[start_row:end_row, i] = 1.0
            
    return J

def dist(a, b):
    return np.linalg.norm(a - b, 2)

class BadInitialGuess(Exception):
    def __init__(self, ot):
        self.ot = ot


class OptimalTransport:
    def __init__(self, positions=None, weights=None, domain=None, masses=None, radial_func=RadialFuncUnit(),
                 obj_max_dw=1e-8, obj_max_dm=0, linear_solver="Petsc",
                 petsc_options=None, verbosity=0):
        """
           stopping criterion = first obj_max_xy that is != 0
             * obj_max_dw => delta weights between two iterations
             * obj_max_dm => 
        """

        self.pd = PowerDiagram(positions, weights, domain, radial_func)
        self.obj_max_dw = obj_max_dw
        self.obj_max_dm = obj_max_dm
        self.masses = masses
        
        self.linear_solver = linear_solver
        self.verbosity = verbosity
        self.max_iter = 1000
        self.delta_m = []
        self.delta_w = []

        # --- MODIFICATION: Set PETSc options upon initialization ---
        if self.linear_solver == "Petsc" and petsc_options is not None:
            set_petsc_options(petsc_options)
        # -----------------------------------------------------------

        self._linear_solver_inst = None
        self._masses_are_new = True

    def get_positions(self):
        return self.pd.positions

    def set_positions(self, new_positions):
        self.pd.set_positions(new_positions)

    def get_masses(self):
        return self.masses

    def set_masses(self, new_masses):
        self._masses_are_new = True
        self.masses = new_masses

    def get_domain(self):
        return self.pd.get_domain()

    def set_domain(self, new_domain):
        self.pd.set_domain(new_domain)

    def get_weights(self):
        return self.pd.weights

    def set_weights(self, new_weights):
        self.pd.set_weights(new_weights)

    def adjust_weights(self, initial_weights=None, ret_if_err=False, relax=1.0):
        assert( self.obj_max_dw or self.obj_max_dm )

        if not ( initial_weights is None ):
            self.set_weights( initial_weights )
            
        if self.pd.domain is None:
            domain = ConvexPolyhedraAssembly()
            if self.pd.positions.shape[ 1 ] == 2:
                domain.add_box([0, 0], [1, 1])
            else:
                domain.add_box([0, 0, 0], [1, 1, 1])
            self.pd.set_domain( domain )

        if self.masses is None:
            N = self.pd.positions.shape[0]
            if isinstance(self.pd.radial_func, RadialFuncUnit):
                self.masses = np.ones(N) * self.pd.domain.measure() / N
            elif isinstance(self.pd.radial_func, RadialFuncInBall):
                self.masses = np.ones(N) * 1e-6
            elif isinstance(self.pd.radial_func, RadialFuncEntropy):
                self.masses = np.ones(N) * self.pd.domain.measure() / N
            else:
                TODO

        if self.pd.weights is None:
            self.pd.weights = np.sqrt( self.masses )

        linear_solver = self._get_linear_solver()
        old_weights = self.pd.weights + 0.0
        for num_iter in range(self.max_iter):
            # derivatives
            mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)
            if mvs.error:
                if num_iter == 0:
                    # print("initial guess for the weight lead to void cells, trying with 0 weights")
                    # self.pd.set_weights( self.pd.weights * 0.0 )
                    # mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)
                    # if mvs.error:
                    #     raise BadInitialGuess( self )
                    raise BadInitialGuess( self )
                else:
                    ratio = 0.5
                    self.pd.set_weights(
                        (1 - ratio) * old_weights + ratio * self.pd.weights
                    )
                    if ret_if_err:
                        return True
                    if (self.verbosity > 1):
                        print("bim (going back)")
                    continue
            old_weights = self.pd.weights

            #
            if self.pd.radial_func.need_rb_corr():
                mvs.m_values[0] *= 2
            mvs.v_values -= self.masses

            # "dm" stopping criterion
            nm = np.max(np.abs(mvs.v_values))
            self.delta_m.append(nm)
            if self.obj_max_dm:
                if self.verbosity > 1:
                    print("max dm:", nm)
                if nm < self.obj_max_dm:
                    break

            # linear system
            A = linear_solver.create_matrix(
                self.pd.weights.shape[0],
                mvs.m_offsets,
                mvs.m_columns,
                mvs.m_values
            )

            b = linear_solver.create_vector(
                mvs.v_values
            )

            x = linear_solver.solve(A, b)

            # update weights
            loc_relax = relax
            cpt_loc = 0
            while True:
                W = self.pd.get_weights() - loc_relax * x
                if self.pd.radial_func.ball_cut() == False or np.all( W >= 0 ): # HUM
                    self.pd.set_weights( W )
                    break
                if self.verbosity > 1:
                    print("negative weight, loc_relax=", loc_relax)
                loc_relax *= 0.75

                cpt_loc += 1
                if cpt_loc == 50:
                    print( "impossible to get positive weights" )
                    return True

            # "dw" stopping criterion
            nw = np.max(np.abs(x))
            self.delta_w.append(nw)
            if self.obj_max_dw:
                if self.verbosity > 1:
                    print("max dw:", nw)
                if nw < self.obj_max_dw:
                    break
                
        return False

    def get_centroids(self):
        return self.pd.centroids()

    def display_vtk(self, filename, points=False, centroids=False):
        self.pd.display_vtk(filename, points, centroids)

    def display_asy(self, filename, preamble="", closing="", output_format="pdf", linewidth=0.02, dotwidth=0.0, values=np.array([]), colormap="inferno", avoid_bounds=False, min_rf=1, max_rf=0):
        self.pd.display_asy(filename, preamble, closing, output_format, linewidth, dotwidth, values, colormap, avoid_bounds, min_rf, max_rf)

    def nb_diracs(self):
        return self.pd.positions.shape[0]

    def dim(self):
        return self.pd.positions.shape[1]

    def set_stopping_criterion(self, value, type="max delta masses"):
        """
            Possible values for type
            * "max delta masses" => max(abs(actual masses - target ones))
            * "max delta weights" => max(abs(weights - weights last iteration))
        """

        self.obj_max_dw = 0
        self.obj_max_dm = 0

        if type == "max delta masses":
            self.obj_max_dm = value
            return

        if type == "max delta weights":
            self.obj_max_dw = value
            return

        raise "'{}' is not a known stopping criterion type".format( type )

    def _get_linear_solver(self):
        if (self._linear_solver_inst is not None):
            return self._linear_solver_inst
        
        sparse_linear_solvers = ('CuPyx', 'Petsc', 'Scipy')
        msg = 'Available solvers are: {}.'.format(', '.join(sparse_linear_solvers))
        assert self.linear_solver in sparse_linear_solvers, msg

        for solver in (self.linear_solver,)+sparse_linear_solvers:
            try:
                module = importlib.import_module('pysdot.solvers.{}'.format(solver))
                break
            except ImportError:
                continue
        else:
            msg='Could not import any of the solver modules.'
            raise ImportError(msg)

        if (self.verbosity > 0):
            print('Sucessfully imported sparse linear solver {}.'.format(solver))

        self._linear_solver_inst = module.Solver()
        return self._linear_solver_inst

    def adjust_weights_periodic(self, Z_initial, f, c_p, Pi_0, Lx, ret_if_err=False, relax=1.0):
        """
        Corrected Newton's method for the periodic case using robust
        block matrix algebra to compute the effective Hessian.
        """
        assert(self.obj_max_dw or self.obj_max_dm)

        num_real = Z_initial.shape[0]
        num_total = self.pd.positions.shape[0]
        num_replicas_per_side = (num_total // num_real - 1) // 2

        # Get the Jacobian that defines the dependency between weights
        Jacobian = get_periodic_jacobian(num_real, num_replicas_per_side)

        old_weights = self.pd.weights.copy()

        for num_iter in range(self.max_iter):
            # 1. GET FULL DERIVATIVES
            mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)
            if mvs.error:
                if num_iter == 0: raise Exception("BadInitialGuess: Void cells on first iteration.")
                # Standard error handling
                ratio = 0.5
                self.pd.set_weights((1 - ratio) * old_weights + ratio * self.pd.weights)
                continue
            old_weights = self.pd.weights.copy()

            # 2. CONSTRUCT THE EFFECTIVE N x N SYSTEM ROBUSTLY
            
            # A. Create the full sparse Hessian H. The library computes -H.
            #    This method is proven to work from your debugging code.
            H_full_sparse = -csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets), shape=(num_total, num_total))
            H_full = H_full_sparse.toarray()

            # B. Extract the required blocks using slicing
            H_rr = H_full[:num_real, :num_real]
            H_rv = H_full[:num_real, num_real:]
            
            # C. Calculate the effective Hessian using the derived formula
            H_eff = H_rr + H_rv @ Jacobian

            # D. Calculate the mass error vector for ONLY the real cells
            #    mvs.v_values contains the current masses of all cells.
            F_r = mvs.v_values[:num_real] - self.masses[:num_real]
            
            # 3. SOLVE THE N x N SYSTEM: H_eff * dw_r = -F_r
            try:
                dw_r = np.linalg.solve(H_eff, -F_r)
            except np.linalg.LinAlgError:
                print("Error: Effective Hessian is singular. Cannot solve.")
                print("H_eff:\n", H_eff)
                print("F_r:", F_r)
                return True

            # 4. UPDATE ALL WEIGHTS
            current_real_weights = self.pd.get_weights()[:num_real]
            updated_real_weights = current_real_weights - relax * dw_r
            
            temp_full_weights = self.pd.get_weights().copy()
            temp_full_weights[:num_real] = updated_real_weights
            
            # This function correctly enforces the dependency rule after the update
            psi_consistent = recalculate_replica_weights(
                Z_initial, temp_full_weights, f, c_p, Pi_0, Lx
            )
            self.pd.set_weights(psi_consistent)

            # 5. STOPPING CRITERIA (based on real quantities)
            nm = np.max(np.abs(F_r))
            self.delta_m.append(nm)
            if self.obj_max_dm and nm < self.obj_max_dm:
                if self.verbosity > 0: print(f"Stopped on max_dm criterion: {nm}")
                break

            nw = np.max(np.abs(dw_r))
            self.delta_w.append(nw)
            if self.obj_max_dw and nw < self.obj_max_dw:
                if self.verbosity > 0: print(f"Stopped on max_dw criterion: {nw}")
                break
                
        return False

    def adjust_weights_hybrid_periodic(self, Z_initial, f, c_p, Pi_0, Lx, relax=1.0, epsilon=1e-8):
        """
        A robust hybrid optimization method. This version uses a corrected
        loop structure and state management to ensure both BFGS and Newton
        phases are functional.
        """
        assert(self.obj_max_dw or self.obj_max_dm)

        num_real = Z_initial.shape[0]
        num_total = self.pd.positions.shape[0]
        Jacobian = get_periodic_jacobian(num_real, (num_total // num_real - 1) // 2)

        # --- CORRECTED INITIALIZATION BLOCK ---
        mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)
        if mvs.error:
            if self.verbosity > 0: print("Info: Initial state has void cells. Starting BFGS with identity Hessian.")
            B_k = np.eye(num_real)
            g_k = np.ones(num_real)
        else:
            if self.verbosity > 0: print("Info: Initial state is valid. Starting BFGS with regularized Hessian.")
            H_full_sparse = -csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets), shape=(num_total, num_total))
            H_eff = H_full_sparse[:num_real, :num_real].toarray() + H_full_sparse[:num_real, num_real:].toarray() @ Jacobian
            try:
                B_k = np.linalg.inv(H_eff + epsilon * np.eye(num_real))
            except np.linalg.LinAlgError:
                B_k = np.eye(num_real)
            g_k = mvs.v_values[:num_real] - self.masses[:num_real]
        
        old_weights = self.pd.weights.copy()
        
        for num_iter in range(self.max_iter):
            # 1. GET DERIVATIVES AND USE RATIO RECOVERY IF A STEP FAILED
            mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)
            if mvs.error:
                if num_iter == 0: raise BadInitialGuess(self)
                ratio = 0.5
                self.pd.set_weights((1 - ratio) * old_weights + ratio * self.pd.weights)
                if self.verbosity > 0: print(f"Info: Method -> {mode}. Iter {num_iter} -> void cells. Readjusting & retrying.")
                continue

            # 2. IF STEP WAS SUCCESSFUL, UPDATE CHECKPOINT AND DECIDE MODE
            old_weights = self.pd.weights.copy()
            current_masses = mvs.v_values
            
            if np.all(current_masses[:num_real] > 1e-12):
                mode = 'Newton'
            else:
                mode = 'BFGS'
                if num_iter == 0 and self.verbosity > 0:
                    print("--- Initial state has zero-mass cells. Starting with BFGS. ---")
            
            # 3. PERFORM A STEP
            if mode == 'Newton':
                if num_iter > 0 and self.delta_m[-1] < 1e-9: break # Already converged
                H_full_sparse = -csr_matrix((mvs.m_values, mvs.m_columns, mvs.m_offsets), shape=(num_total, num_total))
                H_eff = H_full_sparse[:num_real, :num_real].toarray() + H_full_sparse[:num_real, num_real:].toarray() @ Jacobian
                F_r = current_masses[:num_real] - self.masses[:num_real]
                try: dw_r = np.linalg.solve(H_eff, -F_r)
                except np.linalg.LinAlgError: return True
                
                updated_real_weights = old_weights[:num_real] - relax * dw_r
                temp_full_weights = old_weights.copy()
                temp_full_weights[:num_real] = updated_real_weights
                self.pd.set_weights(recalculate_replica_weights(Z_initial, temp_full_weights, f, c_p, Pi_0, Lx))
                dw_r_final = relax * dw_r

            else: # mode == 'BFGS'
                g_k = current_masses[:num_real] - self.masses[:num_real]
                if np.linalg.norm(g_k) < 1e-9: break
                
                p_k = -B_k @ g_k
                alpha = relax
                # Backtracking Line Search
                while True:
                    dw_r = alpha * p_k
                    temp_psi = old_weights.copy()
                    temp_psi[:num_real] += dw_r
                    self.pd.set_weights(recalculate_replica_weights(Z_initial, temp_psi, f, c_p, Pi_0, Lx))
                    mvs_new = self.pd.der_integrals_wrt_weights(stop_if_void=True)
                    if not mvs_new.error: break
                    alpha *= 0.5
                    if alpha < 1e-8: return True

                # Update BFGS approximation
                s_k = alpha * p_k
                g_k_plus_1 = mvs_new.v_values[:num_real] - self.masses[:num_real]
                y_k = g_k_plus_1 - g_k
                if abs(y_k.T @ s_k) > 1e-12:
                    rho_k = 1.0 / (y_k.T @ s_k)
                    I = np.eye(num_real)
                    term1 = I - rho_k * np.outer(s_k, y_k)
                    B_k = term1 @ B_k @ term1.T + rho_k * np.outer(s_k, s_k)
                dw_r_final = s_k
            
            # 4. STOPPING CRITERIA
            nm = np.max(np.abs(self.pd.integrals()[:num_real] - self.masses[:num_real]))
            self.delta_m.append(nm)
            if self.obj_max_dm and nm < self.obj_max_dm: break
            nw = np.max(np.abs(dw_r_final))
            self.delta_w.append(nw)
            if self.obj_max_dw and nw < self.obj_max_dw: break
        
        return False
    
    def adjust_weights_bfgs_to_positive(self, Z_initial, f, c_p, Pi_0, Lx):
        """
        A robust BFGS solver with a Wolfe line search to find a state 
        with all positive masses.
        """
        num_real = Z_initial.shape[0]
        
        # 1. INITIALIZATION
        w_k = self.pd.get_weights()[:num_real].copy()
        B_k = np.eye(num_real) # Initial Hessian approximation

        # Evaluate the initial state
        mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)

        # --- PHASE 1: FIND A VALID STATE ---
        # If the initial guess is bad, we must first find a valid geometry
        # before we can begin proper BFGS updates.
        if mvs.error:
            if self.verbosity > 0: print("Initial guess is bad. Starting search for a valid state...")
            # Use a simple search to escape the void-cell region
            search_dir = np.ones(num_real) # Simple, non-zero direction
            alpha = 1.0
            while alpha > 1e-9:
                temp_w_real = w_k - alpha * search_dir # Move away from initial guess
                
                # Update weights and check validity
                temp_full = self.pd.get_weights().copy(); temp_full[:num_real] = temp_w_real
                psi_consistent = recalculate_replica_weights(Z_initial, temp_full, f, c_p, Pi_0, Lx)
                self.pd.set_weights(psi_consistent)
                mvs = self.pd.der_integrals_wrt_weights(stop_if_void=True)

                if not mvs.error:
                    if self.verbosity > 0: print("Found a valid state. Proceeding with BFGS.")
                    break
                alpha *= 0.5
            
            if mvs.error:
                if self.verbosity > 0: print("BFGS failed: Could not find any valid state from the initial guess.")
                return True # Failed to find a valid state

        # Now we have a valid state, get the initial gradient
        g_k = mvs.v_values[:num_real] - self.masses[:num_real]
        w_k = self.pd.get_weights()[:num_real].copy()

        # --- PHASE 2: BFGS OPTIMIZATION WITH WOLFE SEARCH ---
        for num_iter in range(self.max_iter):
            # Check for success before starting the iteration
            if np.all(mvs.v_values > 1e-12):
                if self.verbosity > 0:
                    print(f"\nSuccess on iter {num_iter}: All cells have positive mass. Terminating BFGS.")
                return False

            # 2. GET SEARCH DIRECTION
            try:
                p_k = np.linalg.solve(B_k, -g_k)
            except np.linalg.LinAlgError:
                if self.verbosity > 0: print("Hessian approximation is singular. Resetting.")
                B_k = np.eye(num_real)
                p_k = -g_k

            # 3. LINE SEARCH
            alpha, g_k_plus_1, mvs_new, success = self._wolfe_line_search(
                w_k, p_k, g_k, Z_initial, f, c_p, Pi_0, Lx
            )
            if not success:
                return True # Line search failed, cannot continue

            # 4. UPDATE STATE
            s_k = alpha * p_k
            w_k += s_k
            y_k = g_k_plus_1 - g_k
            
            # 5. UPDATE HESSIAN APPROXIMATION
            if abs(np.dot(y_k, s_k)) > 1e-12: # Check curvature condition
                rho_k = 1.0 / np.dot(y_k, s_k)
                # Standard Sherman-Morrison update for B_k
                term1 = rho_k * np.outer(s_k, y_k)
                B_k = (np.eye(num_real) - term1) @ B_k @ (np.eye(num_real) - term1.T) + rho_k * np.outer(s_k,s_k)
            
            g_k = g_k_plus_1
            mvs = mvs_new

        if self.verbosity > 0:
            print(f"BFGS failed to achieve positive masses within {self.max_iter} iterations.")
        return True
    
    def adjust_weights_to_convergence_bfgs(self, Z_initial, f, c_p, Pi_0, Lx, tolerance=1e-3, max_iter=200):
        """
        Single‐pass “rescue + BFGS”:
        1) Rescue: as long as any m_i==0, do mini‐GD steps on those faces.
        2) Once all m_i>0, switch to full BFGS w/ backtracking‐Armijo on φ=½‖m*−m‖².
        """
        N      = Z_initial.shape[0]
        total  = self.pd.positions.shape[0]
        nrep   = (total//N - 1)//2
        Jrep   = get_periodic_jacobian(N, nrep)

        def set_w(w):
            tmp      = self.pd.get_weights().copy()
            tmp[:N]  = w
            psi      = recalculate_replica_weights(Z_initial, tmp, f, c_p, Pi_0, Lx)
            self.pd.set_weights(psi)

        def masses():
            return self.pd.integrals()[:N]

        def build_Heff():
            mvs   = self.pd.der_integrals_wrt_weights(stop_if_void=True)
            Hfull = -csr_matrix(
                (mvs.m_values, mvs.m_columns, mvs.m_offsets),
                shape=(total, total)
            ).toarray()
            return Hfull[:N,:N] + Hfull[:N,N:] @ Jrep

        target = self.get_masses()[:N]
        w_k    = self.pd.get_weights()[:N].copy()
        B_k    = np.eye(N)

        # 1) RESCUE subloop: fill any void cells via targeted GD
        print("=== Rescue phase (fill zero‐mass cells) ===")
        for rescue_it in range(20):
            m_k = masses()
            zeros = np.where(m_k <= 0)[0]
            if len(zeros)==0:
                print("All cells have positive mass — switching to BFGS.")
                break

            # gradient only on zero entries
            g_k = target - m_k
            p   = np.zeros_like(g_k)
            p[zeros] = -g_k[zeros]        # push w_i down to expand those regions
            if np.linalg.norm(p) < 1e-12:
                p = -g_k                 # fallback to full GD
            print(f" rescue_it={rescue_it}, zero_count={len(zeros)}, ‖p‖={np.linalg.norm(p):.3e}")

            # simple backtracking to raise the minimum zero cell to >0
            alpha = 1.0
            for bt in range(20):
                w_try = w_k + alpha*p
                set_w(w_try)
                m_try = masses()
                if np.all(m_try>0):
                    print(f"  → rescue success α={alpha:.2e}")
                    w_k = w_try
                    break
                alpha *= 0.5
            else:
                print("  ! rescue FAILED to fill voids")
                return False
        else:
            # if rescue loop never broke
            print("  ✗ rescue phase gave up")
            return False

        # 2) BFGS subloop
        print("=== BFGS phase (converge to tolerance) ===")
        # rebuild initial grad & φ
        m_k      = masses()
        g_k      = target - m_k
        φ_k      = 0.5 * g_k.dot(g_k)
        J_k      = build_Heff()
        gradφ_k  = J_k.T.dot(g_k)

        backtrack_max = 30
        for it in range(1, max_iter+1):
            err     = np.linalg.norm(g_k)
            if err < tolerance:
                print(f"✓ converged in {it-1} iters; ‖mass_err‖={err:.3e}")
                return True

            # quasi‑Newton step
            p_k = -B_k.dot(gradφ_k)
            print(f"\niter {it:3d}: ‖g‖={err:.3e}, φ={φ_k:.3e}, ‖p‖={np.linalg.norm(p_k):.3e}")

            # backtracking‐Armijo
            alpha = 1.0
            φ0    = φ_k
            for bt in range(backtrack_max):
                w_try = w_k + alpha*p_k
                set_w(w_try)
                m_try = masses()
                φ_try = 0.5 * ((target-m_try)**2).sum()
                print(f" bt{bt:2d}: α={alpha:.2e}, φ_try={φ_try:.3e}, min_mass={m_try.min():.3e}")
                # require strict φ drop and no voids
                if φ_try < φ0 and np.all(m_try>0):
                    break
                alpha *= 0.5
            else:
                print("  ! line‐search failed")
                return False

            # accept step
            w_k      = w_try
            m_k      = m_try
            g_try    = target - m_k
            J_try    = build_Heff()
            gradφ_try= J_try.T.dot(g_try)

            # BFGS update
            s = w_k - (self.pd.get_weights()[:N] - alpha*p_k)  # or simply s=alpha*p_k
            y = gradφ_try - gradφ_k
            sy= s.dot(y)
            print(f" s⋅y = {sy:.3e}")
            if abs(sy)>1e-12:
                rho = 1.0/sy
                V   = np.eye(N) - rho*np.outer(s,y)
                B_k = V.dot(B_k).dot(V.T) + rho*np.outer(s,s)
                print("  → updated B_k")
            else:
                B_k[:] = np.eye(N)
                print("  → reset B_k")

            # update for next iter
            g_k     = g_try
            gradφ_k = gradφ_try
            φ_k     = φ_try

        print("✗ failed to converge within max_iter")
        return False