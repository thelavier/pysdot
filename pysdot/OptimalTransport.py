from .domain_types import ConvexPolyhedraAssembly
from .radial_funcs import RadialFuncEntropy
from .radial_funcs import RadialFuncInBall
from .radial_funcs import RadialFuncUnit
from .PowerDiagram import PowerDiagram
import numpy as np
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
