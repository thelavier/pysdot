
#
class CompressibleFunc:
    def __init__(self, kappa = 1, gamma = 1.41, g = 1, f_cor = 1, pi_0 = 1, c_p = 1, w_offset = 0, Int = True, Int_res = 9):
        self.s = "compressible_func({:.16f} {:.16f} {:.16f} {:.16f} {:.16f} {:.16f} {:.16f} {:.16f})".format(kappa, gamma, g, f_cor, pi_0, c_p, w_offset, Int, Int_res)

    def name(self):
        return self.s

    def need_rb_corr(self):
        return False
        
    def ball_cut(self):
        return False
