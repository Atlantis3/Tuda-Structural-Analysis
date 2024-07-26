# Created on 26.07.2024 at Akaflieg Darmstadt, Author @ akram_metar@akaflieg.tu-darmstadt.de
# This is the official file containing all the functions and codes related to the structural analysis

def Q_bars(E11,E22,G12,nu12,theta):
        '''Returns the Stiffness matrix or the Q_bar matrix for an individual lamina
        
        Parameters
        ----------
        E11 : float
            The Young's modulas in fibre direction, E11
        E22 : float
            The young's modulas in perpendicular to fibre direction
        G12 : float
            The Shear modulas of elasticity G12
        nu12 : float
            The poisson's ratio
        theta : float
            The orinetation angle of fibres in degrees

        Returns
        -------
        Q_bar : numpy array 3x3 matrix
            Returns the numpy array in matrix form
            np.array([[Q11_bar,Q12_bar,Q16_bar],[Q12_bar,Q22_bar,Q26_bar],[Q16_bar,Q26_bar,Q66_bar]])
        '''
        import numpy as np
        import math
        # first let us determine the value of nu21
        nu21 = (E22/E11)*nu12

        # now let us determine the values of Q11,Q22,Q12 and Q66 as if were the case for isotropic plate
        Q11 = E11/(1-(nu12*nu21))
        Q22 = E22/(1-(nu12*nu21))
        Q12 = (E22/(1-(nu12*nu21)))*nu12
        Q66 = G12

        # define cos theta and sin theta
        C = math.cos(math.radians(theta))
        S = math.sin(math.radians(theta))

        # now let us determine the values Qbar, for a ply laminate at an angle theta
        Q11_bar = (Q11*C**4)+(2*(Q12+2*Q66)*C**2*S**2)+(Q22*S**4)
        Q22_bar = (Q11*S**4)+(2*(Q12+2*Q66)*C**2*S**2)+(Q22*C**4)
        Q12_bar = ((Q11+Q22-4*Q66)*C**2*S**2)+(Q12*(C**4+S**4))
        Q66_bar = ((Q11+Q22-2*Q12-2*Q66)*C**2*S**2)+(Q66*(C**4+S**4))
        Q16_bar = ((Q11-Q12-2*Q66)*C**3*S)+((Q12-Q22+2*Q66)*C*S**3)
        Q26_bar = ((Q11-Q12-2*Q66)*C*S**3)+((Q12-Q22+2*Q66)*C**3*S)

        # Define the matrix as an numpy array and retun the array
        Q_bar = np.array([[Q11_bar,Q12_bar,Q16_bar],[Q12_bar,Q22_bar,Q26_bar],[Q16_bar,Q26_bar,Q66_bar]])
        return Q_bar