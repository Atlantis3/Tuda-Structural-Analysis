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



def material_type(stiffness_matrix):
    '''Checks the given stiffness matrix for the type of linearity, whether the material is isotropic, orthotropic, monotropic, transversal isotropic or fully anisotropic
    
    Parameters
    ----------
    stiffness_matrix: numpy array, 3x3 or 6x6 matrix
        The Stifness matrix in its standard form, a 3x3 matrix for 2d problems and a 6x6 matrix for 3d problems

    Returns
    -------
    output: string
        Tells the user if the material is isotropic, orthotropic, monotropic, transversal isotropic or fully anisotropic
    '''
    if len(stiffness_matrix[0]) == len(stiffness_matrix[:,0]):
        if len(stiffness_matrix[0]) == 3:
            pass
        elif (len(stiffness_matrix[0]) == 6):
            c11 = stiffness_matrix[0,0]
            c12 = stiffness_matrix[0,1]
            c13 = stiffness_matrix[0,2]
            c14 = stiffness_matrix[0,3]
            c15 = stiffness_matrix[0,4]
            c16 = stiffness_matrix[0,5]
            
            c22 = stiffness_matrix[1,1]
            c23 = stiffness_matrix[1,2]
            c24 = stiffness_matrix[1,3]
            c25 = stiffness_matrix[1,4]
            c26 = stiffness_matrix[1,5]
            
            c33 = stiffness_matrix[2,2]
            c34 = stiffness_matrix[2,3]
            c35 = stiffness_matrix[2,4]
            c36 = stiffness_matrix[2,5]
            
            c44 = stiffness_matrix[3,3]
            c45 = stiffness_matrix[3,4]
            c46 = stiffness_matrix[3,5]
            
            c55 = stiffness_matrix[4,4]
            c56 = stiffness_matrix[4,5]
            
            c66 = stiffness_matrix[5,5]
            
            
            
            
            
            
            if (c14==0)&(c15==0)&(c24==0)&(c25==0)&(c34==0)&(c35==0)&(c46==0)&(c56==0)&(c16!=0)&(c26!=0)&(c36!=0)&(c45!=0):
                return 'The material is Monotropic and it requires 13 independent material properties, it is symmetric with respect to x1-x2 plane. '
            elif (c14==0)&(c15==0)&(c16==0)&(c24==0)&(c25==0)&(c26==0)&(c34==0)&(c35==0)&(c36==0)&(c45==0)&(c46==0)&(c56==0):
                if (c22==c33)&(c12==c13)&(c55==c66)&(c44==(c22-c23)/2.0):
                    return 'The material is Transversal Isotropy and it requires 5 independent material properties'
                
                elif (c11==c22==c33)&(c12==c13==c23)&(c44==c55==c66==(c11-c12)/2.0):
                    return 'The material is Isotropic and it requires 2 Independent material properties'
                
                else:
                    return 'The material is Orthotropic and it requires 9 independent material properties, it is symmetric with respect to three mutually perpendicular planes. '
            else:
                return 'The material is fully anisotropic and requires 21 independent material properties'
            
            
            
        else:
            raise Exception('The given matrix is not of order 3x3 (2d problem) or 6x6 (3d problem)')
        
    else:
        raise Exception('The given matrix is not a square matrix')