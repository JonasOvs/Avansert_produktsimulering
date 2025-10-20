# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:43:51 2018

@author: bjohau
"""

import numpy as np
import math

def rot_matrix(theta):
    """
    Return the 2x2 rotation matrix representing a rotation theta
    :param theta:  rotation angle in radians
    :return: Rotation matrix (or tensor)
    """
    s = math.sin(theta)
    c = math.cos(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    return R

def beam2local_def_disp(ex,ey, disp_global):
    """

    :param ex: element x coordinate [x1, x2] in undeformed position
    :param ey: element y coordinate [y1, y2] in undeformed position
    :param disp_global:  displacement vector [u1, v1, r1, u2, v2, r2] in global directions
    :return: disp_local_def: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = math.sqrt(eVec12 @ eVec12)

    #Ld = L0 #TODO: correct this, fikset jonas
    # Calculate deformed positions
    x1_def = ex[0] + disp_global[0]
    y1_def = ey[0] + disp_global[1]
    x2_def = ex[1] + disp_global[3]
    y2_def = ey[1] + disp_global[4]
    
    # Calculate deformed length
    eVec12_def = np.array([x2_def - x1_def, y2_def - y1_def])
    Ld = math.sqrt(eVec12_def @ eVec12_def)

    # TODO: Quite a bit here, fikset jonas
    theta_def = math.atan2(eVec12_def[1], eVec12_def[0])

    r1 = disp_global[2]
    r2 = disp_global[5]

    theta1_def = r1 - theta_def
    theta2_def = r2 - theta_def

    # theta1_def = 0.0  # TODO: correct this
    # theta2_def = 0.0  # TODO: correct this

    def_disp_local = np.array([ -0.5*(Ld - L0),
                                0.0,
                                theta1_def,
                                0.5 * (Ld - L0),
                                0.0,
                                theta2_def])

    # print("theta_def:", theta_def)
    # print("theta1_def:", theta1_def)
    # print("theta2_def:", theta2_def)
    return def_disp_local

#Gammel:
# def beam2corot_Ke_and_Fe(ex,ey,ep, disp_global):
#     """
#     Compute the stiffness matrix and internal forces for a two dimensional beam element
#     relative to deformed configuration.
    
#     :param list ex: element x coordinates [x1, x2]
#     :param list ey: element y coordinates [y1, y2]
#     :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
#     :param list disp_global displacement vector for the element [tx1,ty1,rz1,tx2,ty2,rz2]


#     :return mat Ke_global: element stiffness matrix [6 x 6]
#     :return mat fe_int_global: element internal forces [6 x 1]
#     """
#     # Undeformed length and unit vector along element
#     eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
#     L0 = math.sqrt(eVec12 @ eVec12)
#     eVec12 /= L0

#     # Deformed position and unit vector along element
#     ex_def = ex + [disp_global[0], disp_global[3]] 
#     ey_def = ey + [disp_global[1], disp_global[4]] 

#     # TODO: Quite a bit here
#     Ke_global = np.zeros((6,6))
#     fe_int_global = np.zeros(6)

#     return Ke_global, fe_int_global
#Ny
# def beam2corot_Ke_and_Fe(ex, ey, ep, disp_global):
#     """
#     Compute the stiffness matrix and internal forces for a 2D beam element 
#     using a simplified corotational approach (without projection matrix).

#     :param ex: [x1, x2] element x coordinates
#     :param ey: [y1, y2] element y coordinates
#     :param ep: [E, A, I] element properties
#     :param disp_global: [tx1, ty1, rz1, tx2, ty2, rz2] displacement vector

#     :return: Ke_global (6x6), fe_int_global (6,)
#     """
#     E, A, I = ep

#     # Original (undeformed) element length
#     dx = ex[1] - ex[0]
#     dy = ey[1] - ey[0]
#     L0 = math.sqrt(dx*dx + dy*dy)
#     l0 = dx / L0
#     m0 = dy / L0

#     # Direction cosines
#     c = l0
#     s = m0

#     # Local linear stiffness matrix for beam (6x6)
#     # Axial + bending (Euler-Bernoulli)
#     L = L0
#     EA_L = E*A / L
#     EI_L3 = E*I / (L**3)
#     EI_L2 = E*I / (L**2)
#     EI_L = E*I / L

#     # 6x6 local stiffness in local coords (axial + bending)
#     Ke_local = np.array([
#         [ EA_L,       0,        0, -EA_L,       0,        0],
#         [  0,  12*EI_L3,  6*EI_L2,   0, -12*EI_L3,  6*EI_L2],
#         [  0,   6*EI_L2,   4*EI_L,   0,  -6*EI_L2,   2*EI_L],
#         [-EA_L,       0,        0,  EA_L,       0,        0],
#         [  0, -12*EI_L3, -6*EI_L2,   0,  12*EI_L3, -6*EI_L2],
#         [  0,   6*EI_L2,   2*EI_L,   0,  -6*EI_L2,   4*EI_L]
#     ])

#     # Simple transformation to global coordinates (rotation only)
#     T = np.array([
#         [ c,  s, 0, 0, 0, 0],
#         [-s,  c, 0, 0, 0, 0],
#         [ 0,  0, 1, 0, 0, 0],
#         [ 0,  0, 0,  c,  s, 0],
#         [ 0,  0, 0, -s,  c, 0],
#         [ 0,  0, 0,  0,  0, 1]
#     ])

#     # Global stiffness
#     Ke_global = T.T @ Ke_local @ T

#     # Internal force vector (linear assumption)
#     fe_int_global = Ke_global @ disp_global

#     return Ke_global, fe_int_global
def beam2corot_Ke_and_Fe(ex, ey, ep, disp_global):
    """
    Compute the stiffness matrix and internal forces for a 2D beam element 
    using a simplified corotational approach (linear Euler-Bernoulli).

    :param ex: [x1, x2] element x coordinates
    :param ey: [y1, y2] element y coordinates
    :param ep: [E, A, I] element properties
    :param disp_global: [tx1, ty1, rz1, tx2, ty2, rz2] displacement vector

    :return: Ke_global (6x6), fe_int_global (6,)
    """

    E, A, I = ep

    # Element length
    dx = ex[1] - ex[0]
    dy = ey[1] - ey[0]
    L0 = np.sqrt(dx**2 + dy**2)

    # Retningcosinus
    c = dx / L0
    s = dy / L0

    # Lokal stivhetsmatrise (Euler-Bernoulli, 2D, 6x6)
    EA_L = E * A / L0
    EI_L3 = E * I / L0**3
    EI_L2 = E * I / L0**2
    EI_L = E * I / L0

    Ke_local = np.array([
        [ EA_L,      0,        0, -EA_L,       0,        0],
        [   0,  12*EI_L3,  6*EI_L2,   0, -12*EI_L3,  6*EI_L2],
        [   0,   6*EI_L2,   4*EI_L,   0,  -6*EI_L2,   2*EI_L],
        [-EA_L,      0,        0,  EA_L,       0,        0],
        [   0, -12*EI_L3, -6*EI_L2,   0,  12*EI_L3, -6*EI_L2],
        [   0,   6*EI_L2,   2*EI_L,   0,  -6*EI_L2,   4*EI_L]
    ])

    # Blokk-diagonal transformasjon for 2D element (3x3 per node)
    R = np.array([
        [ c, s, 0],
        [-s, c, 0],
        [ 0, 0, 1]
    ])
    T = np.zeros((6,6))
    T[0:3,0:3] = R
    T[3:6,3:6] = R

    # Global stivhetsmatrise og interne krefter
    Ke_global = T.T @ Ke_local @ T
    fe_int_global = Ke_global @ disp_global

    return Ke_global, fe_int_global
    
def beam2corot_Te(ex,ey):
    """
    Compute the transformation matrix for an element
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Te: element transformation from global to local
    """

    n = np.array([ex[1]-ex[0],ey[1]-ey[0]])
    L = np.linalg.norm(n)
    n = n / L  
    
    Te=np.array([
        [ n[0], n[1],  0.,    0.,    0.,   0.],
        [-n[1], n[0],  0.,    0.,    0.,   0.],
        [0.,    0.,    1.,    0.,    0.,   0.],
        [0.,    0.,    0.,   n[0],  n[1],  0.],
        [0.,    0.,    0.,  -n[1],  n[0],  0.],
        [0.,    0.,    0.,    0.,    0.,   1.]
    ])
    

    return Te
    
    
def beam2local_stiff(L,ep):
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    :param list L : element length
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :return mat Kle: element stiffness matrix [6 x 6]
    """
    
    E=ep[0]
    A=ep[1]
    I=ep[2]
        
    Kle = np.array([
        [E*A/L,              0.,           0.,    -E*A/L,            0.,           0. ],
        [   0.,    12*E*I/L**3.,  6*E*I/L**2.,        0., -12*E*I/L**3.,  6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      4*E*I/L,        0.,  -6*E*I/L**2.,     2*E*I/L  ],
        [-E*A/L,             0.,           0.,     E*A/L,            0.,           0. ],
        [   0.,   -12*E*I/L**3., -6*E*I/L**2.,        0.,  12*E*I/L**3., -6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      2*E*I/L,        0.,  -6*E*I/L**2.,      4*E*I/L ]
    ])
     
    return Kle


def beam2e(ex, ey, ep, eq=None):
    """
    Compute the linear stiffness matrix for a two dimensional beam element.
    Largely from CALFEM core module

    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: element consistent force for distributed load [6 x 1] (if eq!=None)
    """

    n = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L = np.linalg.norm(n)
    n = n / L

    qx = 0.
    qy = 0.
    if not eq is None:
        qx = eq[0]
        qy = eq[1]

    Kle = beam2local_stiff(L,ep)

    fle = L * np.array([[qx / 2, qy / 2, qy * L / 12, qx / 2, qy / 2, -qy * L / 12]]).T

    Te = beam2corot_Te(ex,ey)

    Ke = Te.T @ Kle @ Te
    fe = Te.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe