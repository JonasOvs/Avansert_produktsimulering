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
import numpy as np
import math

def beam2local_def_disp(ex, ey, disp_global):
    """
    Transform global nodal displacements to local displacements
    relative to the deformed configuration (corotational formulation).

    :param ex: element x coordinates [x1, x2]
    :param ey: element y coordinates [y1, y2]
    :param disp_global: displacement vector [u1, v1, r1, u2, v2, r2]
    :return: def_disp_local: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """

    # --- Undeformed geometry ---
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = math.sqrt(eVec12 @ eVec12)

    # --- Deformed coordinates ---
    x1_def = ex[0] + disp_global[0]
    y1_def = ey[0] + disp_global[1]
    x2_def = ex[1] + disp_global[3]
    y2_def = ey[1] + disp_global[4]

    # --- Deformed length ---
    eVec_def = np.array([x2_def - x1_def, y2_def - y1_def])
    Ld = math.sqrt(eVec_def @ eVec_def)

    # --- Orientation angles (global x-axis to element axis) ---
    alpha0 = math.atan2(eVec12[1], eVec12[0])      # initial
    alphad = math.atan2(eVec_def[1], eVec_def[0])  # deformed
    dAlpha = alphad - alpha0                       # rigid-body rotation of the element

    # --- Compute local displacements ---
    # Transformation matrix from global to local (initial configuration)
    T0 = np.array([
        [ math.cos(alpha0),  math.sin(alpha0), 0, 0, 0, 0],
        [-math.sin(alpha0),  math.cos(alpha0), 0, 0, 0, 0],
        [ 0, 0, 1, 0, 0, 0],
        [ 0, 0, 0, math.cos(alpha0),  math.sin(alpha0), 0],
        [ 0, 0, 0,-math.sin(alpha0),  math.cos(alpha0), 0],
        [ 0, 0, 0, 0, 0, 1]
    ])

    disp_local_init = T0 @ disp_global

    # --- Local rotations relative to the deformed axis ---
    # r1 and r2 are nodal rotations in global coordinates.
    r1_global = disp_global[2]
    r2_global = disp_global[5]

    # Subtract rigid-body rotation (dAlpha) to get relative bending rotations
    theta1_def = r1_global - dAlpha
    theta2_def = r2_global - dAlpha

    # --- Define local displacement vector in deformed configuration ---
    def_disp_local = np.array([
        -0.5 * (Ld - L0),
         0.0,
         theta1_def,
         0.5 * (Ld - L0),
         0.0,
         theta2_def
    ])

    return def_disp_local



def beam2corot_Ke_and_Fe(ex,ey,ep, disp_global):
    """
    Compute the stiffness matrix and internal forces for a two dimensional beam element
    relative to deformed configuration.
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list disp_global displacement vector for the element [tx1,ty1,rz1,tx2,ty2,rz2]


    :return mat Ke_global: element stiffness matrix [6 x 6]
    :return mat fe_int_global: element internal forces [6 x 1]
    """
    # Undeformed length and unit vector along element
    import numpy as np
import math

def beam2corot_Ke_and_Fe(ex, ey, ep, disp_global):
    """
    Compute the stiffness matrix and internal forces for a 2D corotational beam element.
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I]
    :param list disp_global: displacement vector [ux1, uy1, rz1, ux2, uy2, rz2]
    :return: Ke_global (6x6), fe_int_global (6,)
    """
    E, A, I = ep

    # --- Undeformed geometry ---
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = np.sqrt(eVec12 @ eVec12)
    nx0, ny0 = eVec12 / L0

    # --- Deformed geometry ---
    x1_def = ex[0] + disp_global[0]
    y1_def = ey[0] + disp_global[1]
    x2_def = ex[1] + disp_global[3]
    y2_def = ey[1] + disp_global[4]

    eVec_def = np.array([x2_def - x1_def, y2_def - y1_def])
    L = np.sqrt(eVec_def @ eVec_def)
    nx, ny = eVec_def / L

    # --- Rotation from global to local (current configuration) ---
    T_rot = np.array([
        [nx, ny, 0,   0,  0, 0],
        [-ny, nx, 0,  0,  0, 0],
        [0, 0, 1,     0,  0, 0],
        [0, 0, 0,    nx, ny, 0],
        [0, 0, 0,   -ny, nx, 0],
        [0, 0, 0,     0,  0, 1]
    ])

    # --- Local displacements ---
    disp_local = T_rot @ disp_global

    # Local displacements: [u1, v1, θ1, u2, v2, θ2]
    u1, v1, th1, u2, v2, th2 = disp_local

    # --- Axial deformation ---
    du = u2 - u1
    axial_strain = (L - L0) / L0
    N = E * A * axial_strain  # axial force

    # --- Local internal forces (Euler-Bernoulli beam) ---
    kL = E * I / (L ** 3)
    k_axial = E * A / L

    # Local tangent stiffness (linearized about current config)
    Ke_local = np.array([
        [ k_axial,    0,           0,   -k_axial,    0,           0],
        [ 0,     12*kL,   6*L*kL,     0,  -12*kL,   6*L*kL],
        [ 0,     6*L*kL,  4*L*L*kL,   0,  -6*L*kL,  2*L*L*kL],
        [-k_axial,   0,           0,    k_axial,    0,           0],
        [ 0,    -12*kL,  -6*L*kL,    0,   12*kL,  -6*L*kL],
        [ 0,     6*L*kL,  2*L*L*kL,   0,  -6*L*kL,  4*L*L*kL]
    ])

    # Local internal force vector (axial + bending)
    fe_int_local = Ke_local @ disp_local
    fe_int_local[0] += N
    fe_int_local[3] -= N

    # --- Transform back to global coordinates ---
    Ke_global = T_rot.T @ Ke_local @ T_rot
    fe_int_global = T_rot.T @ fe_int_local

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
