
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
    exvec0 = np.array([dx/L0, dy/L0])

    ex_def = ex + [disp_global[0],disp_global[3]]
    ey_def = ey + [disp_global[1],disp_global[4]]
    dx_def = ex_def[1] - ex_def[0]
    dy_def = ey_def[1] - ey_def[0]
    Ld = np.sqrt(dx_def**2 + dy_def**2)   
    exvecD = np.array([dx_def/Ld, dy_def/Ld])
    eyvecD = np.array([-exvecD[1],exvecD[0]])

    displ_local = np.zeros(6)

    for i in range(2):
        R = rot_matrix(disp_global[i*3+2])
        tvec = R @ exvec0
        theta_d = math.asin(eyvecD.dot(tvec))
        displ_local[i*3+2] = theta_d
        
    delta_L = Ld - L0

    displ_local[0] = -delta_L *0.5
    displ_local[3] = delta_L * 0.5

    Ke_local = beam2local_stiff(L0,ep)

    #f_int_local = Ke_local @ displ_local
    # --- Geometrisk stivhet (lokal), basert på DEFORMERT lengde ---
    # Deformert lengde Ld via de deformer­te nodale koordinatene (ex_def, ey_def)
    Ld = float(np.hypot(ex_def[1] - ex_def[0], ey_def[1] - ey_def[0]))
    if Ld == 0.0:
        Ld = L0  # fallback

    E, A, I = ep

    # Aksialkraft i lokal akse (liten tøyning, store rotasjoner – korotasjon)
    N = E * A * (Ld - L0) / L0

    # Konsistent geometrisk stivhet for 2D ramme (Euler–Bernoulli) i lokal basis
    coeff = N / (30.0 * Ld)
    Kg_local = coeff * np.array([
        [0.0,    0.0,        0.0,   0.0,     0.0,         0.0],
        [0.0,   36.0,    3.0*Ld,   0.0,   -36.0,     3.0*Ld],
        [0.0, 3.0*Ld,  4.0*Ld*Ld,  0.0,  -3.0*Ld,  -1.0*Ld*Ld],
        [0.0,    0.0,        0.0,   0.0,     0.0,         0.0],
        [0.0,  -36.0,   -3.0*Ld,   0.0,    36.0,    -3.0*Ld],
        [0.0, 3.0*Ld, -1.0*Ld*Ld,  0.0,  -3.0*Ld,   4.0*Ld*Ld]
    ], dtype=float)

    # Total lokal tangent
    Kt_local = Ke_local + Kg_local 
    f_int_local = Ke_local @ displ_local


    # Blokk-diagonal transformasjon for 2D element (3x3 per node)

    T = beam2corot_Te(ex_def,ey_def)

    # Global stivhetsmatrise og interne krefter
    Ke_global = T.T @ Kt_local @ T
    f_int_local = Ke_local @ displ_local
    fe_int_global = T.T @ f_int_local

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

