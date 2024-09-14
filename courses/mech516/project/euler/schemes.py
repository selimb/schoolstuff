"""
    Grid definition:
       node        edge
    |---o---|---o---|---o ... ---|---o---|
        0   0   1   1   2       J-1  J
"""

import numpy as np
import pdb
from riemann.riemann import Riemann
from euler import limiters


def roe(rho, u, p, dt, dx,  gamma=1.4):
    """Roe scheme for 1D Euler equations (no entropy fix).

    Calculates rho, u, p at next time step.

    Parameters
    ----------
        rho, u, p: 1xN ndarray
            N is the number of NODES.
    """
    # http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-17.pdf
    gam1 = gamma - 1

    ###############################
    ## Initialize Node Variables ##
    ###############################
    rho_e = p/(gamma - 1)  # Internal energy times density
    e = rho_e/rho
    m = u*rho  # Momentum rho*u
    rho_E = rho_e + 0.5*m*u  # Total energy rho*E
    h = e + p/rho
    H = h + u*u*0.5

    ##########################
    ## Calculate Roe Values ##
    ##########################
    # Fluxes are EDGE variables.
    # Pre-calculate Left and right values.
    rho_l = rho[:-1]
    rho_r = rho[1:]
    u_l = u[:-1]
    u_r = u[1:]
    p_l = p[:-1]
    p_r = p[1:]
    m_l = m[:-1]
    m_r = m[1:]
    H_l = H[:-1]
    H_r = H[1:]
    rho_E_l = rho_E[:-1]
    rho_E_r = rho_E[1:]
    sqrt_rho_l = np.sqrt(rho_l)
    sqrt_rho_r = np.sqrt(rho_r)

    # Calculate Roe averages
    u_av = (sqrt_rho_l*u_l + sqrt_rho_r*u_r)/(sqrt_rho_l + sqrt_rho_r)
    H_av = (sqrt_rho_l*H_l + sqrt_rho_r*H_r)/(sqrt_rho_l + sqrt_rho_r)
    c_av2 = gam1*(H_av - 0.5*u_av*u_av)
    c_av = np.sqrt(c_av2)
    f1av = 0.5*(m_l + m_r)
    f2av = 0.5*(m_l*m_l/rho_l + p_l + m_r*m_r/rho_r + p_r)
    f3av = 0.5*(m_l*H_l + m_r*H_r)
    M_av = u_av/c_av  # Mach number

    # Calculate necessary deltas
    delta_rho = rho_r - rho_l
    delta_m = m_r - m_l
    delta_rho_E = rho_E_r - rho_E_l

    alpha1 = (0.25*M_av*(2+gam1*M_av)*delta_rho -
              0.5/c_av*(1 + gam1*M_av)*delta_m +
              0.5*gam1/c_av2*delta_rho_E)
    alpha2 = ((1 - 0.5*gam1*M_av**2)*delta_rho +
              gam1*M_av/c_av*delta_m -
              gam1*M_av/c_av*delta_rho_E)
    alpha3 = (-0.25*M_av*(2-gam1*M_av)*delta_rho +
              0.5/c_av*(1 - gam1*M_av)*delta_m +
              0.5*gam1/c_av2*delta_rho_E)

    # Calculate lambdas (eigenvalues)
    lambda1 = np.abs(u_av - c_av)
    lambda2 = np.abs(u_av)
    lambda3 = np.abs(u_av + c_av)

    ######################
    ## Calculate Fluxes ##
    ######################
    roeflux1 = f1av - 0.5*(lambda1*alpha1 + lambda2*alpha2 + lambda3*alpha3)
    roeflux2 = f2av - 0.5*(lambda1*alpha1*(u_av - c_av) +
                           lambda2*alpha2*(u_av) +
                           lambda3*alpha3*(u_av + c_av))
    roeflux3 = f3av - 0.5*(lambda1*alpha1*(H_av - u_av*c_av) +
                           lambda2*alpha2*u_av*u_av +
                           lambda3*alpha3*(H_av + u_av*c_av))
    # print(roeflux1, roeflux2, roeflux3)

    ###################################
    ## Calculate new node primitives ##
    ###################################
    # Initialization
    rho_new = np.zeros_like(rho)
    m_new = np.zeros_like(m)
    rho_E_new = np.zeros_like(rho_E)

    # Set boundaries to remain the same
    rho_new[0] = rho[0]
    rho_new[-1] = rho[-1]
    m_new[0] = m[0]
    m_new[-1] = m[-1]
    rho_E_new[0] = rho_E[0]
    rho_E_new[-1] = rho_E[-1]

    # Calculate non-boundary primitives
    lamb = dt/dx
    for j in range(1, len(rho_new) - 1):
        rho_new[j] = rho[j] - lamb*(roeflux1[j] - roeflux1[j-1])
        m_new[j] = m[j] - lamb*(roeflux2[j] - roeflux2[j-1])
        rho_E_new[j] = rho_E[j] - lamb*(roeflux3[j] - roeflux3[j-1])

    u_new = m_new/rho_new
    p_new = gam1*(rho_E_new - 0.5*m_new*u_new)
    # pdb.set_trace()
    return rho_new, u_new, p_new


def mccormack(rho, u, p, dt, dx, gamma=1.4):
    # http://sitemaker.umich.edu/anand/files/project_2.pdf
    lamb = dt/dx
    # Boundary values are the same.
    rho_new = rho.copy()
    rho_new[1:-1] = 0.0
    u_new = u.copy()
    u_new[1:-1] = 0.0
    p_new = p.copy()
    p_new[1:-1] = 0.0

    # Update fluxes and conservative
    U = np.empty((len(rho), 3))
    F = U.copy()
    for i in range(len(U)):
        W_i = np.array([rho[i], u[i], p[i]])
        U[i] = W2U(W_i, gamma=gamma)
        F[i] = W2F(W_i, gamma=gamma)

    # Update predictor
    U_bar = U.copy()
    U_bar[1:-1] = U[1:-1] - lamb*(F[2:] - F[1:-1])

    # Update predictor fluxes
    F_bar = F.copy()
    for i in range(len(F_bar)):
        W_i = U2W(U_bar[i], gamma=gamma)
        F_bar[i] = W2F(W_i, gamma=gamma)

    # Final update
    U_new = U_bar.copy()
    assert((U_new[0] == U[0]).all())
    assert((U_new[-1] == U[-1]).all())
    U_new[1:-1] = 0.5*((U[1:-1] + U_bar[1:-1]) -
                       lamb*(F_bar[1:-1] - F_bar[:-2]))

    # Update primitives
    for i in range(1, len(rho_new)):
        U_i = U_new[i]
        rho_new[i], u_new[i], p_new[i] = U2W(U_i)

    # pdb.set_trace()
    return rho_new, u_new, p_new


def godunov(rho, u, p, dt, dx, gamma=1.4):
    # Boundary values are the same.
    rho_new = rho.copy()
    rho_new[1:-1] = 0.0
    u_new = u.copy()
    u_new[1:-1] = 0.0
    p_new = p.copy()
    p_new[1:-1] = 0.0

    for i in range(1, len(rho) - 1):
        # Left flux
        W_i = np.array([rho[i], u[i], p[i]])
        W_l = np.array([rho[i - 1], u[i - 1], p[i - 1]])
        W_r = W_i
        riemann = Riemann(W_l, W_r, gamma, x0=0)
        W = riemann.evaluate_single_state(0, 1)  # dx/dt = 0
        F_left = W2F(W, gamma=gamma)
        # Right flux
        W_l = W_i
        W_r = np.array([rho[i + 1], u[i + 1], p[i + 1]])
        riemann = Riemann(W_l, W_r, gamma, x0=0)
        W = riemann.evaluate_single_state(0, 1)
        F_right = W2F(W, gamma=gamma)
        # Compute at next time step
        U_i = W2U(W_i, gamma=gamma)
        U_new = U_i - dt/dx*(F_right - F_left)
        rho_new[i], u_new[i], p_new[i] = U2W(U_new)

    return rho_new, u_new, p_new


def MUSCL(rho, u, p, dt, dx, gamma=1.4, limiter='zero'):
    # Boundary values are the same.
    rho_new = rho.copy()
    rho_new[1:-1] = 0.0
    u_new = u.copy()
    u_new[1:-1] = 0.0
    p_new = p.copy()
    p_new[1:-1] = 0.0
    # Update fluxes and conservative
    W = np.empty((len(rho), 3))
    U = np.empty_like(W)
    for i in range(len(U)):
        W[i] = np.array([rho[i], u[i], p[i]])
        U[i] = W2U(W[i], gamma=gamma)

    ##################
    # Predictor Step #
    ##################
    # Compute limited slope
    delta_W = limiters.limiter(W, name=limiter)
    F_l = np.empty_like(delta_W)
    F_r = F_l.copy()
    U_tilde = F_l.copy()
    W_tilde = F_l.copy()
    for i in range(len(delta_W)):
        F_l[i] = W2F(W[i] + 0.5*delta_W[i])
        F_r[i] = W2F(W[i] - 0.5*delta_W[i])
        U_tilde[i] = U[i] - dt/dx*(F_l[i] - F_r[i])
        W_tilde[i] = U2W(U_tilde[i])

    ####################
    ## Corrector Step ##
    ####################
    # Godunov-like
    W_l = 0.5*(W + W_tilde + delta_W)
    W_r = 0.5*(W + W_tilde - delta_W)
    W_rs = np.zeros_like(W_l)
    F_rs = W_rs.copy()
    for i in range(0, len(W_l) - 1):
        # Left flux
        riemann = Riemann(W_l[i], W_r[i+1], gamma, x0=0)
        W_rs[i] = riemann.evaluate_single_state(0, 1)
        F_rs[i] = W2F(W_rs[i])

    U_new = U.copy()
    U_new[1:-1] = U[1:-1] - dt/dx*(F_rs[1:-1] - F_rs[0:-2])

    for i in range(len(U_new)):
        rho_new[i], u_new[i], p_new[i] = U2W(U_new[i])

    # pdb.set_trace()
    return rho_new, u_new, p_new


def W2F(W, gamma=1.4):
    """Compute fluxes vector F from primitives (rho, u, P)"""
    e = W[2]/(gamma - 1)/W[0]
    E = e + W[1]**2/2
    return np.array([W[0]*W[1],
                     W[0]*W[1]**2 + W[2],
                     W[0]*W[1]*(E + W[2]/W[0])])


def W2U(W, gamma=1.4):
    """Compute conservative vector U from primitives (rho, u, P)"""
    e = W[2]/(gamma - 1)/W[0]
    E = e + W[1]**2/2
    return np.array([W[0],
                     W[0]*W[1],
                     W[0]*E])


def U2W(U, gamma=1.4):
    """Compute primitive vector W from conservative"""
    rho = U[0]
    u = U[1]/rho
    E = U[2]/rho
    p = (gamma - 1)*rho*(E - 0.5*u**2)
    return np.array([rho, u, p])
