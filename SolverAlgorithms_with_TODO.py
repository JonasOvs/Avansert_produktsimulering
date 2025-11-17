import math
import numpy as np


def solveArchLength(model, archLength=0.02, max_steps=50, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)
    res_Vec = np.zeros(num_dofs)
    Lambda = 0.0

    w_q0_prev = np.zeros(num_dofs)
    d_q_prev = np.zeros(num_dofs)

    for iStep in range(max_steps):

        #TODO: Implement this
        # residual og stivhet i nåværende likevektspunkt
        res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)
        K_mat   = model.get_K_sys(uVec)

        # inkrementell lastvektor q(λ) ≈ dp/dλ via differanse i residual
        # (for proporsjonal last blir dette en konstant lastvektor)
        q_vec = model.get_incremental_load(loadFactor=Lambda)

        # prediktor-retning w_q0 fra K w_q0 = q(λ)  (Haugen eq. 3.2.6)
        w_q0 = np.linalg.solve(K_mat, q_vec)

        f = math.sqrt(1.0 + w_q0.dot(w_q0)) 
        dLambda = archLength / f

        if w_q0.dot(d_q_prev) < 0.0:
            dLambda = -dLambda
  
        du_pred = dLambda * w_q0
        d_q_prev = du_pred

        # gå til prediktor-konfigurasjon
        Lambda += dLambda
        uVec   += du_pred

        # start-residual for korrektor-iterasjonene
        #res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)

        bConverged = False

        for iIter in range(max_iter):

            # TODO: Implement this
            K_mat = model.get_K_sys(uVec)

            # inkrementell lastvektor q(λ)
            res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)
            q_vec = model.get_incremental_load(loadFactor=Lambda)
            #q_vec = model.get_residual(loadFactor=Lambda-1.0, disp_sys=uVec) - res_Vec

            # løs K w_r = r  og  K w_q = q
            w_r = np.linalg.solve(K_mat, res_Vec)
            w_q = np.linalg.solve(K_mat,  q_vec)

            # δλ = - (w_q^T w_r) / (1 + w_q^T w_q)
            num   = np.dot(w_q, w_r)
            denom = 1.0 + np.dot(w_q, w_q)
            deltaLambda = - num / denom

            # δv = w_r + δλ w_q
            du = w_r + deltaLambda * w_q

            # oppdater til neste residual-evaluering
            uVec   += du
            Lambda += deltaLambda

            res_Vec = model.get_residual(disp_sys=uVec, loadFactor=Lambda)
            eps = res_Vec.dot(res_Vec)
            print("     Iter {:}  eps = {:12.3e}".format(iIter, eps))
            if (eps < 1.0e-11):
                bConverged = True
                break

        model.append_solution(Lambda, uVec)
        print("Linear arc step {:}  load_factor= {:12.3e}".format(iStep, Lambda))
        if not bConverged:
            print("  WARNING: No convergence in arc-length step {:}".format(iStep))
            break


def solveNonlinLoadControl(model, load_steps=0.001, max_steps=100, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec   = np.zeros(num_dofs)
    d_uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):
        Lambda = load_steps * iStep
        #TODO: Implement this: jonas fikser
        d_uVec[:] = 0
        # res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)
        # print("res_vec", res_Vec)
        #TODO slutt, jonas fikset
        for iIter in range(max_iter):

            # TODO: Implement this
            print(uVec)
            res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)
            K_mat = model.get_K_sys(uVec)
            d_uVec = np.linalg.solve(K_mat, res_Vec)

            uVec += d_uVec
            # TODO slutt, jonas fikset

            res_Vec = model.get_residual(loadFactor=Lambda, disp_sys=uVec)
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                break

        model.append_solution(Lambda, uVec)
        print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))



def solveLinearSteps(model, load_steps=0.01, max_steps=100):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        q_Vec   = model.get_incremental_load(Lambda)

        K_mat = model.get_K_sys(uVec)

        d_q = np.linalg.solve(K_mat, q_Vec)

        uVec = d_q * Lambda

        model.append_solution(Lambda, uVec)
        print("Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))
