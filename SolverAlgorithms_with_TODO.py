import math
import numpy as np






def solveArchLength(model, archLength=0.02, max_steps=50, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)
    res_Vec = np.zeros(num_dofs)
    Lambda = 0.0

    d_q_prev = np.zeros(num_dofs)

    for iStep in range(max_steps):

        #TODO: Implement this
        print(f"\n=== Arc‑length step {iStep} ===")

        q = model.get_incremental_load(1.0)
        K = model.get_K_sys(uVec)
        dU = np.linalg.solve(K, q)
        dU_norm = np.sqrt(dU.dot(dU))
        dU /= dU_norm
        dLambda = archLength / dU_norm
        u_k = uVec + dU * dLambda
        Lambda_k = Lambda + dLambda

        for iIter in range(max_iter):

            # TODO: Implement this

            res_Vec = model.get_residual(uVec, Lambda)
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                print(f"  Converged in {iIter+1} iterations.")
                break

            K = model.get_K_sys(u_k)
            q = model.get_incremental_load(1.0)

            
            dx = np.linalg.solve(K, -res)
            dv = np.linalg.solve(K, q)

            
            a = dv.dot(dU)
            b = 2.0 * (dx.dot(dU) + (Lambda_k - Lambda) * a)
            c = dx.dot(dx) + (Lambda_k - Lambda)**2 * a**2 - archLength**2
            disc = b**2 - 4*a*c
            dLambda = (-b + math.copysign(math.sqrt(abs(disc)), dLambda)) / (2*a)
            dU_iter = dx + dv * dLambda

            u_k += dU_iter
            Lambda_k += dLambda

    
        uVec = u_k
        Lambda = Lambda_k


        model.append_solution(Lambda, uVec)
        print("Linear arc step {:}  load_factor= {:12.3e}".format(iStep, Lambda))

def solveNonlinLoadControl(model, load_steps=0.01, max_steps=100, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec   = np.zeros(num_dofs)
    d_uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        #TODO: Implement this
        print(f"\n=== Load step {iStep}, Lambda = {Lambda:.3f} ===")

        res_Vec = model.get_residual(uVec, Lambda)

        for iIter in range(max_iter):

            # TODO: Implement this
            K_mat = model.get_K_sys(uVec)

            du = np.linalg.solve(K_mat, -res_Vec)
            uVec += du

            res_Vec = model.get_residual(uVec, Lambda)
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


