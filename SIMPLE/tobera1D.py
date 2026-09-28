import numpy as np
from tdma import TDMA

def tobera1D(n=100, alphap=0.2, alphau=0.8, max_iter=1_000, tol=1e-6):

    L = 1
    deltax = L / (n + 1/2)
    AA = 0.5
    AB = 0.1
    pB = 1.0
    uA = 1.0
    rhoA = 1.0
    rhoB = 1.0
    muA = 1.0
    muB = 1.0
    Su = 0.0

    i = np.arange(n)

    xu = deltax + deltax * i
    xP = deltax / 2 + deltax * i

    Au = AA - (AA - AB) / L * xu
    AP = AA - (AA - AB) / L * xP

    p = np.full(n, pB)
    u = np.full(n, uA)

    rho = np.full(n, rhoA)
    rho_ext = np.append(rho, rhoB)
    rhou = (rho_ext[0:n] + rho_ext[1:n+1]) / 2
    mu = np.full(n, muA)

    DW = mu * AP / deltax
    DE = np.zeros(n)
    DE[0:n-1] = mu[1:n] * AP[1:n] / deltax
    DE[n-1] = 0

    p_corr = 0.5
    index = 0

    res_masa = []

    divergido = False 

    while np.max(np.abs(p_corr)) > tol and index < max_iter:

        FW = np.zeros(n)
        FE = np.zeros(n)
        Sp = np.zeros(n)
        b = np.zeros(n)

        for k in range(n):

            if k == 0:
                FW[k] = rho[k] * AP[k] * (u[k] + uA) / 2
            else:
                FW[k] = rho[k] * AP[k] * (u[k] + u[k-1]) / 2
            if k == n-1:
                r = rhou[k] * Au[k] / rhoB / AB
                FB = rhoB * AB * r * u[k]
                Sp[k] = - FB * r
            else:
                FE[k] = rho[k+1] * AP[k+1] * (u[k+1] + u[k]) / 2

            if k == n-1:
                b[k] = -(pB - p[k]) * Au[k] + Su * Au[k] * deltax + muB * AB * (r-1) * u[k] / (deltax/2)
            elif k == 0:
                b[k] = (1/2 * FW[0] + DW[0]) * uA - (p[1] - p[0]) * Au[0] + Su * Au[0] * deltax
            else:
                b[k] = -(p[k+1] - p[k]) * Au[k] + Su * Au[k] * deltax

        Sp[0] = - (DW[0] + 1/2 * FW[0])

        aW = DW + 1 / 2 * FW
        aW[0] = 0
        aE = DE - 1 / 2 * FE
        aE[n-1] = 0
        aP = aW + aE + (FE - FW) - Sp

        # Subrelajación de las velocidades
        aP = aP / alphau
        b = b + aP * (1-alphau) * u

        A = np.zeros((n,n))

        for k in range(n):
            A[k,k] = aP[k]
            if k-1 >= 0: A[k, k-1] = - aW[k]
            if k+1 < n: A[k, k+1] = - aE[k]

        u_star = TDMA(A, b)

        # ---- Coeficientes de presión ----

        aEp = np.zeros(n)
        aWp = np.zeros(n)
        bp = np.zeros(n)
        Spp = np.zeros(n)
        d = Au / aP

        for k in range(n):

            aEp[k] = 0 if k == n-1 else rhou[k] * d[k] * Au[k]

            aWp[k] = 0 if k == 0 else rhou[k-1] * Au[k-1] * d[k-1]
            if k == 0:
                bp[k] = rhoA * uA * AA - rhou[k] * u_star[k] * Au[k]
            else:
                bp[k] = rhou[k-1] * u_star[k-1] * Au[k-1] - rhou[k] * u_star[k] * Au[k]

            Spp[k] = - rhou[k] * Au[k] * d[k] if k == n-1 else 0

        aPp = aEp + aWp - Spp

        Mp = np.zeros((n,n))

        for k in range(n):
            Mp[k,k] = aPp[k]
            if k-1 >= 0: Mp[k, k-1] = - aWp[k]
            if k+1 < n: Mp[k, k+1] = - aEp[k]

        p_corr = TDMA(Mp, bp)

        p_corr_ext = np.append(p_corr, 0)

        u = u_star + d * (p_corr_ext[0:n] - p_corr_ext[1:n+1])
        # Subrelajación de las presiones
        p = p + alphap * p_corr

        res_masa.append(np.sum(np.abs(bp)))

        index += 1

        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(p)):
            divergido = True
            break

    estado = "divergió" if divergido else "se acabaron las iteraciones" if index >= max_iter else "convergió"

    res_masa = np.array(res_masa) / (rhoA * uA* AA)

    resultado = {
        "u" : u,
        "p" : p,
        "índice" : index,
        "residuos de masa" : res_masa,
        "estado de convergencia": estado,
        "xu" : xu, 
        "xP" : xP,
        "alphap" : alphap,
        "alphau" : alphau
    }

    return resultado

if __name__ == "__main__":
    resultado = tobera1D()
    print(f'Con alphap = {resultado["alphap"]} y alphau = {resultado["alphau"]} {resultado["estado de convergencia"]}')
    print(f"Finaliza en la iteración {resultado["índice"]}")

