from dataclasses import dataclass

@dataclass(frozen=True)
class Contorno:
    tipo: str 
    valor: float

def coeficientes_contorno(contorno: Contorno, delta: float) -> tuple[float, float, float, float]:
    if contorno.tipo == "dirichlet":
        phi_b = contorno.valor
        alfa = 0.0
        beta = phi_b
        gamma = -1.0 / delta
        kappa = phi_b / delta
    elif contorno.tipo == "neumann":
        g = contorno.valor
        alfa = 1.0
        beta = g * delta
        gamma = 0.0
        kappa = g
    else:
        raise ValueError(f"Tipo de contorno desconocido: {contorno.tipo}")
    return alfa, beta, gamma, kappa

if __name__ == "__main__":
    print(coeficientes_contorno(Contorno("dirichlet", 100.0), 0.05))
    print(coeficientes_contorno(Contorno("neumann", 0.0), 0.05))

