import pandas as pd
import numpy as np

#The line data should be inserted in per-unit (p.u.) values of Ohms, 
#while the other data should be in p.u. values of admittance.

class Barra:
    def __init__(self, bar_number, typebar, p, q, v, theta, cshunt):
        self.bar_number = bar_number
        self.typebar = typebar #1 is VTheta, 2 is PV, 3 is PQ
        self.P = p
        self.Q = q
        self.V = v
        self.Theta = theta
        self.Cshunt = cshunt

class System:
    def __init__(self, barras, line_data, admittance_matrix):
        self.barras = barras
        self.line_data = line_data
        self.admittance_matrix = admittance_matrix

    def calculate_line_flux(self):
        n_bars = len(self.barras)
        P = np.zeros((n_bars, n_bars))
        Q = np.zeros((n_bars, n_bars))
        for _, line in self.line_data.iterrows():
            origin = int(line.origin) - 1
            destiny = int(line.destiny) - 1
            bar_a = self.barras[origin]
            bar_b = self.barras[destiny]

            Va = bar_a.V
            Vb = bar_b.V
            b_sh = line.Xshunt
            Theta_a = bar_a.Theta
            Theta_b = bar_b.Theta
            Theta_ab= Theta_a - Theta_b
            Z_line = complex(line.R_line, line.X_line)
            Y_line = 1 / Z_line
            g_ab= Y_line.real
            b_ab = Y_line.imag
            
            P[origin, destiny] = Va ** 2 * g_ab - Va * Vb * g_ab * np.cos(Theta_ab) - Va * Vb * b_ab * np.sin(Theta_ab)
            Q[origin, destiny] = -Va ** 2 * (b_ab + b_sh) + Va * Vb * b_ab * np.cos(Theta_ab) - Va * Vb * g_ab * np.sin(Theta_ab)

        return P, Q
        
    def calculate_submatrices(self):
        n_bars = len(self.barras)
        H = np.zeros((n_bars, n_bars))
        L = np.zeros((n_bars, n_bars))
        M = np.zeros((n_bars, n_bars))
        N = np.zeros((n_bars, n_bars))

        P, Q = calculate_power(self)

        for i, bar_i in enumerate(self.barras):
            for j, bar_j in enumerate(self.barras):
                if i != j and self.admittance_matrix[i, j] != 0:
                    G_ij = self.admittance_matrix[i, j].real
                    B_ij = self.admittance_matrix[i, j].imag

                    theta_ij = bar_i.Theta - bar_j.Theta

                    H[i, j] = bar_i.V * bar_j.V * (G_ij * np.sin(theta_ij) - B_ij * np.cos(theta_ij))
                    N[i, j] = bar_i.V * (G_ij * np.cos(theta_ij) + B_ij * np.sin(theta_ij))
                    M[i, j] = -bar_i.V * bar_j.V * (G_ij * np.cos(theta_ij) + B_ij * np.sin(theta_ij))
                    L[i, j] = bar_i.V * (G_ij * np.sin(theta_ij) - B_ij * np.cos(theta_ij))
                     
        for k in range(n_bars):
            G_kk = self.admittance_matrix[k, k].real
            B_kk = self.admittance_matrix[k, k].imag
            V_k = self.barras[k].V
            if self.barras[k].typebar == 1:  # V Theta bar
                H[k, k] = 1e15
                L[k, k] = 1e15
                N[k, k] = V_k * G_kk + P[k] / V_k
                M[k, k] = -(V_k)**2 * G_kk + P[k]
            elif self.barras[k].typebar == 2:  # PV bar
                L[k, k] = 1e15
                H[k, k] = -V_k**2 * B_kk - Q[k]
                N[k, k] = V_k * G_kk + P[k] / V_k
                M[k, k] = -(V_k)**2 * G_kk + P[k]
            else:
                H[k, k] = -V_k**2 * B_kk - Q[k]
                N[k, k] = V_k * G_kk + P[k] / V_k
                M[k, k] = -(V_k)**2 * G_kk + P[k]
                L[k, k] = -V_k * B_kk + Q[k] / V_k
                
        J = np.block([
            [H[0:, 0:], N[0:, 0:]],
            [M[0:, 0:], L[0:, 0:]]
        ])
        
        # Find the indices of the diagonal elements equal to 1e15
        k = np.where(np.diag(J) == 1e15)[0]
        J_filtered = np.delete(np.delete(J, k, axis=0), k, axis=1)

        return J_filtered

def read_csv_data():
    dadosbarra = pd.read_csv("dadosbarra.csv")
    dadoslinha = pd.read_csv("dadoslinha.csv")
    print("Dadosbarra column names:", dadosbarra.columns)
    return dadosbarra, dadoslinha

def create_system(dadosbarra, dadoslinha):
    barras = []
    #Definição do flat start
    for _, row in dadosbarra.iterrows():
        if row.typebar == 1:  # V Theta [referencia]
            V = row.V
            Theta = row.Theta
        elif row.typebar == 2:  # PV
            V = row.V
            Theta = 0.0
        else:  # PQ
            V = 1.0
            Theta = 0.0

        bar = Barra(row.bar_number, row.typebar, float(row.P), float(row.Q), V, Theta, row.Cshunt)
        barras.append(bar)

    return System(barras, dadoslinha, None)

def calculate_admittance_matrix(system):
    n_bars = len(system.barras)
    admittance_matrix = np.zeros((n_bars, n_bars), dtype = np.complex128)

    for _, line in system.line_data.iterrows():
        origin = int(line.origin) - 1
        destiny = int(line.destiny) - 1
        z_line = complex(line.R_line, line.X_line)
        y_line = 1 / z_line
        tap_ratio = line.TapValue if line.TapValue != 0 else 1  # Get the Tap ratio from line_data

        # Diagonal elements
        admittance_matrix[origin, origin] += y_line + 1j * line.Xshunt + 1j * system.barras[origin].Cshunt
        admittance_matrix[destiny, destiny] += y_line / (tap_ratio**2) + 1j * line.Xshunt + 1j * system.barras[destiny].Cshunt 

        # Off-diagonal elements
        admittance_matrix[origin, destiny] -= y_line / tap_ratio
        admittance_matrix[destiny, origin] -= y_line / tap_ratio
        
    #print(admittance_matrix)

    return admittance_matrix

def calculate_power(system):
    num_bars = len(system.barras)

    P = np.zeros(num_bars, dtype=np.float64)
    Q = np.zeros(num_bars, dtype=np.float64)

    for idx1, bar1 in enumerate(system.barras):
        V_i = np.float64(bar1.V)
        theta_i = np.float64(bar1.Theta)

        # Calculate the diagonal terms (i == j)
        G_ii = np.float64(system.admittance_matrix[idx1, idx1].real)
        B_ii = np.float64(system.admittance_matrix[idx1, idx1].imag)
        P[idx1] += V_i * V_i * G_ii
        Q[idx1] -= V_i * V_i * B_ii

        # Calculate the off-diagonal terms (i != j)
        for idx2, bar2 in enumerate(system.barras):
            if idx1 == idx2:
                continue

            V_j = np.float64(bar2.V)
            theta_j = np.float64(bar2.Theta)
            theta_diff = np.float64(theta_i - theta_j)
            G_ij = np.float64(system.admittance_matrix[idx1, idx2].real)
            B_ij = np.float64(system.admittance_matrix[idx1, idx2].imag)

            P[idx1] += V_i * V_j * (G_ij * np.cos(theta_diff) + B_ij * np.sin(theta_diff))
            Q[idx1] += V_i * V_j * (G_ij * np.sin(theta_diff) - B_ij * np.cos(theta_diff))

    return P, Q

def newton_raphson_method(system, max_iterations=10, tolerance=3e-3):
    print("\n \n \nCondições iniciais: ")
    num_buses = len(system.barras)
    display_results(system, 0)
    for iteration in range(max_iterations):
        P_calc, Q_calc = calculate_power(system)
        #print("Pcalc, Qcalc=", P_calc," ",Q_calc)
        # Create a boolean mask for bar.typebar != 1 and use it to calculate P_mismatch
        mask_P = np.array([bar.typebar != 1 for bar in system.barras])
        P_mismatch = np.zeros(len(system.barras), dtype=np.float64)
        P_mismatch[mask_P] = np.array([bar.P - P_calc[i] for i, bar in enumerate(system.barras) if bar.typebar != 1], dtype=np.float64)
        
        # Create a boolean mask for bar.typebar == 3 and use it to calculate Q_mismatch
        mask_Q = np.array([bar.typebar == 3 for bar in system.barras])
        Q_mismatch = np.zeros(len(system.barras), dtype=np.float64)
        Q_mismatch[mask_Q] = np.array([bar.Q - Q_calc[i] for i, bar in enumerate(system.barras) if bar.typebar == 3], dtype=np.float64)
        
        #Filter the numpy array
        P_mismatch_filtered = P_mismatch[mask_P]
        Q_mismatch_filtered = Q_mismatch[mask_Q]

        mismatch = np.hstack([P_mismatch_filtered[0:], Q_mismatch_filtered[0:]])

        if np.all(np.abs(mismatch)<tolerance):
            print("O sistema convergiu!!!!!!!!!!!!!!!!!!!!")
            P, Q = system.calculate_line_flux()
            print("\nFluxo de linha: \n")
            for i in range(num_buses):
                for j in range(i, num_buses):
                    if i != j:
                        print("Fluxo Linha ", i + 1, " para linha ", j + 1, " P: ", P[i][j], " Q: ", Q[i][j], "\n")
            return system, iteration
            break

        J = system.calculate_submatrices()
        #print("J: \n", J)
        inv=np.linalg.inv(J)
        delta=np.matmul(inv, mismatch)
        
        #Split delta array        
        delta_theta = delta[:len(P_mismatch_filtered)]
        delta_v = delta[len(P_mismatch_filtered):]
        #print(delta_theta)
        # Get indices of PV and PQ buses
        pv_bus_indices = [i for i, bar in enumerate(system.barras) if bar.typebar == 2]
        pq_bus_indices = [i for i, bar in enumerate(system.barras) if bar.typebar == 3]

        # Update angles and voltages for PV and PQ buses
        for i, idx in enumerate(pv_bus_indices):
            system.barras[idx].Theta += delta_theta[i]
            if system.barras[idx].Theta>2*np.pi:
                system.barras[idx].Theta=system.barras[idx].Theta-2*np.pi
        for i, idx in enumerate(pq_bus_indices):
            print(i)
            #print(system.barras[idx].Theta)
            system.barras[idx].Theta += delta_theta[i+len(pv_bus_indices)]
            #system.barras[idx].Theta += delta_theta[i+len(pv_bus_indices)]
            system.barras[idx].V += delta_v[i]
            if system.barras[idx].Theta>2*np.pi:
                system.barras[idx].Theta=system.barras[idx].Theta-2*np.pi
        
        #Details of Magnitude and Angle at each iteration
        print("Number of iterations: {}".format(iteration + 1))
        print("-----------------------------------------------")
        print("Bus Number | Voltage Magnitude | Voltage Angle ")
        print("-----------------------------------------------")
        for bar in system.barras:
            print("{:<11} | {:<17.5f} | {:<13.5f}  ".format(
                bar.bar_number, bar.V, bar.Theta))
        
    return None, max_iterations

def display_results(system, num_iterations):
    print("Number of iterations: {}".format(num_iterations))
    print("--------------------------------------------------------------------------------")
    print("Bus Number | Voltage Magnitude | Voltage Angle(degrees) | Active Power | Reactive Power")
    print("--------------------------------------------------------------------------------")
    for bar in system.barras:
        print("{:<11} | {:<17.5f} | {:<20.5f}° | {:<12.5f} | {:<14.5f}".format(
            bar.bar_number, bar.V, bar.Theta*180/np.pi, bar.P, bar.Q))
        
def update_power_values(system, P_calc, Q_calc):
    for idx, bar in enumerate(system.barras):
        if bar.typebar == 1:
            bar.P = P_calc[idx]
            bar.Q = Q_calc[idx]
        elif bar.typebar == 2:
            bar.Q = Q_calc[idx]
        else:
            continue
            
def gauss_seidel_initial_values(system, tol=1e-6, max_iter=100):
    n_bars = len(system.barras)
    V_prev = np.array([bar.V for bar in system.barras], dtype=complex)

    for iteration in range(max_iter):
        V_new = V_prev.copy()
        for i, bar in enumerate(system.barras):
            if bar.typebar == 3:  # PQ bus
                sum1 = sum(system.admittance_matrix[i, m] * V_new[m] for m in range(i))
                sum2 = sum(system.admittance_matrix[i, m] * V_prev[m] for m in range(i + 1, n_bars))
                V_new[i] = 1 / system.admittance_matrix[i, i] * ((bar.P - 1j * bar.Q) / np.conj(V_prev[i]) - sum1 - sum2)

            elif bar.typebar == 2:  # PV bus
                Saux = np.conj(V_prev[i]) * (np.dot(system.admittance_matrix[i, :], V_prev))
                Q = -(Saux.imag)
                sum1 = sum(system.admittance_matrix[i, m] * V_new[m] for m in range(i))
                sum2 = sum(system.admittance_matrix[i, m] * V_prev[m] for m in range(i + 1, n_bars))
                V_provisorio = 1 / system.admittance_matrix[i, i] * ((bar.P - 1j * Q) / np.conj(V_prev[i]) - sum1 - sum2)
                V_new[i] = np.abs(V_prev[i]) * np.exp(1j * np.angle(V_provisorio))

        if np.allclose(V_new, V_prev, rtol=tol, atol=tol):
            break

        V_prev = V_new.copy()

    for i, bar in enumerate(system.barras):
        bar.V = np.abs(V_new[i])
        bar.Theta = np.angle(V_new[i])

    return system


def main():
    dadosbarra, dadoslinha = read_csv_data()
    system = create_system(dadosbarra, dadoslinha)
    admittance_matrix = calculate_admittance_matrix(system)
    system = System(system.barras, system.line_data, admittance_matrix)
    print("Matriz Y: \n \n", admittance_matrix)
    print("Você quer usar Gauss-Seidel para ter melhores condições iniciais?\n")
    a = int(input("0- Não, 1- Sim: "))
    
    # Apply Gauss-Seidel method for better initial values
    if a==1:
        gauss_seidel_initial_values(system, tol=1e-6, max_iter=1)
    #print(system.admittance_matrix)
    updated_system, num_iterations = newton_raphson_method(system, max_iterations=10, tolerance=3e-3)
    
    if updated_system:
        P_calc, Q_calc = calculate_power(updated_system)
        update_power_values(updated_system, P_calc, Q_calc)
        display_results(updated_system, num_iterations)
        
    else:
        print("O sistema não convergiu após ", num_iterations, "iterações")
        
if __name__ == "__main__":
    main()
