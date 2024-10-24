
#%% Definições

# - Beta_k = ganho médio de potência no canal entre o PB e o k-ésimo dispositivo;
# - s_k = sinal transmitido (variável aleatória gaussiana de média zero e variancia unitária);
# - Psi_k = vetor beamforming analógico (Pertence aos complexos com tamanho 1xN);
# - ||Psi_k||² = Potencia de transmissão (Pt);
# - h_k = Vetor de canal (pertence aos complexos com tamanho 1xN);
# - kappa = fator de rice;

# OBS: k = K. j = 1 e vai até j = k-1


#%%

import numpy as np
import matplotlib.pyplot as plt
import random


K = 50                                  # N° dispositivos
N = 10                                  # N° de antenas
R = 10                                  # Raio [m]
f = 915 *(10**6)                        # Freq. da portadora
kappa = 1.5                             # Fator Rice
mu = 10.73 * 10**(-3)                   # Max pot no IoT
a = 0.2308                              # constante circuito a
b = 5.365                               # constante circuito b  
Pt = 3                                  # Max pot. de transmissão
E_min = 200 * (10**-6)                # Min Energia no IoT [J]
alpha = 2.7                             # Coef. perda de percurso
c = 3*(10**(8))                         # Vel. da luz
Omega = 1/(1+np.exp(a*b))               # Constante resposta 0_in-0_out
Pot_k = 0
tau_k_estrela_k_igual_1_valor = 0
tau_k_valor = 0


def dist_euclid_quad (x, y):
    d_k = np.sqrt(((x**2)+(y**2))) 
    return d_k

def Beta_k(d_k):
    Beta_k = (c**2)/(16 * (np.pi**2) * (f**2) * (d_k**alpha))
    return Beta_k
# Beta_k = Beta_k(1)

def h_barra_k(): # Gerador componente LoS
    h_barra_k = np.random.rand() + 1j*np.random.rand()
    return h_barra_k
# h_barra_k = h_barra_k()

def h_til_k(): # Gerador de componente NLoS
    h_til_k = (np.random.rand() + 1j*np.random.rand())*(1/(np.sqrt(2)))
    return h_til_k
# h_til_k = h_til_k()

def canal_h_k (kappa, h_barra_k, h_til_k): # Gerador de canal
    canal_h_k = (np.sqrt((kappa)/(1+kappa)) * h_barra_k) + (np.sqrt((1)/(1+kappa)) * h_til_k)
    return canal_h_k
# canal_h_k = np.append(canal_h_k, canal_h_k(1.5,h_barra_k, h_til_k))

def P_k(Beta_k, Psi_k, canal_h_k_Hermitiano): # Potencia de RF recebida
    P_k = Beta_k * (np.abs((Psi_k) *(canal_h_k_Hermitiano)))**2
    return P_k

def Gamma_k(P_k): # Parâmetro para calcular a Energia coletada
    Gamma_k = mu / (1 + np.exp(-a*(P_k-b)))
    return Gamma_k

def Phi_k(tau_k, Gamma_k): # Energia coletada
    Phi_k = tau_k * ((Gamma_k - (mu*Omega))/(1-Omega))
    return Phi_k

# def Gamma_k_j(Beta_k, Psi_j, canal_h_k_Hermitiano):
#     Gamma_k_j = mu/(1 + (np.exp(-a*(Beta_k * ((np.abs(Psi_j)*np.abs(canal_h_k_Hermitiano))**2)-b))))
#     return Gamma_k_j

# def Phi_NEIG_j(tau_j, Gamma_k_j):
#     Phi_NEIG_j = tau_j*(Gamma_k_j - (mu*Omega))/(1-Omega)
#     return Phi_NEIG_j

# def Phi_linha_k(Phi_NEIG_j, Phi_k):
#     Phi_linha_k = Phi_NEIG_j + Phi_k
#     return Phi_linha_k

def tau_k_estrela_k_igual_1(Gamma_k):
    tau_k_estrela_k_igual_1 = (E_min * (1 - Omega))/(Gamma_k - (mu*Omega))
    return tau_k_estrela_k_igual_1

# def psi_k_estrela(h_barra_k,):
#     Psi_k_estrela = (np.sqrt(Pt/N))*(h_barra_k/np.abs(h_barra_k))
#     return Psi_k_estrela

# def tau_k_estrela_k_maior_1(Phi_NEIG_j, Gamma_k):
#     for u in range j{
#         Phi_NEIG_j(tau_j)


# def sinal_recebido_y (Beta_k, Psi_k, canal_h_k, sinal_s_k): # Sinal recebido
#     canal_h_k_hermitiano = np.transpose(np.conjugate(canal_h_k))
#     y_k = np.sqrt(Beta_k)*Psi_k * canal_h_k_hermitiano * sinal_s_k
#     return y_k




#%%
canal_h = np.array([])
h_barra = np.array([])
h_til = np.array([])
Beta = np.array([])
Phi = np.array([])
tau = np.array([])
Psi_k_estrela = np.array([])
Energia = np.array([])
h_barra_dividido_pela_norma_vetor = np.array([])
h_barra_vetor = np.array([])
P = np.array([])
Gamma = np.array([])
Gamma_k_j = np.array([])
Phi_NEIG_j = np.array([])
tau_k_estrela_maior_que_1 = np.array([])
tau_k_estrela = np.array([])
Phi_k_vetor = np.array([])
Phi_linha_k = np.array([])
tau_total = np.array([])
tau_total_vetor = np.array([])


seed = np.random.seed(9)

for q in range (0,K):

    for k in range (0, K):
        # Canal h
        h_barra = np.append(h_barra, h_barra_k())
        h_barra_vetor = np.append(h_barra_vetor, h_barra)
        h_til = np.append(h_barra, h_til_k())
        canal_h = np.append(canal_h, canal_h_k(kappa, h_barra_k(), h_til_k()))
        h_k_hermitiano = np.conjugate(np.transpose(canal_h))

        # Beta
        x = np.random.randint(1,10)
        y = np.random.randint(1,10)
        Beta = np.append(Beta, Beta_k(dist_euclid_quad(x,y)))

        #Psi_k_estrela
        h_barra_dividido_pela_norma = h_barra_vetor[k] / np.abs(h_barra_vetor[k])
        h_barra_dividido_pela_norma_vetor = np.append(h_barra_dividido_pela_norma_vetor, h_barra_dividido_pela_norma)
        h_barra_dividido_pela_norma_transposto = np.transpose(h_barra_dividido_pela_norma_vetor)
        Psi_k_estrela = (np.sqrt(Pt/N)) * h_barra_dividido_pela_norma_transposto

    for kk in range (0, K):
        Pot_k = Beta[kk] * (np.abs(Psi_k_estrela[kk] * h_k_hermitiano[kk]))**2
        P = np.append(P, Pot_k)


    # print("Funcionou até aqui :)")

    for k in range (0, K):
        # Gamma_k
        Gamma = np.append(Gamma, Gamma_k(P[k]))

        # Phi_k (Energia coletada)
        # Precisa-se calcular tau_k_estrela que é diferente para k = 1 e k>1
        if k==0:
            tau_k_estrela_k_igual_1_valor = tau_k_estrela_k_igual_1(Gamma[0])
            tau_k_estrela = np.append(tau_k_estrela, tau_k_estrela_k_igual_1_valor)

    for j in range (0, K-1)  :
        # for j in range (0, K-1):
            # Gamma k,j:
            Gamma_k_j_valor = mu / (1 + (np.exp(-a * ((Beta[k] * ((np.absolute(Psi_k_estrela[j-1]) * h_k_hermitiano[j])**2))-b))))
            Gamma_k_j = np.append(Gamma_k_j, Gamma_k_j_valor)
            #Phi_NEIG_j
            Phi_NEIG_j_valor = ((Gamma_k_j[j] - (mu*Omega)) / (1-Omega)) * tau_k_estrela[j]
            Phi_NEIG_j = np.append(Phi_NEIG_j, Phi_NEIG_j_valor)
            #tau_k_estrela
            tau_k_valor = (E_min - ((np.sum(Phi_NEIG_j))*(1-Omega))) / (Gamma[k] - (mu*Omega))
            tau_k_estrela = np.append(tau_k_estrela, tau_k_valor)



    # Cálculo Phi_k
    for k in range (0, K):    
        Phi_k_valor = tau_k_estrela[k] * ((Gamma[k] - (mu*Omega)) / (1 - Omega))
        Phi_k_vetor = np.append(Phi_k_vetor, Phi_k_valor)

        # Cálculo de tau_k_total


    # Cálculo Phi_linha_k
    for j in range(0, K-1):
        Phi_linha_k_valor = Phi_NEIG_j[j] + Phi_k_vetor[k]
        Phi_linha_k = np.append(Phi_linha_k, Phi_linha_k_valor)    

    # tau_t total
    tau_total = np.sum(np.abs(tau_k_estrela))


    tau_total_vetor = np.append(tau_total_vetor, tau_total)


print(tau_total)






    
# tau_total_valor = sum(tau_k_estrela_maior_que_1)
    

            
            
        
print(len(Gamma_k_j))
print(len(Phi_NEIG_j))


print(len(Psi_k_estrela))
print(Psi_k_estrela)
print(len(Beta))
print(len(Psi_k_estrela))
print(len(h_k_hermitiano))
print(len(P))

# P_k


# %% Plot

num_dispositivos = np.arange(0,len(tau_total_vetor),1)
# print(len(num_dispositivos))
# print(len(tau_total_vetor))


plt.semilogy(num_dispositivos, tau_total_vetor, 'o-')

plt.xlabel('Números de dispostivos (K)')
plt.ylabel('tau_T')
plt.grid(True)
plt.legend()
plt.show()

# %%
