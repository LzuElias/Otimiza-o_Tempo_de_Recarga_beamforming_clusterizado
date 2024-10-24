
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


K = 10                                   # N° dispositivos
N = 4                                   # N° de antenas
R = 50                                  # Raio [m]
f = 915 *(10**6)                        # Freq. da portadora
kappa = 1.5                             # Fator Rice
mu = 10.73 * 10**(-3)                   # Max pot no IoT
a = 0.2308                              # constante circuito a
b = 5.365                               # constante circuito b  
Pt = 2                                  # Max pot. de transmissão
E_min = 200 * (10**-6)                  # Min Energia no IoT [J]
alpha = 2.7                             # Coef. perda de percurso
c = 3*(10**(8))                         # Vel. da luz
Omega = 1/(1+(np.exp(a*b)))             # Constante resposta 0_in-0_out
Pot_k = 0


def dist_euclid_quad (x, y):
    d_k = np.sqrt(((x**2)+(y**2))) 
    return d_k

def Beta_k(d_k):
    Beta_k = (c**2)/(16 * (np.pi**2) * (f**2) * (d_k**alpha))
    return Beta_k

def h_barra_k(N): # Gerador componente LoS
    h_barra_k = np.random.rand(1,N) + 1j*np.random.rand(1,N)
    return h_barra_k

def h_til_k(N): # Gerador de componente NLoS
    h_til_k = (np.random.rand(1,N) + 1j*np.random.rand(1,N))*(1/(np.sqrt(2)))
    return h_til_k

def canal_h_k (kappa, h_barra_k, h_til_k): # Gerador de canal
    canal_h_k = (np.sqrt((kappa)/(1+kappa)) * h_barra_k) + (np.sqrt((1)/(1+kappa)) * h_til_k)
    return canal_h_k

# def P_k(Beta_k, Psi_k, canal_h_k_Hermitiano): # Potencia de RF recebida
#     P_k = Beta_k * (np.abs((Psi_k) *(canal_h_k_Hermitiano)))**2
#     return P_k

# def Gamma_k(P_k): # Parâmetro para calcular a Energia coletada
#     Gamma_k = mu / (1 + np.exp(-a*(P_k-b)))
#     return Gamma_k

# def Phi_k(tau_k, Gamma_k): # Energia coletada
#     Phi_k = tau_k * ((Gamma_k - (mu*Omega))/(1-Omega))
#     return Phi_k

# def tau_k_igual_1(Gamma_k):
#     tau_k_estrela_k_igual_1 = (E_min * (1 - Omega))/(Gamma_k - (mu*Omega))
#     return tau_k_estrela_k_igual_1



#%%
#Rodando
canal_h = np.zeros((K,N)).astype(complex)
h_barra = np.zeros((K,N)).astype(complex)
h_til = np.zeros((K,N)).astype(complex)
Beta = np.zeros(K)
Pot = np.zeros(K)
Psi_k_estrela = np.zeros((K,N)).astype(complex)
Gamma = np.array([])
Gamma_k_j = np.array([])
Phi_NEIG_j = np.array([])
tau_k_estrela = np.zeros([K])
tau = np.array([])
tau_total = np.array([])
# H = np.zeros((K,N))
Phi_NEIG_j_vetor = np.array([])
h = np.zeros((N,K)).astype(complex)
h_k = np.zeros((N,K)).astype(complex)
seed = np.random.seed(9)

# Canal h
h_barra = np.random.rand(K,N) + 1j*np.random.rand(K,N)
h_til = np.random.rand(K,N) + 1j*np.random.rand(K,N)
canal_h = (np.sqrt((kappa)/(1+kappa)) * h_barra) + (np.sqrt((1)/(1+kappa)) * h_til)
h_k = np.conjugate(np.transpose(canal_h)) # Hermitiano (4x1)
h = h_k

for k in range (0,K):
    # Beta
    x =  np.random.randint(1,10)
    y =  np.random.randint(1,10)
    Beta[k] = Beta_k(dist_euclid_quad(x,y)) # (Constante)

    # Beamforming S-CSI
    Psi_k_estrela[k][:] = (np.sqrt(Pt/N)) *(h_barra[:][k]/np.abs(h_barra[:][k])) #(1x4)


for k in range (0, K):
  
    # Potência   
    Pot[k] = (Beta[k]*(np.abs(Psi_k_estrela[k][:].dot(h[:][k])))**2) # array(10)

    # Gamma (função logística tradicional)
    # Gamma = Gamma_k(Pot) # array(10)  
    aux_1 = -a*(Pot - b)
    Gamma = mu/(1 + (np.exp(aux_1)))
    
    # Tau_k para k=0
    
    if k == 0:
        tau_0 = (E_min*(1-Omega)) / (Gamma[0] - (mu*Omega)) # Valor fixo
        # tau_k_estrela = np.append(tau_k_estrela, tau_0) # Vetor para calcular Phi_NEIG (array(10))
        tau_k_estrela[0] = tau_0
    else:
        for j in range (0,k):
            # print(j)
            produto_escalar = np.abs(Psi_k_estrela[j][:] @ h[:][k])**2 # Absolut (Psi * H_hermitiano) **2
            Gamma_k_j = mu / (1+(np.exp(-a*((Beta[k]*produto_escalar)-b)))) # Valor
            Phi_NEIG_j = (tau_k_estrela[j]) * ((Gamma_k_j - (mu*Omega))/(1-Omega))
            Phi_NEIG_j_vetor = np.append(Phi_NEIG_j_vetor, Phi_NEIG_j)
            Phi_NEIG_j_sum = np.sum(Phi_NEIG_j_vetor)
            if(Phi_NEIG_j_sum > E_min): 
                tau=0
            else:
                tau = ((E_min - Phi_NEIG_j_sum) * (1-Omega)) / (Gamma[k] - (mu*Omega))
                # if tau < 0: tau = 0
            tau_k_estrela[k] = tau
    # print(f"k = {k}")
tau_total = np.append(tau_total, sum(tau_k_estrela))




