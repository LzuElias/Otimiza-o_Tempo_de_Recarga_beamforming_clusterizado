
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


K = 10                                  # N° dispositivos
N = 4                                   # N° de antenas
R = 10                                  # Raio [m]
f = 915 *(10**6)                        # Freq. da portadora
kappa = 1.5                             # Fator Rice
mu = 10.73 * 10**(-3)                   # Max pot no IoT
a = 0.2308                              # constante circuito a
b = 5.365                               # constante circuito b  
Pt = 3                                  # Max pot. de transmissão
E_min = 200 * (10**-6)                  # Min Energia no IoT [J]
alpha = 2.7                             # Coef. perda de percurso
c = 3*(10**(8))                         # Vel. da luz
Omega = 1/(1+(np.exp(a*b)))             # Constante resposta 0_in-0_out


#%% 
seed = np.random.seed(9)

# Parâmetros
canal = np.zeros((K,N)).astype(complex)
h_LoS = np.zeros((K,N)).astype(complex)
h_NLoS = np.zeros((K,N)).astype(complex)
Beta = np.zeros(K)
Pr = np.zeros(K)
Beamform_SCSI = np.zeros((K,N)).astype(complex)
Gamma = np.zeros(K)
tempo_carregamento = np.zeros(K)
Gamma_kj = np.zeros(K)
Phi_NEIG = np.zeros(K)
tempo_carregamento_total = np.zeros(K)
beam_hH_aux = np.zeros(K).astype(complex)
beam_hH_aux_2 = np.zeros([K])
beam_hH_aux_3 = np.zeros([K])


# %%

# Canal
h_LoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
h_NLoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
canal = (np.sqrt((kappa)/(1+kappa)) * h_LoS) + (np.sqrt((1)/(1+kappa)) * h_NLoS)
hH = np.conjugate(np.transpose(canal)) # Hermitiano (NxK)

for k in range (0,K):
    # Beta:
    x = np.random.randint(1,R)
    y = np.random.randint(1,R)
    d = np.sqrt(x**2 + y**2)
    Beta[k] = (c**2) / ((16*((np.pi)**2)) * (f**2) * (d**alpha))

for k in range (0,K): 
    # Beamforming S-CSI (h_LoS)
    Beamform_SCSI[k][:] = (np.sqrt(Pt/N) * (h_LoS[k][:] / np.abs(h_LoS[k][:])))

Pr = Beta @ np.abs((Beamform_SCSI @ hH))

for k in range (0,K): 
    # Função Logística Tradicional (Gamma):
    Gamma[k] = mu / (1 + (np.exp(-a*(Pr[k] - b))))

for k in range(0, K):
    if k == 0:
        beam_hH_aux[k] = np.dot(Beamform_SCSI[k][:], hH[:, k])
    else:
        for j in range(0, k):
            beam_hH_aux[k] += np.dot(Beamform_SCSI[j][:], hH[:, k]) 

for k in range (0,K):
    beam_hH_aux_2[k] =  Beta[k] * (((np.abs(beam_hH_aux[k]))**2) -b)

beam_hH_aux_3 = -a*beam_hH_aux_2


#%%

for k in range (0,K):     
    # Tempo de carregamento para o primeiro dispositivo:
    if k==0:
        tempo_carregamento[0] = (E_min * (1-Omega)) / (Gamma[0] - (mu*Omega))
    else:
        for j in range (0, k):
            # print(j)
            # Energia total coletada carregando o vizinho (Phi_NEIG):
            #Gamma_kj[k] = mu / (1 + (np.exp(-a*(Beta[k] * (((np.abs(np.dot(Beamform_SCSI[k][:], hH[:, k])))**2)-b)))))
            Gamma_kj[j] = mu / (1 + (np.exp(beam_hH_aux_3[k])))
            Phi_NEIG = tempo_carregamento[j] * ((Gamma_kj - (mu*Omega)) / (1 - Omega))
            if np.sum(Phi_NEIG) > E_min:
                tempo_carregamento[k] = 0
            else:
                tempo_carregamento[k] = ((E_min - np.sum(Phi_NEIG))*(1-Omega)) / (Gamma[k] - (mu*Omega))
    # print(f"k = {k}")

tempo_carregamento_total = np.sum(tempo_carregamento)

        
np.shape(Beamform_SCSI)
np.shape(Gamma_kj)
np.shape(Beta)


print("Terminou de rodar\n")