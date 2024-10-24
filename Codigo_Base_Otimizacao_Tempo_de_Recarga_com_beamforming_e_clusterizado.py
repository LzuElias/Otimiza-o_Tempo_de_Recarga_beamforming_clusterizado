
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
E_min = 1* (10**-6)                  # Min Energia no IoT [J]
alpha = 2.7                             # Coef. perda de percurso
c = 3*(10**(8))                         # Vel. da luz
Omega = 1/(1+(np.exp(a*b)))             # Constante resposta 0_in-0_out
Pot_k = 0


#%% Rodando
seed = np.random.seed(9)

# Parâmetros
hH = np.zeros((K,N)).astype(complex)
h_LoS = np.zeros((K,N)).astype(complex)
h_NLoS = np.zeros((K,N)).astype(complex)
Beta = np.zeros(K)
Pr = np.zeros(K)
Beamform_SCSI = np.zeros((K,N)).astype(complex)
Gamma = np.zeros(K)
tempo_carregamento = np.zeros(K)
Gamma_kj = np.zeros(K)
tempo_carregamento_total = np.zeros(K)


# Canal
h_LoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
h_NLoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
hH = (np.sqrt((kappa)/(1+kappa)) * h_LoS) + (np.sqrt((1)/(1+kappa)) * h_NLoS)

# Beamforming S-CSI (h_LoS)
for k in range (0,K):
    Beamform_SCSI[k] = np.transpose(np.sqrt(Pt/N) * (h_LoS[k]/np.abs(h_LoS[k])))

#hH = np.conjugate(np.transpose(canal)) # Hermitiano (NxK)

for k in range (0,K):
    Phi_NEIG = np.zeros(K)
    
    # Beta:
    x = np.random.randint(1,R)
    y = np.random.randint(1,R)
    d = np.sqrt(x**2 + y**2)
    
    Beta[k] = (c**2) / ((16*((np.pi)**2)) * (f**2) * (d**alpha))
        

    # Potência
    Pr[k] = Beta[k]*(np.abs((Beamform_SCSI[k].dot(hH[k]))))**2
    
    # Função Logística Tradicional (Gamma):
    Gamma[k] = mu / (1 + (np.exp(-a*(Pr[k] - b))))
    
    # Tempo de carregamento para o primeiro dispositivo:
    if k==0:
        tempo_carregamento[0] = (E_min * (1-Omega)) / (Gamma[0] - (mu*Omega))
    else:
        for j in range (0, k):
            print(j)
            # Energia total coletada carregando o vizinho (Phi_NEIG):
            Gamma_kj[j] = mu / (1 + (np.exp(-a*((Beta[k]*(np.abs(Beamform_SCSI[j].dot(hH[k]))**2))-b))))
            Phi_NEIG[j] = tempo_carregamento[j] * ((Gamma_kj[j] - (mu*Omega)) / (1 - Omega))
            
        if np.sum(Phi_NEIG) > E_min:
            tempo_carregamento[k] = 0
        else:
            tempo_carregamento[k] = ((E_min - np.sum(Phi_NEIG))*(1-Omega)) / (Gamma[k] - (mu*Omega))
                
                
    
tempo_carregamento_total = np.sum(tempo_carregamento)
                




        
np.shape(Beamform_SCSI)
np.shape(Gamma_kj)
np.shape(Beta)


print("Terminou de rodar\n")