
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
L = 1                                   # N° de RIS
f = 915 *(10**6)                        # Freq. da portadora
kappa = 1.5                             # Fator Rice
mu = 10.73 * 10**(-3)                   # Max pot no IoT
a = 0.2308                              # constante circuito a
b = 5.365                               # constante circuito b  
Pt = 3                                  # Max pot. de transmissão
E_min = 1* (10**-6)                     # Min Energia no IoT [J]
alpha = 3.5                             # Coef. perda de percurso sem visada direta
alpha2 = 2                              # Coef. perda de percurso com visada direta
c = 3*(10**(8))                         # Vel. da luz
Omega = 1/(1+(np.exp(a*b)))             # Constante resposta 0_in-0_out
Pot_k = 0
# Parâmetros Ris
M = 5                                   # Linha e Coluna das células da Ris
a = 1                                   # Valor absoluto Ris


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

#  para acréscimo da RIS
theta = np.zeros((M, M)).astype(complex)                    # theta é a matriz da RISe Beta_RIS_Disp
h_LoS_PB_RIS = np.zeros((M,N)).astype(complex)
hH_Ris = np.zeros((M,N)).astype(complex)
Beta2 = np.zeros((L))

# Posições
RIS_position = np.array([R/2, R/2])
PB_position = np.array([0,0])

#%%

# Montagem da RIS
for linha in range (0,M):
    for coluna in range(0,M):
        ris_phase = np.random.random()
        theta[linha][coluna] = a + np.exp(-1j*(ris_phase*2*np.pi))
        

#%%

# Canal
# PB - Dispositivo
h_LoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
h_NLoS = np.random.rand(K,N) + 1j*np.random.rand(K,N) # Matriz (KxN)
hH = (np.sqrt((kappa)/(1+kappa)) * h_LoS) + (np.sqrt((1)/(1+kappa)) * h_NLoS) # Canal PB-Device

# PB-RIS
h_LoS_PB_RIS = np.random.rand(M,N) + 1j*np.random.rand(M,N) # Matriz (MxN)
hH_Ris = np.transpose(h_LoS_PB_RIS) * theta # erro aqui

np.shape(h_LoS_PB_RIS)
np.shape(theta)

# Localização dos dispositivos em torno da RIS
x_c = R/2
y_c = R/2
    # Número de pontos
n = 20
    # Desvio padrão (controla o quão dispersos os pontos estarão do centro)
std_dev = 2
    # Gerando coordenadas x e y aleatórias com distribuição normal em torno do ponto central
x_disp = np.random.normal(loc=x_c, scale=std_dev, size=n)
y_disp = np.random.normal(loc=y_c, scale=std_dev, size=n)
loc_dispositivos = np.array((x_disp, y_disp))

# Exibindo os pontos em um gráfico
plt.scatter(x_disp, y_disp)
plt.plot(x_c, y_c, 'ro')  # Ponto central marcado em vermelho
plt.xlabel('X_Dispositivos')
plt.ylabel('Y_Dispositivos')
plt.title('Localização dos dispositivos em torno da RIS')
plt.show()

#%%

# Beamforming S-CSI (h_LoS)
for k in range (0,K):
    Beamform_SCSI[k] = np.transpose(np.sqrt(Pt/N) * (h_LoS[k]/np.abs(h_LoS[k])))

#hH = np.conjugate(np.transpose(canal)) # Hermitiano (NxK)

#%%
for k in range (0,K):
    Phi_NEIG = np.zeros(K)
        
    # Beta:
    #x = np.random.randint(1,R)
    #y = np.random.randint(1,R)
    d = np.sqrt(x_disp[k]**2 + y_disp[k])
    # Beta PB_Disp
    Beta[k] = (c**2) / ((16*((np.pi)**2)) * (f**2) * (d**alpha))
    # Beta PB_RIS e Beta_RIS_Disp
    Beta2[k] = (((c/f)**2)/((4*(np.pi))**2)) * d**(-alpha2)
    
    
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
# %%
