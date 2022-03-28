num_nu = 3  # quantos neutrinos, 2 ou 3
t_ini = 0.0
t_fin = 1.0

# dados do NuFit
dm2_21 = 34.0
dm2_32 = 170.5
theta12 = 33.44
theta23 = 49.2  # todos os angulos em graus
theta13 = 8.57
deltacp = 195.0

# parametros da EDO
passo = 1e-7
eps_abs = 1e-8
eps_rel = 1e-8

# condicao inicial
re1 = 1.0
re2 = 0.0
re3 = 0.0
im1 = 0.0
im2 = 0.0
im3 = 0.0

# constantes
G_F = 3.0
N_A = 6.02214076 * 10e23
energ = 1.0
metodo = 'DOP853'   # qual metodo utilizar

# energia em MeV, distancia em R_solar, velocidade em c=1
