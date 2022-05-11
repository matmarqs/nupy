# simulacao
metodo = 'DOP853'   # qual metodo utilizar

num_nu = 3  # quantos neutrinos, 2 ou 3
t_fin = 1.0

# dados do NuFit
dm2_21 = 7.42e-5    # eV^2
dm2_32 = 2.514e-3   # eV^2
theta12 = 33.44
theta23 = 49.2  # todos os angulos em graus
theta13 = 8.57
deltacp = 194.0

# parametros da EDO
eps_abs = 1e-5
eps_rel = 1e-5

# condicao inicial
re1 = 1.0
re2 = 0.0
re3 = 0.0
im1 = 0.0
im2 = 0.0
im3 = 0.0

# constantes
#G_F = 1.1663787e-11    # em MeV
G_F = 1.1663787e-23     # eV^{-2}
N_A = 6.02214076e+23    # Avogadro

# energia em eV, distancia em R_solar, velocidade em c=1
