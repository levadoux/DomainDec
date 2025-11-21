import numpy as np
import matplotlib.pyplot as plt

# Liste des valeurs de h
h_values = [0.025, 0.01, 0.005, 0.0025]

plt.figure(figsize=(7,5))

filename1 = "dat/iter.dat"
data1 = np.loadtxt(filename1, delimiter=';')
h = data1[:, 0]
niter1 = data1[:, 1]

filename2 = "dat/iter_precond.dat"
data2 = np.loadtxt(filename2, delimiter=';')
niter2 = data2[:, 1]

plt.plot(h, niter1, linestyle = '-', marker = 'o', label="Sans préconditionneur")
plt.plot(h, niter2, linestyle = '-', marker = 'o', label="Avec préconditionneur")

plt.xlabel("Pas du maillage h")
plt.ylabel("Nombre d'itérations pour une erreur H1 < 10e-6")
plt.title("Convergence du gradient conjugué selon le pas de maillage")
plt.grid(True, which='both', ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()


output_file = "iterations.png"
plt.savefig(output_file, dpi=300)
print(f"Figure enregistrée sous : {output_file}")