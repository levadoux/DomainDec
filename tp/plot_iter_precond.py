import numpy as np
import matplotlib.pyplot as plt

# Liste des valeurs de h
h_values = [0.025, 0.01, 0.005, 0.0025]

plt.figure(figsize=(7,5))

filename = "iter_precond.dat"
data = np.loadtxt(filename, delimiter=';')
h = data[:, 0]
niter = data[:, 1]
plt.plot(h, niter, linestyle = '-', marker = 'o', label="Nombre d'itérations")

plt.xlabel("Pas du maillage h")
plt.ylabel("Nombre d'itérations pour une erreur H1 < 10e-6")
plt.title("Convergence du gradient conjugué préconditionné selon le pas de maillage")
plt.grid(True, which='both', ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()


output_file = "iterations_precond.png"
plt.savefig(output_file, dpi=300)
print(f"Figure enregistrée sous : {output_file}")