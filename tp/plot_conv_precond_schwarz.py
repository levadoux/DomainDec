import numpy as np
import matplotlib.pyplot as plt

# Liste des valeurs de h
h_values = [0.025, 0.01, 0.005]

plt.figure(figsize=(7,5))

# Boucle sur chaque fichier
for h in h_values:
    filename = f"dat/errH1precond_{h}.dat"
    data = np.loadtxt(filename, delimiter=';')
    niter = data[:, 0]
    errH1 = data[:, 1]
    plt.semilogy(niter, errH1, '--', label=f"h = {h}")

plt.xlabel("Nombre d'itérations")
plt.ylabel("Erreur relative H1")
plt.title("Convergence de l'erreur H1 avec préconditionnement pour différentes valeurs de h")
plt.grid(True, which='both', ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()


output_file = "errH1precond_schwarz.png"
plt.savefig(output_file, dpi=300)
print(f"Figure enregistrée sous : {output_file}")