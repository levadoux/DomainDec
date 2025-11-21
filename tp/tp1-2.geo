//+ Utilisation du kernel OpenCASCADE
SetFactory("OpenCASCADE");

//+ Taille caractéristique de la maille
lc = 0.02;

//+ Définition du cube : x0,y0,z0, dx,dy,dz
Box(1) = {0, 0, 0, 1, 1, 1};

//+ Définition du volume physique
Physical Volume("Omega") = {1};

//+ Définir une maille sur toutes les faces et le volume
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Mesh 3;
