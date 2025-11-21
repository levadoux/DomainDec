//+ Utilisation du kernel OpenCASCADE
SetFactory("OpenCASCADE");

h = 0.0025;
//+ Définition des points (coins du carré)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

//+ Définition des lignes du contour
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+ Boucle fermée pour la surface
Curve Loop(1) = {1, 2, 3, 4};

//+ Définition de la surface plane
Plane Surface(1) = {1};

//+ Définition de la surface physique
Physical Surface("Omega") = {1};
