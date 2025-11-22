#include <cmath>
#include <femtool.hpp>

int main(int argc, char** argv) {
  std::string h = argv[1];  
  Mesh2D Omega;
  Read(Omega,"tp1-3_"+h+".mesh");
  auto Vh   = FeSpace(Omega);
  auto A    = Stiffness(Vh)+Mass(Vh);
  
  // Partition en 4 sous-domaines avec nl=2 couches
  auto partition = Partition4(Vh, 1);
  
  // Extraire les matrices R_j (prolongement = transposée de la restriction)
  std::vector<CooMatrix<double>> R(4);
  for(std::size_t j = 0; j < 4; ++j){
    R[j] = partition[j].second.T();  // P^T pour avoir le prolongement
  }
  
  // Créer le préconditionneur de Schwarz additif
  auto Q = AdditiveSchwarz(A, R);
  double k = 10*M_PI; 
  auto F = [&k](const R3& x){return std::cos(k*(x[0]));};
  auto ue = Vh(F);  
  auto b    = A(ue);
  
  // Numerical solution obtained by conjugate gradient
  auto [uh, niter]   = PCGSolver_tp(A,b,Q,ue,"dat/errH1precond_AdditiveSchwarz"+h+".dat"); 
  Plot(Vh,uh,"output/output_tp3_exo3_"+h+".");
  std::cout<<"Nombre d'itérations :"<<std::endl;
  std::cout<<niter<<std::endl;

  return 0;

}