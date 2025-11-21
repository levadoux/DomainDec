#include <cmath>
#include <femtool.hpp>

int main(int argc, char** argv) {
  std::string h = argv[1];  
  Mesh2D Omega;
  Read(Omega,"tp1-3_"+h+".mesh");
  auto Vh   = FeSpace(Omega);
  // Assembly of finite element matrix of the operator
  // -Delta + 1 with Neumann BC
  auto A    = Stiffness(Vh)+Mass(Vh);
  double k = 10*M_PI; 
  auto F = [&k](const R3& x){return std::cos(k*(x[0]));};
  auto ue = Vh(F);  
  auto b    = A(ue);
  // Numerical solution obtained by conjugate gradient
  auto [uh, niter]   = cgsolve_tp(A,b,ue,"dat/errH1_"+h+".dat"); 
  Plot(Vh,uh,"output/output_tp1_exo3_"+h+".");
  std::cout<<"Nombre d'itérations :"<<std::endl;
  std::cout<<niter<<std::endl;

  return 0;

}