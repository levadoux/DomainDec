#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh3D Omega;
  Read(Omega,"tp1-2.mesh");
  auto [Bd, index] = Boundary(Omega);
  auto VhB = FeSpace(Bd); 
  double k = 5*M_PI; 
  auto F    = [&k](const R3& x){return std::cos(k*(x[0]+x[1]+x[2]));}; 
  auto u    = VhB(F); 
  Plot(VhB,u,"output/output_tp1_exo2"); 

  return 0;
}