#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh2D Omega;
  Read(Omega,"tp1-1.mesh");
  auto Vh   = FeSpace(Omega); 
  double k = 10*M_PI; 
  auto F    = [&k](const R3& x){return std::cos(k*(x[0]+x[1]));}; 
  auto u    = Vh(F); 
  Plot(Vh,u,"output/output_tp1_exo1"); 

  return 0;
}