#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh2D Omega;
  Read(Omega,"tp1-3_0.01.mesh");
  
  auto partition = Partition4(Omega); 
  std::vector<Mesh2D> Sigma = partition.first; 
  CooMatrix<double> Q = partition.second;
  Plot(Sigma,"partition/partition4.mesh"); 
}