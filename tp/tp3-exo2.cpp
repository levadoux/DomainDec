#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh2D Omega;
  Read(Omega,"tp1-3_0.01.mesh");
  auto partition4 = Partition4(Omega,10); 
  std::vector<Mesh2D> Gamma4 = partition4.first; 
  CooMatrix<double> R4 = partition4.second;
  Write(Gamma4[0],"partition/part4_q0");
  Write(Gamma4[1],"partition/part4_q1");
  Write(Gamma4[2],"partition/part4_q2");
  Write(Gamma4[3],"partition/part4_q3");
}