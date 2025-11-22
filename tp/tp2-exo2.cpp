#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh2D Omega;
  Read(Omega,"tp1-3_0.01.mesh");
  
  auto partition4 = Partition4(Omega); 
  std::vector<Mesh2D> Sigma = partition4.first; 
  CooMatrix<double> Q = partition4.second;
  Plot(Sigma,"partition/partition4.mesh");
  
  auto partition4r = Partition4(Omega,2); 
  std::vector<Mesh2D> Gamma4 = partition4r.first; 
  CooMatrix<double> R4 = partition4r.second;
  Write(Gamma4[0],"partition/part4_q0");
  // Write(Gamma4[1],"partition/part4_q1");
  // Write(Gamma4[2],"partition/part4_q2");
  // Write(Gamma4[3],"partition/part4_q3");

  // auto partition16 = Partition16(Omega); 
  // std::vector<Mesh2D> Sigma16 = partition16.first; 
  // CooMatrix<double> Q16 = partition16.second;
  // Plot(Sigma16,"partition/partition16.mesh");
  
  auto partition16r = Partition16(Omega,2); 
  // std::vector<Mesh2D> Gamma16 = partition16r.first; 
  // CooMatrix<double> R16 = partition16r.second;
  // Write(Gamma16[0],"partition/part16_q0");
  // Write(Gamma16[1],"partition/part16_q1");
  // Write(Gamma16[2],"partition/part16_q2");
  // Write(Gamma16[3],"partition/part16_q3");
  // Write(Gamma16[4],"partition/part16_q4");
  // Write(Gamma16[5],"partition/part16_q5");
  // Write(Gamma16[6],"partition/part16_q6");
  // Write(Gamma16[7],"partition/part16_q7");
  // Write(Gamma16[8],"partition/part16_q8");
  // Write(Gamma16[9],"partition/part16_q9");
  // Write(Gamma16[10],"partition/part16_q10");
  // Write(Gamma16[11],"partition/part16_q11");
  // Write(Gamma16[12],"partition/part16_q12");
  // Write(Gamma16[13],"partition/part16_q13");
  // Write(Gamma16[14],"partition/part16_q14");
  // Write(Gamma16[15],"partition/part16_q15");
  return 0;
}