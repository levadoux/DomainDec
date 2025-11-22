#include <cmath>
#include <femtool.hpp>

int main(){
  Mesh2D Omega;
  Read(Omega,"tp1-3_0.01.mesh");
  auto Vh   = FeSpace(Omega); 
  double k = 7*M_PI; 
  auto F    = [&k](const R3& x){return std::cos(k*(x[0]+x[1]));}; 
  auto fh    = Vh(F);

  auto partition4 = Partition4(Omega,1); 
  std::vector<Mesh2D> Gamma4 = partition4.first; 
  CooMatrix<double> R4 = partition4.second;

  std::vector<std::size_t> tbl_0;
  for (const auto& [i,j,val]:R4){
    if (j==0 && val==1){
      tbl_0.push_back(i);
    }
  }
  auto [Uh_0,P_0] = Restrict(Vh,Gamma4[0],tbl_0);
  
  std::cout << "dim(Uh_0) = " << dim(Uh_0) << std::endl;
  auto fh_0 = P_0*fh;  
  //std::cout<<min(fh.begin(),fh.end())<<std::endl;
  Plot(Uh_0,fh_0,"output/output_tp3_exo1"); 

  
  auto pair_list = Partition4(Vh,1);
  auto pair_0 = pair_list[0];
  auto Uh_0_bis = pair_list[0].first;
  auto P_0_bis = pair_list[0].second;
  auto fh_h_bis = P_0_bis*fh;
  Plot(Uh_0_bis,fh_h_bis,"output/output_tp3_exo1_bis");

  return 0;
}