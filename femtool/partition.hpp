#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <vector>
#include <filesystem>
using Mesh2DPart = std::vector<Mesh2D>;


std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega){
   
    std::size_t ne = Omega.size();    
    auto Sigma = Mesh2DPart(4);    
    CooMatrix<double> Q(ne,4);    
    for (int p = 0; p < 4; p++){
        Sigma[p] = Mesh2D(Omega.nodes());
    } 
    
    double x;
    double y;
    int p;
    for (std::size_t e = 0; e < ne; e++){
        auto elt = Omega[e];
        R3 bary = Ctr(elt);
        x = bary[0];
        y = bary[1];

        if (x<0.5){
            if (y<0.5){
                p=0;
            }
            else{
                p=1;
            }
        }
        else {
            if (y<0.5){
                p=2;
            }
            else{
                p=3;
            }
        }
        Sigma[p].push_back(elt);
        Q.push_back(e,p,1);
    
    }

    return std::make_pair(Sigma,Q);
}


bool AreNeighbors(const Element<2>& e1, const Element<2>& e2){
    int common = 0;
    for (std::size_t i = 0; i < 3; ++i){
        for (std::size_t j = 0; j < i; ++j){
            if (Norm(e1[i]-e2[j]) < 10e-8){
                ++common;
            }
        }    
    }
    return (common ==2);
}



std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega, const std::size_t& nl){
    std::size_t ne = Omega.size();
    auto partition = Partition4(Omega);
    Mesh2DPart Gamma = partition.first;
    CooMatrix<double> R = partition.second;
    for (std::size_t n = 0; n<nl; n++){
        for (int p = 0; p<4; p++){
            auto m = Gamma[p];
            for(const auto& elt1:m){
                for (std::size_t e = 0; e < ne; e++){
                    auto elt2 = Omega[e]; 
                    if (AreNeighbors(elt1,elt2)){
                        Gamma[p].push_back(elt2);
                        R.push_back(e,p,1);
                    }
                }
            }
        }
    }
    return std::make_pair(Gamma,R);
}

void
Plot(const std::vector<Mesh2D>& Sigma, const std::string& filename){
    std::size_t DIM = 2;
    std::vector<std::string> tag = {"Vertices", "Edges", "Triangles", "Tetrahedra"};

    std::size_t ne = 0;
    for (const auto& m:Sigma){
        ne += m.size();
    }

    // Ouverture fichier
    std::ofstream f;
    f.open(filename.c_str());

    // Preambule
    f << "MeshVersionFormatted 3\n\n";
    f << "Dimension\n3\n\n";

    // Section noeuds
    auto v = Sigma[0].nodes();
    const auto& v0 = v[0];
    f << "Vertices\n";
    f << v.size() << "\n";
    for(const auto& x:v){
    f << x << "\t1\n";}
    f << "\n";

    f << tag[DIM]   << "\n"; // maillage 2D
    f << ne << "\n";

    // Section elements
    int p=0;
    for (const auto& m:Sigma){
        p +=1;
        for(const auto& e:m){
            for(std::size_t j=0; j<DIM+1; ++j){
                f << 1+int(&e[j]-&v0) << "\t";}
                f << p+1 <<"\n";}
    }
    
  

  // Fermeture
  f << "\nEnd";
  f.close();
  
}

#endif