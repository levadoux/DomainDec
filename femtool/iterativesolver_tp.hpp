#ifndef ITERATIVE_SOLVER_TP_HPP
#define ITERATIVE_SOLVER_TP_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>

std::pair<std::vector<double>,int>
cgsolve_tp(const CooMatrix<double>&   A,
	const std::vector<double>& b, const std::vector<double>& ue, std::string filename) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = (1e-8)*r2;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  std::size_t niter = 0;    
  double errH1 = 1.0;
  std::ofstream datafile(filename);

  while (errH1> 1e-6 && niter++<2000){
   
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;

    auto err = ue-x;
    errH1 = std::real((A(err)|err)/(A(ue)|ue));
    errH1 = std::sqrt(errH1);
    datafile<<niter<<" ; "<<errH1<<std::endl;

    
      
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
    
      
  }
  datafile.close();
    
  return std::make_pair(x,niter);
    
}

std::pair<std::vector<double>,int>
PCGSolver_tp(const CooMatrix<double>&   A,
	const std::vector<double>& b, const CholeskyPrec& Q, const std::vector<double>& ue, std::string filename) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto z   = Q * r; // Ici on utilise le préconditionneur
  auto    p   = z;
  auto   Ap   = A*p; 
  double rz   = (r|z);  
  double eps2 = (1e-8)*rz;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  
  std::size_t niter = 0;    
  double errH1 = 1.0;
  std::ofstream datafile(filename);  

  while( errH1>1e-6 && niter++<2000 ){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = rz/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;  
    z     = Q*r;  
    rz    = (r|z);
    beta  = rz/(alpha*pAp);
    p     = beta*p+z;

    auto err = ue-x;
    errH1 = std::real((A(err)|err)/(A(ue)|ue));
    errH1 = std::sqrt(errH1);
    datafile<<niter<<" ; "<<errH1<<std::endl;
      
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
      
  }
  datafile.close();
    
  return std::make_pair(x,niter);
    
}



#endif
