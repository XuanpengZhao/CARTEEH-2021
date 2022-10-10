function x=Tridiag_Solver_ADM(e,f,g,h)
%----------------------------------------------
% This program solves a tridiagonal matrix
% Author: Akula Venkatram  Date: 6/14/2022
%----------------------------------------------

    N=length(f);
    
    for k=2:N                       
        
        mult=e(k)/f(k-1);
        
        f(k)=f(k)-mult*g(k-1);
        
        h(k)=h(k)-mult*h(k-1);
        
    end
    
    x(N)=h(N)/f(N);
    
    for k=N-1:-1:1
        
        x(k)=(h(k)-g(k)*x(k+1))/f(k);
        
    end
        
  % End of program