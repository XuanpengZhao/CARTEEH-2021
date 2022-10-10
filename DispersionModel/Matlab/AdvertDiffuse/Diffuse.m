function Vn=Diffuse(Vn,X,K,delt)

    NX=length(X); V(1:NX)=Vn;

    for ix=2:NX-1
        
        dxp=(X(ix+1)-X(ix)); dxm=(X(ix)-X(ix-1)); 
        
        Kxp=(K(1)+K(1))/2; Kxm=(K(1)+K(1))/2;
        
        dx=(dxm+dxp)/2;
        
        alphap=Kxp*delt/(dxp*dx); alpham=Kxm*delt/(dxm*dx); 
        
        e(ix)=-alpham; g(ix)=-alphap;
        
        f(ix)=(1+alpham+alphap);
        
        h(ix)=V(ix);
        
    end
    
    alpha1=(2*K(1))*delt/(X(2)-X(1))^2;
    
    e(1)=0.0; f(1)=(1+alpha1); g(1)=-alpha1; h(1)=V(1);
    
    alpha2=(K(1)*2)*delt/(X(NX)-X(NX-1))^2;
    
    e(NX)=-alpha2; f(NX)=(1+alpha2); g(NX)=0; h(NX)=V(NX);
    
    Vn=Tridiag_Solver_ADM(e,f,g,h);
        
    
   
    