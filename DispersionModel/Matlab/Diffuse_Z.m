function Vn=Diffuse_Z(V,Z,K,delt)

    NZ=length(Z); 

    for iz=2:NZ-1
        
        dzp=Z(iz+1)-Z(iz); dzm=Z(iz)-Z(iz-1); 
        
        dz=(dzp+dzm)/2;
        
        Kzp=(K(iz+1)+K(iz))/2;  Kzm=(K(iz)+K(iz-1))/2;
        
        alphap=Kzp*delt/(dzp*dzm); alpham=Kzm*delt/(dz*dzm); 
        
        e(iz)=-alpham; g(iz)=-alphap;
        
        f(iz)=(1+alpham+alphap);
        
        h(iz)=V(iz);
        
    end
    
    e(1)=0.0; f(1)=1.0; g(1)=0; h(1)=V(1);
    
    e(NZ)=1.0; f(NZ)=-1; g(NZ)=0; h(NZ)=0;
    
    Vn=Tridiag_Solver_ADM(e,f,g,h);
        
    
   
    