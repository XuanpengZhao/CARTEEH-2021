function Kz=Businger_ADM(z,Met)

    kappa=0.35; 
    
    zi=Met.zi; turbz=Met.turbz; L=Met.Lmon; ustar=Met.ustar;
    
    z=min(z,0.1*zi); 
    
    if z<=turbz;z=turbz;end
    
    psi=z/L; 
    
    if psi>0.0
        
        phih=0.74*(1+6.3*psi);
        
    else
        
        phih=0.74/(1-9*psi)^0.5;
        
    end
    
    if z<=zi
    
        Kz=kappa*ustar*z/phih;
        
    else
        
        Kz=0.01;
        
    end
        
        
        

    
     
    
 
 %  End of program