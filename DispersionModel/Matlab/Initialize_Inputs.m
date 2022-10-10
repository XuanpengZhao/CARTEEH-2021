function [X,Y,Z,Emis,V]=Initialize_Inputs(Xmax,Ymax,hmix,NX,NY,NZ,fz)

    X=linspace(0,Xmax,NX);  Y=linspace(0,Ymax,NY); 
    
    Z(1)=hmix; 
    
    for iz=2:NZ
        
        Z(iz)=Z(iz-1)*fz;
        
    end
    
    Emis(NX,NY,NZ)=0.0;  
    
    V(NX,NY,NZ)=0.0;
    
    
    
    
    
    