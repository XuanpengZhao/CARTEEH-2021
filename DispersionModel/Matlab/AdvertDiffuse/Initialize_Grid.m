function [X,Y,Z,Emis,V]=Initialize_Grid

    G=readtable('Grid_Information.xlsx');
    
    Xmax=G.Xmax; Ymax=G.Ymax; ztop=G.ztop; hmix=G.hmix;
    
    NX=G.NX; NY=G.NY;  NZ=G.NZ; 

    X=linspace(0,Xmax,NX);  Y=linspace(0,Ymax,NY); 
    
    %%%%%%%%%%%%%%%%%%%%% XZ
    NZ = 6;
    Z(1) = hmix;
    for i = 2:NZ
        Z(i) = Z(i-1)*4;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Z=logspace(log10(hmix),log10(ztop),NZ-1);
    %logspace(log10(0.3),log10(76.8),5)
    %logspace(log10(0.3),log10(19.2),7)
    %logspace(log10(0.3),log10(153.6),10)
    Z=[0,Z]; NZ=length(Z);
    
    Emis(NX,NY)=0.0;  
    
    V(NX,NY,NZ)=0.0;
    
    
    
    
    
    