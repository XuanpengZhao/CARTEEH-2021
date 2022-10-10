function [Kz,u,v]=Read_Met(Z)
%     Z = [0.3, 0.6, 1.2];
    M=readtable('Met_Information.xlsx');
    
    Met.wspd=M.wspd; Met.wdir=M.wdir; Met.Lmon=M.Lmon;
    
    Met.zi=M.zi; Met.z0=M.z0; Met.zref=M.zref; Met.turbz=M.turbz;
    
    Met=Compute_Ustar_ADM(Met);
    
    wd=(270-Met.wdir)*pi/180; cost=cos(wd); sint=sin(wd);
     
    NZ=length(Z);
    
    for iz=1:NZ
        
        Kz(iz)=Businger_ADM(Z(iz),Met);
        
        U(iz)=Similarity_Wind_Airport_ADM(Z(iz),Met);
        
        u(iz)=U(iz)*cost; v(iz)=U(iz)*sint;
       
    end
 a = 1   
   