function V=Run_Dispersion_Model(X,Y,Z,V,Emis,Xe,Ye,Erate,wspd,wdir,sigw,delt)
%-----------------------------------------------------
% This program updates concentrations over a time delt
% Author: Akula Venkatram  Date: 6/14/2022
%-----------------------------------------------------
% V(NX,NY,NZ)=Two dimensional array of concentrations at t
% Vn(NX,NY,NZ)=Two dimensional array of concentrations at t+delt
% delt=time step
% Emis(NX,NY)=Array of emission rates
% u(NZ)=Vector of velocities along the z-axis
% v(NZ)=Vector of velocities along the z-axis
% sigw=standard deviation of vertical velocity fluctuations
% X=Vector of uniformly spaced grid points along x-axis
% Y=Vector of uniformly spaced grid points along y-axis
% Z=Vector of grid points along the z-axis
% hmix=Height of first Z level
%------------------------------------------------  

    Nw=0.3; NZ=length(Z);
    
    delx=X(2)-X(1); dely=Y(2)-Y(1);
    
    wdir=(270-wdir); 
    
    u(1)=wspd*cosd(wdir); v(1)=wspd*sind(wdir);
    
    Kx=wspd*delx; Ky=wspd*dely; % Tentative Kx Ky
    
    for iz=1:NZ
        
        fw=(Z(iz)/Z(1))^Nw;
        
        u(iz)=u(1)*fw; v(iz)=v(1)*fw;
        
        Kz(iz)=sigw*Z(iz);
        
    end
    
    NXe=Locate_ADM(X,Xe); NYe=Locate_ADM(Y,Ye);
    
    Emis(NXe,NYe,1)=Erate; 
    
    V=Run_Dispersion_Model4(V,X,Y,Z,Emis,u,v,Kz,Kx,Ky,delt);
    
    
    
    