% function Simulate_Vehicle_Dispersion(Xe,Ye,Erate,wspd,wdir,sigw,delt,Tsim)
function Simulate_Vehicle_Dispersion(delt,Tsim,Kxy)
%-----------------------------------------------------
% This program updates concentrations over a time delt
% Author: Akula Venkatram  Date: 6/14/2022
%-----------------------------------------------------
% Xe,Ye=Position of emitting vehicle
% Erate= Emission rate of vehicle
% delt=time step for updating concentration
% Tsim=Total time for simulation
% wspd= Wind speed during delta
% wdir= Wind direction using meteorological convention
% Emis(NX,NY)=Array of emission rates
% u(NZ)=Vector of velocities along the z-axis
% v(NZ)=Vector of velocities along the z-axis
% sigw=standard deviation of vertical velocity fluctuations
% X=Vector of uniformly spaced grid points along x-axis
% Y=Vector of uniformly spaced grid points along y-axis
% Z=Vector of grid points along the z-axis
% hmix=Height of first Z level
% V(NX,NY,NZ)=Two dimensional array of concentrations at t
% Vn(NX,NY,NZ)=Two dimensional array of concentrations at t+delt
%------------------------------------------------  
    
%     G=readtable('Grid_Information.xlsx');

%     [X,Y,Z,Emis,V]=Initialize_Inputs(G.Xmax,G.Ymax,G.hmix,G.NX,G.NY,G.NZ,1.5); 
%     sigw = 5;
%     wspd = 2;
%     wdir = 225;
%     Erate = 20;
%     Xe = 500;
%     Ye = 500;

% 
    [X,Y,Z,Emis,V]=Initialize_Grid; 
    
    [Kz,u,v]=Read_Met(Z);  NX=length(X); NY=length(Y);
%     
    Kx(1:NX)=Kxy; Ky(1:NY)=Kxy;
%     u = [14.1421, 15.9714, 18.0372, 20.3703, 23.0051];
%     v = u;
    for it=1:ceil(Tsim/delt)
        
        Emis=Read_Emissions(X,Y,Emis,'Emis_1.xlsx');
        V=Run_Dispersion_Model4(V,X,Y,Z,Emis,u,v,Kz,Kx,Ky,delt);
%         V=Run_Dispersion_Model...
%         (X,Y,Z,V,Emis,Xe,Ye,Erate,wspd,wdir,sigw,delt);
    
        %h(it)=figure;
        center = V(50,50,3)
        area = sum(sum(V(:,:,3)))
        %pcolor(V(:,:,3)'); shading interp
        
    end
    a=1
 