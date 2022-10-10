function Vn=Spread_Emissions(V,X,Y,Z,Emis,u,v,Kx,Ky,Kz,delt)
%-----------------------------------------------------
% This program updates concentrations over a time delt
% Author: Akula Venkatram  Date: 6/14/2022
%-----------------------------------------------------
% V(NX,NY,NZ)=Two dimensional array of concentrations at t
% Vn(NX,NY,NZ)=Two dimensional array of concentrations at t+delt
% delt=time step
% Emis(NX,NY)=Array of emission rates
% u(NZ)=Array of velocities along the z-axis
% v(NZ)=Array of velocities along the z-axis
% sigv=standard deviation of horizontal velocity fluctuations
% sigw=standard deviation of vertical velocity fluctuations
% X=Vector of uniformly spaced grid points along x-axis
% Y=Vector of uniformly spaced grid points along y-axis
% Z=Vector of grid points along the z-axis
%----------------------------------------------

    NX=length(X); NY=length(Y); NZ=length(Z);
    
    delx=X(2)-X(1); dely=Y(2)-Y(1); 
    
    %------------------------------------------------------
    
    % Compute change associated with emissions
    
    for iz=1:NZ
    
        for ix=1:NX

            for iy=1:NY
                % Why times Z(1) ?
                V(ix,iy,1)=V(ix,iy,1)+Emis(ix,iy)*delt/(delx*dely*Z(1)); 

            end

        end
    
    end
    
    %-----------------------------------------------------
    
     % Change associated with transport by mean wind (u,v)
    
    for iz=1:NZ
    
        V(:,:,iz)=Advect(u(iz),v(iz),delt,V(:,:,iz),X,Y,NX,NY);
       
    %------------------------------------------------------
    
    % Change associated with diffusion in the x-direction
     
        for iy=1:NY
         
            V(:,iy,iz)=Diffuse(V(:,iy,iz),X,Kx,delt);
         
        end
     
     %------------------------------------------------------
     
     % Change associated with diffusion in the y-direction
         
            for ix=1:NX

                V(ix,:,iz)=Diffuse(V(ix,:,iz),Y,Ky,delt);

            end
            
    end
    
    %------------------------------------------------------
    
    % Change associated with diffusion in the z-direction
      
    for ix=1:NX

        for iy=1:NY
             
            V(ix,iy,:)=Diffuse_Z(V(ix,iy,:),Z,Kz,delt);

        end

    end
     
     Vn=V;
    
             
           

    