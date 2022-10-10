function V=Run_Dispersion_Model4(V,X,Y,Z,Emis,u,v,Kz,Kx,Ky,delt)
%-----------------------------------------------------
% This program updates concentrations over a time delt
% Author: Akula Venkatram  Date: 7/9/2022
%-----------------------------------------------------
% V(NX-1,NY-1,NZ-1)=Three dimensional array of concentrations at t
% Vn(NX-1,NY,NZ)=Three dimensional array of concentrations at t+delt
% delt=time step
% Emis(NX-,NY-1)=Array of surface emission rates
% u(NZ)=Array of velocities along the z-axis
% v(NZ)=Array of velocities along the z-axis
% sigv=standard deviation of horizontal velocity fluctuations
% sigw=standard deviation of vertical velocity fluctuations
% X=Vector of uniformly spaced grid points along x-axis
% Y=Vector of uniformly spaced grid points along y-axis
% Z=Vector of grid points along the z-axis
%----------------------------------------------

    NX=length(X); NY=length(Y); NZ=length(Z);
    
    %------------------------------------------------------
    
    % Compute change associated with emissions
    
    delz=Z(2)-Z(1); delx=X(2)-X(1); dely=Y(2)-Y(1);
 
    for ix=1:NX

        for iy=1:NY

            V(ix,iy,1)=V(ix,iy,1)+...
                Emis(ix,iy)*delt/(delx*dely*delz); 

        end

    end
    
    %---------------------------------------------------------------------
    
   for iz=1:NZ
    
        V(:,:,iz)=Advect(u(iz),v(iz),delt,V(:,:,iz),X,Y,NX,NY);
        
    end
 %----------------------------------------------------------------
    
   % Change associated with diffusion in the z-direction
    
    for ix=1:NX

        for iy=1:NY

            V(ix,iy,:)=Diffuse(V(ix,iy,:),Z,Kz,delt);

        end

    end
%     
% %--------------------------------------------------------------
    
%    Change associated with diffusion in the x-direction

     for iz=1:NZ
     
        for iy=1:NY
         
            V(:,iy,iz)=Diffuse(V(:,iy,iz),X,Kx,delt);
         
        end
     
%     % ------------------------------------------------------
%      
     % Change associated with diffusion in the y-direction
         
            for ix=1:NX

                V(ix,:,iz)=Diffuse(V(ix,:,iz),Y,Ky,delt);

            end
            
    end
    
    %---------------------------------------------------- 
    
             
           

    