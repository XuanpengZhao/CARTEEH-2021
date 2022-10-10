function U=Similarity_Wind_Airport_ADM(z,Met)

    kappa=0.4;
    
    z0=Met.z0; zi=Met.zi; L=Met.Lmon; dh=5*z0; ust=Met.ustar; 
    
    z=min(z,0.1*zi);  z=max(z,1.1*(z0+dh));

    if L>0.0
        
       psi1=-17*(1-exp(-0.29*(z-dh)/L));
       
       psi2=-17*(1-exp(-0.29*z0/L));
       
   else
       
       x1=(1-16*(z-dh)/L).^0.25; x2=(1-16*z0/L)^0.25;
       
       psi1=2*log((1+x1)/2)+log((1+x1.*x1)/2)-2*atan(x1)+pi/2;
       
       psi2=2*log((1+x2)/2)+log((1+x2.*x2)/2)-2*atan(x2)+pi/2;
       
    end
  
     U=ust*(log((z-dh)/z0)-psi1+psi2)/kappa;
             