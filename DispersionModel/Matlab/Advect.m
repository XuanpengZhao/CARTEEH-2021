function Vn=Advect(u,v,delta,V,X,Y,NI,NJ)

    for ix=1:NI
        
        for iy=1:NJ
        
            Xn=X(ix)-u*delta;
      
            in=Locate_ADM(X,Xn);
            Yn=Y(iy)-v*delta;
          
            jn=Locate_ADM(Y,Yn);
            Vn(ix,iy)=Interpolate_2D(X,Y,V,Xn,Yn,in,jn);%inter_2D_ADM(V,X,Y,in,jn,Xn,Yn);
            
          
        end
        
    end
        
        
            
            
            
            

    