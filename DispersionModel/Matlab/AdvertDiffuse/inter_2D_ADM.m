function Vint=inter_2D_ADM(V,x,z,jx,jz,xd,Zr)

        p=max(0,(xd-x(jx))/(x(jx+1)-x(jx)));
        
%         p=(xd-x(jx))/(x(jx+1)-x(jx));
        
        q=max(0,(Zr-z(jz))/(z(jz+1)-z(jz)));
        
%         q=(Zr-z(jz))/(z(jz+1)-z(jz));
        
    

        
        V1= V(jx,jz); V2=V(jx+1,jz);
        
        V3= V(jx+1,jz+1); V4=V(jx,jz+1);
        
        Vint=(V1*(1-p)+V2*p)*(1-q)+(V4*(1-p)+V3*p)*q;