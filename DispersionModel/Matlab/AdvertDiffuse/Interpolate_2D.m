function Vint=Interpolate_2D(X,Y,V,Xn,Yn,i,j)

    imax=length(X)-1; jmax=length(Y)-1;

    p=(Xn-X(i))/(X(i+1)-X(i));
    
    q=(Yn-Y(j))/(Y(j+1)-Y(j));
    
    ip1=min(imax,i+1); jp1=min(jmax,j+1); 
    
    Vb=p*V(ip1,j)+(1-p)*V(i,j);
    
    Vt=p*V(ip1,jp1)+(1-p)*V(i,jp1);
    
    Vint=q*Vt+(1-q)*Vb;

    
