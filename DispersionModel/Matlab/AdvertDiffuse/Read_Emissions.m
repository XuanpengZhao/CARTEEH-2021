function Emis=Read_Emissions(X,Y,Emis,flname)

    Em=readmatrix(flname);
    
    Xv=Em(:,1); Yv=Em(:,2); Erate=Em(:,3);
    
    for iv=1:length(Xv)
        
        ix=Locate_ADM(X,Xv(iv)); iy=Locate_ADM(Y,Yv(iv));
        
        w1=1/((Xv-X(ix))^2+(Yv-Y(iy))^2); w2=1/((Xv-X(ix+1))^2+(Yv-Y(iy))^2);
        
        w3=1/((Xv-X(ix+1))^2+(Yv-Y(iy+1))^2); w4=1/((Xv-X(ix))^2+(Yv-Y(iy+1))^2);
        
        W=w1+w2+w3+w4;
        
        Emis(ix,iy)=Erate(iv)*w1/W; Emis(ix+1,iy)=Erate(iv)*w2/W;
        
        Emis(ix+1,iy+1)=Erate(iv)*w3/W; Emis(ix,iy+1)=Erate(iv)*w4/W;
        
    end