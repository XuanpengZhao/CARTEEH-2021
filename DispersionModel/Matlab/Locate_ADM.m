function jb=Locate_ADM(v,vs)

    N=length(v);
    
    if vs<=v(1); jb=1; return; end
    
    if vs>=v(N); jb=N-1; return; end
    
    jb=1; jt=N;
    
    jm=floor((jb+jt)/2);
    
    while (jt-jb)>1
    
        mid=v(jm);
    
        if vs>=mid
        
            jb=jm;
        
        else
        
            jt=jm;
        
        end
        
        jm=floor((jb+jt)/2);
        
    end
    
 %  End of program
        
        