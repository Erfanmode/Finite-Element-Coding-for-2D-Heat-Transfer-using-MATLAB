function T_array=init_T(T1,T2,dx,dy,h1,w1,h2,w2,xc,yc)
N1x=w1/dx;
N1y=h1/dy;
N2x=w2/dx;
N2y=h2/dy;
Ncx=xc/dx;
Ncy=yc/dy;
T_array=T2*ones(N2x,N2y);
for i=1:N2x
    for j=1:N2y
        if ( (j> Ncy-N1y/2) && (j< Ncy+N1y/2) ) && ( (i> Ncx-N1x/2) && (i< Ncx+N1x/2) )
                T_array(i,j)=T1;  
        end
    end

end