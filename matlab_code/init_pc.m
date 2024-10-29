function pc_array=init_pc(p1,c1,p2,c2,dx,dy,h1,w1,h2,w2,xc,yc)
N1x=w1/dx;
N1y=h1/dy;
N2x=w2/dx;
N2y=h2/dy;
Ncx=xc/dx;
Ncy=yc/dy;
pc_array=p2*c2*ones(N2x,N2y);
for i=0:N2x
    for j=0:N2y
        if ( (j> Ncy-N1y/2) && (j< Ncy+N1y/2) ) && ( (i> Ncx-N1x/2) && (i< Ncx+N1x/2) )
                pc_array(i,j)=p1*c1;  
        end
    end


end