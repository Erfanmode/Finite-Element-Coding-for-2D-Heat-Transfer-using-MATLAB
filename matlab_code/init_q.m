function q_array=init_q(q_dot1,q_dot2,T2,dx,dy,h1,w1,h2,w2,xc,yc)
N1x=w1/dx;
N1y=h1/dy;
N2x=w2/dx;
N2y=h2/dy;
Ncx=xc/dx;
Ncy=yc/dy;
q_array=q_dot2*ones(N2x,N2y);
for i=0:N2x
    for j=0:N2y
        if ( (j> Ncy-N1y/2) && (j< Ncy+N1y/2) ) && ( (i> Ncx-N1x/2) && (i< Ncx+N1x/2) )
                q_array(i,j)=q_dot1;  
        end
    end

end