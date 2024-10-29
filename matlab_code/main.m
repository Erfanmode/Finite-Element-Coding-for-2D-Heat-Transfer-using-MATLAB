clear
clc
q_dot1=0;
q_dot2=1e5;
T_inf=20; %Ambient temperature
T1=0; % initial Temperature
T2=0; % initial Temperature
sim_time=10000;%maximum simulation run time in second
dt=0.02;% in second
H=40;  % convection Coefficient
h1= 0.10; % height in meter
w1= 0.10; % width in meter
xc=0.10;  % center of inside rectangle
yc=0.10;  % center of inside rectangle
h2=0.20; % height in meter
w2=0.20 ;% width in meter
dx=0.001; % in meter
dy=0.001;% in meter
p1=2000; %density in SI
c1=400 ; %Thermal Capacity in SI
p2=3000; %density in SI
c2=500 ; %Thermal Capacity in SI
k1=1;
k2=2;
T_center=[];
T_contour=[];
k=init_k(k1,k2,dx,dy,h1,w1,h2,w2,xc,yc);
T=init_T(T1,T2,dx,dy,h1,w1,h2,w2,xc,yc);
pc=init_pc(p1,c1,p2,c2,dx,dy,h1,w1,h2,w2,xc,yc);
q_dot=init_q(q_dot1,q_dot2,T2,dx,dy,h1,w1,h2,w2,xc,yc);
N1x=w1/dx; %it should be even
N1y=h1/dy; %it should be even
N2x=w2/dx;
N2y=h2/dy;
Ncx=xc/dx; %it should be even
Ncy=yc/dy; %it should be even
Nt=sim_time/dt;
T_next=T;
count=0;
for t=1:Nt
    for i=1:N2x
        for j=1:N2y
            
            if(i==N2x)
                q1=dy*H*(T_inf - T(i,j));
            else
                Tm= ( k(i,j)*T(i,j)+k(i+1,j)*T(i+1,j) )/(k(i,j)+k(i+1,j));
                q1=dy*k(i+1,j)*( T(i+1,j)-Tm )/(dx) + dy*k(i,j)*(- T(i,j)+Tm )/(dx);
            end
           
            
            if(j==N2y)
                q2=dx*H*(T_inf - T(i,j));
            else
                Tm= ( k(i,j)*T(i,j)+k(i,j+1)*T(i,j+1) )/(k(i,j)+k(i,j+1));
                q2=dx*k(i,j+1)*( T(i,j+1)-Tm )/(dy) + dx*k(i,j)*(- T(i,j)+Tm )/(dy);
            end
            
            
            if(i==1)
                q3=dy*H*(T_inf - T(i,j));
            else
                Tm= ( k(i,j)*T(i,j)+k(i-1,j)*T(i-1,j) )/(k(i,j)+k(i-1,j));
                q3=dy*k(i-1,j)*( T(i-1,j)-Tm )/(dx) + dy*k(i,j)*( -T(i,j)+Tm )/(dx);
            end
            
            
            if(j==1)
                q4=dx*H*(T_inf - T(i,j));
            else
                Tm= ( k(i,j)*T(i,j)+k(i,j-1)*T(i,j-1) )/(k(i,j)+k(i,j-1));
                q4=dx*k(i,j-1)*( T(i,j-1)-Tm )/(dy) + dx*k(i,j)*(- T(i,j)+Tm )/(dy);
            end
           
            
            T_next(i,j)= ((q1+q2+q3+q4)+ q_dot(i,j)*dx*dy)*dt / (pc(i,j)*dx*dy) + T(i,j);
            
        end
    end
    clc
    fprintf("%.2f percent of simulation is done",t/Nt*100);
    T=T_next;
    T_center(t)=T(Ncx,Ncy);
    if( t==round(count*(Nt/7))+1 )
        count=count+1;
        T_contour(:,:,count)=T;
    end

end
%contour((1:N2x)*dx,(1:N2y)*dy,T_contour(:,:,5));
mesh((1:N2x)*dx,(1:N2y)*dy,T);
figure
plot(dt*(1:Nt),T_center);
xlabel("time_{s}");
ylabel("TEmperature of center_{c}");


