clear all
close all
clc
format long

%----------------Constants and Initial Conditions--------------------------

Dx=0.001; %Length of the Cells
Dt=1*10^(-7); %Timestep 
N=1000; %Number of cells = 2N 
T=3000; %Number of Timesteps 
R=287.05; %Specific Gas Constant for air

P_A=zeros(2*N,1); %Air Pressure of shock tube 
T_A=zeros(2*N,1); %Air Temperature of shock tube 
rho_A=zeros(2*N,1); %Air density of shock tube 
u_A=zeros(2*N,1); %Air velocity profile of shock tube 
c_A=zeros(2*N,1); %Air sound velocity profile of shock tube

Ta=20; %Atmospheric Temperature in Celcius 
Pa=101325; %Atmospheric Pressure in Pascals 
rho_a=(Pa)/(R*(Ta+273.15));%density of air (Ideal Gas Law)

u_a=0; %initial shock tube air velocity 
c_a=331.3*sqrt(1+Ta/273.15); %Air sound speed 
e_a=5/2*R*(Ta+273.15); %Air Internal Energy

PD=50000; % Pressure difference in the tube

%-----------Setting the initial condition in the tube----------------

for i=1:N
    P_A(i)=Pa;
    T_A(i)=Ta;
    rho_A(i)=rho_a;
    u_A(i)=u_a;
    c_A(i)=c_a;
    P_A(N+i)=Pa-PD;
    T_A(N+i)=Ta;
    rho_A(N+i)=P_A(N+i)/(R*(T_A(N+1)+273.15));
    u_A(N+i)=u_a;
    c_A(N+i)=c_a;
end

e_A=5/2.*R.*(T_A+273.15);

%----------------Creating initial set of conserved variables---------------
for j=1:2*N
    Q_A(1,j)=rho_A(j);
    Q_A(2,j)=rho_A(j)*u_A(j);
    Q_A(3,j)=rho_A(j)*(e_A(j)+0.5*u_A(j)^2);
end
    Q_OLD=Q_A; AIRPRESSURE=P_A;

%--------------------HLL RIEMANN SOLVER----------------------------------
    for timestep=1:T
        time=T-timestep;
        Q_NEW(:,1)=Q_OLD(:,1);
        Q_NEW(:,2*N)=Q_OLD(:,2*N);
        
        for i=2:2*N-1
            %--------------Calculating F(i+1/2)--------------------------
            SL=min([u_A(i)-c_A(i),u_A(i+1)-c_A(i+1)]);
            SR=max([u_A(i)+c_A(i),u_A(i+1)+c_A(i+1)]);
            RHO_L=Q_A(1,i);
            U_L=Q_A(2,i)/RHO_L;
            E_L=Q_A(3,i)/RHO_L;
            P_L=AIRPRESSURE(i);
            RHO_R=Q_A(1,i+1);
            U_R=Q_A(2,i+1)/RHO_R;
            E_R=Q_A(3,i+1)/RHO_R;
            P_R=AIRPRESSURE(i+1);

            if SL>=0
                F_A_PLUS=[RHO_L*U_L;RHO_L*U_L^2+P_L;U_L*RHO_L*E_L+U_L*P_L];
                marker=1;
            elseif SR<=0
                F_A_PLUS=[RHO_R*U_R;RHO_R*U_R^2+P_R;U_R*RHO_R*E_R+U_R*P_R];
                marker=2;
            else
                F_A_PLUS=(SR.*[RHO_L.*U_L;RHO_L.*U_L.^2+P_L;RHO_L.*E_L.*U_L+U_L.* P_L]-SL.*[RHO_R.*U_R;RHO_R.*U_R.^2+P_R;RHO_R.*E_R.*U_R+U_R.* P_R]+SL.*SR.*([RHO_R;RHO_R.*U_R;RHO_R.*E_R]-[RHO_L;RHO_L.*U_L;RHO_L.*E_L]))./(SR-SL);
            end

        F(:,i)=F_A_PLUS;
        end %end cell iteration

    F(:,1)=F(:,2);
    F(:,2*N)=F(:,2*N-1);

    %----------------Updating Vector of Conserved Variables---------------
    for i=2:2*N-1
        Q_NEW(:,i)=Q_OLD(:,i)+(Dt/Dx).*(F(:,i-1)-F(:,i));
    end

    %------------------Updating Primitive variables-----------------------
    rho_A=Q_NEW(1,:);
    u_A=Q_NEW(2,:)./rho_A;
    e_A=Q_NEW(3,:)./rho_A-0.5.*u_A.^2;
    T_A=(2/5).*e_A./R-273.15;
    c_a=331.3.*sqrt(1+Ta./273.15);
    AIRPRESSURE=(2/5).*rho_A.*e_A;

    %---------------------------------------------------------------------
    clear Q_OLD
    Q_OLD=Q_NEW;
    clear Q_NEW
end %end timestep
    plot(1:2000,rho_A)