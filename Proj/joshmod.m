clear;
clc;
close all;

%Model parameters
%------------------------------------------------------------------------%
n=10001;                  %Number of grid points
L=1;                    %Length of domain
h=L/(n-1);              %Spatial step size
CFL=0.3;               %CFL number for stability
t_final=500e-8;         %Final time
x=0:h:L;

%Material params
%------------------------------------------------------------------------%
G= 76e9;                %Shear Modulus
rho0= 2785;               %Density
Grun= 2.0;              %Gruneisen parameter
Sigma= 1.338;           %related to dK/dp
specheat= 897;           %J/kg K

%Initial conditions
%------------------------------------------------------------------------%
p_l=8e9;            
p_r=0;            
rho_l=4000;            
rho_r=rho0;        
u_l=0;              
u_r=-2000;

p(1:1:(n+1)/2)=p_l;
p((n+3)/2:1:n)=p_r;
rho(1:1:(n+1)/2)=rho_l;
rho((n+3)/2:1:n)=rho_r;
u(1:1:(n+1)/2)=u_l;
u((n+3)/2:1:n)=u_r;
T(1:n)= 298;
rho0v(1:n)=rho0;

e0v=specheat*T;                                      %Undeformed internal energy
E=e0v+0.5*p.*(1./rho0v-1./rho)+0.5*u.^2;              %Total Energy
a=sqrt(G./rho);                                     %Speed of sound
a0=sqrt(G/rho0);
dt=CFL*h/max(abs(u+a));
step=0;

%Time iteration
%------------------------------------------------------------------------%
for t=dt:dt:t_final
    %Define q & F matrix
    q=[rho; rho.*u; rho.*E];
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)];
    %Update q matrix and flow parameters
    q(1:3,2:n-1)=0.5*(q(1:3,3:n)+q(1:3,1:n-2))-dt/(2*h)*(F(1:3,3:n)-F(1:3,1:n-2));
    rho=q(1,1:n);
    u=q(2,1:n)./rho(1:n);
    E=q(3,1:n)./rho;
    %Using linear shock speed model:
    %Calculate pressure and energy along shock Hugoniot
    xi= 1./rho0v-1./rho;
    pH= a0^2*xi./(1./rho0v-Sigma*xi).^2;
    eH= e0v+(pH/2).*xi;
    %Calculate pressure based on first-order Mie-Gruneisen approx
    p=pH+Grun*rho.*(E-eH);
    step=step+1;
end
%calculation of flow parameters
a=sqrt(G./rho);
M=u./a;
p_ref=101325;           %Reference air pressure (N/m^2)
rho_ref=1.225;          %Reference air density (kg/m^3)
%s=1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho));   %Entropy w.r.t reference condition
Q=rho.*u;                                               %Mass Flow rate per unit area
%------------------------------------------------------------------------%
offset=0.05;
subplot(231);plot(x,p,'k');xlabel('X-Coordinate (m)');ylabel('Pressure (Pa)');ylim([min(p)-(offset)*max(p) (1+offset)*max(p)]);
%subplot(232);plot(x,s,'k');xlabel('X-Coordinate (m)');ylabel('Entropy/R gas');ylim([min(s)-(offset)*max(s) (1+offset)*max(s)]);
subplot(233);plot(x,u,'k');xlabel('X-Coordinate (m)');ylabel('Velocity (m/s)');ylim([min(u)-(offset)*max(u) (1+offset)*max(u)]);
subplot(234);plot(x,M,'k');xlabel('X-Coordinate (m)');ylabel('Mach number');ylim([min(M)-(offset)*max(M) (1+offset)*max(M)]);
subplot(235);plot(x,rho,'k');xlabel('X-Coordinate (m)');ylabel('Density (kg/m^3)');ylim([min(rho)-(offset)*max(rho) (1+offset)*max(rho)]);
subplot(236);plot(x,Q,'k');xlabel('X-Coordinate (m)');ylabel('Mass Flow (kg/m^2s)');ylim([min(Q)-(offset)*max(Q) (1+offset)*max(Q)]);
%------------------------------------------------------------------------%
step