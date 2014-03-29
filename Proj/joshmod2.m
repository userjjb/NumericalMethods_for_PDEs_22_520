%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Josh Bevan 2013
%Num. Methods for PDEs, 22.520
%
%Solves shock wave propogation through a 1D solid
%
%Uses Mie-Gruneisen constitutive EoS to model compressive effects
%Uses Rankine-Hugoniot jump conditions to bridge shock discontinuity
%Uses Lax-Friedrich flux scheme to time-iterate conserved quantities
%Based on code by Anand Dhariya, UMich
%
%Theory heavily based on work by Geoffrey Ward, Caltech
%http://thesis.library.caltech.edu/6211/1/ward_thesis.pdf
%
%Limitations:
%-First-order approx for EoS limited to compressive states where rho>rho0
%If tensile effects need to be resolved try second-order Murnaghan isentrope
%(see Ward pg 8)
%-Artificial max density due to singulatiry (Ward pg 9)
%-Artificial minimum pressure (Ward pg 9)

clear;
clc;
close all;

%Model parameters
%------------------------------------------------------------------------%
n=10001;            %Number of grid points
L=1;                %Length of domain
h=L/(n-1);          %Spatial step size
CFL=0.57;           %CFL number for stability
t_final=8000e-8;    %Final time in sim (secs)
x=0:h:L;

%Material params
%------------------------------------------------------------------------%
Material='Aluminum'
G= 76e9;            %Shear Modulus (Pa)
rho0= 2785;         %Density (kg/m^3)
Grun= 2.0;          %Gruneisen parameter
Sigma= 1.338;       %related to dK/dp uses linear assumption for Us to Up
specheat= 897;      %Cv (J/kg K)

%Initial conditions
%------------------------------------------------------------------------%
%Separated into left and right domains with x=0.5 dividing
p_l=-6.4e+10;            
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
T(1:n)= 0;
rho0v(1:n)=rho0;

e0v=specheat*T;                                     %Undeformed internal energy
%Total internal Energy
E=e0v+0.5*p.*(1./rho0v-1./rho);%+0.5*u.^2; Kinetic energy should not be included, it is not internal energy
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
    q(1:3,2:n-1)=0.5*(q(1:3,1+2:n)+q(1:3,1:n-2))-dt/(2*h)*(F(1:3,1+2:n)-F(1:3,1:n-2));
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

%Info for plotting
a=sqrt(G./rho);
M=u./a;
s=specheat*log(rho./rho0);   %Entropy w.r.t reference condition, this is likely incorrect
Q=rho.*u;                    %Mass Flow rate per unit area
%------------------------------------------------------------------------%
offset=0.05;
subplot(231);plot(x,p,'k');xlabel('X-loc (m)');ylabel('Pressure (Pa)');ylim([min(p)-(offset)*max(p) (1+offset)*max(p)]);
%subplot(232);plot(x,s,'k');xlabel('X-loc (m)');ylabel('Entropy');ylim([min(s)-(offset)*max(s) (1+offset)*max(s)]);
subplot(233);plot(x,u,'k');xlabel('X-loc (m)');ylabel('Velocity (m/s)');ylim([min(u)-(offset)*max(u) (1+offset)*max(u)]);
subplot(234);plot(x,M,'k');xlabel('X-loc (m)');ylabel('Mach number');ylim([min(M)-(offset)*max(M) (1+offset)*max(M)]);
subplot(235);plot(x,rho,'k');xlabel('X-loc (m)');ylabel('Density (kg/m^3)');ylim([min(rho)-(offset)*max(rho) (1+offset)*max(rho)]);
subplot(236);plot(x,Q,'k');xlabel('X-loc (m)');ylabel('Mass Flow (kg/m^2s)');ylim([min(Q)-(offset)*max(Q) (1+offset)*max(Q)]);
%------------------------------------------------------------------------%
step