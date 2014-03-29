clear all;
close all;
clc;

R=287;
%Cv=717;
gamma=1.4;
rhol=1;
Ul=0;
Pl=100000;

rhor=0.125;
Ur=0;
Pr=10000;

c=1000;
delx=20/c;
delt=.0004278;
lambda=.001069;


for i=1:1:c+1
x(i)=-10+(i-1)*delx;
end

for i=1:1:c+1

if (x(i)<0)

rho(1,i)=rhol;
u(1,i)=Ul;
p(1,i)=Pl;

else
rho(1,i)=rhor;
u(1,i)=Ur;
p(1,i)=Pr;
end
end

for i=1:1:c+1

mom(1,i)=rho(1,i)*u(1,i);
%t(1,i)=p(1,i)/(rho(1,i)*R);
e(1,i)=p(1,i)/((gamma-1)*rho(1,i));
energ(1,i)=(e(1,i)+.5*u(1,i)^2);
energrho(1,i)=energ(1,i)*rho(1,i);
end

for n=2:1:30

for i=2:1:c

rho(n,i)= .5*(rho(n-1,i-1)+ rho(n-1,i+1))-.5*lambda*(u(n-1,i+1)*rho(n-1,i+1)-u(n-1,i-1)*rho(n-1,i-1));
mom(n,i)= .5*(mom(n-1,i-1)+ mom(n-1,i+1))-.5*lambda*(u(n-1,i+1)*mom(n-1,i+1)+p(n-1,i+1)-p(n-1,i-1)-u(n-1,i-1)*mom(n-1,i-1));
energrho(n,i)= .5*(energrho(n-1,i-1)+energrho(n-1,i+1))-.5*lambda*((u(n-1,i+1)*energrho(n-1,i+1)+p(n-1,i+1)*u(n-1,i+1) - u(n-1,i-1)*energrho(n-1,i-1)+p(n-1,i-1)*u(n-1,i-1))); 

energ(n,i)=energrho(n,i)/rho(n,i);
u(n,i)=mom(n,i)/rho(n,i);
e(n,i)=energ(n,i)-.5*u(n,i)^2;
%t(n,i)=e(n,i)/Cv;
p(n,i)=(gamma-1)*rho(n,i)*e(n,i);

end


rho(n,1)=rho(1,1);
rho(n,c+1)=rho(1,c+1);
mom(n,1)=mom(1,1);
mom(n,c+1)=mom(1,c+1); 
e(n,1)=e(1,1);
e(n,c+1)=e(1,c+1);
%t(n,1)=t(n,2);
%t(n,1)=e(n,1)/Cv;
%t(n,c+1)=t(n,c);
%t(n,c+1)=e(n,c+1)/Cv;
p(n,1)=p(1,1);
%p(n,1)=rho(n,1)*R*t(n,1);
p(n,c+1)=p(1,c+1);
%p(n,c+1)=rho(n,c+1)*R*t(n,c+1);
u(n,1)=u(1,1);
%u(n,1)=mom(n,1)/rho(n,1);
u(n,c+1)=u(1,c+1);
%u(n,c+1)=mom(n,c+1)/rho(n,c+1);

end


figure(1)
t=plot(x,p(30,: ),'r');
axis([-10 10 0 100000])

figure(2)
t=plot(x,u(30,: ),'r');
axis([-10 10 -150 450])

figure(3)
t=plot(x,rho(25,: ),'r');
axis([-10 10 0 1.1])