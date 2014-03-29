%Nodes are numbered as follows:
%    1  2  3  4
% O==O==O==O==O
%The element is associated with the node to its right, they share the same
%'i'

close all
clear all
clc

rhoMax= 1.0;
rhoL= 0.8;
rhoR= 0.0;
uMax= 1.0;
dx= 1/100;
dt= 0.8*dx/uMax;
TimeIncrements= 9/dt; %Mimimum time run to ensure steady state at all points for moving average of flux
xLeft= -1;
xRight= 1;

%Allow user to select flux function
[FluxFunction,ok] = listdlg('PromptString','Select flux function:',...
                            'SelectionMode','single',...
                            'ListString',{'Godunov';'Average';'Lax-Friedrich';'Richtmyer'},...
                            'ListSize',[160 100]);
                        
if FluxFunction==2 %Set finer timestep for average to see errors propogate
    warndlg('Setting smaller dt to show how errors propogate','Attention')
    key=2;
    pause(1)
    dt=dt/10;
    TimeIncrements= 2/dt;
    new=0;
end
                       
% Spatial Discretization
x=xLeft:dx:xRight;

% Initial Conditions
mid=(length(x)-1)/2;
rho(1:mid+1)=rhoL;
rho(mid+2:2*mid+1)=rhoR;

%Pre-define matrices for speed
saved= zeros(1,2/dt);
FluxNodal= zeros(1,2*mid);

%Flux triggers
lights(1,:)= [2/dt, 1/dt, mid+1, 0]; %Sets flux=0 with [freq, duty cycle, x_loc, offset]
freq=1;
trig=2;
loc=3;
offset=4;

counter=0;
    %To optimize light2 offset uncomment this
for iter=0:.005:2
    total=0;
    counter=counter+1;
    %To optimize light2 offset uncomment this
lights(2,:)= [2/dt, 1/dt, mid+1+(.25/dx), iter/dt];

for t=0:TimeIncrements;
    % Calculate flux inside elements
	FluxElemental = uMax.*rho.*(1-rho./rhoMax);		

	% Spatial loop
	for i= 1:length(x)-1;
%         switch FluxFunction
%             case 1 %Godunov
                if rho(i)<rho(i+1) %Choked downstream
                    FluxNodal(i) = min(FluxElemental(i), FluxElemental(i+1));
                elseif rho(i)>=rho(i+1) %Unchoked
                    if rho(i)>rhoMax/2 && rho(i+1)<rhoMax/2 %Check to see if quadratic peak bounded between
                        FluxNodal(i) = rhoMax*uMax/4;
                    else
                        FluxNodal(i) = max(FluxElemental(i), FluxElemental(i+1));
                    end
                end
%             case 2 %Average
%                 FluxNodal(i)=mean([FluxElemental(i), FluxElemental(i+1)]);
%             case 3 %Lax-Friedrich
%                 FluxNodal(i)=0.5*(FluxElemental(i)+FluxElemental(i+1))-0.5*((dx/dt)*(rho(i+1)-rho(i)));
%             case 4 %Richtmyer
%                 rhoRicht= ((dx/dt)*(FluxElemental(i)-FluxElemental(i+1))+rho(i+1)+rho(i))/2;
%                 FluxNodal(i)= rhoRicht*(1-rhoRicht);
%         end
    end
    
    %Trigger for each of the "lights" with a particular frequency, at a
    %specified time and place
    for l=1:size(lights,1)
        if mod(t-lights(l,offset),lights(l,freq)) >= lights(l,trig)
            FluxNodal(lights(l,loc)) = 0;
        end
    end
    
%     Calculate moving average of flux through desired element
    %Continual moving average method
    saved(floor(mod(t,2/dt)+1))= FluxElemental(50);
    FluxAvg= sum(saved)*(dt/2);
    %Average at "end" of time period
%     if t>=375
%         total= total+FluxElemental(50);
%     end
	
	% Spatial loop for new Rho vals
	for i = 2:length(x)-1
		rho(i) = rho(i) - dt/dx*(FluxNodal(i) - FluxNodal(i-1));
    end
% 	plot(x, rho,'o-')
% 	hold on
% 	plot([-1,1],[1.2,-.2],'w.')
% 	hold off
%     text(1,1,num2str(t*dt));
%     text(1,.8,num2str(FluxAvg));
%     pause(0.001)
end
    %To optimize light2 offset uncomment this
FluxAvg2=total*(dt/2);
savedAvg(counter)= FluxAvg;
end
plot(0:.005:2,savedAvg./max(savedAvg))