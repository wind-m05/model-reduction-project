%% Project model reduction
clear all; close all; clc

% Properties
Lx  = 0.2;
Ly  = 0.3;
rho = [2328 2300]; % Material densities (yellow,blue respectively)
c   = [700 680]; % Heat capacities
kappa = [148 148]; % Thermal conductivities
X1 = Lx/4 ; Y1 = Ly/2; % Heat source 1 location
X2 = (3*Lx)/4; Y2 = Ly/2; % Heat source 2 location
W = 0.05; % Actuator width [m]
Tamb = 309; % Ambient temp Kelvin
TK2C = 273;

% Simulation parameters
show_visuals = false; % Will show all the visualization plots subsequently
Nx = 100; % Number of samples x,y,t
Ny = 100;
Nt = 1000;
tend = 600; % 10 minutes
tstep = tend/(Nt-1);
xstep = Lx/(Nx-1);
ystep = Ly/(Ny-1);
X = 0:xstep:Lx;
Y = 0:ystep:Ly;
time = 0:tstep:tend;
K = 5; 
L = 5;

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(K+1,L+1);
a = zeros(K+1,L+1,length(time));

% Initial temperature
k=2; % Frequency of basis in x
l=2; % Frequency of basis in y
[T0,T0dx,T0dy] = initialTemp(X,Y,k,l,'fourier',true);

%% Initial conditions
for k = 0:K 
    for l = 0:L
        for x = 1:length(X)
            for y = 1:length(Y)
             a0(k+1,l+1) = a0(k+1,l+1) + T0(x,y)*basisx(X(x),k,Lx)*basisy(Y(y),l,Ly);
            end
        end
        a0(k+1,l+1) = a0(k+1,l+1)*xstep*ystep;
        a(k+1,l+1,:) = aODE(time,a0(k+1,l+1),kappa(1),rho(1),c(1),Lx,Ly,k,l);
    end
end

%% Simulation 
for t = 1:length(time)
    for x = 1:length(X)
        for y = 1:length(Y)
            sumT = 0;
            for k = 0:K
                for l = 0:L 
                sumT = sumT + (a(k+1,l+1,t)*basisxy(X(x),Y(y),k,l,Lx,Ly));
                end
            end
            T(x,y,t) = sumT;
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_visuals

%% full simulation
[X_mesh,Y_mesh] = meshgrid(X,Y);
TaxisMin = min(min(T0))-1;
TaxisMax = max(max(T0))+1;
font = 15;
figure()
for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T(:,:,t));  
    axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    pause(0.001)
end

%% Initial behaviour
[X_mesh,Y_mesh] = meshgrid(X,Y);
time_redux = 0.05; % percentage of shown time instances
TaxisMin = min(min(T0))-0.2;
TaxisMax = max(max(T0))+0.2;
font = 15;
figure()
for t = 1:round((length(time)*time_redux))
    mesh(X_mesh,Y_mesh,T(:,:,t));  
    axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    colorbar
    pause(0.01)
end

%% End behaviour
[X_mesh,Y_mesh] = meshgrid(X,Y);
time_redux = 0.2; % percentage of shown time instances
TaxisMin = min(min(T(:,:,end)))-min(min(T(:,:,end)))*0.25;
TaxisMax = max(max(T(:,:,end)))+max(max(T(:,:,end)))*0.25;
font = 15;
figure()
for t = time(end)-round((length(time)*time_redux)):time(end)
    mesh(X_mesh,Y_mesh,T(:,:,t));  
    axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    pause(0.01)
end

%% Sample code for fixed colorbar
[X_mesh,Y_mesh] = meshgrid(X,Y);
TaxisMin = min(min(T0));
TaxisMax = max(max(T0));
font = 15;
figure()

for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T(:,:,t));  
    axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    caxis manual;          % allow subsequent plots to use the same color limits
    caxis([TaxisMin TaxisMax]); 
    colorbar;
    pause(0.001)
end
caxis auto; 
%% Solutions of ODE solver for a_k,l over time
figure()
for k = 0:K
    for l = 0:L
       plot(squeeze(a(k+1,l+1,:)))
       hold on
    end
end
xlabel('Time [t]')
ylabel('[-]')
title('Time solutions of a_k,l')
grid on

%% Plot initial conditions a_kl
figure()
k = 0:K
l = 0:L
surf(k,l,a0)
colormap default
colorbar
xlabel('k')
ylabel('l')
title('Initial conditions for a, (dot product of T0 wrt phi_k,l)')

%% Plot basis 1D
for x = 1:length(X)
    for k = 0:K
        phix(x,k+1) = basisx(X(x),k,Lx);
    end
end

figure()
hold on
for k = 0:K
    txt = ['K = ',num2str(k)];
    plot(phix(:,k+1),'DisplayName',txt,'Linewidth',1.5)
end
grid on
xlabel('x [m]')
ylabel('[-]')
title('2D basis functions on spatial domain x looped over K')
hold off
legend show

%% plot basis 2D
% To test and include quadrants Q2,Q3 and Q4
% Xtest = -0.3:0.01:0.3;
% Ytest = -0.3:0.01:0.3;
% X = Xtest;
% Y = Ytest;

for x = 1:length(X)
    for y = 1:length(Y)
        for k = 0:2
            for l = 0:2
                phi_kl(y,x,k+1,l+1) = basisxy(X(x),Y(y),k,l,Lx,Ly);
            end
        end
    end
end

figure()
[X_mesh,Y_mesh] = meshgrid(X,Y);
for k = 0:2
    for l = 0:2
        mesh(X_mesh,Y_mesh,phi_kl(:,:,k+1,l+1))
        hold on
    end
end
xlabel('x [m]')
ylabel('y [m]')
zlabel('[-]')
title('Basis functions \phi_{k,l} on spatial domain x,y and looped over k and l')

%% Show intital temperature of the plate
figure()
surf(X,Y,T0)
xlabel('x')
ylabel('y')
zlabel('Temperature [-]')
title('Intial temperature')

%% Show partial derivatives of the intital temperature 
figure()
subplot(1,2,1)
surf(X,Y,T0dx)
xlabel('x')
ylabel('y')
title('dT0/dx')
subplot(1,2,2)
surf(X,Y,T0dy)
xlabel('x')
ylabel('y')
title('dT0/dy')

end
