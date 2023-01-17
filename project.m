%% Project model reduction Michiel Wind, Nick Cruijsen
clear all; close all; clc
% Parameters
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
show_visuals = true;
Nx = 30; % Number of samples x,y,t
Ny = 30;
Nt = 100;
tend = 10;
tstep = tend/(Nt-1);
xstep = Lx/(Nx-1);
ystep = Ly/(Ny-1);
X = 0:xstep:Lx;
Y = 0:ystep:Ly;
time = 0:tstep:tend;
K = 2; 
L = 2;

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(K+1,L+1);
a = zeros(K+1,L+1,length(time));

% Initial temperature
% T0 = ones(length(X),length(Y));

k=2;
l=2;
[T0,T0dx,T0dy] = initialTemp(X,Y,k,l,'gauss',true);

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


%% Simulation (must go wrong here somewhere
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

[X_mesh,Y_mesh] = meshgrid(X,Y);
for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T(:,:,t));  
%     axis([0 Lx 0 Ly 1 1.2]);
    title(sprintf('Time = %f seconds', t));
    pause(0.01)
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations for understanding %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_visuals

% Solutions of aSol over time in 1 plot
figure()
for k = 0:K
    for l = 0:L
       plot(squeeze(a(k+1,l+1,:)))
       hold on
    end
end

% Plot initial conditions a_kl
figure()
surf(a0)
xlabel('k')
ylabel('l')

% Plot basis 1D
for x = 1:length(X)
    for k = 0:K
        phix(x,k+1) = basisx(X(x),k,Lx);
    end
end

figure()
hold on
for k = 0:K
    txt = ['K = ',num2str(k)];
    plot(phix(:,k+1),'DisplayName',txt)
end
hold off
legend show

% plot basis 2D
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

%% Show intital temperature of the plate
figure()
subplot(1,3,1)
surf(X,Y,T0)
xlabel('x')
ylabel('y')
title('Intial temperature')
subplot(1,3,2)
surf(X,Y,T0dx)
xlabel('x')
ylabel('y')
title('dT0/dx')
subplot(1,3,3)
surf(X,Y,T0dy)
xlabel('x')
ylabel('y')
title('dT0/dy')

end
