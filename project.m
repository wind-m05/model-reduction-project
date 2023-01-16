%% Project model reduction Michiel Wind, Nick Cruijsen
close all; clear all; clc
%% Parameters
Lx  = 0.2;
Ly  = 0.3;
rho = [2328 2300]; % Material densities (yellow,blue respectively)
c   = [700 680]; % Heat capacities
kappa = [148 148]; % Thermal conductivities
X1 = Lx/4 ; Y1 = Ly/2; % Heat source 1 location
X2 = (3*Lx)/4; Y2 = Ly/2; % Heat source 2 location
W = 0.05; % Actuator width [m]
T_amb = 309; % Ambient temp Kelvin
T_K2C = 273;

%% Exercise 5
Nx = 30; % Number of samples x,y,t
Ny = 30;
Nt = 100;
t_end = 50;
t_step = t_end/(Nt-1);
x_step = Lx/(Nx-1);
y_step = Ly/(Ny-1);
X = 0:x_step:Lx;
Y = 0:y_step:Ly;
time = 0:t_step:tend;
K = 10; 
L = 10;
T = zeros(length(X),length(Y),length(time));

aSol = zeros(K+1,L+1,length(time));
T0 = ones(length(X),length(Y));

% Initial conditions
a_kl = zeros(L+1,K+1);
for x = 1:length(X)
    for y = 1:length(Y)
        for k = 0:K
            for l = 0:L
             a_kl(k+1,l+1) = a_kl(k+1,l+1) + T0(x,y)*basisx(X(x),k,Lx)*basisy(Y(y),l,Ly);
            end
        end
    end
end

% a_kl = a_kl/a_kl(1,1);
% a_kl(k+1,l+1)
%Compute a_k,l for till end time for each k and l
for k = 0:K
    for l = 0:L
        a = aODE(time,a_kl(k+1,l+1),kappa(1),rho(1),c(1),Lx,Ly,k,l);
        aSol(k+1,l+1,:) = a;
    end
end



% Simulation
for t = 1:length(time)
    for x = 1:length(X)
        for y = 1:length(Y)
            sum = 0;
            for l = 1:L 
                for k = 1:K
                sum = sum + (aSol(k,l,t)*basisx(X(x),k,Lx)*basisy(Y(y),l,Ly));
                end
            end
            T(x,y,t) = sum + T_amb - T_K2C;
        end
    end
end

[X_mesh,Y_mesh] = meshgrid(Y,X);
for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T(:,:,t))
    title(sprintf('Time = %f seconds', t));
    pause(0.01)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations for understanding %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solutions of aSol over time in 1 plot

for k = 0:K
    for l = 0:L
       plot(squeeze(aSol(k+1,l+1,:)))
       hold on
    end
end

%% Plot initial conditions a_kl
surf(a_kl)
%% Plot basis 1D

for x = 1:length(X)
    for k = 0:K
        phix(x,k+1) = basisx(X(x),k,Lx);
    end
end

figure()
for k = 0:K
    plot(phix(:,k+1))
    hold on
end

%% plot basis 2D
% Xtest = -0.3:0.01:0.3;
% Ytest = -0.3:0.01:0.3;
% X = Xtest;
% Y = Ytest;

for x = 1:length(X)
    for y = 1:length(Y)
        for k = 0:2
            for l = 0:2
             phi_kl(y,x,k+1,l+1) = basisx(X(x),k,Lx)*basisy(Y(y),l,Ly);
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

