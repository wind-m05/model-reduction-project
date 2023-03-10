%% Main Fourier
% Authors: Michiel Wind and Jelle Cruijsen - TU/e 2023
clear all; close all; clc

% Properties
Lx  = 0.2; % Length of plate in x direction
Ly  = 0.3; % Length of plate in y direction
rho = [2328 2300]; % Material densities (yellow,blue respectively)
c   = [700 680]; % Heat capacities 
kappa = [148 148]; % Thermal conductivities 
X1 = Lx/4 ; Y1 = Ly/2; % Heat source 1 location
X2 = (3*Lx)/4; Y2 = Ly/2; % Heat source 2 location
W = 0.05; % Actuator width [m]
Tamb = 309; % Ambient temp Kelvin
TK2C = 273; % Kelvin to Celcius offset

% Simulation parameters
tend = 600; % User input 
Nt = tend*2;
tstep = tend/(Nt-1);
wstep = 10; % How many discrete steps does the actuator contain 
xstep = W/wstep;
ystep = W/wstep;
Nx = Lx/(xstep)+1; 
Ny = Ly/(ystep)+1;
time = 0:tstep:tend;
X = 0:xstep:Lx;
Y = 0:ystep:Ly;
K = 10; % User input
L = 10; % User input

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(K+1,L+1);
a = zeros(K+1,L+1,Nt);

% User parameters
input.switch = false; % User input - Turn the input source (on/off)-(true/false)
input.par.type = 'sine'; % User input {const,sine} constant/sinusoidal input
input.par.freq = 0.01; % [Hz] % User input 
input.par.tstart = 0; % [s] % User input 
input.par.tend = tend; % [s] % User input 
input.par.amp1 = 0.4; % [-] % User input 
input.par.amp2 = 0.4; % User input 

% Initial temperature
kinit=2; % Frequency of basis in x
linit=2; % Frequency of basis in y
[T0,T0dx,T0dy] = initialTemp(X,Y,kinit,linit,'blockup',true); % User input choose ...
% ({gauss,blockup,blockdown,fourier},normalization = true/false)

%% Calculate phi_kl for x,y positions overlapping with u

% Actuator positions
yindex = (Ly/2-W/2):ystep:(Ly/2+W/2);
xindex1 = (Lx/4-W/2):xstep:(Lx/4+W/2);
xindex2 = ((3*Lx/4)-W/2):xstep:((3*Lx/4)+W/2);

indx1 = discretize(xindex1,X);
indx2 = discretize(xindex2,X);
indy  = discretize(yindex,Y); 
input.u1 = zeros(length(X),length(Y));
input.u2 = zeros(length(X),length(Y));
input.u1(indx1,indy) = 1;
input.u2(indx2,indy) = 1;

%% Define basis
phi_kl = zeros(length(X),length(Y),K+1,L+1);
for x = 1:length(X)
    for y = 1:length(Y)
        for k = 0:K
            for l = 0:L
                phi_kl(x,y,k+1,l+1) = basisxy(X(x),Y(y),k,l,Lx,Ly);
            end
        end
    end
end

%% Initial conditions and ODE solver for a
for k = 0:K 
    for l = 0:L
        a0(k+1,l+1) = sum(T0.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
    end
end

%% With input
A = zeros(K,L);
for k = 0:K
    for l = 0:L
    A(k+1,l+1) = ((k^2*pi^2/Lx^2)+(l^2*pi^2/Ly^2));
    end
end
A = -kappa(1)/(rho(1)*c(1)).*A;
A = reshape(A,(K+1)*(L+1),1);
a0 = reshape(a0,(K+1)*(L+1),1);
A_ = diag(A);

source = zeros(length(time),2);
for k = 0:K
    for l = 0:L
        B_1(k+1,l+1) = sum(input.u1.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
        B_2(k+1,l+1) = sum(input.u2.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
    end
end
B_1 = reshape(B_1,(K+1)*(L+1),1);
B_2 = reshape(B_2,(K+1)*(L+1),1);
B = (1/rho(1)*c(1)).*[B_1 B_2];

if input.switch
[u1,u2] = heatInput(time,input.par);
source = [u1' u2'];
end

C = eye(1,(K+1)*(L+1))*0;
sys = ss(A_,B,C,0);
[y,~,x] = lsim(sys,source,time,a0);
x = reshape(x,length(time),K+1,L+1);

%% Temperature over time
for t = 1:length(time)
    sumT = 0;
    for k = 0:K
        for l = 0:L
        sumT = sumT + x(t,k+1,l+1)*phi_kl(:,:,k+1,l+1); 
        end
    end
    T(:,:,t) = sumT;
end
save('T_snap.mat','T')
