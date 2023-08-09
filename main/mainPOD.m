%% Main POD
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

%% Define basis
T_snap = load('T_snap.mat').T; % Load data for POD basis
reduced_energy_remaint = 0.99; % User input - Ratio of information included into the reduced order model
[phiPOD,diagn] = PODbasis(T_snap,Nx,Ny,xstep,ystep,Nt,reduced_energy_remaint);

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(diagn.R,1);
a = zeros(diagn.R,Nt);

% User parameters
input.switch = false; % User input Turn the input source on or off
input.par.type = 'const'; % User input - {const,sine} constant/sinusoidal input
input.par.freq = 0.01; % User input [Hz]
input.par.tstart = 0; % User input [s]
input.par.tend = tend; % User input [s]
input.par.amp1 = 0.4; % User input [-]
input.par.amp2 = 0.4; % User input


% Initial temperature
kinit=2; % Frequency of basis in x
linit=2; % Frequency of basis in y
[T0,T0dx,T0dy] = initialTemp(X,Y,kinit,linit,'gauss',true);  % User input choose ...
% ({gauss,blockup,blockdown,fourier},normalization = true/false)

%% Calculate phiPOD for x,y positions overlapping with u

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


%% Initial conditions and ODE solver for a
for r = 1:diagn.R 
        a0(r) = sum(T0.*phiPOD.xy(:,:,r),'all')*xstep*ystep;
end

%% With input
A = kappa(1)/(rho(1)*c(1))*phiPOD.dotp;
source = zeros(length(time),2);
for r = 1:diagn.R 
    B_1(r) = sum(input.u1.*phiPOD.xy(:,:,r),'all')*xstep*ystep;
    B_2(r) = sum(input.u2.*phiPOD.xy(:,:,r),'all')*xstep*ystep;
end
B_1 = reshape(B_1,diagn.R,1);
B_2 = reshape(B_2,diagn.R,1);
B = (1/rho(1)*c(1)).*[B_1 B_2];
if input.switch
[u1,u2] = heatInput(time,input.par);
source = [u1' u2'];
end
C = eye(1,diagn.R)*0;
sys = ss(A,B,C,0);
[y,~,x] = lsim(sys,source,time,a0);

%% Temperature over time
for t = 1:length(time)
    sumT = 0;
    for r = 1:diagn.R
        sumT = sumT + x(t,r)*phiPOD.xy(:,:,r); 
    end
    T(:,:,t) = sumT;
end