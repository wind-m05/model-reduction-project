%% POD basis 
clear all, close all, clc

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
tend = 600; % 10 minutes
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

% User parameters
show_visuals = false; % Will show all the visualization plots subsequently
input.switch = false; % Turn the input source on or off
input.par.type = 'sinusoid'; % {const,sinusoid} What type of input
input.par.freq = 0.1; % [Hz]
input.par.tstart = 5; % [s]
input.par.tend = tend; % [s]
input.par.amp1 = 10; % [-]
input.par.amp2 = 10;


% Initial temperature
kinit=2; % Frequency of basis in x
linit=2; % Frequency of basis in y
[T0,T0dx,T0dy] = initialTemp(X,Y,kinit,linit,'gauss',true);

T_snap = load('T_snap.mat').T;

reduced_energy_remaint = 0.9; % Ratio of information included into the reduced order model

[phiPOD,diagn] = PODbasis(T_snap,Nx,Ny,xstep,ystep,Nt,reduced_energy_remaint);

% Calculate phiPOD for x,y positions overlapping with u
phi_u1 = zeros(diagn.R);
phi_u2 = zeros(diagn.R);
a0 = zeros(diagn.R,1);
a = zeros(diagn.R,Nt);

yindex = (Ly/2-W/2):ystep:(Ly/2+W/2);
xindex1 = (Lx/4-W/2):xstep:(Lx/4+W/2);
xindex2 = ((3*Lx/4)-W/2):xstep:((3*Lx/4)+W/2);

indx1 = discretize(xindex1,X);
indx2 = discretize(xindex2,X);
indy  = discretize(yindex,Y);  

% for r = 1:diagn.R
%     phi_u1(r) = sum(phiPOD(indx1,indy,r),'all')*xstep*ystep;
%     phi_u2(r) = sum(phiPOD(indx2,indy,r),'all')*xstep*ystep;
% end



% Initial conditions and ODE solver for a
for r = 1:diagn.R 
        a0(r) = sum(T0.*phiPOD.xy(:,:,r),'all')*xstep*ystep;
%         a(r,:) = aODE_POD(time,a0(r),kappa(1),rho(1),c(1),Lx,Ly,diagn.R,phiPOD,phi_u1,phi_u2,input);
end
a= aODE_POD(time,a0,kappa(1),rho(1),c(1),Lx,Ly,diagn.R,phiPOD,phi_u1,phi_u2,input);

T = zeros(length(X),length(Y),length(time));
% Temperature over time
for t = 1:length(time)
    sumT = 0;
    for r = 1:diagn.R
        sumT = sumT + a(r,t)*phiPOD.xy(:,:,r);
    end
    T(:,:,t) = sumT;
end
