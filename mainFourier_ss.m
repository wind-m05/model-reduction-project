%% Project model reduction
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
K = 10; 
L = 10;

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(K+1,L+1);
a = zeros(K+1,L+1,length(time));


% User parameters
show_visuals = false; % Will show all the visualization plots subsequently
input.switch = false; % Turn the input source on or off
input.par.type = 'const'; % {const,sinusoid} What type of input
input.par.freq = 0.01; % [Hz]
input.par.tstart = 5; % [s]
input.par.tend = tend; % [s]
input.par.amp1 = 0.4; % [-]
input.par.amp2 = 0.4;


% Initial temperature
kinit=2; % Frequency of basis in x
linit=2; % Frequency of basis in y
[T0,T0dx,T0dy] = initialTemp(X,Y,kinit,linit,'gauss',true);

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

% Make phi_kl 4D matrix
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


% Initial conditions and ODE solver for a
for k = 0:K 
    for l = 0:L
        a0(k+1,l+1) = sum(T0.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
        a(k+1,l+1,:) = aODE(time,a0(k+1,l+1),kappa(1),rho(1),c(1),Lx,Ly,xstep,ystep,k+1,l+1,phi_kl,input);
    end
end

% Temperature over time
for t = 1:length(time)
    sumT = 0;
    for k = 0:K
        for l = 0:L 
        sumT = sumT + a(k+1,l+1,t)*phi_kl(:,:,k+1,l+1);
        end
    end
    T(:,:,t) = sumT;
end
% save('T_snap.mat','T')

