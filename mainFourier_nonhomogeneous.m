%% Project model reduction
clear all; close all; clc

% Properties
lx = 0.04;
ly = 0.08;
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
Xqq = 0:xstep/10:Lx;
Yqq = 0:ystep/10:Ly;

% Preallocation
T = zeros(length(X),length(Y),length(time));
a0 = zeros(K+1,L+1);
e = zeros(K+1,L+1);
a = zeros(K+1,L+1,length(time));


% User parameters
show_visuals = false; % Will show all the visualization plots subsequently
input.switch = false; % Turn the input source on or off
input.par.type = 'sine'; % {const,sine} What type of input
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
X1 = Lx/4 ; Y1 = Ly/2; % Heat source 1 location
X2 = (3*Lx)/4; Y2 = Ly/2; % Heat source 2 location

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


%% Define nonhomegenous rho and c
lxindex = (Lx/2-lx/2):xstep:(Lx/2+lx/2);
lyindex = (Ly/2-ly/2):xstep:(Ly/2+ly/2);
indlx = discretize(lxindex,X);
indly = discretize(lyindex,Y);
M = ones(length(X),length(Y))*rho(1)*c(1);
M(17:25,22:39) = rho(2)*c(2);
M_test = reshape(M,length(X)*length(Y),1);
% M = M./(max(M));

% Create interpolated cubic spline representation

[X_mesh,Y_mesh] = meshgrid(X,Y);
[Xq,Yq] = meshgrid(Xqq,Yqq);
Vq = interp2(X_mesh,Y_mesh,M',Xq,Yq,'cubic');

% Make phi_kl 4D matrix for interpolation
phi_kl_intep = zeros(length(Xqq),length(Yqq),K+1,L+1);
for x = 1:length(Xqq)
    for y = 1:length(Yqq)
        for k = 0:K
            for l = 0:L
                phi_kl_interp(x,y,k+1,l+1) = basisxy(Xqq(x),Yqq(y),k,l,Lx,Ly);
            end
        end
    end
end

for k = 0:K 
    for l = 0:L
    e(k+1,l+1) = sum(Vq'.*phi_kl_interp(:,:,k+1,l+1),'all')*xstep*ystep/100; 
    end
end
einv = pinv(e);
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

% Test plate
sumE = 0;
for k = 0:K
    for l = 0:L
    sumE = sumE + e(k+1,l+1)*phi_kl(:,:,k+1,l+1); 
    end
end
figure()
surf(sumE)
hold on
surf(M)

%% Initial conditions and ODE solver for a
for k = 0:K 
    for l = 0:L
        a0(k+1,l+1) = sum(T0.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
    end
end

% for k = 1:K
%     for l = 1:L
%     [phi_kl_dx, phi_kl_dy] = gradient(phi_kl(:,:,k,l),xstep,ystep);
%     [phi_kl_ddx(:,:,k,l), phi_kl_dxdy(:,:,k,l)] = gradient(phi_kl_dx,xstep,ystep);
%     [phi_kl_dydx(:,:,k,l), phi_kl_ddy(:,:,k,l)] = gradient(phi_kl_dy,xstep,ystep);
%     end
% end
% rhoc=0;
% 
% for k = 0:K
%     for l=0:L
%         rhoc = rhoc + sum(e(k+1,l+1).*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
%     end
% end

A = zeros(K,L);
for k = 0:K
    for l = 0:L
    A(k+1,l+1) = -kappa(1)*((k^2*pi/Lx^2)+(l^2*pi^2/Ly^2));
    end
end
A = einv*A;
A = reshape(A,(K+1)*(L+1),1);
a0 = reshape(a0,(K+1)*(L+1),1);
A_ = diag(A);

%% Without input
% B = eye((K+1)*(L+1),1)*0;
% C = eye(1,(K+1)*(L+1))*0;
% sys = ss(A_,B,C,0);
% [y,~,x] = initial(sys,a0,time);
% x = reshape(x,length(time),K+1,L+1);
%% With input
source = zeros(length(time),2);
for k = 0:K
    for l = 0:L
        B_1(k+1,l+1) = sum(input.u1.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
        B_2(k+1,l+1) = sum(input.u2.*phi_kl(:,:,k+1,l+1),'all')*xstep*ystep;
    end
end
B_1 = einv*B_1; 
B_2 = einv*B_2;
B_1 = reshape(B_1,(K+1)*(L+1),1);
B_2 = reshape(B_2,(K+1)*(L+1),1);
B = [B_1 B_2];

if input.switch
[u1,u2] = heatInput(time,input.par);
source = [u1' u2'];
end

C = eye(1,(K+1)*(L+1))*0;
sys = ss(A_,B,C,0);
[y,~,x] = lsim(sys,source,time,a0,'zoh');
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
