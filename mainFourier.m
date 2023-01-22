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
K = 20; 
L = 20;

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
save('T_snap.mat','T')

















%% I seperated the scripts below into new scripts to keep it organised, keep it for now for a backup
% %% POD basis (without input) 
% 
% 
% reduced_energy_remaint = 0.9; % Ratio of information included into the reduced order model
% 
% 
% [phiPOD,diagn] = PODbasis(T,Nx,Ny,Nt,reduced_energy_remaint);
% 
% % Calculate phiPOD for x,y positions overlapping with u
% phi_u1 = zeros(diagn.R);
% phi_u2 = zeros(diagn.R);
% a0 = zeros(diagn.R,1);
% a = zeros(diagn.R,Nt);
% 
% yindex = (Ly/2-W/2):ystep:(Ly/2+W/2);
% xindex1 = (Lx/4-W/2):xstep:(Lx/4+W/2);
% xindex2 = ((3*Lx/4)-W/2):xstep:((3*Lx/4)+W/2);
% 
% indx1 = discretize(xindex1,X);
% indx2 = discretize(xindex2,X);
% indy  = discretize(yindex,Y);  
% 
% for r = 1:diagn.R
%     phi_u1(r) = sum(phiPOD(indx1,indy,r),'all')*xstep*ystep;
%     phi_u2(r) = sum(phiPOD(indx2,indy,r),'all')*xstep*ystep;
% end
% 
% 
% 
% % Initial conditions and ODE solver for a
% for r = 1:diagn.R 
%         a0(r) = sum(T0.*phiPOD(:,:,r),'all');
% %         a0(r) = sum(T0.*phiPOD(:,:,r),'all')*xstep*ystep;
%         a(r,:) = aODE_POD(time,a0(r),kappa(1),rho(1),c(1),Lx,Ly,diagn.R,phi_u1,phi_u2,input);
% end
% 
% TPOD = zeros(length(X),length(Y),length(time));
% % Temperature over time
% for t = 1:length(time)
%     sumT = 0;
%     for r = 1:diagn.R
%         sumT = sumT + a(r,t)*phiPOD(:,:,r);
%     end
%     TPOD(:,:,t) = sumT;
% end
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if show_visuals
% %% full simulation
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% TaxisMin = min(min(T0));
% TaxisMax = max(max(T0));
% TaxisMin_switch = min(min(T0));
% TaxisMax_switch = max(max(T(:,:,end)));
% font = 15;
% fps = 60;
% figure()
% for t = 1:length(time)
%     mesh(X_mesh,Y_mesh,T(:,:,t));  
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
% 
%     title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
%     xlabel('x [m]',Interpreter='latex',FontSize=font); 
%     ylabel('y [m]',Interpreter='latex',FontSize=font); 
%     zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
%     pause(1/fps)
% end
% 
% %% Initial behaviour
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% time_redux = 0.1; % percentage of shown time instances
% TaxisMin = min(min(T0))-0.2;
% TaxisMax = max(max(T0))+0.2;
% TaxisMax_switch = max(max(T(:,:,round(Nt*time_redux)+1)))
% font = 15;
% fps = 60;
% figure()
% for t = 1:round((Nt*time_redux))
%     mesh(X_mesh,Y_mesh,T(:,:,t));
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
%     title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
%     xlabel('x [m]',Interpreter='latex',FontSize=font); 
%     ylabel('y [m]',Interpreter='latex',FontSize=font); 
%     zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
%     colorbar
%     pause(1/fps)
% end
% 
% %% End behaviour
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% time_redux = 0.2; % percentage of shown time instances
% TaxisMin = min(min(T(:,:,end)))-min(min(T(:,:,end)))*0.025;
% TaxisMax = max(max(T(:,:,end)))+max(max(T(:,:,end)))*0.025;
% TaxisMin_switch = min(min(T(:,:,Nt-round((Nt*time_redux)))))-min(min(T(:,:,Nt-round((Nt*time_redux)))))*0.025;
% font = 15;
% fps = 60;
% figure()
% for t = Nt-round((Nt*time_redux)):Nt
%     mesh(X_mesh,Y_mesh,T(:,:,t));
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
%     title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
%     xlabel('x [m]',Interpreter='latex',FontSize=font); 
%     ylabel('y [m]',Interpreter='latex',FontSize=font); 
%     zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
%     pause(1/fps)
% end
% 
% %% Sample code for fixed colorbar
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% TaxisMin = min(min(T0));
% TaxisMax = max(max(T0));
% TaxisMax_switch = max(max(T(:,:,end)));
% font = 15;
% fps = 60;
% figure()
% for t = 1:Nt
%     mesh(X_mesh,Y_mesh,T(:,:,t));  
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
%     title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
%     xlabel('x [m]',Interpreter='latex',FontSize=font); 
%     ylabel('y [m]',Interpreter='latex',FontSize=font); 
%     zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
%     caxis manual;          % allow subsequent plots to use the same color limits
%     caxis([TaxisMin TaxisMax]); 
%     colorbar;
%     pause(1/fps)
% end
% caxis auto; 
% %% Solutions of ODE solver for a_k,l over time
% figure()
% for k = 0:K
%     for l = 0:L
%        plot(squeeze(a(k+1,l+1,:)))
%        hold on
%     end
% end
% xlabel('Time [t]')
% ylabel('[-]')
% title('Time solutions of a_k,l')
% grid on
% 
% %% Plot initial conditions a_kl
% figure()
% k = 0:K
% l = 0:L
% surf(k,l,a0)
% colormap default
% colorbar
% xlabel('k')
% ylabel('l')
% title('Initial conditions for a, (dot product of T0 wrt phi_k,l)')
% 
% %% Plot basis 1D
% for x = 1:length(X)
%     for k = 0:K
%         phix(x,k+1) = basisx(X(x),k,Lx);
%     end
% end
% 
% figure()
% hold on
% for k = 0:K
%     txt = ['K = ',num2str(k)];
%     plot(phix(:,k+1),'DisplayName',txt,'Linewidth',1.5)
% end
% grid on
% xlabel('x [m]')
% ylabel('[-]')
% title('2D basis functions on spatial domain x looped over K')
% hold off
% legend show
% 
% %% plot basis 2D
% % To test and include quadrants Q2,Q3 and Q4
% % Xtest = -0.3:0.01:0.3;
% % Ytest = -0.3:0.01:0.3;
% % X = Xtest;
% % Y = Ytest;
% 
% figure()
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% for k = 0:2
%     for l = 0:2
%         mesh(X_mesh,Y_mesh,phi_kl(:,:,k+1,l+1))
%         hold on
%     end
% end
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('[-]')
% title('Basis functions \phi_{k,l} on spatial domain x,y and looped over k and l')
% 
% %% Show intital temperature of the plate
% figure()
% surf(X,Y,T0')
% xlabel('x')
% ylabel('y')
% zlabel('Temperature [-]')
% title('Intial temperature')
% 
% %% Show partial derivatives of the intital temperature 
% figure()
% subplot(1,2,1)
% surf(X,Y,T0dx')
% xlabel('x')
% ylabel('y')
% title('dT0/dx')
% subplot(1,2,2)
% surf(X,Y,T0dy')
% xlabel('x')
% ylabel('y')
% title('dT0/dy')
% 
% %% Runtime simulation visual
% psi = [1,2,5,7,10,20]; % Increasing order of K and L
% t = [1,2.11,9.98,18.35,34.6,142]; % Increasing runtime
% p = polyfit(psi,t,5)
% figure()
% plot(psi,t)
% psi_des = 1:100;
% f = polyval(p,psi_des);
% figure()
% plot(psi_des,f./(3600))
% xlabel('order (K and L)')
% ylabel('time [h]')
% 
% %% Show singular values of POD basis
% figure()
% plot(diagn.svp,'ko','Linewidth',2)
% figure()
% plot(real(diagn.U(:,1:3)))
% title('The first 3 dominant modes')
% 
% %% Check if the POD basis is orthogonal in its columns
% % test = sum(phiPOD(:,:,1)'*phiPOD(:,:,2),'all')*xstep*ystep
% 
% end
