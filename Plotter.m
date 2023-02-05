%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_snap = load('T_snap_backup.mat').T;
%% Plot critical point
figure()
surf(T(:,:,35))
title('POD')
figure()
surf(T_snap(:,:,35))
title('Fourier')
%% full simulation
[X_mesh,Y_mesh] = ndgrid(X,Y);
TaxisMin = min(min(T0));
TaxisMax = max(max(T0));
TaxisMin_switch = min(min(T0));
TaxisMax_switch = max(max(T(:,:,end)));
font = 15;
fps = 60;
figure()
for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T(:,:,t));  
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
    axis([0 Lx 0 Ly -0.1 +1]) 
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    pause(1/fps)
end

%% Residual simulation
[X_mesh,Y_mesh] = ndgrid(X,Y);
for t = 1:length(time)
T_res(:,:,t) = T(:,:,t)-T_snap(:,:,t);
end
font = 15;
fps = 60;
TaxisMin = min(min(T_res(:,:,end)));
TaxisMax = max(max(T_res(:,:,end)));
figure()
for t = 1:length(time)
    mesh(X_mesh,Y_mesh,T_res(:,:,t));  
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
    axis([0 Lx 0 Ly -1 1])  
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    pause(1/fps)
end

%% H inf norm over time
for t = 1:length(time)
[U,S,V] = svd(abs(T_res(:,:,t)));
T_maxsvd(t) = max(diag(S));
end
figure()
plot(time,T_maxsvd)
%% Initial behaviour
[X_mesh,Y_mesh] = ndgrid(X,Y);
time_redux = 0.1; % percentage of shown time instances
TaxisMin = min(min(T0))-0.2;
TaxisMax = max(max(T0))+0.2;
TaxisMax_switch = max(max(T(:,:,round(Nt*time_redux)+1)));
font = 15;
fps = 60;
figure()
for t = 1:round((Nt*time_redux))
    mesh(X_mesh,Y_mesh,T(:,:,t));
%     if input.switch
%         axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
%     else
%         axis([0 Lx 0 Ly TaxisMin TaxisMax]);
%     end
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    colorbar
    pause(1/fps)
end

%% End behaviour
[X_mesh,Y_mesh] = ndgrid(X,Y);
time_redux = 0.2; % percentage of shown time instances
TaxisMin = min(min(T(:,:,end)))-min(min(T(:,:,end)))*0.025;
TaxisMax = max(max(T(:,:,end)))+max(max(T(:,:,end)))*0.025;
TaxisMin_switch = min(min(T(:,:,Nt-round((Nt*time_redux)))))-min(min(T(:,:,Nt-round((Nt*time_redux)))))*0.025;
font = 15;
fps = 60;
figure()
for t = Nt-round((Nt*time_redux)):Nt
    mesh(X_mesh,Y_mesh,T(:,:,t));
    if input.switch
        axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    else
        axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    end
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    pause(1/fps)
end

%% Sample code for fixed colorbar
[X_mesh,Y_mesh] = ndgrid(X,Y);
TaxisMin = min(min(T0));
TaxisMax = max(max(T0));
TaxisMax_switch = max(max(T(:,:,end)));
font = 15;
fps = 60;
figure()
for t = 1:Nt
    mesh(X_mesh,Y_mesh,T(:,:,t));  
    if input.switch
        axis([0 Lx 0 Ly TaxisMin TaxisMax_switch]);
    else
        axis([0 Lx 0 Ly TaxisMin TaxisMax]);
    end
    title(sprintf('Plate temperature for time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
    xlabel('x [m]',Interpreter='latex',FontSize=font); 
    ylabel('y [m]',Interpreter='latex',FontSize=font); 
    zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
    caxis manual;          % allow subsequent plots to use the same color limits
    caxis([TaxisMin TaxisMax]); 
    colorbar;
    pause(1/fps)
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
%% Solutions of ODE solver for a_r over time
figure()
for r= 1:diagn.R
   plot(squeeze(a(:,r)))
   hold on
end
xlabel('Time [t]')
ylabel('[-]')
title('Time solutions of a_r')
grid on

%% Plot initial conditions a_kl
figure()
k = 0:K;
l = 0:L;
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

figure()
[X_mesh,Y_mesh] = ndgrid(X,Y);
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
surf(X,Y,T0')
xlabel('x')
ylabel('y')
zlabel('Temperature [-]')
title('Intial temperature')

%% Show partial derivatives of the intital temperature 
figure()
subplot(1,2,1)
surf(X,Y,T0dx')
xlabel('x')
ylabel('y')
title('dT0/dx')
subplot(1,2,2)
surf(X,Y,T0dy')
xlabel('x')
ylabel('y')
title('dT0/dy')

%% Runtime simulation visual
psi = [1,2,5,7,10,20]; % Increasing order of K and L
t = [1,2.11,9.98,18.35,34.6,142]; % Increasing runtime
p = polyfit(psi,t,5);
figure()
plot(psi,t)
psi_des = 1:100;
f = polyval(p,psi_des);
figure()
plot(psi_des,f./(3600))
xlabel('order (K and L)')
ylabel('time [h]')

%% Show singular values of POD basis
figure()
plot(diagn.svp,'ko','Linewidth',2)
figure()
plot(real(diagn.U(:,1:3)))
title('The first 3 dominant modes')

%% Check if the POD basis gradients are actually what you want
figure()
surf(phiPOD.xy(:,:,1))
figure()
surf(phiPOD.ddx(:,:,1))
figure()
surf(phiPOD.ddy(:,:,1))
figure()
% test = sum(phiPOD(:,:,1)'*phiPOD(:,:,2),'all')*xstep*ystep
