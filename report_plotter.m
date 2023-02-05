[X_mesh,Y_mesh] = ndgrid(X,Y);
xindex = 1:length(X);
yindex = 1:length(Y);
font = 18;
TaxisMin = min(min(min(T)))+Tamb;
TaxisMax = max(max(max(T)))+Tamb;

s = get(0, 'ScreenSize');
figure('Position', [10 s(4)/4 1500 400]);

tindex = 1;
TaxisMin = min(min(min(T(:,:,tindex))))+Tamb-0.1;
TaxisMax = max(max(max(T(:,:,tindex))))+Tamb+0.1;
subplot(1,3,1);
Z = T(xindex,yindex,tindex)+Tamb;
s=mesh(X_mesh,Y_mesh,Z); 
axis([0 Lx 0 Ly TaxisMin TaxisMax]);
title(sprintf('Plate temperature at time t = %g[s]', round(time(tindex))),Interpreter='latex',FontSize=font);
xlabel('x [m]',Interpreter='latex',FontSize=font); 
ylabel('y [m]',Interpreter='latex',FontSize=font); 
zlabel('T(x,y,t) [K]',Interpreter='latex',FontSize=font);
s.FaceColor = 'flat';
colorbar
caxis([TaxisMin TaxisMax]);

tindex = 20; %length(t)/10;
TaxisMin = min(min(min(T(:,:,tindex))))+Tamb-0.1;
TaxisMax = max(max(max(T(:,:,tindex))))+Tamb+0.1;
subplot(1,3,2);
Z = T(xindex,yindex,tindex)+Tamb;
s=mesh(X_mesh,Y_mesh,Z); 
axis([0 Lx 0 Ly TaxisMin TaxisMax]);
title(sprintf('Plate temperature at time t = %g[s]', round(time(tindex))),Interpreter='latex',FontSize=font);
xlabel('x [m]',Interpreter='latex',FontSize=font); 
ylabel('y [m]',Interpreter='latex',FontSize=font); 
zlabel('T(x,y,t) [K]',Interpreter='latex',FontSize=font);
s.FaceColor = 'flat';
colorbar
caxis([TaxisMin TaxisMax]);

tindex = length(time);
subplot(1,3,3);
TaxisMin = min(min(min(T(:,:,tindex))))+Tamb-0.1;
TaxisMax = max(max(max(T(:,:,tindex))))+Tamb+0.1;
Z = T(xindex,yindex,tindex)+Tamb;
s=mesh(X_mesh,Y_mesh,Z); 
axis([0 Lx 0 Ly TaxisMin TaxisMax]);
title(sprintf('Plate temperature at time t = %g[s]', round(time(tindex))),Interpreter='latex',FontSize=font);
xlabel('x [m]',Interpreter='latex',FontSize=font); 
ylabel('y [m]',Interpreter='latex',FontSize=font); 
zlabel('T(x,y,t) [K]',Interpreter='latex',FontSize=font);
s.FaceColor = 'flat';
colorbar
caxis([TaxisMin TaxisMax]);

set(gcf,'Renderer','Painter')
hgexport(gcf,'figure.eps');







% %% Plot the three most dominant basis functions
% 
% [X_mesh,Y_mesh] = ndgrid(X,Y);
% xindex = 1:length(X);
% yindex = 1:length(Y);
% font = 18;
% TaxisMin = min(min(min(phiPOD.xy(:,:,1))));
% TaxisMax = max(max(max(phiPOD.xy(:,:,1))));
% 
% s = get(0, 'ScreenSize');
% figure('Position', [10 s(4)/4 1500 400]);
% 
% % tindex = 1;
% % TaxisMin = min(min(min(T(:,:,tindex))))+Tamb-0.1;
% % TaxisMax = max(max(max(T(:,:,tindex))))+Tamb+0.1;
% subplot(1,3,1);
% Z = phiPOD.xy(:,:,1);
% s=mesh(X_mesh,Y_mesh,Z); 
% % axis([0 Lx 0 Ly TaxisMin TaxisMax]);
% title('1st dominant basis function',Interpreter='latex',FontSize=font);
% xlabel('x [m]',Interpreter='latex',FontSize=font); 
% ylabel('y [m]',Interpreter='latex',FontSize=font); 
% zlabel('Amplitude [-]',Interpreter='latex',FontSize=font);
% s.FaceColor = 'flat';
% colorbar
% caxis([TaxisMin TaxisMax]);
% 
% % tindex = 20; %length(t)/10;
% TaxisMin = min(min(min(phiPOD.xy(:,:,2))));
% TaxisMax = max(max(max(phiPOD.xy(:,:,2))));
% subplot(1,3,2);
% Z = phiPOD.xy(:,:,2);
% s=mesh(X_mesh,Y_mesh,Z); 
% % axis([0 Lx 0 Ly TaxisMin TaxisMax]);
% title('2nd dominant basis function',Interpreter='latex',FontSize=font);
% xlabel('x [m]',Interpreter='latex',FontSize=font); 
% ylabel('y [m]',Interpreter='latex',FontSize=font); 
% zlabel('Amplitude [-]',Interpreter='latex',FontSize=font);
% s.FaceColor = 'flat';
% colorbar
% caxis([TaxisMin TaxisMax]);
% 
% subplot(1,3,3);
% TaxisMin = min(min(min(phiPOD.xy(:,:,3))));
% TaxisMax = max(max(max(phiPOD.xy(:,:,3))));
% Z = phiPOD.xy(:,:,3);
% s=mesh(X_mesh,Y_mesh,Z); 
% % axis([0 Lx 0 Ly TaxisMin TaxisMax]);
% title('3rd dominant basis function',Interpreter='latex',FontSize=font);
% xlabel('x [m]',Interpreter='latex',FontSize=font); 
% ylabel('y [m]',Interpreter='latex',FontSize=font); 
% zlabel('Amplitude [-]',Interpreter='latex',FontSize=font);
% s.FaceColor = 'flat';
% colorbar
% caxis([TaxisMin TaxisMax]);
% 
% 
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'figure.eps');

%% plot for e(k,l)
font = 18;
[X_mesh,Y_mesh] = ndgrid(X,Y);
sumE = 0;
for k = 0:K
    for l = 0:L
    sumE = sumE + e(k+1,l+1)*phi_kl(:,:,k+1,l+1); 
    end
end
Z = sumE;

figure()
mesh(X_mesh,Y_mesh,Z)
title('Smooth approximation $\rho c(x,y)$',Interpreter='latex',FontSize=font);
xlabel('x $[m]$',Interpreter='latex',FontSize=font); 
ylabel('y $[m]$',Interpreter='latex',FontSize=font); 
zlabel('e(x,y) $[J/m^3K]$',Interpreter='latex',FontSize=font)
% hold on
% surf(M)

set(gcf,'Renderer','Painter')
hgexport(gcf,'figure.eps');


