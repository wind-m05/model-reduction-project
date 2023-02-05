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

