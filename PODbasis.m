function [phiPOD,diagn] = PODbasis(T,Nx,Ny,xstep,ystep,Nt,reduced_energy_remaint)
%PODBASIS Summary of this function goes here
T_stacked = reshape(T,[Nx*Ny,Nt]);
[U,S,V] = svd(T_stacked,'econ');
svp = diag(S)/sum(diag(S)); % singular value percentiles

indeces = find(svp<=1-reduced_energy_remaint);
R = indeces(1);
fprintf('POD basis is of order R =%d \n',R);

phiPOD.xy = zeros(Nx,Ny,R);
phiPOD.ddx = zeros(Nx,Ny,R);
phiPOD.ddy = zeros(Nx,Ny,R);
phiPOD.xy = U(:,1:R);


phiPOD.xy = reshape(phiPOD.xy,[Nx,Ny,R]);

% G = eye(Nx)*1/sqrt((xstep*ystep)); % This is G^(-1/2) THIS IS ASSUMED!
% for r = 1:R
% phiPOD.xy(:,:,r) = G * phiPOD.xy(:,:,r);
% end

for r = 1:R
phiPOD.xy(:,:,r) = (phiPOD.xy(:,:,r)).*sqrt(1/(xstep*ystep));
end


for r = 1:R
[X, Y] = gradient(phiPOD.xy(:,:,r),xstep,ystep);
[phiPOD.ddx(:,:,r), phiPOD.dxdy(:,:,r)] = gradient(X,xstep,ystep);
[phiPOD.dydx(:,:,r), phiPOD.ddy(:,:,r)] = gradient(Y,xstep,ystep);
end

phiPOD.grad = phiPOD.ddx+phiPOD.ddy; %phiPOD.dxdy+phiPOD.dydx;
phiPOD.dotp = zeros(R,R);

% for i = 1 :R
%     for j = 1:R
%     phiPOD.dotp(i,j) = (sum(phiPOD.xy(:,:,i).*phiPOD.ddx(:,:,j),'all')+ ...
%     sum(phiPOD.xy(:,:,i).*phiPOD.ddy(:,:,j),'all'))*xstep*ystep;
%     end
% end

% Redefinition of the gradient trial
for i = 1 :R
    for j = 1:R
    phiPOD.dotp(i,j) = sum(phiPOD.xy(:,:,i).*phiPOD.grad(:,:,j),'all')*xstep*ystep;
    end
end

% phiPOD.ddx = reshape(phiPOD.ddx,[Nx,Ny,R]);
% phiPOD.ddy = reshape(phiPOD.ddy,[Nx,Ny,R]);

diagn.R = R; % Diagnostics
diagn.U = U;
diagn.S = S;
diagn.V = V;
diagn.svp = svp;
end




