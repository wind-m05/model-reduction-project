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

G = eye(Nx)*sqrt(1/(xstep*ystep)); % This is G^(-1/2) THIS IS ASSUMED!
G = eye(Nx); % This is G^(-1/2) THIS IS ASSUMED!
phiPOD.xy = reshape(phiPOD.xy,[Nx,Ny,R]);

for r = 1:R
phiPOD.xy(:,:,r) = G * phiPOD.xy(:,:,r);
end

for r = 1:R
[X, Y] = gradient(phiPOD.xy(:,:,r),xstep,ystep);
[phiPOD.ddx(:,:,r), ~] = gradient(X,xstep,ystep);
[~, phiPOD.ddy(:,:,r)] = gradient(Y,xstep,ystep);
end

phiPOD.dotp = zeros(R,R);
for i = 1 :R
    for j = 1:R
    phiPOD.dotp(i,j) = (sum(phiPOD.xy(:,:,i).*phiPOD.ddx(:,:,j),'all')+ ...
    sum(phiPOD.xy(:,:,i).*phiPOD.ddy(:,:,j),'all'))*xstep*ystep;
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




