function dxdt = afun_POD(t,a,kappa,rho,c,Lx,Ly,xstep,ystep,r,phiPOD,input)
    source = 0;
        if input.switch  
        [u1, u2] = heatInput(t,input.par); 
        source = sum(u1*input.u1.*phiPOD.xy(:,:,r)+u2*input.u2.*phiPOD.xy(:,:,r),'all')*xstep*ystep; % !!! I want to solve this for every r instantly
        end
    
     dxdt = -kappa/(rho*c)*phiPOD.dotp*a;  % 1/(rho*c)*source  ; % This should be a column vector of 9 if R = 3
end
