function a = aODE_POD(t,a0,kappa,rho,c,Lx,Ly,xstep,ystep,R,phiPOD,input)
    [t,a] = ode45(@afun_POD,t,a0,[],kappa,rho,c,Lx,Ly,xstep,ystep,R,phiPOD,input);
end
