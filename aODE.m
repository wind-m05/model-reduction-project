function a = aODE(t,a0,kappa,rho,c,Lx,Ly,xstep,ystep,k,l,phi_kl,input)
    [t,a] = ode45(@afun,t,a0,[],kappa,rho,c,Lx,Ly,xstep,ystep,k,l,phi_kl,input);
end
