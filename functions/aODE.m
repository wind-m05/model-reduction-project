function a = aODE(t,a0,kappa,rho,c,Lx,Ly,k,l)
    [t,a] = ode45(@afun,t,a0,[],kappa,rho,c,Lx,Ly,k,l);
end

