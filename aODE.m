function a = aODE(t,a0,kappa,rho,c,Lx,Ly,k,l,phi_u1,phi_u2,input)
    [t,a] = ode45(@afun,t,a0,[],kappa,rho,c,Lx,Ly,k,l,phi_u1,phi_u2,input);
end

