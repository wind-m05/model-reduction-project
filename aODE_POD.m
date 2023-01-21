function a = aODE(t,a0,kappa,rho,c,Lx,Ly,R,phiPOD,phi_u1,phi_u2,input)
    [t,a] = ode45(@afun_POD,t,a0,[],kappa,rho,c,Lx,Ly,R,phiPOD,phi_u1,phi_u2,input);
end
