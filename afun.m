function dxdt = afun(t,a,kappa,rho,c,Lx,Ly,k,l)
    dxdt = -(kappa/(rho*c))*((k^2*pi^2/Lx^2)+(l^2*pi^2/Ly^2))*a;   %+1/(rho*c)*innerProduct(u,phi)
end