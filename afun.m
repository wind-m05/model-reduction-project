function dxdt = afun(t,a,kappa,rho,c,Lx,Ly,k,l,phi_u1,phi_u2,input)
    source = 0;
        if input.switch  
        [u1, u2] = heatInput(t,input.par); 
        source = u1*phi_u1(k+1,l+1)+u2*phi_u2(k+1,l+1);
        end
    dxdt = -kappa/(rho*c)*(k^2*pi^2/Lx^2+l^2*pi^2/Ly^2)*a + 1/(rho*c)*source; 
end