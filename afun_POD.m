function dxdt = afun_POD(t,a,kappa,rho,c,Lx,Ly,r,phiPOD,phi_u1,phi_u2,input)
    source = 0;
        if input.switch  
        [u1, u2] = heatInput(t,input.par); 
        source = u1*phi_u1(r)+u2*phi_u2(r);
        end
    dxdt = -kappa/(rho*c)*(r^2/(Lx^2+Ly^2))*a + source; 
end
