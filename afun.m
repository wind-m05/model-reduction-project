function dxdt = afun(t,a,kappa,rho,c,Lx,Ly,xstep,ystep,k,l,phi_kl,input)
    source = 0;
        if input.switch  
        [u1, u2] = heatInput(t,input.par); 
        source = sum(u1*input.u1.*phi_kl(:,:,k,l)+u2*input.u2.*phi_kl(:,:,k,l),'all')*xstep*ystep;
        end
    dxdt = -kappa/(rho*c)*(k^2*pi^2/Lx^2+l^2*pi^2/Ly^2)*a + source; % source /rho c?
end