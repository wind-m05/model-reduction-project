function phi_kl = basisxy(x,y,k,l,Lx,Ly)
    if((k == 0)&&(l == 0))
        phi_kl = 1/sqrt(Lx*Ly);
    end
    if((k > 0)&&(l == 0))
        phi_kl = sqrt(2/(Lx*Ly))*cos(k*pi*x/Lx);
    end
    if((k == 0)&&(l > 0))
        phi_kl = sqrt(2/(Lx*Ly))*cos(l*pi*y/Ly);
    end
    if((k > 0)&&(l > 0))
        phi_kl = sqrt(4/(Lx*Ly))*cos(k*pi*x/Lx)*cos(l*pi*y/Ly);    
    end
end

