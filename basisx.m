function phix = basisx(x,k,Lx)
    if k == 0
        phix = 1/sqrt(Lx);
    elseif k>0
        phix = sqrt(2/Lx)*cos(k*pi*x/Lx);
    end
end

