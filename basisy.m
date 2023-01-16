function phiy = basisy(y,l,Ly)
    if l == 0
        phiy = 1/sqrt(Ly);
    elseif l>0
        phiy = sqrt(2/Ly)*cos(l*pi*y/Ly);
    end
end

