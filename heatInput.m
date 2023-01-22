function [u1,u2] = heatInput(t,par)
    if(t>=par.tstart && t<=par.tend)

        if strcmp(par.type, 'const')
        u1 = par.amp1;
        u2 = par.amp2;
        
        elseif strcmp(par.type, 'sinusoid')
        u1 = par.amp1*cos(2*pi*par.freq*t);
        u2 = par.amp2*sin(2*pi*par.freq*t);
        end
    else
        u1 = 0;
        u2 = 0;
    end
end