function [u1,u2] = heatInput(time,par)
for t = 1 : length(time)
    if(time(t)>=par.tstart && time(t)<=par.tend)
        if strcmp(par.type, 'const')
        u1(t) = par.amp1;
        u2(t) = par.amp2;
        
        elseif strcmp(par.type, 'sine')
        u1(t) = par.amp1*abs(cos(2*pi*par.freq*t));
        u2(t) = par.amp2*abs(sin(2*pi*par.freq*t));
        end
    else
        u1(t) = 0;
        u2(t) = 0;
    end
 
end

% u1 = input.par.amp1*cos(2*pi*input.par.freq*time);
% u2 = input.par.amp2*sin(2*pi*input.par.freq*time);
end