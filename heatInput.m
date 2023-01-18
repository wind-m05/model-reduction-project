function [u1,u2] = heatInput(t,selector)
    u1 = 1;%abs(cos(0.2*pi*t));
    u2 = 2;%abs(sin(0.2*pi*t));
    if strcmp(selector, 'const')
        u1 = 1;
        u2 = 2;
    end
    if strcmp(selector, 'tempconst')
        if(t>15)
            u1 = 0;
            u2 = 0;
        elseif(t>5)
            u1 = 1;
            u2 = 1;
        else
            u1 = 0;
            u2 = 0;
        end
    end
    if strcmp(selector, 'sinusoid')
        u1 = abs(cos(0.1*pi*t));
        u2 = abs(sin(0.1*pi*t));
    end
end