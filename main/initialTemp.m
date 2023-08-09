function [T0,TX,TY] = initialTemp(X,Y,k,l,selector,normalize)
%initialTemp Calculates a 2D discrete gaussian for an intitial condition of
% [T0] = initialTemp(X,Y,Lx,Ly,k,l,selector,normalize)
% X,Y = discretized spatial domain
% k,l = fourier number for x,y direction fourier basis respectively
% selector = {'gauss','fourier','blockup','blockdown'}
% Normalize = true/false
% Output: T0     = is temperature of a Gaussian shape at every x,y
% element
stdx = 0.03;
stdy=stdx;
pitsize = 1;

T0 = zeros(length(X),length(Y));
if strcmp(selector, 'gauss')
center_index_x = length(X)/2;
center_index_y = length(Y)/2;
meanx = round(center_index_x);
meany = round(center_index_y);
T0x = 1/(stdx*sqrt(2*pi))*(exp((-(X-X(meanx)).^2)/(2*stdx^2)));
T0x = T0x/max(T0x);
T0y = 1/(stdy*sqrt(2*pi))*(exp((-(Y-Y(meany)).^2)/(2*stdy^2)));
T0y = T0y/max(T0y);
T0 = T0x'*T0y;
elseif strcmp(selector, 'fourier')
    if((k == 0)&&(l == 0))
        T0  = 1/sqrt(X(end)*Y(end))*X'*Y;
    end
    if((k > 0)&&(l == 0))
        T0  = sqrt(2/(X(end)*Y(end)))*cos(k*pi*X'/X(end))*Y;
    end
    if((k == 0)&&(l > 0))
        T0  = sqrt(2/(X(end)*Y(end)))*X'*cos(l*pi*Y/Y(end));
    end
    if((k > 0)&&(l > 0))
        T0 = sqrt(4/(X(end)*Y(end)))*cos(k*pi*X'/X(end))*cos(l*pi*Y/Y(end));    
    end
elseif strcmp(selector, 'blockup')
T0((round(length(X)/4): round(3*length(X)/4)), (round(length(Y)/4):round(3*length(Y)/4))) = pitsize;
elseif strcmp(selector, 'blockdown')
T0((round(length(X)/4): round(3*length(X)/4)), (round(length(Y)/4):round(3*length(Y)/4))) = -pitsize;
end
if normalize == true
T0 = T0./(max(abs(T0),[],'all'));
end
[TX,TY] = gradient(T0,Y,X);
end

