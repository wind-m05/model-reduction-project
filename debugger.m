% Check the source over time
source = 0;
% User parameters
show_visuals = false; % Will show all the visualization plots subsequently
input.switch = true; % Turn the input source on or off
input.par.type = 'sinusoid'; % {const,sinusoid} What type of input
input.par.freq = 0.1; % [Hz]
input.par.tstart = 5; % [s]
input.par.tend = tend; % [s]
input.par.amp1 = 1; % [-]
input.par.amp2 = 1;
source_dotproduct_profile_matrix = cell(K,L,Nt)
for k = 1:K 
    for l = 1:L
        for t = 1:length(time)
            if input.switch  
            [u1, u2] = heatInput(time(t),input.par); 
            source = sum(u1*input.u1.*phi_kl(:,:,k,l)+u2*input.u2.*phi_kl(:,:,k,l),'all')*xstep*ystep;
            source_matrix(k,l,t) = source;
            source_dotproduct_profile_matrix{k,l,t} = u1*input.u1.*phi_kl(:,:,k,l)+u2*input.u2.*phi_kl(:,:,k,l);
            end
        end
    end
end
   
% figure()
% for k = 1:K 
%     for l = 1:L
%     plot(time,squeeze(source_matrix(k,l,:)))   
%     hold on
%     end
% end
% 
% figure()
% plot(time,squeeze(source_matrix(5,5,:)))   
figure()
fps = 60;
font = 15;
tcounter=0
changes = K*L;
for t = 1:length(time)

for k = 1:K
    for l = 1:L
surf(source_dotproduct_profile_matrix{k,l,t})
title(sprintf('Input projected onto basis in time = %g [s]', round(time(t))),Interpreter='latex',FontSize=font);
xlabel('x [m]',Interpreter='latex',FontSize=font); 
ylabel('y [m]',Interpreter='latex',FontSize=font); 
zlabel('T(x,y,t) [$^\circ \mathrm{C}]$',Interpreter='latex',FontSize=font);
pause(1/fps)
tcounter = tcounter+1;
    end
end
end

%% Test input 2D in time 

        u1 = abs(cos(2*pi*par.freq*t));
        u2 = abs(sin(2*pi*par.freq*t));
