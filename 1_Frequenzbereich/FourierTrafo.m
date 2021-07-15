function [x_Fourier, u_Fourier, G, f, Y_x, Y_u] = FourierTrafo(x, u, t, smooth)

tdiff = zeros(length(t)-1,1);
for i=1:length(t)-1
    tdiff(i)=t(i+1,1)-t(i,1);
end

T = mean(tdiff);            % Sampling period 
Fs = 1/T;                   % Sampling frequency                          
L = length(t);              % Length of signal
f = 2*Fs*(0:(L/2))/L;

% Fourier-Trafo Zustand
Y_x=fft(x);
P2_x = Y_x/L;               % 2-sided spectrum
P1_x = P2_x(1:L/2+1,:);     % 1-sided spectrum
P1_x(2:end-1,:) = 2*P1_x(2:end-1,:);
x_Fourier = P1_x(1:L/2+1,:);

% Fourier-Trafo Steuergrößen
Y_u=fft(u);
P2_u = Y_u/L;
P1_u = P2_u(1:L/2+1,:);
P1_u(2:end-1,:) = 2*P1_u(2:end-1,:);
u_Fourier = P1_u(1:L/2+1,:);

% Glättung
if smooth == 1
    x_Fourier = smoothdata(x_Fourier,'movmedian',15);
    u_Fourier = smoothdata(u_Fourier,'movmedian',15);
end
% Berechnung der Übertragungsmatrix für jede Frequenz
G = zeros(4,2,L/2+1);
% for k=1:L/2+1
%     G(:,:,k) = x_Fourier(k,:)'/u_Fourier(k,:)';
% end

for k=1:L/2+1
    for i=1:4
        for j=1:2
            G(i,j,k) = x_Fourier(k,i)/u_Fourier(k,j);
        end
    end
end

end