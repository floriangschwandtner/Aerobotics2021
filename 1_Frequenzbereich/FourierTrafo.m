function [x_Fourier, u_Fourier, G, f] = FourierTrafo(x, u, t)

tdiff = zeros(length(t)-1,1);
for i=1:length(t)-1
    tdiff(i)=t(i+1,1)-t(i,1);
end

T = mean(tdiff);            % Sampling period 
Fs = 1/T;                   % Sampling frequency                          
L = length(t);              % Length of signal
f = Fs*(0:(L/2))/L;

% Fourier-Trafo Zustand
Y=fft(x);
P2 = Y/L;
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
x_Fourier = Y(1:L/2+1,:);

% Fourier-Trafo Steuergrößen
Y=fft(u);
P2_u = Y/L;
P1_u = P2_u(1:L/2+1,:);
P1_u(2:end-1,:) = 2*P1_u(2:end-1,:);
u_Fourier = Y(1:L/2+1,:);

% Berechnung der Übertragungsmatrix für jede Frequenz
G = zeros(4,2,L/2+1);
for k=1:L/2+1
    G(:,:,k) = x_Fourier(k,:)'/u_Fourier(k,:)';
end

% for k=1:L/2+1
%     for i=1:4
%         for j=1:2
%             G(i,j,k) = x_Fourier(k,i)/u_Fourier(k,j);
%         end
%     end
% end

end