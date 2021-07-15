clear 
close all
T= 0.11;
Fs = 1/T;                   % Sampling frequency                          
t = 0:T:5.01;
x=sin(2*pi*t);%+sin(4*pi*t+0.5);
L = length(t);              % Length of signal
f = 2*Fs*(0:(L/2))/L;

% Fourier-Trafo Zustand
Y_x=fft(x);
P2_x = Y_x/L;               % 2-sided spectrum
P1_x = P2_x(1:L/2+1);     % 1-sided spectrum
P1_x(2:end-1) = 2*P1_x(2:end-1);
x_Fourier = P1_x(1:L/2+1);

Y_x_back = x_Fourier;
Y_x_back(2:end-1) = 0.5*(Y_x_back(2:end-1));
Y_x_back_ = [Y_x_back(1:end),conj(Y_x_back(end-1:-1:2))];
Y_x_back = L*Y_x_back_;
 
x_back = ifft(Y_x_back);

figure
plot(t,x,'b')
title('Signal')
hold on
plot(t,x_back,'r--')
legend('Originalsignal', 'rücktransformiertes Signal')

figure
plot(f,abs(x_Fourier))
title('Fourier-Trafo (1-seitig)')
hold on
legend('Originalsignal', 'rücktransformiertes Signal')