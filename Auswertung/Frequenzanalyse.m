clc
clear
close all

filename = "data_11_14_29";
load(filename+'.mat');
load('data_11_14_29_model_5b_2_648_668.mat');

% Vorbereiten
yo = -(data.vehicle_attitude_0.phi);
state = -(data.vehicle_attitude_0.phi);
roll_stick = data.ad_conv_heli_ifr_0.d_alpha;
stepsize = 0.02;
simtime = linspace(0,data.vehicle_attitude_0.timeSec(end)-data.vehicle_attitude_0.timeSec(1), length(state));

% KP = s.KP;
% T_L = s.T_L;
% T_l = s.T_l;
% tau_e = s.tau_e;
% 
% T_p = tf([KP*T_L KP],[T_l 1], 'InputDelay', tau_e);
% outputBestFit = lsim(T_p, state, simtime);

Fs = 50;                 % Sampling frequency                    
T = 1/Fs;                % Sampling period       
L1 = length(yo);          % Length of signal
L2 = length(roll_stick);
t = (0:L-1)*T;           % Time vector

% Zustand
YO = fft(yo);                       % FFT
P2O = abs(YO/L);                    % Zweiseitiges Leistungsdichtespektrum
P1O = P2O(1:L/2+1);                 % Einseitig
P1O(2:end-1) = 2*P1O(2:end-1);
fO = Fs*(0:(L/2))/L;                % Frequenzen

% Input
YI = fft(roll_stick);
P2I = abs(YI/L);
P1I = P2I(1:L/2+1);
P1I(2:end-1) = 2*P1I(2:end-1);
fI = Fs*(0:(L/2))/L;

% % Modell
% YM = fft(outputBestFit)
% P2M = abs(YM/L);
% P1M = P2M(1:L/2+1);
% P1M(2:end-1) = 2*P1M(2:end-1);
% fM = Fs*(0:(L/2))/L;

% Plot

figure
hold on
plot(fO,P1O)
plot(fI,P1I)
%plot(fM,P1M)
title('Einseitiges Leistungsdichtespektrum')
legend('Zustand', 'Stickwinkel Gemessen', 'Stickwinkel Simulation');
xlabel('f (Hz)')
ylabel('|P1(f)|')

print('Leistungsdichtespektrum, '-depsc')