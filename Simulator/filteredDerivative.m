function [out, diffOut] = filteredDerivative(in,dt,fEck)
%FILTEREDDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

%% Filter Erstellung
s = tf([1 0],1);
% fEck = 4; % Filter Eckfrequenz in Hz
omegaFilt = 2*pi*fEck;
dFilt = 1/sqrt(2); % Filter Dämpfung (1/sqrt(2)= kritische Dämpfung) 
filter = omegaFilt^2/(s^2 + 2* dFilt*omegaFilt *s + omegaFilt^2); % Tiefpass 2. Ordnung
filter = filter ^2; % Dann z.B. Teifpass 4. Ordnung
% figure;
% h = bodeplot(filter);
% setoptions(h,'FreqUnits','Hz');

%% Doppelte Anwendung des Filters auf Signal:
% matlab Variante als diskretes Filter:
% besser, weil initial conditions sinnvoll gewählt werden:
filterD = c2d(filter,dt); 
y = filtfilt(filterD.Numerator{:},filterD.Denominator{:},in);

out = y(1:end-1);
diffOut = (y(1:end-1)-y(2:end))./dt;
end