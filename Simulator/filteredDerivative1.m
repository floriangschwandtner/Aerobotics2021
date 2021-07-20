function [out, diffOut] = filteredDerivative1(in,dt,fEck,dFilt)
%FILTEREDDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

%% Filter Erstellung
s = tf([1 0],1);
% fEck = 4; % Filter Eckfrequenz in Hz
omegaFilt = 2*pi*fEck;
%dFilt =200;% 1/sqrt(2); % Filter Dämpfung (1/sqrt(2)= kritische Dämpfung) 
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
diffOut = -(y(1:end-1)-y(2:end))./dt;
%%
% out=smoothdata(in,'movmean',100);
% diffOut(1)=(out(2)-out(1))/(dt);
% diffOut(end)=(out(end)-out(end-1))/(dt);
% for i=2:length(out)-1
%     diffOut(i) = (out(i+1)-out(i-1))/(2*dt);
% end
end