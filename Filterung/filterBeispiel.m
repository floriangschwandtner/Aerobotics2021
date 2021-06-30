close all 
clear
%% Beispiel Anwendung Vorwärt-Rückwärts Filterung (zero-phase filtering)
% Universität Stuttgart
% iFR - Institut für Flugmechanik und Flugregelung
% Tobias Richter, Benjamin Rothaupt
% 07.08.2019
%% Beispielhaftes Zeitsignal erstellen
dt = 0.01;
t = 0:dt:5;
t = t(:);
omega = 2*pi* [0.5  2 6 20];
omega = omega(:);
phase = rand(size(omega))*2*pi*0;
a = [1 0.8 0.5 0.2]
a = a(:);
z = sum(ones(size(t))*a' .* sin(t*omega' + ones(size(t))*phase'),2);
%% Filter Erstellung
s = tf([1 0],1);
fEck = 4; % Filter Eckfrequenz in Hz
omegaFilt = 2*pi*fEck;
dFilt = 1/sqrt(2); % Filter Dämpfung (1/sqrt(2)= kritische Dämpfung) 
filter = omegaFilt^2/(s^2 + 2* dFilt*omegaFilt *s + omegaFilt^2); % Tiefpass 2. Ordnung
filter = filter ^2; % Dann z.B. Teifpass 4. Ordnung
figure;
h = bodeplot(filter);
setoptions(h,'FreqUnits','Hz');

%% Doppelte Anwendung des Filters auf Signal:
% selbst gebaut:
y1 = lsim(filter,z,t);
y1 = lsim(filter,y1(end:-1:1),t);
y1 = y1((end:-1:1));
% matlab Variante als diskretes Filter:
% besser, weil initial conditions sinnvoll gewählt werden:
filterD = c2d(filter,dt); 
y2 = filtfilt(filterD.Numerator{:},filterD.Denominator{:},z);

% plot data und Amplitudengang des doppel-Filters
figure;
h1 = bodeplot(filter^2);
setoptions(h1,'FreqUnits','Hz','PhaseVisible','off');
figure,
plot(t,[z,y1,y2]);
legend({'data','without initial conditions','matlab: filtfilt()'})