%% Parameterschätzer für Maja Modell
%2021-07-04     Florian Gschwandtner Filterung+Schätzung

clc
clear
close all

%% Daten einlesen

load("data_even.mat")

t0 = 0;
tend = 595;
[t0err, t0ind] = min(abs(t-t0));
[tenderr, tendind] = min(abs(t-tend));

fEck = 10; %Hz
dt = t(2,1)-t(1,1)
[x1filt, x1_dot] = filteredDerivative(x(t0ind:tendind,1), dt, fEck);
[x2filt, x2_dot] = filteredDerivative(x(t0ind:tendind,2), dt, fEck);
[x3filt, x3_dot] = filteredDerivative(x(t0ind:tendind,3), dt, fEck);
[x4filt, x4_dot] = filteredDerivative(x(t0ind:tendind,4), dt, fEck);

[u1filt, u1_dot] = filteredDerivative(u(t0ind:tendind,1), dt, fEck);
[u2filt, u2_dot] = filteredDerivative(u(t0ind:tendind,2), dt, fEck);


figure
plot(t(t0ind:tendind-1), x1filt)
hold on
plot(t(t0ind:tendind), x(t0ind:tendind,1))

%% Schätzer

z = [x1_dot, x2_dot, x3_dot, x4_dot];

H = [x1filt, x2filt, x3filt, x4filt, u1filt, u2filt];

xhat = (H'*H)\H'*z

A = xhat(1:4, 1:4)'
B = xhat(5:6, 1:4)'

save('ABMatrix.mat', 'A', 'B');