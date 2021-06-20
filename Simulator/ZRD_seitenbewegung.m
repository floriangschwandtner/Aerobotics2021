% 2021-06-19 Florian Gschwandtner: Grundstruktur, Modell

clc
clear all
close all

V0 = 30;

% Beiwerte für A
Nr = 1;
Nb = 1;
Np = 1;
Yz = 1;
Lr = 1;
g = 9.81;
Lb = 1;
Lp = 1;

%Beiwerte für B
Nc = 1;
Nz = 0;
Yz = 0;
Le = 1;
Lz = 1;

%Erster längerer Geradeausflug
t0 = 906;
tend = 1119;

% Erstellen der Matrixen für Zustand und Eingang
[X, U, t_vec] = createStateAndInput(t0, tend);

t = t_vec-t_vec(1);

U1sim = [t, U(:,1)]
U2sim = [t, U(:,2)]

output = sim('ZRD_seitw_Simulator.slx', t_vec-t_vec(1))

figure(1)
subplot(2,2,1)
plot(output.tout, output.ysim(:,1))
hold on
plot(t,X(:,1))
xlabel('t[s]')
ylabel('r')

subplot(2,2,2)
plot(output.tout, output.ysim(:,2))
hold on
plot(t,X(:,2))
xlabel('t[s]')
ylabel('\beta')


subplot(2,2,3)
plot(output.tout, output.ysim(:,3))
hold on
plot(t,X(:,3))
xlabel('t[s]')
ylabel('p')


subplot(2,2,4)
plot(output.tout, output.ysim(:,4))
hold on
plot(t,X(:,4))
xlabel('t[s]')
ylabel('\Phi')


