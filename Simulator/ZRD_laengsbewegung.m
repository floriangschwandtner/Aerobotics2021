% 2021-06-19 Florian Gschwandtner: Grundstruktur, Modell

clc
clear all
close all

V0 = 30;

% Beiwerte für A
Ma = 1;
Mq = 1;
Mv = 1;
Xa = 1;
Xv = 1;
g = 9.81;
Za = 1;
Zv = 1;

%Beiwerte für B
Xdf = 1;
a0 = 0;
i_f = 0;
Meta = 1;
Mdf = 1;
Xeta = 1;
Zeta = 1;

%Erster längerer Geradeausflug
t0 = 906;
tend = 1119;

% Erstellen der Matrixen für Zustand und Eingang
[X, U, t_vec] = createStateAndInput(t0, tend);

t = t_vec-t_vec(1);

U1sim = [t, U(:,1)]
U2sim = [t, U(:,2)]

output = sim('ZRD_laengs_Simulator.slx', t_vec-t_vec(1))

figure(1)
subplot(2,2,1)
plot(output.tout, output.ysim(:,1))
hold on
plot(t,X(:,1))
xlabel('t[s]')
ylabel('\Delta \alpha')

subplot(2,2,2)
plot(output.tout, output.ysim(:,2))
hold on
plot(t,X(:,2))
xlabel('t[s]')
ylabel('q')


subplot(2,2,3)
plot(output.tout, output.ysim(:,3))
hold on
plot(t,X(:,3))
xlabel('t[s]')
ylabel('\Delta V_A')


subplot(2,2,4)
plot(output.tout, output.ysim(:,4))
hold on
plot(t,X(:,4))
xlabel('t[s]')
ylabel('\Delta \gamma')


