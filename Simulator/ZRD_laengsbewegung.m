% 2021-06-19 Florian Gschwandtner: Grundstruktur, Modell

clc
%clear all
close all

V0 = 30;

% Beiwerte für A
Ma = -0.3215;
Mq = 3.3984;
Mv = 5.2744;
Xa = -3.6481;
Xv = 6.7118;
g = 9.81;
Za = 9.9196;
Zv = -6.3420;

%Beiwerte für B
Xdf = 2.8712;
a0 = 0;
i_f = 0;
Meta = -0.0967;
Mdf = -1.9727;
Xeta = -0.5580;
Zeta = 0.8461;

%Erster längerer Geradeausflug
t0 = 906;
tend = 1119;

% Erstellen der Matrixen für Zustand und Eingang
%[X, U, t_vec] = createStateAndInput(t0, tend);

%t = t_vec-t_vec(1);
X = x;
U = u;
t_vec = t;
U1sim = [t, u(:,1)];
U2sim = [t, u(:,2)];

output = sim('ZRD_laengs_Simulator.slx', t_vec);

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


