close all
clear
load ('result_202107202005.mat');

theta= [];

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

%% Trimmwerte
%a0 = 0.0061;
V0 = 26.9962;
eta0 = -0.1326;
deltaF0 = 0.4244;

x(:,1) = x(:,1)-a0;
x(:,3) = x(:,3)-V0;
u(:,1) = u(:,1)-eta0;
u(:,2) = u(:,2)-deltaF0;

%% 
x = x(1:7500,:);
u = u(1:7500,:);
t = t(1:7500);
%% Simulation
x0 = [a0 0 V0 0];
U1sim = [t, u(:,1)];
U2sim = [t, u(:,2)];
output = sim('ZRD_laengs_Simulator.slx',t);

%% Visualization
figure(1)
subplot(2,2,1)
plot(output.tout, output.ysim(:,1))
hold on
plot(t,x(:,1))
xlabel('t[s]')
ylabel('\Delta \alpha')

subplot(2,2,2)
plot(output.tout, output.ysim(:,2))
hold on
plot(t,x(:,2))
xlabel('t[s]')
ylabel('q')


subplot(2,2,3)
plot(output.tout, output.ysim(:,3))
hold on
plot(t,x(:,3))
xlabel('t[s]')
ylabel('\Delta V_A')


subplot(2,2,4)
plot(output.tout, output.ysim(:,4))
hold on
plot(t,x(:,4))
xlabel('t[s]')
ylabel('\Delta \gamma')