% 2021-06-19 Florian Gschwandtner: Grundstruktur, Modell

clc
clear all
close all

%% Daten-Vorbereitung %%%%%
load("data_even.mat");

% Parameter 
n1 = 1;     %Startindex
n2 = 1200;  %Endindex
fEck =1;    %Hz
dt = 0.01;  % Zeitschritt
g    = 9.81;

% Anfangswerte
a0   = 0;   
i_f  = 0;
V0   = 26.92;
gamma0 = -0.019  ;
q0    = 0.083;

% Interpolation 

x1=interp1(t(n1:n2),x(n1:n2,1)-a0,t(n1):dt:t(n2));
x2=interp1(t(n1:n2),x(n1:n2,2),t(n1):dt:t(n2));
x3=interp1(t(n1:n2),x(n1:n2,3)-V0,t(n1):dt:t(n2));
x4=interp1(t(n1:n2),x(n1:n2,4)-gamma0,t(n1):dt:t(n2));
u1=interp1(t(n1:n2),u(n1:n2,1)+ 0.1326,t(n1):dt:t(n2));
u2=interp1(t(n1:n2),u(n1:n2,2)-0.4244,t(n1):dt:t(n2));

time = (t(n1):dt:t(n2))';
X    = [x1',x2',x3',x4'];
U    = [u1',u2'];

%% Filterung %%%%%%


[x1filt, x1_dot] = filteredDerivative(x1', dt, fEck);
[x2filt, x2_dot] = filteredDerivative(x2', dt, fEck);
[x3filt, x3_dot] = filteredDerivative(x3', dt, fEck);
[x4filt, x4_dot] = filteredDerivative(x4', dt, fEck);

[u1filt, u1_dot] = filteredDerivative(u1', dt, fEck);
[u2filt, u2_dot] = filteredDerivative(u2', dt, fEck);

Xfilt=[x1filt,x2filt,x3filt,x4filt];
Ufilt=[u1filt,u2filt];





%% LSQ- Schätzer 
% Beiwerte 
xhat = LSQ(Xfilt,Ufilt,dt,V0,a0, i_f,g ,'L');

% Beiwerte für A 
Ma = xhat(8);
Mq = xhat(9);
Mv = xhat(10);
Xa = xhat(1);
Xv = xhat(2);
Za = xhat(5);
Zv = xhat(6);

%Beiwerte für B
Xdf  = xhat(4);
Meta = xhat(11);
Mdf  = xhat(12);
Xeta = xhat(3);
Zeta = xhat(7);


%% Matrix-Schätzer
z = [x1_dot, x2_dot, x3_dot, x4_dot];

H = [x1filt, x2filt, x3filt, x4filt, u1filt, u2filt];

xhat1 = (H'*H)\H'*z

A  = [Za/V0, 1, Zv/V0, 0; Ma, Mq, Mv, 0; Xa, 0, Xv, -g; -Za/V0, 0, -Zv/V0, 0]
A1 = xhat1(1:4, 1:4)'
B  = [Zeta/V0, -Xdf/V0*sin(a0+i_f); Meta, Mdf; Xeta, Xdf*cos(a0+i_f); -Zeta/V0, Xdf/V0*sin(a0+i_f)]
B1 = xhat1(5:6, 1:4)'


%% Simulation %%%%%%
U1sim = [time, U(:,1)];
U2sim = [time, U(:,2)];
x0    = [0;0;0;0];
output = sim('ZRD_laengs_Simulator.slx', time);


%% LSQ Schätzer %%%%%%
figure(1)
subplot(2,2,1)

plot(output.tout, output.ysim(:,1))
hold on
plot(time(1:length(time)-1),x1filt)
%plot(time,X(:,1))
xlabel('t[s]')
ylabel('\Delta \alpha')

subplot(2,2,2)

plot(output.tout, output.ysim(:,2))
hold on
plot(time(1:length(time)-1),x2filt)
%plot(time,X(:,2))
xlabel('t[s]')
ylabel('q')


subplot(2,2,3)

plot(output.tout, output.ysim(:,3))
hold on
plot(time(1:length(time)-1),x3filt)
%plot(time,X(:,3))
xlabel('t[s]')
ylabel('\Delta V_A')


subplot(2,2,4)

plot(output.tout, output.ysim(:,4))
hold on
plot(time(1:length(time)-1),x4filt)
%plot(time,X(:,4))
xlabel('t[s]')
ylabel('\Delta \gamma')

legend('sim','true')


%% A&B Schätzer %%%%%%%%

figure(2)
subplot(2,2,1)

plot(output.tout, output.ysim1(:,1))
hold on
plot(time(1:length(time)-1),x1filt)
%plot(time,X(:,1))
xlabel('t[s]')
ylabel('\Delta \alpha')

subplot(2,2,2)

plot(output.tout, output.ysim1(:,2))
hold on
plot(time(1:length(time)-1),x2filt)
%plot(time,X(:,2))
xlabel('t[s]')
ylabel('q')


subplot(2,2,3)

plot(output.tout, output.ysim1(:,3))
hold on
plot(time(1:length(time)-1),x3filt)
%plot(time,X(:,3))
xlabel('t[s]')
ylabel('\Delta V_A')


subplot(2,2,4)

plot(output.tout, output.ysim1(:,4))
hold on
plot(time(1:length(time)-1),x4filt)
%plot(time,X(:,4))
xlabel('t[s]')
ylabel('\Delta \gamma')
legend('sim','true')
