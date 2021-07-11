%Initialisierung der Testdaten
clear all
close all

time   = linspace(0,500,5001); %Zeit in 
time   = time';
u_test = zeros(length(time),2);
x_test = zeros(length(time),4);

%gewünschte Frequenzen
f_1 = 0.5;
f_2 = 1;
f_3 = 2;

%u_1 Höhenruder

eta0 = -0.1;
u_test(:,1) = eta0;
u_test(1:601,1) = 0.05 * sin(time(1:601)*f_1*2*pi)+eta0;
u_test(702:1301,1) = 0.05 *sin(time(702:1301)*f_2*2*pi)+eta0;
u_test(1402:2001,1) = 0.05 * sin(time(1402:2001)*f_3*2*pi)+eta0;

%u_2 Schub
delta0 = 0.4;
u_test(:,2) = delta0;
u_test(2100:2701,2) = 0.1 * sin(time(2100:2701)*f_1*2*pi)+delta0;
u_test(2802:3401,2) = 0.1 *sin(time(2702:3301)*f_2*2*pi)+delta0;
u_test(3502:4101,2) = 0.1 * sin(time(3402:4001)*f_3*2*pi)+delta0;

%x_1 Anstellwinkel alpha

alpha0 = 0;
x_test(:,1) = alpha0;
x_test(1:601,1) = 0.03 * sin((time(1:601)+0.5*(1/f_1*ones(601,1)))*f_1*2*pi)+alpha0;
x_test(702:1301,1) = 0.02 *sin((time(702:1301)+0.5*(1/f_2*ones(600,1)))*f_2*2*pi)+alpha0;
x_test(1402:2001,1) = 0.0001 * sin((time(1402:2001)+0.5*(1/f_3*ones(600,1)))*f_3*2*pi)+alpha0;
x_test(2100:2701,1) = 0.005 * sin(time(2100:2701)*f_1*2*pi)+alpha0;
x_test(2802:3401,1) = 0.0035 *sin(time(2802:3401)*f_2*2*pi)+alpha0;
x_test(3502:4101,1) = 0.00001 * sin(time(3502:4101)*f_3*2*pi)+alpha0;

%x_2 Drehrate

q0 = 0;
x_test(:,2) = q0;
x_test(1:601,2) = 0.4 * sin((time(1:601)+0.5*(1/f_1*ones(601,1)))*f_1*2*pi)+q0;
x_test(702:1301,2) = 0.3 *sin((time(702:1301)+0.5*(1/f_2*ones(600,1)))*f_2*2*pi)+q0;
x_test(1402:2001,2) = 0.01 * sin((time(1402:2001)+0.5*(1/f_3*ones(600,1)))*f_3*2*pi)+q0;
x_test(2100:2701,2) = 0.06 * sin(time(2100:2701)*f_1*2*pi)+q0;
x_test(2802:3401,2) = 0.0035 *sin(time(2802:3401)*f_2*2*pi)+q0;
x_test(3502:4101,2) = 0.00001 * sin(time(3502:4101)*f_3*2*pi)+q0;

%x_3 Geschwindigkeit VA

V0 = 27;
x_test(:,3) = V0;
x_test(1:601,3) = 5 * sin((time(1:601)+0.5*(1/f_1*ones(601,1)))*f_1*2*pi)+V0;
x_test(702:1301,3) = 3 *sin((time(702:1301)+0.5*(f_2*ones(600,1)))*f_2*2*pi)+V0;
x_test(1402:2001,3) = 0.1 * sin((time(1402:2001)+0.5*(1/f_3*ones(600,1)))*f_3*2*pi)+V0;
x_test(2100:2701,3) = 6 * sin(time(2100:2701)*f_1*2*pi)+V0;
x_test(2802:3401,3) = 3 *sin(time(2802:3401)*f_2*2*pi)+V0;
x_test(3502:4101,3) = 0.2 * sin(time(3502:4101)*f_3*2*pi)+V0;

%X_4 Bahnwinkel

gamma0 = 0;
x_test(:,4) = gamma0;
x_test(1:601,4) = 0.03 * sin((time(1:601)+0.5*(1/f_1*ones(601,1)))*f_1*2*pi)+gamma0;
x_test(702:1301,4) = 0.02 *sin((time(702:1301)+0.5*(1/f_2*ones(600,1)))*f_2*2*pi)+gamma0;
x_test(1402:2001,4) = 0.0001 * sin((time(1402:2001)+0.5*(1/f_3*ones(600,1)))*f_3*2*pi)+gamma0;
x_test(2100:2701,4) = 0.005 * sin(time(2100:2701)*f_1*2*pi)+gamma0;
x_test(2802:3401,4) = 0.0035 *sin(time(2802:3401)*f_2*2*pi)+gamma0;
x_test(3502:4101,4) = 0.00001 * sin(time(3502:4101)*f_3*2*pi)+gamma0;

[x_Fourier_test, u_Fourier_test, G, f] = FourierTrafo(x_test, u_test, time(1:5000));

figure
plot(time,u_test(:,1));
hold on
plot(time,u_test(:,2));
xlabel('time [s]');
ylabel('u');
figure 
plot(f,abs(u_Fourier_test(:,1)));
hold on
plot(f,abs(u_Fourier_test(:,2)));
xlim([0.1 4]);
figure
plot(time,x_test(:,1));
hold on
plot(time,x_test(:,2));
hold on
plot(time,x_test(:,3));
hold on
plot(time,x_test(:,4));
xlabel('time [s]');
ylabel('x');
figure 
plot(f,abs(x_Fourier_test(:,1)));
xlim([0.1 4]);
save('testdata.mat','u_test','x_test','time');