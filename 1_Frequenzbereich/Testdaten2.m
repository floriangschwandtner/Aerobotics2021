addpath('F:\Dokumente\Studium\Master\Semester 3\Aerobotics-Seminar\Aerobotics2021\Simulator')
addpath('F:\Dokumente\Studium\Master\Semester 3\Aerobotics-Seminar\Aerobotics2021\1_Frequenzbereich\OutputErrorMethod')
addpath('F:\Dokumente\Studium\Master\Semester 3\Aerobotics-Seminar\Aerobotics2021\1_Frequenzbereich')
load('Matrices.mat')

Z_alpha = -70;
Z_V = 50; 
M_alpha = 20;
M_q = 10;
M_V = 1;
X_alpha = 2;
X_V = 3;
Z_eta = 60;
X_deltaF = 7;
M_eta = 2;
M_deltaF = 2;
X_eta = 3;
g = 9.81;
theta = [Z_alpha Z_V M_alpha M_q M_V X_alpha X_V Z_eta X_deltaF M_eta M_deltaF X_eta];

t   = linspace(0,500,5000); %Zeit in 
t   = t';
u_test = zeros(length(t),2);
x_test = zeros(length(t),4);

%gewünschte Frequenzen
f_1 = 0.5;
f_2 = 1;
f_3 = 2;

%u_1 Höhenruder

eta0 = -0.1;
u_test(:,1) = eta0;
u_test(1:5000,1) = 0.05 * sin(t(1:5000)*f_1*2*pi)+eta0;
%u_test(602:1201,1) = 0.05 *sin(t(602:1201)*f_2*2*pi)+eta0;
%u_test(1202:1801,1) = 0.05 * sin(t(1202:1801)*f_3*2*pi)+eta0;

%u_2 Schub
delta0 = 0.4;
u_test(:,2) = delta0;
%u_test(2100:2701,2) = 0.1 * sin(t(2100:2701)*f_1*2*pi)+delta0;
% u_test(2802:3401,2) = 0.1 *sin(time(2702:3301)*f_2*2*pi)+delta0;
% u_test(3502:4101,2) = 0.1 * sin(time(3402:4001)*f_3*2*pi)+delta0;

tdiff = zeros(length(t)-1,1);
for i=1:length(t)-1
    tdiff(i)=t(i+1,1)-t(i,1);
end

T = mean(tdiff);            % Sampling period 
Fs = 1/T;                   % Sampling frequency                          
L = length(t);              % Length of signal
f = 2*Fs*(0:(L/2))/L;

Y_u=fft(u_test);
P2_u = Y_u/L;
P1_u = P2_u(1:L/2+1,:);
P1_u(2:end-1,:) = 2*P1_u(2:end-1,:);
u_Fourier = P1_u(1:L/2+1,:);

G_sub  = subs(G,theta_sym(1:15),[theta V0 g alpha0]);
GF     = matlabFunction(G_sub);

x_hat = zeros(4,length(t)/2+1);
for k=1:length(t)/2+1
    x_hat(:,k) = GF(f(k))*u_Fourier(k,:)';
end
x_hat = x_hat';

x = InvFourierTrafo(x_hat,length(t));