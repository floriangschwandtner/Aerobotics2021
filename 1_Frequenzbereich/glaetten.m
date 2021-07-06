clear
close all
load('data_even.mat')

f_end = 200;

alpha0 = 0.0061;
V0 = 26.9962;
eta0=-0.1326;
deltaF0 = 0.4244;

x=x(1:7500,:);
u=u(1:7500,:);
t=t(1:7500);

x(:,1) = x(:,1)-alpha0;
x(:,3) = x(:,3)-V0;
u(:,1) = u(:,1)-eta0;
u(:,2) = u(:,2)-deltaF0;

[x_Fourier, u_Fourier, G_exp, f] = FourierTrafo(x, u, t);

x_Fourier=x_Fourier(2:f_end,:);
u_Fourier=u_Fourier(2:f_end,:);
G_exp=G_exp(:,:,2:f_end);
f=f(2:f_end);
N = length(f);

semilogx(f,20*log10(abs(x_Fourier(:,1))))
hold on
x_Fourier_smooth(:,1) = smoothdata(20*log10(abs(x_Fourier(:,1))),'sgolay',15);

semilogx(f,x_Fourier_smooth(:,1),'--')

x_Fourier_backtrafo = 10^(x_Fourier_smooth/20);