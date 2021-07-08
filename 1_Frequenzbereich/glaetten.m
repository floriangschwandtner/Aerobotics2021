load('data_even.mat')

x=x(1:7500,:);
u=u(1:7500,:);
t=t(1:7500);

x(:,1) = x(:,1)-alpha0;
x(:,3) = x(:,3)-V0;
u(:,1) = u(:,1)-eta0;
u(:,2) = u(:,2)-deltaF0;

[x_Fourier, u_Fourier, G_exp, f] = FourierTrafo(x, u, t);