clear
load('data_even.mat')
x=x(1:49254,:);
u=u(1:49254,:);
t=t(1:49254);
T = t(2)-t(1);

[x_Fourier, u_Fourier, G, f, Y_x, Y_u] = FourierTrafo(x, u, t);

k_max = length(f);

G_xx=zeros(2,2,k_max);
G_yy=zeros(4,4,k_max);
G_xy=zeros(2,4,k_max);
H_xy=zeros(2,4,k_max);

for k=1:k_max
    G_xx(:,:,k)=2/T*(conj(u_Fourier(k,:)')*u_Fourier(k,:));
    G_yy(:,:,k)=2/T*(conj(x_Fourier(k,:)')*x_Fourier(k,:));
    G_xy(:,:,k)=2/T*(conj(u_Fourier(k,:)')*x_Fourier(k,:));
    H_xy(:,:,k)=inv(G_xx(:,:,k))*G_xy(:,:,k);
    gamma()=det(G_xy)/
end