function [x] = InvFourierTrafo(x_Fourier, u_Fourier, f, L)

x_Fourier = 0.5*L*[flip(x_Fourier(2:end,:)); x_Fourier(1:end,:)];

x = real(ifft(x_Fourier));


end