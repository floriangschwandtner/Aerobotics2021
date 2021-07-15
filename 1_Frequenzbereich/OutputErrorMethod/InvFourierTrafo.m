function [x] = InvFourierTrafo(x_Fourier, L)

Y_x = x_Fourier;
Y_x(2:end-1,:) = 0.5*(Y_x(2:end-1,:));
Y_x = L*[Y_x(1:end,:);conj(Y_x(end-1:-1:2,:))];
 
x = ifft(Y_x);

end