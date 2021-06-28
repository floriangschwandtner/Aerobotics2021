function res = compute_residual(z, u, theta, N)
    for i = 1:N
       res(i) = z(i) - G(:,:,i)*u(i,:)'; 
    end
end