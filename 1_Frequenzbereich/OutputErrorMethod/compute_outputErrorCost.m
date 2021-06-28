function J = compute_outputErrorCost(z, u, theta, N)
    for k = 1:N
       y = compute_G(theta)
       Svv = Svv + (z(k)-y)*(z(k)-y);
    end
    
    J = N*log(abs(Svv));
    for i = 1:N
       J = J + ; 
    end
    
    function y = y(k,theta)
        y = G*u
    end
end