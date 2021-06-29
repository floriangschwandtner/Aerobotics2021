function [G_k] = compute_G(N, f, GF)  
 G_k = zeros(4,2,N);
     for k = 2:N
        G_k(:,:,k) = GF(1j*f(k)); 
     end
end
% 
% function [G] = compute_G(theta, V0, alpha0, N, f)
%     [A, B] = compute_AB(theta, V0, alpha0);
%     
%     G = zeros(4,2,N);
%     for k = 2:N
%         G(:,:,k) = (1j*f(k)-A)\B;
%     end
% end