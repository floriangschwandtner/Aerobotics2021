function [dJ] = compute_dJdtheta(N,f,G_k,dGF,Svv,Suu,Szu)
dJ = zeros(length(dGF),1);
inv_Svv = inv(Svv);
for i=1:length(dGF)
    for k = 1:N
       dJ(i) = dJ(i) + trace(conj(dGF{i}(1j*f(k)))'*inv_Svv*(Szu(:,:,k)-G_k(:,:,k)*Suu(:,:,k))); 
    end
end
dJ = -2 * real(dJ);
end