function [dJ] = compute_dJdtheta(N,f,G_k,dGF_conj,Svv,Suu,Szu)
dJ = zeros(length(dGF_conj),1);
inv_Svv = inv(Svv);
for i=1:length(dGF_conj)
    for k = 1:N
       dJ(i) = dJ(i) + 1/f(k)*trace((dGF_conj{i}(1j*f(k)))*inv_Svv*(Szu(:,:,k)-G_k(:,:,k)*Suu(:,:,k))); 
    end
end
dJ = -2 * real(dJ);
end