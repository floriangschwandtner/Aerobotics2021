function [M] = compute_M(N,f,dGF,Svv,Suu)
M = zeros(length(dGF),length(dGF)); 
inv_Svv = inv(Svv);
for ai=1:length(dGF)
    for aj=1:length(dGF)
        for k = 1:N
            M(ai,aj) = M(ai,aj) + trace(conj(dGF{aj}(1j*f(k)))'* inv_Svv* dGF{ai}(1j*f(k))*Suu(:,:,k));
        end
    end
end
M = 2*real(M);
end