function [M] = compute_M(N,f,dGF,Svv,Suu)
M = zeros(length(dGF),length(dGF)); 
for a=1:length(dGF)
    for b=1:length(dGF)
        for k = 1:N
            M(a,b) = M(a,b) + trace(conj(dGF{a}(1j*f(k)))'* inv(Svv)* dGF{b}(1j*f(k))*Suu(k));
        end
        M(a,b) = 2*real(M(a,b));
    end
end
end