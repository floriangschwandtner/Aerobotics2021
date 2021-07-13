function [xhat]= LSQ(X,U,dt,V0,a0, i_f,g ,C )
%% Estimation Problem
% Längs-& Seitenbewegung 
% z = H*x + e
%dt = (t(length(t))-t(1))/length(t) ;
switch C 

    case 'L'
alpha = X(:,1);
q     = X(:,2);
V     = X(:,3);
gamma = X(:,4);
nu    = U(:,1);
deltaF= U(:,2);
N     = length(alpha);
%dt    = (tend-t0)/N;
jj    = (a0+i_f)*pi/180; 

z=reshape([((alpha(2:N)-alpha(1:N-1))./(dt))- q(2:N);
     (q(2:N) - q(1:N-1))./(dt) ;
         ((V(2:N)-V(1:N-1))./(dt))- g*gamma(2:N)]',[(N-1)*3,1]);
     
H1 = [zeros(N-1,3), -sin(jj)*deltaF(2:N), alpha(2:length(alpha))/V0, V(2:length(V))/V0, nu(2:length(nu))/V0,zeros(N-1,5)];
H2 = [zeros(N-1,7), alpha(2:length(alpha)), q(2:length(q)), V(2:length(V)),nu(2:length(nu)),deltaF(2:length(deltaF))];
H3 = [alpha(2:length(alpha)), V(2:length(V)),nu(2:length(nu)), cos(jj)*deltaF(2:length(deltaF)), zeros(N-1,8)];    
H  = zeros (3*(N-1),12);

s=1;
for i=1: N-1  
H(s:3*i,:)=[H1(i,:);H2(i,:);H3(i,:)];
s=s+3;
end

    case 'S'
r    = X(:,1);
beta = X(:,2);
p    = X(:,3);
phi  = X(:,4);
xii  = U(:,1);
xio  = U(:,2);
ceta = U(:,3);
N    = length(r);
%dt    = (tend-t0)/N;

 z=reshape([((r(2:N)-r(1:N-1))./(dt));
     (beta(2:N) - beta(1:N-1))./(dt) + r(2:N) - g*phi(2:N)/V0 ;
         (p(2:N)-p(1:N-1))./(dt)]',[(N-1)*3,1]); 
H1 = [zeros(N-1,7), r(2:N), beta(2:N),p(2:N),xii(2:N),xio(2:N),ceta(2:N)];   
H2 = [(beta(2:N)+ ceta(2:N))/V0, zeros(N-1,12)];
H3 = [zeros(N-1,1),  r(2:N), beta(2:N),p(2:N),xii(2:N),xio(2:N),ceta(2:N),zeros(N-1,6)];
H  = zeros(3*(N-1),13);
s=1;
for i=1: N-1  
H(s:3*i,:)=[H1(i,:);H2(i,:);H3(i,:)];
s=s+3;
end

end


xhat = (pinv(H)*z);
end
