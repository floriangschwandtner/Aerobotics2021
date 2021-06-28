%% Pfad
cd ..
addpath(cd)
cd OutputErrorMethod

clear
load ('data_even.mat');
x=x(1:49254,:);
u=u(1:49254,:);
t=t(1:49254,:);

%% 
threshold = 1e-3;
iter_max = 100;

V0 = 25;
alpha0 = 0;

%% Fourier-Trafos
[x_Fourier, u_Fourier, G_exp, f] = FourierTrafo(x, u, t);
N = length(x_Fourier);


%% Initialisierung des Parametervektors
Z_alpha = 1;
Z_V = 1; 
M_alpha = 1;
M_q = 1;
M_V = 1;
X_alpha = 1;
X_V = 1;
Z_eta = 1;
X_deltaF = 1;
M_eta = 1;
M_deltaF = 1;
X_eta = 1;
theta = [Z_alpha; Z_V; M_alpha; M_q; M_V; X_alpha; X_V; Z_eta; X_deltaF; M_eta; M_deltaF; X_eta];

%% Newton-Raphson-Algorithmus
iter = 0;
while iter <= iter_max
    % Kostenfunktion
    y = zeros(4,N);
    Svv = 0;
    G = compute_G(theta, V0, alpha0, N, f);
    for k = 1:N
       y(:,k) = G(:,:,k)*u(k,:)';
       Svv = Svv + (x(k,:)'-y(:,k))*conj((x(k,:)'-y(:,k))');
    end
    
    J = N*log(abs(Svv));
    for k = 1:N
       J = J + N*( conj((x(k,:)'-y(:,k))')/Svv*(x(k,:)'-y(:,k)) ); 
    end
    
    % Gradient
    Suu = zeros(2,2,N);
    
    
    % Update des Parametervektors
    dtheta = -inv(M)*compute_dJdtheta(theta);
    theta = theta + dtheta;
    iter = iter + 1;
end
