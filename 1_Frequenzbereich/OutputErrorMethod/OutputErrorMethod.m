%% Pfad

cd .. 
addpath(cd)
cd OutputErrorMethod


clear
load ('data_even.mat');
load ('Matrices.mat');
x=x(1:49254,:);
u=u(1:49254,:);
t=t(1:49254,:);

%% 
threshold = 1;
iter_max = 15;

V0 = 26.992;
alpha0 = 0;
eta0 = -0.1326;
deltaF0 = 0.4244;
g = 9.81;

%% 
x(:,1) = x(:,1)-alpha0;
x(:,3) = x(:,3)-V0;
u(:,1) = u(:,1)-eta0;
u(:,2) = u(:,2)-deltaF0;
%% Fourier-Trafos
t_end = 7500;
f_end = 200;

[x_Fourier, u_Fourier, G_exp, f] = FourierTrafo(x(1:t_end,:), u(1:t_end,:), t(1:t_end));
x_Fourier=x_Fourier(2:f_end,:);
u_Fourier=u_Fourier(2:f_end,:);
G_exp=G_exp(:,:,2:f_end);
f=f(2:f_end);
N = length(f);

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
theta = [Z_alpha Z_V M_alpha M_q M_V X_alpha X_V Z_eta X_deltaF M_eta M_deltaF X_eta];


%% Newton-Raphson-Algorithmus
iter = 1;
dtheta = 0;
J = zeros(iter_max,1);
J(1) = 5;
while iter <= iter_max && abs(J(iter)) > threshold
    % Update Parametervektor
    theta  = theta + dtheta';
    
    % Kostenfunktion
    
    y = zeros(4,N);
    Svv = zeros(4,4);
    Suu = zeros(2,2,N);
    Szu = zeros(4,2,N);
    
    %Besetzung von G und deren Ableitung mit aktuellem theta
    
    G_sub    = subs(G,theta_sym(1:14),[theta V0 g]);
    dG1_sub   = subs(dG1,theta_sym(1:14),[theta V0 g]);
    dG2_sub   = subs(dG2,theta_sym(1:14),[theta V0 g]);
    dG3_sub   = subs(dG3,theta_sym(1:14),[theta V0 g]);
    dG4_sub   = subs(dG4,theta_sym(1:14),[theta V0 g]);
    dG5_sub   = subs(dG5,theta_sym(1:14),[theta V0 g]);
    dG6_sub   = subs(dG6,theta_sym(1:14),[theta V0 g]);
    dG7_sub   = subs(dG7,theta_sym(1:14),[theta V0 g]);
    dG8_sub   = subs(dG8,theta_sym(1:14),[theta V0 g]);
    dG9_sub   = subs(dG9,theta_sym(1:14),[theta V0 g]);
    dG10_sub  = subs(dG10,theta_sym(1:14),[theta V0 g]);
    dG11_sub  = subs(dG11,theta_sym(1:14),[theta V0 g]);
    dG12_sub  = subs(dG12,theta_sym(1:14),[theta V0 g]);
    
    %Umwandeln in Frequenzabhängige Funktion
    
    GF     = matlabFunction(G_sub);
    dGF{1} = matlabFunction(dG1_sub);
    dGF{2} = matlabFunction(dG2_sub);
    dGF{3} = matlabFunction(dG3_sub);
    dGF{4} = matlabFunction(dG4_sub);
    dGF{5} = matlabFunction(dG5_sub);
    dGF{6} = matlabFunction(dG6_sub);
    dGF{7} = matlabFunction(dG7_sub);
    dGF{8} = matlabFunction(dG8_sub);
    dGF{9} = matlabFunction(dG9_sub);
    dGF{10}= matlabFunction(dG10_sub);
    dGF{11}= matlabFunction(dG11_sub);
    dGF{12}= matlabFunction(dG12_sub);
    
    G_k   = compute_G(N, f, GF);
    
    for k = 1:N
       y(:,k)     = G_k(:,:,k)*u_Fourier(k,:)';
       %Svv        = Svv + (x_Fourier(k,:)'-y(:,k)) * conj((x_Fourier(k,:)'-y(:,k))');
       Svv = eye(4,4);
       Suu(:,:,k) = u_Fourier(k,:)'*conj(u_Fourier(k,:)')';
       Szu(:,:,k) = x_Fourier(k,:)'*conj(u_Fourier(k,:)')';
    end
    
    J(iter+1) = N*log(norm(Svv));
    inv_Svv = inv(Svv);
    for k = 1:N
       J(iter+1) = J(iter+1) + N*( conj((x_Fourier(k,:)'-y(:,k))')*inv_Svv*(x_Fourier(k,:)'-y(:,k)) ); %#ok<MINV>
    end
    
    
    % Update des Parametervektors
    dJ     = compute_dJdtheta(N,f,G_k,dGF,Svv,Suu,Szu);
    M      = compute_M(N,f,dGF,Svv,Suu);
    dtheta = - inv(M)*dJ;
    iter   = iter + 1;
end
