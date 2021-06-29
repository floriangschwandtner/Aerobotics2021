%% Pfad
cd ..
cd OutputErrorMethod

clear
load ('data_even.mat');
load ('Matrices.mat');
x=x(1:49254,:);
u=u(1:49254,:);
t=t(1:49254,:);

%% 
threshold = 1e-3;
iter_max = 1;

V0 = 25;
alpha0 = 0;
g = 9.81;

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
theta = [Z_alpha Z_V M_alpha M_q M_V X_alpha X_V Z_eta X_deltaF M_eta M_deltaF X_eta];


%% Newton-Raphson-Algorithmus
iter = 0;
while iter <= iter_max
    % Kostenfunktion
    
    y = zeros(4,N);
    Svv = zeros(4,4);
    Suu = zeros(length(u_Fourier),1);
    Szu = zeros(4,2,N);
    
    %Besetzung von G und deren Ableitung mit aktuellem theta
    
    G    = subs(G,theta_sym(1:14),[theta V0 g]);
    dG1  = subs(dG1,theta_sym(1:14),[theta V0 g]);
    dG2  = subs(dG2,theta_sym(1:14),[theta V0 g]);
    dG3  = subs(dG3,theta_sym(1:14),[theta V0 g]);
    dG4  = subs(dG4,theta_sym(1:14),[theta V0 g]);
    dG5  = subs(dG5,theta_sym(1:14),[theta V0 g]);
    dG6  = subs(dG6,theta_sym(1:14),[theta V0 g]);
    dG7  = subs(dG7,theta_sym(1:14),[theta V0 g]);
    dG8  = subs(dG8,theta_sym(1:14),[theta V0 g]);
    dG9  = subs(dG9,theta_sym(1:14),[theta V0 g]);
    dG10 = subs(dG10,theta_sym(1:14),[theta V0 g]);
    dG11 = subs(dG11,theta_sym(1:14),[theta V0 g]);
    dG12 = subs(dG12,theta_sym(1:14),[theta V0 g]);
    
    %Umwandeln in Frequenzabhängige Funktion
    
    GF     = matlabFunction(G);
    dGF{1} = matlabFunction(dG1);
    dGF{2} = matlabFunction(dG2);
    dGF{3} = matlabFunction(dG3);
    dGF{4} = matlabFunction(dG4);
    dGF{5} = matlabFunction(dG5);
    dGF{6} = matlabFunction(dG6);
    dGF{7} = matlabFunction(dG7);
    dGF{8} = matlabFunction(dG8);
    dGF{9} = matlabFunction(dG9);
    dGF{10}= matlabFunction(dG10);
    dGF{11}= matlabFunction(dG11);
    dGF{12}= matlabFunction(dG12);
    
    G_k   = compute_G(N, f, GF);
    
    for k = 1:N
       y(:,k)     = G_k(:,:,k)*u_Fourier(k,:)';
       Svv        = Svv + (x_Fourier(k,:)'-y(:,k))*conj((x_Fourier(k,:)'-y(:,k))');
       Suu(k)     = u_Fourier(k,:)*conj(u_Fourier(k,:))';
       Szu(:,:,k) = x_Fourier(k,:)'*conj(u_Fourier(k,:));
    end
    
    J = N*log(abs(Svv));
    for k = 1:N
       J = J + N*( conj((x(k,:)'-y(:,k))')/Svv*(x(k,:)'-y(:,k)) );
    end
    
    
    % Update des Parametervektors
    dJ     = compute_dJdtheta(N,f,G_k,dGF,Svv,Suu,Szu);
    M      = compute_M(N,f,dGF,Svv,Suu);
    dtheta = - inv(M)*dJ;
    theta  = theta + dtheta;
    iter   = iter + 1;
end
