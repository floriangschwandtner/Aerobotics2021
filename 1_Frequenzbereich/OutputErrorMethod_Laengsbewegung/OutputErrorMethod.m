%% Pfad
cd .. 
addpath(cd)
cd OutputErrorMethod_Laengsbewegung
clear
close all
%load ('data_laengs_even.mat');
load ('Matrices.mat');
load ('testdata2.mat');
x=x_test(1:5000,:);
u=u_test(1:5000,:);
t=t(1:5000,:);

%% Datenbereich
t_start = 1;
t_end   = 550;
f_start = 2;
f_end   = 275;%500;
normierung = false;

%% Trimmzustand
g = 9.8067;
% TP1
alpha0 = 0;
V0 = 26.9996;
gamma0 = 0;
eta0 = -0.1;
deltaF0 = 0.4;
 
%% Algorithmus-Optionen
iter_max = 10;
nugget = 0.01;
threshold = 1e-5;

%% Berechnung der Delta-Werte und Normierung
x = x(t_start:t_end,:);
u = u(t_start:t_end,:);
t = t(t_start:t_end);

% Abweichungen vom Trimmpunkt
deltax(:,1) = x(:,1)-alpha0;
deltax(:,2) = x(:,2);
deltax(:,3) = x(:,3)-V0;
deltax(:,4) = x(:,4)-gamma0;
deltau(:,1) = u(:,1)-eta0;
deltau(:,2) = u(:,2)-deltaF0;
if (normierung)
    % Normierung
    norm_x = abs(max(deltax));
    deltax(:,1) = deltax(:,1)/norm_x(1);
    deltax(:,2) = deltax(:,2)/norm_x(2);
    deltax(:,3) = deltax(:,3)/norm_x(3);
    deltax(:,4) = deltax(:,4)/norm_x(4);
end
%% Fourier-Trafos
[x_Fourier_orig, u_Fourier_orig, G_exp_orig, f_orig] = FourierTrafo(deltax(:,:), deltau(:,:), t, 0);
x_Fourier=x_Fourier_orig(f_start:f_end,:);
u_Fourier=u_Fourier_orig(f_start:f_end,:);
G_exp=G_exp_orig(:,:,f_start:f_end);
f=f_orig(f_start:f_end);
N = length(f);

%% Initialisierung des Parametervektors
Z_alpha = [-0.473186580900319]*V0;
Z_V = [-0.000763325096957492]*V0; 
M_alpha = [-87.6450530495459];
M_q = [-58.2128494877069];
M_V = [1.66437872449425];
X_alpha = [-6.26938451560701];
X_V = [-0.0741111475602237];
Z_eta = [-0.0634198470314440]*V0;
X_deltaF = [-0.0226287671386684];%[-0.0226287671386684]*(-1)*V0/(sin(alpha0));
M_eta = [-97.6955223213428];
M_deltaF = [-53.7881728478177];
X_eta = [-0.351042475926444];
theta = zeros(iter_max,12);
theta(1,:) = [Z_alpha Z_V M_alpha M_q M_V X_alpha X_V Z_eta X_deltaF M_eta M_deltaF X_eta];

%% Newton-Raphson-Algorithmus
fprintf('==================================\n')
fprintf('Newton-Raphson-Algorithm started.\n')
fprintf('==================================\n')
dtheta = 100*ones(length(theta),1);
J = zeros(iter_max,1);
konv = zeros(iter_max,1);
M_pd = zeros(iter_max,1);       % M(iter) is pd -> M_pd(iter) = 1 (else 0)
iter = 1;
while iter <= iter_max && norm(dtheta)/norm(theta(iter,:)) > threshold
    % Update Parametervektor
    if iter > 1
        theta(iter,:)  = theta(iter-1,:) + dtheta';
    end
    
    % Kostenfunktion
    
    y = zeros(4,N);
    Svv = zeros(4,4);
    Suu = zeros(2,2,N);
    Szu = zeros(4,2,N);
    
    %Besetzung von G und deren Ableitung mit aktuellem theta
    
    G_sub    = subs(G,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG1_sub   = subs(dG1,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG2_sub   = subs(dG2,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG3_sub   = subs(dG3,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG4_sub   = subs(dG4,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG5_sub   = subs(dG5,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG6_sub   = subs(dG6,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG7_sub   = subs(dG7,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG8_sub   = subs(dG8,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG9_sub   = subs(dG9,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG10_sub  = subs(dG10,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG11_sub  = subs(dG11,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG12_sub  = subs(dG12,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    
    dG1_conj_sub   = subs(dG1_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG2_conj_sub   = subs(dG2_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG3_conj_sub   = subs(dG3_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG4_conj_sub   = subs(dG4_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG5_conj_sub   = subs(dG5_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG6_conj_sub   = subs(dG6_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG7_conj_sub   = subs(dG7_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG8_conj_sub   = subs(dG8_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG9_conj_sub   = subs(dG9_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG10_conj_sub  = subs(dG10_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG11_conj_sub  = subs(dG11_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    dG12_conj_sub  = subs(dG12_conj,theta_sym(1:15),[theta(iter,:) V0 g alpha0]);
    
    %Umwandeln in frequenzabhängige Funktion

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

    dGF_conj{1} = matlabFunction(dG1_conj_sub);
    dGF_conj{2} = matlabFunction(dG2_conj_sub);
    dGF_conj{3} = matlabFunction(dG3_conj_sub);
    dGF_conj{4} = matlabFunction(dG4_conj_sub);
    dGF_conj{5} = matlabFunction(dG5_conj_sub);
    dGF_conj{6} = matlabFunction(dG6_conj_sub);
    dGF_conj{7} = matlabFunction(dG7_conj_sub);
    dGF_conj{8} = matlabFunction(dG8_conj_sub);
    dGF_conj{9} = matlabFunction(dG9_conj_sub);
    dGF_conj{10}= matlabFunction(dG10_conj_sub);
    dGF_conj{11}= matlabFunction(dG11_conj_sub);
    dGF_conj{12}= matlabFunction(dG12_conj_sub);
    
    G_k   = compute_G(N, f, GF);
    
    for k = 1:N
       y(:,k)     = G_k(:,:,k)*u_Fourier(k,:)';
       Svv        = Svv + (x_Fourier(k,:)'-y(:,k)) * conj((x_Fourier(k,:)'-y(:,k))');
       %Svv = eye(4,4);
       Suu(:,:,k) = u_Fourier(k,:)'*conj(u_Fourier(k,:)')';
       Szu(:,:,k) = x_Fourier(k,:)'*conj(u_Fourier(k,:)')';
    end
    
    J(iter) = 1/N*log(norm(Svv));
    inv_Svv = inv(Svv);
    for k = 1:N
       J(iter) = J(iter) + 1/(N*f(k))*( conj((x_Fourier(k,:)'-y(:,k))')*inv_Svv*(x_Fourier(k,:)'-y(:,k)) ); %#ok<MINV>
    end
    
   
    % Update des Parametervektors
    dJ     = compute_dJdtheta(N,f,G_k,dGF_conj,Svv,Suu,Szu);
    M      = compute_M(N,f,dGF,dGF_conj,Svv,Suu);
    M = M + nugget*eye(size(M));
%     lambda=eig(M); 
%     flag = 0;
%     for i = 1:rank(M)
%         if lambda(i)<=0
%             flag = 1;
%         end
%     end
%     if flag == 0   % flag==0 --> M ist positiv definit
%         dtheta = - inv(M)*dJ;
%         M_pd(iter) = 1;
%     else
%         dtheta = inv(M)*dJ;
%     end
    dtheta = - inv(M)*dJ;
    konv(iter) = norm(dtheta)/norm(theta(iter,:));
    
    fprintf('Iteration #%d done. | J = %.5f\n',iter, abs(J(iter)))
    
    iter   = iter + 1;
end

[minJ, minIndex] = min(abs(J));
theta_best = theta(minIndex,:);

G_sub  = subs(G,theta_sym(1:15),[theta_best V0 g alpha0]);
GF     = matlabFunction(G_sub);

fprintf('==================================\n')
fprintf('Algorithm finished.\n')
fprintf('==================================\n')

x_hat = zeros(4,length(t)/2+1);
G_k = zeros(4,2,length(t)/2+1);
for k=2:51%length(t)/2+1
    x_hat(:,k) = GF(f_orig(k))*u_Fourier_orig(k,:)';
    G_k(:,:,k) = GF(f_orig(k));
end
x_hat = x_hat';

%% Plotting
figure
subplot(2,2,1);
semilogx(f_orig,20*log10(abs(x_hat(:,1))))
hold on
semilogx(f_orig,20*log10(abs(x_Fourier_orig(:,1))),'--')
ylabel('alpha')
subplot(2,2,2);
semilogx(f_orig,20*log10(abs(x_hat(:,2))))
hold on
semilogx(f_orig,20*log10(abs(x_Fourier_orig(:,2))),'--')
ylabel('q')
subplot(2,2,3);
semilogx(f_orig,20*log10(abs(x_hat(:,3))))
hold on
semilogx(f_orig,20*log10(abs(x_Fourier_orig(:,3))),'--')
ylabel('V_A')
subplot(2,2,4);
semilogx(f_orig,20*log10(abs(x_hat(:,4))))
hold on
semilogx(f_orig,20*log10(abs(x_Fourier_orig(:,4))),'--')
ylabel('gamma')

figure
subplot(2,2,1)
semilogx(f_orig,20*log10(abs(squeeze(G_exp_orig(1,1,:)))))
hold on
semilogx(f_orig,20*log10(abs(squeeze(G_k(1,1,:)))))
title('G(1,1)');

subplot(2,2,2)
semilogx(f_orig,20*log10(abs(squeeze(G_exp_orig(1,2,:)))))
hold on
semilogx(f_orig,20*log10(abs(squeeze(G_k(1,2,:)))))
title('G(1,2)');

subplot(2,2,3)
semilogx(f_orig,20*log10(abs(squeeze(G_exp_orig(3,1,:)))))
hold on
semilogx(f_orig,20*log10(abs(squeeze(G_k(3,1,:)))))
title('G(3,1)');

subplot(2,2,4)
semilogx(f_orig,20*log10(abs(squeeze(G_exp_orig(3,2,:)))))
hold on
semilogx(f_orig,20*log10(abs(squeeze(G_k(3,2,:)))))
title('G(3,2)');

[deltax_time_hat] = InvFourierTrafo(x_hat,length(t));

if (normierung)
    % Rück-Normierung
    deltax(:,1) = deltax(:,1)/norm_x(1);
    deltax(:,2) = deltax(:,2)/norm_x(2);
    deltax(:,3) = deltax(:,3)/norm_x(3);
    deltax(:,4) = deltax(:,4)/norm_x(4);
end

x_time_hat(:,1) = deltax_time_hat(:,1)+alpha0;
x_time_hat(:,2) = deltax_time_hat(:,2);
x_time_hat(:,3) = deltax_time_hat(:,3)+V0;
x_time_hat(:,4) = deltax_time_hat(:,4)+gamma0;

titles = {'\Delta\alpha', 'q', '\Delta V_A', '\Delta \gamma'};
for i=1:4
    figure
    plot(t,x_time_hat(:,i),'--' )
    hold on
    plot(t,x(:,i),'r--' )
    plot(t,u(:,1),'g')
    plot(t,u(:,2),'y')
    legend(strcat(titles{i},'_{Näherung}'),titles{i},'\Delta\eta', '\Delta\delta_{F}')
end