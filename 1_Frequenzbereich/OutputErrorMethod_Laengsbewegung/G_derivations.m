% smbolische Inverse von A, Übertragungsmatrix G und deren Ableitungen

theta_sym = sym('a',[1 16]);
% Z_alpha Z_V M_alpha M_q M_V X_alpha X_V Z_eta X_deltaF M_eta M_deltaF X_eta V0  g   alpha0  
% a1      a2  a3      a4  a5  a6      a7  a8    a9       a10   a11      a12   a13 a14 a15   
% f
% a16

A = [theta_sym(1)/theta_sym(13) 1 theta_sym(2)/theta_sym(13) 0;...
        theta_sym(3) theta_sym(4) theta_sym(5) 0;...
        theta_sym(6) 0 theta_sym(7) -theta_sym(14);...
        -theta_sym(1)/theta_sym(13) 0 -theta_sym(2)/theta_sym(13) 0];
    
B = [theta_sym(8)/theta_sym(13) -theta_sym(9)/theta_sym(13)*sin(theta_sym(15));...
        theta_sym(10) theta_sym(11);...
        theta_sym(12) theta_sym(9)*cos(theta_sym(15));...
        -theta_sym(8)/theta_sym(13) theta_sym(9)/theta_sym(13)*sin(theta_sym(15))];
    
G = inv(1j*2*pi*theta_sym(16)*eye(size(A))-A)*B;

dG1  = diff(G,theta_sym(1));
dG2  = diff(G,theta_sym(2));
dG3  = diff(G,theta_sym(3));
dG4  = diff(G,theta_sym(4));
dG5  = diff(G,theta_sym(5));
dG6  = diff(G,theta_sym(6));
dG7  = diff(G,theta_sym(7));
dG8  = diff(G,theta_sym(8));
dG9  = diff(G,theta_sym(9));
dG10 = diff(G,theta_sym(10));
dG11 = diff(G,theta_sym(11));
dG12 = diff(G,theta_sym(12));

dG1_conj  = diff(conj(G)',theta_sym(1));
dG2_conj  = diff(conj(G)',theta_sym(2));
dG3_conj  = diff(conj(G)',theta_sym(3));
dG4_conj  = diff(conj(G)',theta_sym(4));
dG5_conj  = diff(conj(G)',theta_sym(5));
dG6_conj  = diff(conj(G)',theta_sym(6));
dG7_conj  = diff(conj(G)',theta_sym(7));
dG8_conj  = diff(conj(G)',theta_sym(8));
dG9_conj  = diff(conj(G)',theta_sym(9));
dG10_conj = diff(conj(G)',theta_sym(10));
dG11_conj = diff(conj(G)',theta_sym(11));
dG12_conj = diff(conj(G)',theta_sym(12));

filename = 'Matrices.mat';
save Matrices.mat