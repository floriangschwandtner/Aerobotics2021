function [A,B] = compute_AB(theta, V0, alpha0)
% theta = [1      ; 2  ; 3      ; 4  ; 5  ; 6      ; 7  ; 8    ; 9       ; 10   ; 11      ; 12   ];
% theta = [Z_alpha; Z_V; M_alpha; M_q; M_V; X_alpha; X_V; Z_eta; X_deltaF; M_eta; M_deltaF; X_eta];
    
    iF = 0;
    g = 9.81;
    
    A = [theta(1)/V0 1 theta(2)/V0 0;...
        theta(3) theta(4) theta(5) 0;...
        theta(6) 0 theta(7) -g;...
        -theta(1)/V0 0 -theta(2)/V0 0];
    
    B = [theta(8)/V0 -theta(9)/V0*sin(alpha0+iF);...
        theta(10) theta(11);...
        theta(12) theta(9)*cos(alpha0+iF);...
        -theta(8)/V0 theta(9)/V0*sin(alpha0+iF)];
%     load 'AB.mat'
%     A=A1;
%     B=B1;
end