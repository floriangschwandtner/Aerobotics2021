%% 
t_start = 915;
t_end = 1115;

%%
[X, U, t_vec] = createStateAndInput(t_start, t_end);


tdiff = zeros(length(t_vec)-1,1);
for i=1:length(t_vec)-1
    tdiff(i)=t_vec(i+1,1)-t_vec(i,1);
end

Fs = 1/mean(tdiff);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(t_vec);             % Length of signal
f = Fs*(0:(L/2))/L;

% Fourier-Trafo Zustand
Y=fft(X);
P2 = Y/L;
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
X_Fourier = Y(1:L/2+1,:);

% Fourier-Trafo Steuergrößen
Y=fft(U);
P2_U = Y/L;
P1_U = P2_U(1:L/2+1,:);
P1_U(2:end-1,:) = 2*P1_U(2:end-1,:);
U_Fourier = Y(1:L/2+1,:);

%% Berechnung der Übertragungsmatrix für jede Frequenz
G = zeros(4,2,L/2+1);
for k=1:L/2+1
    for i=1:4
        for j=1:2
            G(i,j,k) = X_Fourier(k,i)/U_Fourier(k,j);
        end
    end
end