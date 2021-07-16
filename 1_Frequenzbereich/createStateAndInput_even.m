%% Funktion um X und U aus den Messdaten zu bestimmen
%2021-06-19 Florian Gschwandtner: Funktion erstellt
%2021-06-24 FT: Interpolation für x_4 hinzugefügt
%2021-06-24 FT: Interpolation für alle Variablen hinzugefügt, sodass
%Abtastrate gleichmäßig

clear
load("x_state_original.mat")
load("u_control_original.mat")

% Bestimmung der mittleren Abtastrate
tdiff = zeros(length(x_1(:,1))-1,1);
for i=1:length(x_1(:,1))-1
    tdiff(i)=t(i+1)-t(i);
end
dt = mean(tdiff);

t_vec = t(1):dt:t(end);

x_1Inter = interp1(x_1(:,1), x_1(:,2),(t_vec));
x_2Inter = interp1(x_2(:,1), x_2(:,2),(t_vec));
x_3Inter = interp1(x_3(:,1), x_3(:,2),(t_vec));
x_4Inter = interp1(x_4(:,1), x_4(:,2),(t_vec));
u_1Inter = interp1(u_1(:,1), u_1(:,2),(t_vec));
u_2Inter = interp1(u_2(:,1), u_2(:,2),(t_vec));

X = [x_1Inter, x_2Inter, x_3Inter, x_4Inter];
U = [u_1Inter, u_2Inter];