%% Funktion um X und U aus den Messdaten zu bestimmen
%2021-06-19 Florian Gschwandtner: Funktion erstellt
%2021-06-24 FT: Interpolation für x_4 hinzugefügt
%2021-06-24 FT: Interpolation für alle Variablen hinzugefügt, sodass
%Abtastrate gleichmäßig


load("data.mat")

% Bestimmung der mittleren Abtastrate
tdiff = zeros(length(t)-1,1);
for i=1:length(t)-1
    tdiff(i)=t(i+1)-t(i);
end
dt = mean(tdiff);

t_vec = t(1):dt:t(end);

x_1Inter = interp1(t, x(:,1),(t_vec))';
x_2Inter = interp1(t, x(:,2),(t_vec))';
x_3Inter = interp1(t, x(:,3),(t_vec))';
x_4Inter = interp1(t, x(:,4),(t_vec))';
u_1Inter = interp1(t, u(:,1),(t_vec))';
u_2Inter = interp1(t, u(:,2),(t_vec))';

X = [x_1Inter, x_2Inter, x_3Inter, x_4Inter];
U = [u_1Inter, u_2Inter];