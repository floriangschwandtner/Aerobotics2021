%% Funktion um X und U aus den Messdaten zu bestimmen
%2021-06-19 Florian Gschwandtner: Funktion erstellt
%2021-06-24 FT: Interpolation für x_4 hinzugefügt

function [X, U, t_vec] = createStateAndInput(t0, tend)
load("x_state_original.mat")
load("u_control_original.mat")

[minValue0 , i0_x2] = min(abs(x_2(:,1)-t0));
[minValueEnd, iend_x2] = min(abs(x_2(:,1)-tend));
[minValue0 , i0_x4] = min(abs(x_4(:,1)-t0));
[minValueEnd, iend_x4] = min(abs(x_4(:,1)-tend));

if(minValue0>1.0 || minValueEnd>1.0)
    error("t0 or tend not in x_2");
end


%Frequenz x1,x3,x4 gleich
[minValue0 , i0_x] = min(abs(x_1(:,1)-t0));
[minValueEnd, iend_x] = min(abs(x_1(:,1)-tend));

if(minValue0>1.0 || minValueEnd>1.0)
    error("t0 or tend not in x_1");
end


% Input
[minValue0 , i0_u] = min(abs(u_1(:,1)-t0));
[minValueEnd, iend_u] = min(abs(u_1(:,1)-tend));

if(minValue0>1.0 || minValueEnd>1.0)
    error("t0 or tend not in x_1");
end

t_vec = x_1(i0_x:iend_x,1);

x_2Inter = interp1(x_2(i0_x2:iend_x2,1), x_2(i0_x2:iend_x2,2),(t_vec));
x_4Inter = interp1(x_4(i0_x4:iend_x4,1), x_4(i0_x4:iend_x4,2),(t_vec));

X = [x_1(i0_x:iend_x, 2), x_2Inter, x_3(i0_x:iend_x, 2), x_4Inter];
% dirty FIX
X(1,2)=X(2,2);
X(1,4)=X(2,4);

u_1Inter = interp1(u_1(i0_u:iend_u,1), u_1(i0_u:iend_u,2),(t_vec));
u_2Inter = interp1(u_2(i0_u:iend_u,1), u_2(i0_u:iend_u,2),(t_vec));

U = [u_1Inter, u_2Inter];

end