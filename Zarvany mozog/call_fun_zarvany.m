clearvars
close all
clc

pl = 0;
u0 = 1;
dx_vec = -0.008:0.0001:0.008; %[mm] pozicio x
dy_vec = -0.09:0.0001:0.09; %[mm] pozicio y

d_vec = dx_vec;
% d_vec = dy_vec;

R = zeros(length(d_vec));

for i = 1:length(d_vec)
    [R,P] = fun_zarvany(d_vec(i),u0, pl);              % x irányban
    % [R,P] = fun_zarvany_xy(0.07,d_vec(i),u0, pl);      % y irányban
    R_P(i,1) = R;
    R_P(i,2) = P;
end

%% Plot
figure()
plot(d_vec,R_P(:,1));
xlabel('d (mm)')
ylabel('R (mOhm)')
grid on
