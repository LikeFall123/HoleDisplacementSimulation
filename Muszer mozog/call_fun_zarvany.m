clearvars
close all
clc

pl = 0;
u0 = 1;

mw = 0.00005; % műszer érintkezés szélesség
zx = 0.075;   % zarvany x poz
zy = -0.03;    % zarvany y poz

d_muszer_1 = -0.1+2*mw:0.0001:0.1-2*mw; %[mm] pozicio muszer 1
d_muszer_2 = -0.1+2*mw:0.0001:0.1-2*mw; %[mm] pozicio muszer 2
% Ezeken végig lehet menni egyszerre mozgatva, vagy külön-külön, vagy
% ellentétesen.



%% Kettő műszer egyszerre
R_vec_1 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_vec_1(i),~] = fun_zarvany_m12(d_muszer_1(i),d_muszer_1(i),mw,zx,zy,u0, pl);
end

% Plot
figure()
title("Műszer 2 csatlakozása egyszerre")
% [X,Y] = meshgrid(dx_vec,dy_vec);
% surf(X,Y,R_mat)
plot(d_muszer_1,R_vec_1);
xlabel('d (mm)')
ylabel('R (mOhm)')
grid on

%% Bal műszer
R_vec_2 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_vec_2(i),~] = fun_zarvany_m12(d_muszer_1(i),0, mw,zx,zy,u0, pl);
end

% Plot
figure()
title("Műszer bal csatlakozása")
plot(d_muszer_1,R_vec_2);
xlabel('d_{bal} (mm)')
ylabel('R (mOhm)')

%% Jobb műszer
R_vec_3 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_vec_3(i),~] = fun_zarvany_m12(0, d_muszer_1(i), mw,zx,zy,u0, pl);
end

% Plot
figure()
title("Műszer jobb csatlakozása")
plot(d_muszer_1,R_vec_3);
xlabel('d bal (mm)')
ylabel('R (mOhm)')

%% Ellentétes irányban
R_vec_4 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_vec_4(i),~] = fun_zarvany_m12(d_muszer_1(i), d_muszer_1(end+1-i), mw,zx,zy,u0, pl);
end

% Plot
figure()
title("Műszer 2 csatlakozása ellentétesen")
plot(d_muszer_1,R_vec_4);
xlabel('d bal (mm)')
ylabel('R (mOhm)')

%% Minden pozíció
R_mat = zeros(length(d_muszer_1), length(d_muszer_2));

for i = 1:length(d_muszer_1)
    for j = 1:length(d_muszer_2)
        [R,~] = fun_zarvany_m12(d_muszer_1(i), d_muszer_2(j), mw,zx,zy,u0, pl);
        R_mat(i,j) = R;
    end
end

%% Plot
figure()
[X,Y] = meshgrid(d_muszer_1,d_muszer_2);
surf(X,Y,R_mat)
xlabel('d_1 (mm)')
ylabel('d_2 (mm)')
zlabel('R (mOhm)')
title("Műszer csatlakozása")


%% Plot the measurements

figure()
title("4 kind of measurements")

subplot(2,2,1)
plot(d_muszer_1,R_vec_1);
xlabel('d (mm)')
ylabel('R (mOhm)')
legend("Műszer 2 csatlakozása egyszerre")
grid on

subplot(2,2,2)
plot(d_muszer_1,R_vec_2);
xlabel('d_{bal} (mm)')
ylabel('R (mOhm)')
legend("Műszer bal csatlakozása")
grid on

subplot(2,2,3)
plot(d_muszer_1,R_vec_3);
xlabel('d_{jobb} bal (mm)')
ylabel('R (mOhm)')
legend("Műszer jobb csatlakozása")
grid on

subplot(2,2,4)
plot(d_muszer_1,R_vec_4);
xlabel('d_{ellen} (mm)')
ylabel('R (mOhm)')
legend("Műszer 2 csatlakozása ellentétesen")
grid on

%% Plot E vector and U
[R,P,u,p,e,t] = fun_zarvany_plot(0,0,1, 0);
[Er, Ez] = pdegrad(p, t, -u);
figure()
pdeplot(p,[],t,'XYData',u,'ZData',u)
title('U zarvany')

figure()
pdeplot(p,[],t,'XYData',Ez,'ZData',Ez)
title('Ez zarvany')

figure()
pdeplot(p,[],t,'XYData',Er,'ZData',Er)
title('Er zarvany')

