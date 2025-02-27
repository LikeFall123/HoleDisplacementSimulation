clearvars
close all
clc

pl = 0;
u0 = 1;

mw = 0.001; % műszer érintkezés szélesség
zx = 0.077;   % zarvany x poz
zy = 0.07;    % zarvany y poz

d_muszer_1 = -0.075:0.001:0.075; %[mm] pozicio muszer 1
d_muszer_2 = -0.075:0.001:0.075; %[mm] pozicio muszer 2
% Ezeken végig lehet menni egyszerre mozgatva, vagy külön-külön, vagy
% ellentétesen.



%% Kettő műszer egyszerre
R_zar_1 = zeros(1, length(d_muszer_1));
R_1 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
    [R_zar_1(i),~] = fun_zarvany_m12(d_muszer_1(i),d_muszer_1(i),mw,zx,zy,u0, pl);
    [R_1(i),~] = fun_m12(d_muszer_1(i),d_muszer_1(i),mw,u0, pl);
end

%% Plot
figure()
title("Műszer 2 csatlakozása egyszerre")
plot(d_muszer_1,R_zar_1);
hold on
plot(d_muszer_1,R_1);
hold off
xlim([-0.08,0.08])
xlabel('d (m)')
ylabel('R (mOhm)')
grid on
legend('Zarvannyal','Zarvany nelkul','Location','northwest')

figure()
title("Inverz optimalizáció")
plot(d_muszer_1,abs(R_1-R_zar_1));
xlim([-0.08,0.08])
xlabel('d (m)')
ylabel('R (mOhm)')
grid on

%% Bal műszer
R_zar_2 = zeros(1, length(d_muszer_1));
R_2 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_zar_2(i),~] = fun_zarvany_m12(d_muszer_1(i),0, mw,zx,zy,u0, pl);
     [R_2(i),~] = fun_m12(d_muszer_1(i),0,mw,u0, pl);
end

%% Plot
figure()
title("Műszer bal csatlakozása")
plot(d_muszer_1,R_zar_2);
hold on
plot(d_muszer_1,R_2);
hold off
xlim([-0.08,0.08])
xlabel('d_{bal} (m)')
ylabel('R (mOhm)')
grid on
legend('Zarvannyal','Zarvany nelkul','Location','northwest')

figure()
title("Inverz optimalizáció")
plot(d_muszer_1,abs(R_2-R_zar_2));
xlim([-0.08,0.08])
xlabel('d_{bal} (m)')
ylabel('R (mOhm)')
grid on

%% Jobb műszer
R_zar_3 = zeros(1, length(d_muszer_1));
R_3 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_zar_3(i),~] = fun_zarvany_m12(0, d_muszer_1(i), mw,zx,zy,u0, pl);
     [R_3(i),~] = fun_m12(0,d_muszer_1(i),mw,u0, pl);
end

%% Plot
figure()
title("Műszer jobb csatlakozása")
plot(d_muszer_1,R_zar_3);
hold on
plot(d_muszer_1,R_3);
hold off
xlim([-0.08,0.08])
xlabel('d_{jobb} (m)')
ylabel('R (mOhm)')
grid on
legend('Zarvannyal','Zarvany nelkul','Location','northwest')

figure()
title("Inverz optimalizáció")
plot(d_muszer_1,abs(R_3-R_zar_3));
xlim([-0.08,0.08])
xlabel('d_{jobb} (m)')
ylabel('R (mOhm)')
grid on

%% Ellentétes irányban
R_zar_4 = zeros(1, length(d_muszer_1));
R_4 = zeros(1, length(d_muszer_1));

for i = 1:length(d_muszer_1)
     [R_zar_4(i),~] = fun_zarvany_m12(d_muszer_1(i), d_muszer_1(end+1-i), mw,zx,zy,u0, pl);
     [R_4(i),~] = fun_m12(d_muszer_1(i),d_muszer_1(end+1-i),mw,u0, pl);
end

%% Plot
figure()
title("Műszer csatlakozás ellentétes irány")
plot(d_muszer_1,R_zar_4);
hold on
plot(d_muszer_1,R_4);
hold off
xlim([-0.08,0.08])
xlabel('d_{ell} (m)')
ylabel('R (mOhm)')
grid on
legend('Zarvannyal','Zarvany nelkul','Location','northwest')

figure()
title("Inverz optimalizáció")
plot(d_muszer_1,abs(R_4-R_zar_4));
xlim([-0.08,0.08])
xlabel('d_{ell} (m)')
ylabel('R (mOhm)')
grid on

%% Minden pozíció
R_mat_zarv = zeros(length(d_muszer_1), length(d_muszer_2));
R_mat = zeros(length(d_muszer_1), length(d_muszer_2));

for i = 1:length(d_muszer_1)
    for j = 1:length(d_muszer_2)
        try
            [R,~] = fun_zarvany_m12(d_muszer_1(i), d_muszer_2(j), mw,zx,zy,u0, pl);
            R_mat_zarv(i,j) = R;
            [R,~] = fun_m12(d_muszer_1(i), d_muszer_2(j), mw,u0, pl);
            R_mat(i,j) = R;
        catch mex
            R_mat_zarv(i,j) = R_mat_zarv(i,j-1);
            R_mat(i,j) = R_mat(i,j-1);
        end

    end
end

%% Plot
[X,Y] = meshgrid(d_muszer_1,d_muszer_2);

figure()
title("Műszer csatlakozás")
surf(X,Y,R_mat_zarv)
hold on
surf(X,Y,R_mat)
hold off
xlim([-0.08,0.08])
xlabel('d (m)')
ylim([-0.08,0.08])
ylabel('d (m)')
zlabel('R (mOhm)')
grid on
legend('Zarvannyal','Zarvany nelkul','Location','northwest')

figure()
title("Inverz optimalizáció")
surf(X,Y,abs(R_mat_zarv-R_mat));
xlim([-0.08,0.08])
xlabel('d_1 (m)')
ylim([-0.08,0.08])
ylabel('d_2 (m)')
zlabel('R (mOhm)')
grid on