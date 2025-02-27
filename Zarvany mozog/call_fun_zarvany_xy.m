clearvars
close all
clc

pl = 0;
u0 = 1;
dx_vec = 0.07+(-0.008:0.001:0.008); %[mm] pozicio x
dy_vec = -0.09:0.001:0.09; %[mm] pozicio y


R_mat = zeros(length(dy_vec), length(dx_vec));

for i = 1:length(dy_vec)
    for j = 1:length(dx_vec)
        [R,~] = fun_zarvany_xy(dx_vec(j),dy_vec(i),u0, pl);
        R_mat(i,j) = R;
        % disp("R="+num2str(R)+"(i="+num2str(i)+",j="+num2str(j))
        % pause
    end
end

%% Plot
figure()
[X,Y] = meshgrid(dx_vec,dy_vec);
surf(X,Y,R_mat)
xlabel('d (mm)')
ylabel('R (mOhm)')

% Plot y
plot(dy_vec,R_mat(:,10));
%%
[R,P,u,p,e,t] = fun_zarvany_xy_plot(0.07,0,1, 0);
[Er, Ez] = pdegrad(p, t, -u);
figure()
pdeplot(p,[],t,'XYData',u,'ZData',u)
figure()
pdeplot(p,[],t,'XYData',Ez,'ZData',Ez)
%% 

d = 0.08;
h = 0.02;
l = 0.2;
s = 2000;
Ra = 1/(2*pi*s*l)*(log(d/(d-h)))*1000
d = 0.08;
h = 0.05;
Ra = 1/(2*pi*s*l)*(log(d/(d-h)))*1000
%Ra1 = 1/(2*pi*s*l)*(log(((d-h)+1/3*(d-h))/(d-h))+log(d/(d-h+2/3*(d-h))))*1000 %?
