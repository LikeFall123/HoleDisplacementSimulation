function [R,P] = fun_zarvany(d,u0, pl)


sig = 2000;

% N = 500; % tekercs menetszam
% A = 90; %[mm2] tekercseles kersztmetszete
% I = 180; %[mA] tekercs aram
% J = N*I/A; % aramsuruseg

geom = [...
        2	2	2	2	4	4	4	4
        0.0799555555555556	0.0799555555555556	0.0599555555555556	0.0599555555555556	0.069+d	    0.070+d	    0.071+d	    0.070+d
        0.0799555555555556	0.0599555555555556	0.0599555555555556	0.0799555555555556	0.070+d	    0.071+d	    0.070+d	    0.069+d
        0.100000000000000	-0.100000000000000	-0.100000000000000	0.100000000000000	0	        -0.001	    0	        0.001
        -0.100000000000000	-0.100000000000000	0.100000000000000	0.100000000000000	-0.001	    0	        0.001	    0
        0	0	0	0	0	0	0	0
        1	1	1	1	1	1	1	1
        0	0	0	0	0.0700000000000000+d	0.0700000000000000+d	0.0700000000000000+d	0.0700000000000000+d
        0	0	0	0	0	0	0	0
        0	0	0	0	0.00100000000000000	0.00100000000000000	0.00100000000000000	0.00100000000000000
        0	0	0	0	0.00100000000000000	0.00100000000000000	0.00100000000000000	0.00100000000000000
        0	0	0	0	0	0	0	0];


bnd = [...
        1	1	1	1	1	1	1	1
        1	0	1	0	0	0	0	0
        1	1	1	1	1	1	1	1
        1	1	1	1	1	1	1	1
        1	48	1	48	48	48	48	48
        1	48	1	48	48	48	48	48
        48	48	48	48	48	48	48	48
        48	48	48	48	48	48	48	48
        49	49	49	49	49	49	49	49
        48	48	49	48	48	48	48	48];



c_coeff = 'x*2000';
a_coeff = '0.0'; 
f_coeff = '0';    

[p, e, t] = initmesh(geom);
[p, e, t] = refinemesh(geom,p,e,t);
[p, e, t] = refinemesh(geom,p,e,t);
[p, e, t] = refinemesh(geom,p,e,t);

u = assempde(bnd,p,e,t,c_coeff,a_coeff,f_coeff);

if pl == 1
    pdegplot(geom)
    hold on 
    pdeplot(p,[],t,'XYData',u,'Contour','on','Levels',20,'ColorBar','off')
    hold off
end

r_tri = pdeintrp(p,t,p(1,:)'); % haromszog kozeppontok r koordinataja

[Er, Ez] = pdegrad(p, t, -u);

n_tri = size(t,2); % haromszog elemek szama

sig_tri = zeros(1, n_tri);
sig_tri(t(4,:) == 1) = sig;

area = pdetrg(p,t);

% joule-ho a teljes modellben 
P = sum(sig_tri.*(Er.^2 + Ez.^2).*(2*pi*r_tri).*area);

% ellenallas
R = (u0^2)/P*10e3; % mOhm