function [R,P] = fun_zarvany_xy(d,l,u0, pl)


sig = 2000;


geom = [...
        2	    2	    2	    2	    4	        4	        4	        4
        0.08	0.08	0.06	0.06	d-0.001	    d	        d+0.001	    d
        0.08	0.06	0.06	0.08	d	        d+0.001	    d	        d-0.001
        0.1	    -0.1	-0.1	0.1	    l	        l-0.001	    l	        l+0.001
        -0.1	-0.1	0.1	    0.1	    l-0.001	    l	        l+0.001	    l
        0	0	0	0	0	0	0	0
        1	1	1	1	1	1	1	1
        0	0	0	0	d	d	d	d
        0	0	0	0	l	l	l	l
        0	0	0	0	0.001	0.001	0.001	0.001
        0	0	0	0	0.001	0.001	0.001	0.001
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