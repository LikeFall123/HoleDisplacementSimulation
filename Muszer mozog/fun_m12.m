function [R,P] = fun_m12(d_m1, d_m2, mw, u0, pl)

% Ellenallas szomolas zarvany nelkul inverz optimalizáláshoz

% d_m1 - x = 0.06-ban található műszer y pozíció
% d_m2 - x = 0.08-ban található műszer y pozíció
% mw - műszer érintkezési szélesség
% zx - zarvany x pozíció  - nincs = 0
% zy - zarvany y pozíció  - nincs = 0
% u0 - Dirichlet peremen a gerjesztő feszültség
% pl - plotoljunk-e

sig = 2000;


geom = [...
        2	2	2	2	2	2	2	2
        0.06	0.08	0.06	0.06	0.06	0.08	0.08	0.08
        0.08	0.06	0.06	0.06	0.06	0.08	0.08	0.08
        -0.10	0.1	    -0.1	d_m1-mw	d_m1+mw	-0.1	d_m2-mw	d_m2+mw
        -0.1	0.1	    d_m1-mw	d_m1+mw	0.1	    d_m2-mw	d_m2+mw	0.1
        1	1	0	0	0	1	1	1
        0	0	1	1	1	0	0	0];


bnd = [...
        1	1	1	1	1	1	1	1
        0	0	0	1	0	0	1	0
        1	1	1	1	1	1	1	1
        1	1	1	1	1	1	1	1
        48	48	48	1	48	48	1	48
        48	48	48	1	48	48	1	48
        48	48	48	48	48	48	48	48
        48	48	48	48	48	48	48	48
        49	49	49	49	49	49	49	49
        48	48	48	49	48	48	48	48];


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

r_tri = pdeintrp(p,t,p(1,:)');

[Er, Ez] = pdegrad(p, t, -u);

n_tri = size(t,2);

sig_tri = zeros(1, n_tri);
sig_tri(t(4,:) == 1) = sig;

area = pdetrg(p,t);

P = sum(sig_tri.*(Er.^2 + Ez.^2).*(2*pi*r_tri).*area);
R = (u0^2)/P*10e3; % mOhm