% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 5.5000000000000009 50.000000000000007]);
set(ax,'PlotBoxAspectRatio',[1 1 1]);
set(ax,'XLim',[0.050000000000000003 0.089999999999999997]);
set(ax,'YLim',[-0.11 0.11]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pdepoly([ 0.059999999999999998,...
 0.059999999999999998,...
 0.059999999999999998,...
 0.059999999999999998,...
 0.080000000000000002,...
 0.080000000000000002,...
 0.080000000000000002,...
 0.080000000000000002,...
],...
[ 0.10000000000000001,...
 0.0050000000000000001,...
 -0.0050000000000000001,...
 -0.10000000000000001,...
 -0.10000000000000001,...
 -0.0050000000000000001,...
 0.0050000000000000001,...
 0.10000000000000001,...
],...
 'P1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','P1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(8,...
'neu',...
1,...
'0',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'neu',...
1,...
'0',...
'0')
pdesetbd(5,...
'neu',...
1,...
'0',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'1')
pdesetbd(3,...
'neu',...
1,...
'0',...
'0')
pdesetbd(2,...
'neu',...
1,...
'0',...
'0')
pdesetbd(1,...
'neu',...
1,...
'0',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')
pdetool('refine')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'x*2000',...
'0.0',...
'0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['x*2000';...
'0.0   ';...
'0     ';...
'1.0   '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','5760','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 2 1 1 1 2 10 1 0 1 0 1 1 0 0 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
