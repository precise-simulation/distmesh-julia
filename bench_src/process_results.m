% cl;

addpath(genpath('C:\Users\jshys\Dropbox\code\featool-repo\featool'))
rmpath('C:\Users\jshys\Dropbox\code\featool-repo\featool\geom1')

lw = 5;
ms = 15;
fs = 20;

env = {'jl','ml','ot','tr'};
lng = {'Julia', 'Matlab', 'Octave','Triangle'};
col = {'r','g','b','c','m'};


f1 = figure; hold on
f2 = figure; hold on
for i=1:length(env)
  v = load( ['circ_',env{i},'.txt'] );
  x = v(:,5);
  y = v(:,1);
  figure(f1)
  plot( x, y, '.-', 'linewidth', lw, 'markersize', ms, 'color', col{mod(i-1,length(col))+1} )


  v = load( ['poly_',env{i},'.txt'] );
  x = v(:,5);
  y = v(:,1);
  figure(f2)
  plot( x, y, '.-', 'linewidth', lw, 'markersize', ms, 'color', col{mod(i-1,length(col))+1} )
end


figure(1)
h = legend( lng );
set( h, 'location', 'northwest', 'fontsize', fs )
title( 'Grid generation timing (Circle)', 'fontsize', fs )
ylabel( 'CPU time [s]', 'fontsize', fs )
xlabel( 'Number of grid cells', 'fontsize', fs )
set( gca, 'fontsize', round(fs*0.8) )
grid on
print -r300 -djpeg circle_distmesh_benchmark_timing.jpg


figure(2)
h = legend( lng );
set( h, 'location', 'northwest', 'fontsize', fs )
title( 'Grid generation timing (Polygon)', 'fontsize', fs )
ylabel( 'CPU time [s]', 'fontsize', fs )
xlabel( 'Number of grid cells', 'fontsize', fs )
set( gca, 'fontsize', round(fs*0.8) )
grid on
print -r300 -djpeg polygon_distmesh_benchmark_timing.jpg




% Create grid visualizations.
close all
figure
subplot(2,2,1)
p = load('p_circ.txt');
t = load('t_circ.txt');
grid1.p = p';
grid1.c = t';
grid1.a = gridadj(grid1.c,2);
grid1.b = gridbdr(grid1.p,grid1.c,grid1.a);
grid1.s = ones(1,size(t,1));
plotgrid(grid1)
axis off
title('DistMesh')
ylabel('Circle')
% print -r300 -djpeg julia_circ_grid.jpg


subplot(2,2,3)
p = load('p_poly.txt');
t = load('t_poly.txt');
grid1.p = p';
grid1.c = t';
grid1.a = gridadj(grid1.c,2);
grid1.b = gridbdr(grid1.p,grid1.c,grid1.a);
grid1.s = ones(1,size(t,1));
plotgrid(grid1)
axis off
ylabel('Polygon')
% print -r300 -djpeg julia_poly_grid.jpg


hmax = 0.2;
pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
      1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];

% fdp = @(p) dpolygon( p, pv );
% fd  = @(p) sqrt(sum(p.^2,2)) - 1;
% fh  = @(p) ones(size(p,1),1);
% [p1,t1] = distmesh( fd,  fh, hmax, [-1 -1;1 1], [], [], 1000, 1 );
% grid1.p = p1';
% grid1.c = t1';
% grid1.a = gridadj(grid1.c,2);
% grid1.b = gridbdr(grid1.p,grid1.c,grid1.a);
% grid1.s = ones(1,size(t1,1));
% clf
% plotgrid(grid1)
% axis off
% print -r300 -djpeg mloct_circ_grid.jpg


% [p2,t2] = distmesh( fdp, fh, hmax, [-1 -1;2 1], pv, [], 1000, 1 );
% grid1.p = p2';
% grid1.c = t2';
% grid1.a = gridadj(grid1.c,2);
% grid1.b = gridbdr(grid1.p,grid1.c,grid1.a);
% grid1.s = ones(1,size(t2,1));
% clf
% plotgrid(grid1)
% axis off
% print -r300 -djpeg mloct_poly_grid.jpg

subplot(2,2,2)
geom.objects = { gobj_circle };
grid1 = gridgen( geom, 'hmax', hmax, 'meshgen', 'triangle' );
plotgrid(grid1)
title('Triangle')
axis off


subplot(2,2,4)
geom.objects = { gobj_polygon(pv) };
grid1 = gridgen( geom, 'hmax', hmax, 'meshgen', 'triangle' );
plotgrid(grid1)
axis off

print -r300 -djpeg distmesh_benchmark_grids.jpg

% exit
