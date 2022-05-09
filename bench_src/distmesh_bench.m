% function distmesh_bench()
addpath('C:\Users\jshys\Dropbox\code\featool-repo\featool\lib\distmesh')

pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
      1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
fdp = @(p) dpolygon( p, pv );
fd  = @(p) sqrt(sum(p.^2,2)) - 1;
fh  = @(p) ones(size(p,1),1);

fid  = 0;
hmax = [ 0.4 0.2 0.1 0.05 0.025 0.0125 0.00625 ];
% hmax = [ 0.4 0.2 0.1 ];



% Warmup
for iwarmup=1:5
  fprintf(1,'  warmup-%i\n',iwarmup);
  h = 0.4;
  [p, t] = distmesh( fd, fh, h, [-1 -1;1 1], [], [], 1000, fid );

  [p, t] = distmesh( fdp, fh, h, [-1 -1;2 1], pv, [], 1000, fid );
end


% Circle
t1 = zeros(length(hmax),5);
for i=1:length(hmax)
  fprintf(1,'  circ-%i\n',i)
  tic()
  [p, t, stat] = distmesh( fd, fh, hmax(i), [-1 -1;1 1], [], [], 1000, fid );
  t1(i,1) = toc();
  t1(i,2) = stat.ttri;
  t1(i,3) = stat.nit;
  t1(i,4) = stat.ntri;
  t1(i,5) = size(t,1);
end
if( exist('OCTAVE_VERSION','builtin') )
  save circ_ot.txt t1 -ascii
else
  save circ_ml.txt t1 -ascii
end


% Polygon
t2 = zeros(length(hmax),5);
for i=1:length(hmax)
  fprintf(1,'  poly-%i\n',i)
  tic()
  [p, t, stat] = distmesh( fdp, fh, hmax(i), [-1 -1;2 1], pv, [], 1000, fid );
  t2(i,1) = toc();
  t2(i,2) = stat.ttri;
  t2(i,3) = stat.nit;
  t2(i,4) = stat.ntri;
  t2(i,5) = size(t,1);
end
if( exist('OCTAVE_VERSION','builtin') )
  save poly_ot.txt t1 -ascii
else
  save poly_ml.txt t1 -ascii
end



exit
