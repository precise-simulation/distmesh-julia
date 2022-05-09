include( "DistMesh.jl" )
include( "dpolygon.jl" )

fid  = 0
hmax = [ 0.4 0.2 0.1 0.05 0.025 0.0125 0.00625 ]
# hmax = [ 0.4 0.2 0.1 ]

pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
      1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]
fh = p -> ones(size(p,1))


# Warmup
for iwarmup=1:5
    print(@sprintf("  warmup-%i\n",iwarmup))
    h = 0.4
    fd = p -> sqrt.(sum(p.^2,2)) - 1
    (p, t) = distmesh( fd, fh, h, [-1 -1;1 1], [], [], 1000, fid )

    fd = p -> dpolygon( p, pv )
    (p, t) = distmesh( fd, fh, h, [-1 -1;2 1], pv, [], 1000, fid )
end


# Circle
fd = p -> sqrt.(sum(p.^2,2)) - 1
t1 = Array{Float64}(length(hmax),5)
for i=1:length(hmax)
    print(@sprintf("  circ-%i\n",i))
    tic()
    (p, t, stat) = distmesh( fd, fh, hmax[i], [-1 -1;1 1], [], [], 1000, fid )
    t1[i,1] = toq()
    t1[i,2] = stat.ttri
    t1[i,3] = stat.nit
    t1[i,4] = stat.ntri
    t1[i,5] = size(t,1)

    if( hmax[i]==0.2 )
        writedlm( "p_circ.txt", p, " " )
        writedlm( "t_circ.txt", t, " " )
    end
end
writedlm( "circ_jl.txt", t1, " " )


# Polygon
fd = p -> dpolygon( p, pv )
t2 = Array{Float64}(length(hmax),5)
for i=1:length(hmax)
    print(@sprintf("  poly-%i\n",i))
    tic()
    (p, t, stat) = distmesh( fd, fh, hmax[i], [-1 -1;2 1], pv, [], 1000, fid )
    t2[i,1] = toq()
    t2[i,2] = stat.ttri
    t2[i,3] = stat.nit
    t2[i,4] = stat.ntri
    t2[i,5] = size(t,1)

    if( hmax[i]==0.2 )
        writedlm( "p_poly.txt", p, " " )
        writedlm( "t_poly.txt", t, " " )
    end
end
writedlm( "poly_jl.txt", t2, " " )
