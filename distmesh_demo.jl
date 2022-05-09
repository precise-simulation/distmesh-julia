# distmesh_demo

# Simple mesh generation in Julia with the DistMesh algorithm where
# geometries/domains are defined by (level set) distance functions.
# Julia implementation of the Matlab distmesh demos
# (https://github.com/precisesimulation/distmesh)
#
# distmesh-julia is a collection of Julia functions for generation and
# manipulation of unstructured meshes. distmesh-julia is Copyright (C)
# 2018 J.S. Hysing and Precise Simulation Limited.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


# using DistMesh
include( "DistMesh.jl" )
using Plots


"""
(P, T) = distmesh_demo( CASES, DO_PLOT ) Run distmesh example
demos. CASES may be used to specify the distmesh examples to run
(default 1:3). DO_PLOT specifies whether to plot or not (default
true).

Example 1: (Uniform mesh on unit circle)
    fd = p -> sqrt.(sum(p.^2,2)) - 1
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.1, [-1 -1;1 1] )
    plotgrid( p, t )

Example 2: (Uniform mesh on ellipse)
    fd = p -> p[:,1].^2/2^2 + p[:,2].^2/1^2 - 1
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.2, [-2 -1;2 1], [], [], 200 )
    plotgrid( p, t )

Example 3: (Uniform mesh on unit square)
    fd = p -> -minimum([minimum([minimum([p[:,2] 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.2, [-1 -1;1 1], [-1 -1;-1 1;1 -1;1 1] )
    plotgrid( p, t )

Example 4: (Uniform mesh on complex polygon)
    pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
          1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]
    fd = p -> dpolygon( p, pv )
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.1, [-1 -1;2 1], pv )
    plotgrid( p, t )

Example 5: (Rectangle with circular hole, refined at circle boundary)
    fd = p -> maximum( [dpolygon(p,[-1 -1;1 -1;1 1;-1 1;-1 -1]) -(sqrt.(sum(p.^2,2))-0.5)], 2 )
    fh = p -> 0.05 + 0.3*(sqrt.(sum(p.^2,2))-0.5)
    (p, t) = distmesh( fd, fh, 0.05, [-1 -1;1 1], [-1 -1;-1 1;1 -1;1 1] )
    plotgrid( p, t )

Example 6: (Square, with size function point and line sources)
    dcircle = (p,xc,yc,r) -> sqrt.((p[:,1]-xc).^2+(p[:,2]-yc).^2)-r
    fd = p -> -minimum([minimum([minimum([p[:,2] 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> minimum([minimum(hcat(0.01+0.3*abs.(dcircle(p,0,0,0)),
                  0.025+0.3*abs.(dpolygon(p,[0.3 0.7;0.7 0.5;0.3 0.7]))),2) 0.15*ones(size(p,1),1)],2)
    (p, t) = distmesh( fd, fh, 0.01, [0 0;1 1], [0 0;1 0;0 1;1 1] )
    plotgrid( p, t )

Example 7: (NACA0012 airfoil)
    hlead = 0.01; htrail = 0.04; hmax = 2; circx = 2; circr = 4
    a = 0.12/0.2*[0.2969 -0.126 -0.3516 0.2843 -0.1036]
    dcircle = (p,xc,yc,r) -> sqrt.((p[:,1]-xc).^2+(p[:,2]-yc).^2)-r
    fd = p -> maximum(hcat( dcircle(p,circx,0,circr),
               -((abs.(p[:,2])-polyval([a[5:-1:2];0],p[:,1])).^2-a[1]^2*p[:,1]) ), 2 )
    fh = p -> minimum([minimum([hlead+0.3*dcircle(p,0,0,0) htrail+0.3*dcircle(p,1,0,0)],2) hmax*ones(size(p,1),1)],2)

    fixx = 1 - htrail*cumsum(1.3.^(0:4)')
    fixy = a[1]*sqrt(fixx) + polyval([a[5:-1:2];0],fixx)
    pfix = [[circx+[-1 1 0 0]*circr; 0 0 circr*[-1 1]]'; 0 0; 1 0; fixx' fixy'; fixx' -fixy']
    bbox = [circx-circr -circr; circx+circr circr]
    h0   = minimum([hlead htrail hmax])

    (p, t) = distmesh( fd, fh, h0, bbox, pfix )
    plotgrid( p, t )

Example 8: (Uniform mesh on unit sphere)
    fd = p -> sqrt.(sum(p.^2,2)) - 1
    fh = p -> ones(size(p,1),1)
    (p,t) = distmesh( fd, fh, 0.2, [-1 -1 -1;1 1 1] )
    plotgrid( p, t )

Example 9: (Uniform mesh on unit cube)
    fd = p -> -minimum([minimum([minimum([minimum([minimum([p[:,3] 1-p[:,3]],2) p[:,2]],2) 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> ones(size(p,1),1)
    pfix = [-1 -1 -1;-1 1 -1;1 -1 -1;1 1 -1; -1 -1 1;-1 1 1;1 -1 1;1 1 1]
    (p,t) = distmesh( fd, fh, 0.2, [-1 -1 -1;1 1 1], pfix )
    plotgrid( p, t )

Example 10: (Uniform mesh on cylinder)
    fd = p -> -minimum([minimum([p[:,3] 4-p[:,3]],2) 1-sqrt.(sum(p[:,1:2].^2,2))],2)
    fh = p -> ones(size(p,1),1)
    pfix = [-1 -1 -1;-1 1 -1;1 -1 -1;1 1 -1; -1 -1 1;-1 1 1;1 -1 1;1 1 1]
    (p,t) = distmesh( fd, fh, 0.5, [-1 -1 0;1 1 4] )
    plotgrid( p, t )

See also distmesh.
"""
function distmesh_demo( cases=1:7, do_plot=true )

    feval(fn_str, args...) = eval(parse(fn_str))(args...)

    p = []
    t = []
    for i=cases

        fn = "test_case" * @sprintf("%i",i)
        (p, t) = feval( fn, do_plot )

        if( do_plot )
            plotgrid( p, t )
        end
        if( do_plot && i!=cases[end] )
            sleep(2000)
            # pause()
        end

        print("\n\n")
    end

    return (p, t)
end


function test_case1( do_plot )
    s = "Example 1: (Uniform mesh on unit circle)\n"
    print(s)
    fd = p -> sqrt.(sum(p.^2,2)) - 1
    fh = p -> ones(size(p,1))
    # fd(p) = sqrt.(sum(p.^2,2)) - 1
    # fh(p) = ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.025, [-1 -1;1 1] )

    return (p, t)
end


function test_case2( do_plot )
    s = "Example 2: (Uniform mesh on ellipse)\n"
    print(s)
    fd = p -> p[:,1].^2/2^2 + p[:,2].^2/1^2 - 1
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.2, [-2 -1;2 1] )

    return (p, t)
end


function test_case3( do_plot )
    s = "Example 3: (Uniform mesh on unit square)\n"
    print(s)
    fd = p -> -minimum([minimum([minimum([p[:,2] 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.2, [-1 -1;1 1], [-1 -1;-1 1;1 -1;1 1] )

    return (p, t)
end


function test_case4( do_plot )
    s = "Example 4: (Uniform mesh on complex polygon)\n"
    print(s)
    pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
          1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]
    fd = p -> dpolygon( p, pv )
    fh = p -> ones(size(p,1))
    (p, t) = distmesh( fd, fh, 0.1, [-1 -1;2 1], pv )

    return (p, t)
end


function test_case5( do_plot )
    s = "Example 5: (Rectangle with circular hole, refined at circle boundary)\n"
    print(s)
    fd = p -> maximum( [dpolygon(p,[-1 -1;1 -1;1 1;-1 1;-1 -1]) -(sqrt.(sum(p.^2,2))-0.5)], 2 )
    fh = p -> 0.05 + 0.3*(sqrt.(sum(p.^2,2))-0.5)
    (p, t) = distmesh( fd, fh, 0.05, [-1 -1;1 1], [-1 -1;-1 1;1 -1;1 1] )

    return (p, t)
end


function test_case6( do_plot )
    s = "Example 6: (Square, with size function point and line sources)\n"
    print(s)
    dcircle = (p,xc,yc,r) -> sqrt.((p[:,1]-xc).^2+(p[:,2]-yc).^2)-r
    fd = p -> -minimum([minimum([minimum([p[:,2] 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> minimum([minimum(hcat(0.01+0.3*abs.(dcircle(p,0,0,0)),
                  0.025+0.3*abs.(dpolygon(p,[0.3 0.7;0.7 0.5;0.3 0.7]))),2) 0.15*ones(size(p,1),1)],2)
    (p, t) = distmesh( fd, fh, 0.01, [0 0;1 1], [0 0;1 0;0 1;1 1] )

    return (p, t)
end


function test_case7( do_plot )
    s = "Example 7: (NACA0012 airfoil)\n"
    print(s)
    hlead = 0.01; htrail = 0.04; hmax = 2; circx = 2; circr = 4
    a = 0.12/0.2*[0.2969 -0.126 -0.3516 0.2843 -0.1036]
    dcircle = (p,xc,yc,r) -> sqrt.((p[:,1]-xc).^2+(p[:,2]-yc).^2)-r
    fd = p -> maximum(hcat( dcircle(p,circx,0,circr),
               -((abs.(p[:,2])-polyval([a[5:-1:2];0],p[:,1])).^2-a[1]^2*p[:,1]) ), 2 )
    fh = p -> minimum([minimum([hlead+0.3*dcircle(p,0,0,0) htrail+0.3*dcircle(p,1,0,0)],2) hmax*ones(size(p,1),1)],2)

    fixx = 1 - htrail*cumsum(1.3.^(0:4)')
    fixy = a[1]*sqrt.(fixx) + polyval([a[5:-1:2];0],fixx)
    pfix = [[circx+[-1 1 0 0]*circr; 0 0 circr*[-1 1]]'; 0 0; 1 0; fixx' fixy'; fixx' -fixy']
    bbox = [circx-circr -circr; circx+circr circr]
    h0   = minimum([hlead htrail hmax])

    (p, t) = distmesh( fd, fh, h0, bbox, pfix )

    return (p, t)
end

function test_case8( do_plot )
    s = "Example 8: Uniform mesh on unit sphere\n"
    print(s)
    fd = p -> sqrt.(sum(p.^2,2)) - 1
    fh = p -> ones(size(p,1),1)
    (p,t) = distmesh( fd, fh, 0.2, [-1 -1 -1;1 1 1] )

    return (p, t)
end

function test_case9( do_plot )
    s = "Example 9 Uniform mesh on unit cube: \n"
    print(s)
    fd = p -> -minimum([minimum([minimum([minimum([minimum([p[:,3] 1-p[:,3]],2) p[:,2]],2) 1-p[:,2]],2) p[:,1]],2) 1-p[:,1]],2)
    fh = p -> ones(size(p,1),1)
    pfix = [-1 -1 -1;-1 1 -1;1 -1 -1;1 1 -1; -1 -1 1;-1 1 1;1 -1 1;1 1 1]
    (p,t) = distmesh( fd, fh, 0.2, [-1 -1 -1;1 1 1], pfix )

    return (p, t)
end

function test_case10( do_plot )
    s = "Example 10: Uniform mesh on cylinder\n"
    print(s)
    fd = p -> -minimum([minimum([p[:,3] 4-p[:,3]],2) 1-sqrt.(sum(p[:,1:2].^2,2))],2)
    fh = p -> ones(size(p,1),1)
    pfix = [-1 -1 -1;-1 1 -1;1 -1 -1;1 1 -1; -1 -1 1;-1 1 1;1 -1 1;1 1 1]
    (p,t) = distmesh( fd, fh, 0.5, [-1 -1 0;1 1 4] )

    return (p, t)
end



#------------------------------------------------------------------------------#
function dpolygon( p, v )

    n_p = size(p,1)
    n_s = size(v,1)-1
    s = Vector{Float64}(n_p)
    dist = Array{Float64}(n_p,n_s)
    for i_s=1:n_s
        v_i = v[[i_s,i_s+1],:]
        w   = v_i[2,:]-v_i[1,:]
        vp  = repmat(v_i[1,:]',n_p,1)-p
        w1  = repmat(w',n_p,1)
        s   = sum(w1.*vp,2)
        u   = -s./(w[1]^2+w[2]^2)
        u[u.<0.0] = 0.0
        u[u.>1.0] = 1.0
        h = w1.*[u u]+vp
        dist[:,i_s] = sqrt.(sum(h.^2,2))
    end
    # dist = (-1.0).^(inpolygonv(p[:,1],p[:,2],v[:,1],v[:,2])).*minimum(dist,2)
    dist = minimum(dist,2)
    ind  = inpolygonv( p[:,1], p[:,2], v[:,1], v[:,2] )
    dist[ind] = -dist[ind]

    return dist
end

#------------------------------------------------------------------------------#
"return true for Points inside the Polygon"
function inpolygonv( x, y, xv, yv )
    # Octave inpolygon implementation.

    npol = length(xv)
    in = fill(false,length(x))
    on = fill(false,length(x))

    j = npol
    for i = 1:npol
        delta_xv = xv[j] - xv[i]
        delta_yv = yv[j] - yv[i]

        ## distance = [distance from (x,y) to edge] * length(edge)
        distance = delta_xv * (y - yv[i]) - (x - xv[i]) * delta_yv

        ## is y between the y-values of edge i,j AND (x,y) on the left of the edge?
        idx1 = ( (( (yv[i] .<= y) .& (y .< yv[j]) ) .| ( (yv[j] .<= y) .& (y .< yv[i]) ))
                 .& (0 .< distance*delta_yv) )
        in[idx1] = .! in[idx1]

        ## Check if (x,y) are actually on the boundary of the polygon.
        idx2 = ((( (yv[i] .<= y) .& (y .<= yv[j]) ) .| ( (yv[j] .<= y) .& (y .<= yv[i]) ))
                .& (( (xv[i] .<= x) .& (x .<= xv[j]) ) .| ( (xv[j] .<= x) .& (x .<= xv[i]) ))
                .& ( (0 .== distance) .| (delta_xv .== 0)))
        on[idx2] = true

        j = i
    end

    return in .| on
    # return (in .| on, on)
end

#------------------------------------------------------------------------------#
function polyval( p, x )
    y = p[1]*ones(size(x))
    for i=2:length(p)
        y = y.*x + p[i]
    end

    return y
end



(p, t) = distmesh_demo(1)
# plotgrid( p, t )
nothing
