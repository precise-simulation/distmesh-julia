# module DistMesh

# Simple mesh generation in Julia with the DistMesh algorithm where
# geometries/domains are defined by (level set) distance functions.
# Julia implementation of the Matlab DistMesh function
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

# export distmesh, plotgrid, dpolygon


using  VoronoiDelaunay
import VoronoiDelaunay: getx, gety #, getz
using  Plots


type IndexedPoint2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _idx::Int64
    IndexedPoint2D(x, y, idx) = new(x, y, idx)
    IndexedPoint2D(x, y) = new(x, y, 0)
end

getx(p::IndexedPoint2D) = p._x
gety(p::IndexedPoint2D) = p._y
getidx(p::IndexedPoint2D) = p._idx

# type IndexedPoint3D <: AbstractPoint3D
#     _x::Float64
#     _y::Float64
#     _z::Float64
#     _idx::Int64
#     IndexedPoint3D(x, y, z, idx) = new(x, y, z, idx)
#     IndexedPoint3D(x, y, z) = new(x, y, z, 0)
# end

# getx(p::IndexedPoint3D) = p._x
# gety(p::IndexedPoint3D) = p._y
# gety(p::IndexedPoint3D) = p._z
# getidx(p::IndexedPoint3D) = p._idx

mutable struct statistics
    conv; nit; ntri; ndcs; dpmx; dpmn; ttot; ttri
end


#------------------------------------------------------------------------------#
"""
    distmesh(fd, fh, f0, bbox, p_fix, e_fix, it_max, fid)

DistMesh is a simple surface (in 2D) and volume (in 3D) mesh gen-
eration algorithm using distance functions to define geometries.

FD is a function handle to the geometry description that should take
evaluation coordinates and points as input. For example fd =
p->sqrt.(sum(p.^2,2))-1 specifies the distance function for a unit
circle (both function handles, string function names, and anonymous
functions are supported). Similar to FD, FH a function describing the
desired relative mesh size distribution. For example fh =
p->ones(size(p,1),1) specifies a uniform distribution where FH
evaluates to 1 at all points. H0 is a numeric scalar specifying the
initial edge lengths, and BBOX is a 2 by 2 in 2D (or 2 by 3 in 3D)
bounding box of the domain (enclosing the zero contour/level set of
FD). P_FIX optionally specifies a number of points that should always
be present (fixed) in the resulting mesh. E_FIX can be sets of edge
vertex indices to constrain, or alternatively a cell array with
function handle to call. IT_MAX sets the maximum number of grid
generation iterations allowed (default 1000). Finally, FID specifies a
file identifies for output (default 1 = terminal output).

The distmesh function returns the grid point vertices in P,
triangulated simplices in T, as well as an optional statistics struct
STAT including timings and convergence information.

   Input:

      FD:        Distance function d(x,y,(z))
      FH:        Scaled edge length function h(x,y,(z))
      H0:        Initial edge length
      BBOX:      Bounding box [xmin ymin (zmin); xmax ymax (zmax)]
      P_FIX:     Fixed node positions [N_P_FIX x 2/3]
      E_FIX:     Constrained edges [N_E_FIX x 2]
      IT_MAX:    Maximum number of iterations
      FID:       Output file id number (default 1 = terminal)

   Output:

      P:         Grid vertex/node coordinates [N_P x 2/3]
      T:         Triangle indices [N_T x 3]
      STAT:      Mesh generation statistics (struct)


Example: (Uniform mesh on unit circle)

   fd = p -> sqrt.(sum(p.^2,2)) - 1;
   fh = p -> ones(size(p,1),1);
   (p, t) = distmesh( fd, fh, 0.2, [-1 -1;1 1] );
   plotgrid( p, t )


See also: distmesh_demo, plotgrid.
"""
function distmesh( fd=p->sqrt.(sum(p.^2,2))-1, fh=p->ones(size(p,1)), h0=0.1,
                   bbox=[-1 -1;1 1], p_fix=[], e_fix=[], it_max=1000, fid=1 )

    t0 = time_ns()
    #------------------------------------------------------------------------------#
    # Initialization and meshing parameters.
    #------------------------------------------------------------------------------#
    IT_MIN  = 20                # Minimum number of iterations.
    IT_MINC = 50                # Minimum number of iter. after which to call constraint function.
    IT_PRT  = 25                # Output every IT_PRT iterations.

    N_RECV  = 2                 # Number of recovery iteration steps to move points outside back to boundary.
    N_DCF   = 3000                # Frequency of density control checks.
    n_sdim  = size(bbox,2)
    if( n_sdim==2 )
        dp_tol   = -0.001*h0    # Abs point rejection tol (p(dist(p)>=dp0_tol) are rejected).
        dtrm_tol = -0.001*h0    # Abs dist tol for tri rejection (t(dist(p_tcent)>=dtrm_tol) are rejected).
        rt_tol   =  0.3         # Rel fraction of h0 to trigger retriangulation.
        F_scale  =  1.2         # Rel force scaling factor.
        F_DCF    =  2.0         # Fraction of L to L_target to allow.
        dp_scale =  0.2         # Rel fraction of computed new distance to move points in update step.
    else
        dp_tol   = -0.1*h0
        dtrm_tol = -0.1*h0
        rt_tol   =  0.1
        F_scale  =  1.1
        F_DCF    =  2.1
        dp_scale =  0.1
    end
    dpc_tol = 0.001*h0          # Abs tol for grid point movements during convergence check.
    gradeps = sqrt(eps())*h0    # Gradient computation offset.
    #------------------------------------------------------------------------------#

    # Initial grid point distribution, p, confined to the bounding box.
    p = l_initpoints( h0, bbox )

    # Remove points outside the region and apply the rejection method.
    p = p[ find(fd(p).<-dp_tol), : ]
    if( isempty(p) )
        return (p, [], [])
    end

    r0 = fh( p )   # Probability to keep point.
    p  = p[ find( rand(size(p,1)) .< minimum(r0)^n_sdim./r0.^n_sdim ), : ]
    p_fix,  = l_deduplicate( p_fix )
    n_p_fix = size(p_fix,1)
    if( !isempty(p_fix) )
        p, = l_deduplicate( [ p_fix; p ] )
    end
    n_p = size( p, 1 )

    l_message( fid, "Grid generation (DistMesh):" )
    tic()
    if( it_max<=0 )
        t = l_delaunay_triangulation( p, e_fix )
        return (p, t, [])
    end
    t_tri = toq()
    it    = 0
    p0    = Inf
    n_tri = 0
    n_dcs = 0
    dist  = []
    delta_p = []
    do_break = false
    is_converged = false
    while( it<it_max )
        it = it + 1

        # Retriangulate, if grid points have moved significantly.
        delta_p_max = maximum( sqrt.(sum((p-p0).^2,2)) )
        if( rt_tol*h0<delta_p_max )
            n_tri = n_tri + 1
            (p, t, td) = l_triangulate( p, fd, e_fix, dtrm_tol )
            # if( !isempty(e_fix) && it>IT_MINC )
            #     [p,t] = l_call_function( e_fix, p, t, n_sdim, 1:n_p_fix )
            # end
            p0  = p
            n_p = size(p,1)
            t_tri = t_tri + td

            # Describe each edge by a unique edge_pairs of nodes.
            e = [ t[:,[1,2]]; t[:,[2,3]]; t[:,[3,1]] ]
            if( n_sdim==3 )
                e = [ e; t[:,[1,4]]; t[:,[2,4]]; t[:,[3,4]] ]
            end
            e = sort(e,2)
            e_max = maximum(e[:])
            if( e_max*(e_max+1)<realmax() )
                ecomp = (e_max+1)*e[:,1] + e[:,2]
                ind, = l_uniqueperm( ecomp )
                edge_pairs = e[ind,:]
            else
                error( "row wise unique not supported" )
                # edge_pairs = unique( e, "rows" )
            end
        end


        # Move mesh points based on edge lengths L and forces F.
        p1 = p[edge_pairs[:,1],:]
        p2 = p[edge_pairs[:,2],:]
        bars = p1 - p2                  # Bar vectors.
        L = sqrt.(sum(bars.^2,2))       # Bar lengths.
        hbars = fh( 0.5*( p1 + p2 ) )   # Rel bar mid point sizes.
        L_target = hbars*F_scale*(sum(L.^n_sdim)/sum(hbars.^n_sdim))^(1/n_sdim)   # Bar target lengths.

        # Density control, remove points that are too close to each other.
        if( mod(it,N_DCF)==0 && any(L_target.>F_DCF*L) )
            n_dcs = n_dcs + 1
            ind_del  = find( L_target .> F_DCF*L )
            ind_del  = setdiff(reshape(edge_pairs[ind,:],2*length(ind),1),1:n_p_fix)
            ind_keep = setdiff(1:n_p,ind_del)
            p   = p[ind_keep,:]
            n_p = size(p,1)
            p0  = Inf
            continue
        end

        # Compute grid point movements.
        F = L_target - L   # Scalar bar forces.
        F[find(F.<0.0)] = 0.0
        F_bar = F./L*ones(1,n_sdim).*bars
        delta_p = zeros(n_p,n_sdim)
        for i=1:n_sdim
            for j=1:size(edge_pairs,1)
                delta_p[edge_pairs[j,1],i] += F_bar[j,i]
                delta_p[edge_pairs[j,2],i] -= F_bar[j,i]
            end
        end
        delta_p[1:n_p_fix,:] = 0.0
        delta_p = dp_scale * delta_p
        p = p + delta_p


        # Move grid points outside geometry back to the boundary.
        for jt=1:N_RECV

            dist = fd( p )
            ix = dist .> 0
            ix[1:n_p_fix] = 0
            if( any(ix) )
                ix = find(ix)
                grad_dist = zeros(length(ix),n_sdim)
                for i=1:n_sdim
                    doff = zeros(1,n_sdim)
                    doff[i] = gradeps
                    dist_offset_i  = fd( p[ix,:]+ones(length(ix),1)*doff )
                    grad_dist[:,i] = ( dist_offset_i - dist[ix] )/gradeps
                end
                gradnm = sum( grad_dist.^2, 2 )
                p[ix,:] -= dist[ix]./gradnm*ones(1,n_sdim) .* grad_dist
            end
        end


        # Statistics/output.
        delta_p_max = abs( maximum( [sqrt.(sum(delta_p[find(dist.<dp_tol),:].^2,2)); -Inf] ) )
        if( mod(it,IT_PRT)==0 )
            s = @sprintf( "Iteration %4i: %i vertices, %i cells, max(delta_p) = %g\n", it, size(p,1), size(t,1), delta_p_max )
            l_message( fid, s )
        end


        # Check for convergence.
        if( (it>IT_MIN && delta_p_max<dpc_tol) || size(t,1)<=2 || it>it_max || do_break )
            if( delta_p_max<dpc_tol )
                is_converged = true
            end
            break
        end

    end


    # Clean up and check final mesh.
    (p, t, td) = l_fixmesh( p, t, fd, e_fix, dtrm_tol )
    t_tri = t_tri + td;


    # Statistics.
    t_tot = ( time_ns() - t0 )/1e9
    st = statistics( is_converged, it, n_tri, n_dcs,
               maximum(sqrt.(sum(delta_p[find(dist.<-dp_tol),:].^2,2))),
               mean(sqrt.(sum(delta_p[find(dist.<-dp_tol),:].^2,2))), t_tot, t_tri )

    s = "done"
    if( do_break )
        s = "stopped"
    end
    s = @sprintf( "Mesh generation %s: t_tot = %f s, %i iterations, %i (re-)triangulations\n",
                 s, t_tot, it, n_tri )

    l_message( fid, s )

    return (p, t, st)
end


#------------------------------------------------------------------------------#
function l_triangulate( p, fd::Function, e_fix=[], dtrm_tol=1e-4 )
    # Generate triangulation for vertices p.

    AV_TOL = eps()*1e1   # Minimum accepted absolute area/volume.

    ind_keep = find( .!isnan.(p) )
    rows, = ind2sub( size(p), ind_keep )
    p  = p[unique(rows),:]
    p, = l_deduplicate( p )

    # Generate triangulation for grid points p.
    tic()
    t = l_delaunay_triangulation( p, e_fix )
    td = toq()

    # Calculate simplex centers.
    n_sdim = size(p,2)
    pc = zeros(size(t,1),n_sdim)
    for i=1:n_sdim
        pc[:,i] = mean(reshape(p[t,i],size(t)),2)
    end

    # Remove simplices with center outside region.
    dist = fd( pc )
    t = t[find(dist.<dtrm_tol),:]

    # # Reorient simplices.
    av = l_simpvol( p, t )
    ix_flip = find(av.<0)
    t[ix_flip,[1 2]] = t[ix_flip,[2 1]]

    # Remove simplices with volume < AV_TOL.
    t = t[find(abs.(av).>=AV_TOL),:]

    if( isempty(t) )
        tic()
        t = l_delaunay_triangulation( p, e_fix )
        td = td + toq()
    end

    return (p, t, td)
end

#------------------------------------------------------------------------------#
function l_delaunay_triangulation( p, c )
    # Generate triangulation for vertices p.

    # Coordinate scaling.
    n_sdim = size(p,2)
    pp = Array{Any}(n_sdim)
    for i=1:n_sdim
        pmin = minimum(p[:,i])
        pmax = maximum(p[:,i])
        pp[i] = ((p[:,i] - pmin)/(pmax - pmin)*(max_coord - min_coord)) + min_coord
    end

    n_p = size(p,1)
    if( n_sdim==2 )
        ip = IndexedPoint2D[ IndexedPoint2D(pp[1][i],pp[2][i],i) for i in 1:n_p]
        tess = DelaunayTessellation2D{IndexedPoint2D}(n_p)
    else
        # ip = IndexedPoint3D[ IndexedPoint3D(pp[1][i],pp[2][i],pp[3][i],i) for i in 1:n_p]
        # tess = DelaunayTessellation3D{IndexedPoint3D}(n_p)
    end

    push!( tess, ip )

    # Extract triangle indices.
    n_t = 0
    ind_t = Vector{Int64}(length(tess._trigs))
    for i in 2:tess._last_trig_index
        isexternal(tess._trigs[i]) && continue
        n_t += 1
        ind_t[n_t] = i
    end
    t = Array{Int64}(n_t,n_sdim+1)
    # if( n_sdim==2 )
        for i=1:n_t
            tri = tess._trigs[ind_t[i]]
            t[i,1] = tri._a._idx
            t[i,2] = tri._b._idx
            t[i,3] = tri._c._idx
        end
    # else
    #     for i=1:n_t
    #         tri = tess._trigs[ind_t[i]]
    #         t[i,1] = tri._a._idx
    #         t[i,2] = tri._b._idx
    #         t[i,3] = tri._c._idx
    #         t[i,4] = tri._d._idx
    #     end
    # end

    return t
end

#------------------------------------------------------------------------------#
function l_fixmesh( p, t, fd::Function, e_fix, dtrm_tol )
    # FIXMESH Remove duplicated/unused nodes and fix element orientation.

    if( isempty(p) || isempty(t) )
        return (p, t, 0.0)
    end

    P_TOL = eps()*1024.0
    (p, ix, ind_p_orig) = l_deduplicate( p, P_TOL )


    # if( !isempty(fd) )

        t = ind_p_orig[t]

        # Final triangulation.
        (p, t, td) = l_triangulate( p, fd, e_fix, dtrm_tol )

        # Calculate simplex centers.
        n_sdim = size(p,2)
        pc = zeros(size(t,1),n_sdim)
        for i=1:n_sdim
            pc[:,i] = mean(reshape(p[t,i],size(t)),2)
        end

        # Remove simplices with center outside region.
        dist = fd( pc )
        ind  = find(dist.<dtrm_tol)
        if( length(ind)<size(t,1) )
            t = t[find(dist.<dtrm_tol),:]
        end

        # Remove unused nodes.
        ix1,jx1 = l_uniqueperm( t )
        ind_p = t[ix1]
        t = reshape( jx1, size(t) )
        p = p[ ind_p, : ]
        ind_p = ix[ ind_p ]

    # end

    return (p, t, td)
end

#------------------------------------------------------------------------------#
function l_initpoints( h0, bbox )
    # Create an initial grid point/vertex distribution.

    n_sdim = size(bbox,2)
    pinit  = Array{Any}(n_sdim)
    for i=1:n_sdim
        if( n_sdim==2 && i==2 )
            pinit[i] = bbox[1,i]:h0*sqrt(3)/2:bbox[2,i]
        else
            pinit[i] = bbox[1,i]:h0:bbox[2,i]
        end
    end
    pp = Array{Any}(n_sdim)
    n  = length(pinit[1])
    m  = length(pinit[2])
    xh = reshape(pinit[1], n, 1)
    yh = reshape(pinit[2], 1, m)
    if( n_sdim==2 )
        pp[1] = repmat(xh, 1, m)
        pp[2] = repmat(yh, n, 1)
    else
        o = length(pinit[3])
        zh = reshape(pinit[3], 1, o)
        pp[1] = repmat(xh, 1, m*o)
        pp[2] = repmat(yh, n, 1)
        pp[2] = repmat(pp[2][:], o, 1)
        pp[3] = repmat(zh, n*m, 1)
    end
    if( n_sdim==2 )
        pp[1][:,2:2:end] = pp[1][:,2:2:end] + h0/2
    end
    p = Array{Float64}(prod(size(pp[1])),n_sdim)
    for i=1:n_sdim
        p[:,i] = pp[i][:]
    end

    return p
end

#------------------------------------------------------------------------------#
function l_simpvol( p, t )
    # SIMPVOL Triangle and tetrahedron volumes.

    n_sdim = size(p,2)
    if( n_sdim==2 )
        d12 = p[t[:,2],:] - p[t[:,1],:]
        d13 = p[t[:,3],:] - p[t[:,1],:]
        v = 0.5*( d12[:,1].*d13[:,2] - d12[:,2].*d13[:,1] )
    elseif( n_sdim==3 )
        d12 = p[t[:,2],:] - p[t[:,1],:]
        d13 = p[t[:,3],:] - p[t[:,1],:]
        d14 = p[t[:,4],:] - p[t[:,1],:]
        v = dot( cross(d12,d13,2), d14, 2 )/6
    else
        v = []
    end

    return v
end

#------------------------------------------------------------------------------#
function l_uniqueperm( a )
    # Unique index permutation vector.

    ind = sortperm( a[:] )
    tmp = a[ind]
    msk = [ true; diff(tmp).!=0 ]

    ind1 = ind[ find(msk) ]
    ind2 = Vector{Int}(length(ind))
    ind2[ind] = cumsum(msk)

    return (ind1, ind2)
end

#------------------------------------------------------------------------------#
function l_deduplicate( a, s=[] )
    # Deduplicate array rows.

    if( isempty(a) )
        return (a, 0, 0)
    end
    if( isempty(s) )
        s = 1e-6*maximum(maximum(a)-minimum(a))
    end

    c = s*round.(a/s)

    ix = 1:size(c,1)
    c = [c ix[:]]
    c = sortrows( c )
    k = trunc.(Int,c[:,end])
    c = c[:,1:end-1]

    ix = any( c[1:size(c,1)-1,:] .!= c[2:size(c,1),:], 2 )
    j  = Vector{Int}(length(k))
    j[k] = cumsum( [true;ix] )
    i  = k[ [1;find(ix)+1] ]

    jj = sortperm(i)
    i  = i[jj]
    jjinv = Vector{Int}(length(jj))
    jjinv[jj] = 1:length(jj)
    j  = vec(jjinv[j])

    b = a[i,:]

    return (b, i, j)
end

#------------------------------------------------------------------------------#
function l_message( fid, s )
    # Message output.

    if( length(fid)==1 && isinteger(fid) && fid>0 )
        if( !any(Int64(s[end]).==[10,13]) )
            s = s * string(Char(10))
        end
        if( isa(fid,typeof(IOBuffer())) )
            print( fid, s )
        else
            print( s )
        end
    else#if( isa(fid,"Function") )
        # fid( s )
    end

end

#------------------------------------------------------------------------------#
function plotgrid( p, t )
    # Plot triangulation.

    n_t = size(t,1)
    x = Array{Real}(0);
    y = Array{Real}(0);
    for i=1:n_t
        x = [x; NaN; p[t[i,:],1]]
        y = [y; NaN; p[t[i,:],2]]
    end

    plot( Shape(x,y), opacity=0.5 )
end


# end # module DistMesh
