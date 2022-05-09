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
end
#------------------------------------------------------------------------------#
