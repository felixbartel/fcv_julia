using VoronoiCells

function voronoiArea(nodes)
  d = size(nodes, 2)
  M = size(nodes, 1)
  if d == 1
    p = sortperm(nodes)
    nodes_sorted = nodes[p]
    w = ([nodes_sorted[2:end]; nodes_sorted[1]+1]-[nodes_sorted[end]-1; nodes_sorted[1:end-1]])./2
    return w[invperm(p)]
  elseif d == 2
    shift = -ones(d)
    nodes_p = Array{Float64}(undef, M*3^d, d)
    for m = 1:3^d
      nodes_p[(m-1)*M+1:m*M,:] = [ nodes[m2, i]+shift[i] for m2 in 1:size(nodes, 1), i in 1:d ]
      shift[end]+=1
      for idx in d:-1:2
        if shift[idx] == 2
          shift[idx-1] += 1
          shift[idx] = -1
        end
      end
    end
    w = voronoiarea(nodes_p[:, 1], nodes_p[:, 2], [-3/2, 3/2, -3/2, 3/2])
    m = Int((3^d+1)/2)
    return w[(m-1)*M+1:m*M]
  else
    println("Voronoi decomposition not implemented for this dimension")
  end
end
