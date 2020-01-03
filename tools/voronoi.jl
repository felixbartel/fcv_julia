function voronoiArea(nodes)
  if size(nodes, 2) == 1
    p = sortperm(nodes)
    nodes_sorted = nodes[p]
    w = ([nodes_sorted[2:end]; nodes[1]+1]-[nodes_sorted[end]-1; nodes_sorted[1:end-1]])./2
    return w[invperm(p)]
  else
    println("Voronoi decomposition not implemented for this dimension")
  end
end
