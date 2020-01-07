export fcv_t_appr, compute

using NFFT, LinearAlgebra, IterativeSolvers, LinearOperators


struct fcv_t_appr
  d::Int
  nodes::Array{Float64}
  f::Vector{Complex{Float64}}
  W::Vector{Complex{Float64}}
  p # NFFT plan
  M # number of nodes
  N # number of frequencies
  x0::Vector{Complex{Float64}} # inital solution for lsqr

  function fcv_t_appr(nodes::Array{T}, f::Vector{T}, W::Union{Vector{T}, Nothing}, N::Integer) where T <: Number
    d = size(nodes, 2)
    p = Plan(ntuple(x -> Int(N), d), size(nodes, 1))
    if d == 1
      p.x = nodes
    else
      p.x = Matrix(transpose(nodes))
    end
    # p.num_threads = 4 # currently not supported :-(
    isnothing(W) && ( W = voronoiArea(nodes) )
    this = new(d, nodes, f, W, p, length(f), N, zeros(N^d))
  end
end


function H(fcv::fcv_t_appr, What::Array{Float64})
# computes the matrix-vector product with the hat matrix F*inv(F'*W*F+What)*F'*W
  M = LinearOperator(fcv.M+fcv.N^fcv.d, fcv.N^fcv.d, false, false,
    fhat -> (fcv.p.fhat = fhat; NFFT.trafo(fcv.p); [sqrt.(fcv.W).*fcv.p.f; sqrt.(What).*fhat]),
    nothing,
    f -> (fcv.p.f = sqrt.(conj(fcv.W)).*f[1:fcv.M]; NFFT.adjoint(fcv.p); fcv.p.fhat+sqrt.(conj(What)).*f[fcv.M+1:end])
    )
  lsqr!(fcv.x0, M, [sqrt.(fcv.W).*fcv.f; zeros(fcv.N^fcv.d)], maxiter = 20)
  fcv.p.fhat = fcv.x0
  NFFT.trafo(fcv.p)
  return fcv.p.f, fcv.x0
end


function diagonals(fcv::fcv_t_appr, What::Array{Float64})
# this is approximative for non-quadrature rules!
  h = fcv.W*(sum(1 ./(1 .+What)))
end


function compute(fcv::fcv_t_appr, What::Array{Float64})
  M = length(fcv.f)
  h = diagonals(fcv, What)
  f_r, fhat_r = H(fcv, What)
  cv = 1/M*norm((fcv.f-f_r)./(1 .- h))^2
  return (gcv = cv, ocv = cv, f_r = f_r, fhat_r = fhat_r)
end
