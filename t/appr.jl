export fcv_t_appr, compute

#push!(LOAD_PATH, "../../nfft/julia/nfft/")
using NFFT, LinearAlgebra, IterativeSolvers, LinearOperators


struct fcv_t_appr
  nodes::Array{Float64}
  f::Vector{Complex{Float64}}
  W::Vector{Complex{Float64}}
  p # NFFT plan
  M # number of nodes
  N # number of frequencies

  function fcv_t_appr(nodes::Array{T}, f::Vector{T}, W::Vector{T}, N::Integer) where T <: Number
    p = Plan(ntuple(x -> Int32(N), size(nodes, 2)), size(nodes, 1))
    p.x = nodes
    this = new(nodes, f, W, p, length(f), N)
  end
end


function H(fcv::fcv_t_appr, What::Array{Float64})
# computes the matrix-vector product with the hat matrix F*inv(F'*W*F+What)*F'*W

  M = LinearOperator(fcv.M+fcv.N, fcv.N, false, false,
    fhat -> (fcv.p.fhat = fhat; NFFT.trafo(fcv.p); [sqrt.(fcv.W).*fcv.p.f; sqrt.(What).*fhat]),
    nothing,
    f -> (fcv.p.f = f[1:fcv.M]; NFFT.adjoint(fcv.p); sqrt.(conj(fcv.W)).*fcv.p.fhat+sqrt.(conj(What)).*f[fcv.M+1:end])
    )
  fhat_r = lsqr(M, [sqrt.(fcv.W).*fcv.f; zeros(fcv.N)])
  fcv.p.fhat = fhat_r
  NFFT.trafo(fcv.p)
  return fcv.p.f, fhat_r
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
