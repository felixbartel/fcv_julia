export fcv_t_equispaced, compute

using FFTW, LinearAlgebra


struct fcv_t_equispaced
  d::Int
  f::Vector{Complex{Float64}}
  W::Real
  F

  function fcv_t_equispaced(d::Int, f::Vector{T}) where T <: Number
    FFTW.set_num_threads(4)
    M = length(f)
    tmp = ntuple(x -> Int(M^(1/d)), d)
    F = plan_fft(similar(reshape(f, tmp)), flags = FFTW.MEASURE)
    this = new(d, f, 1/length(f), F)
  end
end


function H(fcv::fcv_t_equispaced, What::Array{Float64})
# computes the matrix-vector product with the hat matrix F*inv(F'*W*F+What)*F'*W
  M = length(fcv.f)
  tmp = ntuple(x -> Int(M^(1/fcv.d)), fcv.d)
  What = reshape(What, tmp)
  f_matrix = reshape(fcv.f, tmp)
  fhat_r = M*(fcv.W.*f_matrix |> fftshift |> ifft) |> fftshift
  fhat_r ./= (What.+1)
  f_r =  fhat_r |> fftshift |> fft |> fftshift
  return vec(f_r), vec(fhat_r)
end


function diagonals(fcv::fcv_t_equispaced, What::Array{Float64})
  h = fcv.W*(sum(1 ./(1 .+What)))
end


function compute(fcv::fcv_t_equispaced, What::Array{Float64})
  M = length(fcv.f)
  h = diagonals(fcv, What)
  f_r, fhat_r = H(fcv, What)
  cv = 1/M*norm((fcv.f-f_r)./(1 .- h))^2
  return (gcv = cv, ocv = cv, f_r = f_r, fhat_r = fhat_r)
end
