using LinearAlgebra
using ProgressMeter
using Distributed
using Random
import Base:run
import QuantumInformation: KrausOperators, HaarKet, HilbertSchmidtStates
import QuantumInformation: rand, fidelity, trace_distance, diamond_distance
import QuantumInformation: ptrace

⊗ = kron

# Random pure state (density matrix)
function rand_pure_state(d)
  h = HaarKet{2}(d)
  ψ = rand(h)
  ρS = ψ * ψ'
end

# Random mixed state (density matrix)
function rand_mixed_state(d)
  h = HilbertSchmidtStates(d)
  ρS = rand(h)
end

# Ground state density matrix
function ground_state_dm(ham)
  @assert size(ham)[1] < 1000 "Hilbert space restricted to 1000"
  ψ = eigvecs(ham)[:,1]
  ρ = ψ * ψ'
end

# Random unit vector on hypersphere
function rand_unit_vector(n)
  unit_vec = randn(n)
  unit_vec/norm(unit_vec)
end

function get_pauli()
  id2 = Matrix{ComplexF64}(I,2,2)
  σx = Matrix{ComplexF64}([0 1;1 0])
  σy = Matrix{ComplexF64}([0 -1im; 1im 0])
  σz = Matrix{ComplexF64}([1 0; 0 -1])
  id2, σx, σy, σz
end

function LinearAlgebra.normalize!(op::Matrix,norm_type=opnorm)
  op .= op / norm_type(op)
  return op
end

function LinearAlgebra.normalize(op::Matrix,norm_type=opnorm)
  op = op / norm_type(op)
  return op
end

# Model Hamiltonian with random seed
function model(α, β, seed::Int)
  m = MersenneTwister(seed)
  σ = get_pauli()
  #= σ = Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1) =# 
  normalize!.(σ, norm)

  # Bath operators
  B = Vector{Matrix}(undef,4)
  B[1] = sum(rand_unit_vector(4) .* σ[1:4])
  normalize!(B[1], norm)
  for i = 2:4
    B[i] = sum(rand_unit_vector(4) .* σ[1:4])
  end

  H_B = σ[1] ⊗ B[1] 
  H_SB = sum(σ[i] ⊗ B[i]  for i = 2:4)
  #= normalize!(H_B, norm) # FB-norm =#  
  normalize!(H_SB, norm)
  H_tot = α * H_B + β * H_SB
end


# Model Hamiltonian
function model(α, β)
  σ = get_pauli()
  #= σ = Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1) =# 
  normalize!.(σ, norm)

  # Bath operators
  B = Vector{Matrix}(undef,4)
  B[1] = sum( rand_unit_vector(4) .* σ[1:4])
  normalize!(B[1], norm)
  for i = 2:4
    B[i] = sum( rand_unit_vector(4) .* σ[1:4])
  end

  H_B = σ[1] ⊗ B[1] 
  H_SB = sum(σ[i] ⊗ B[i]  for i = 2:4)
  #= normalize!(H_B, norm) # FB-norm =#  
  normalize!(H_SB, norm)
  H_tot = α * H_B + β * H_SB
end


# Model Hamiltonian with given bath operators
function model(α, β, B)
  σ = get_pauli()
  #= σ = Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1) =# 
  normalize!.(σ, norm)
  H_B = σ[1] ⊗ B[1] 
  H_SB = sum(σ[i] ⊗ B[i]  for i = 2:4)
  #= normalize!(H_B, norm) =# 
  normalize!(H_SB, norm)
  H_tot = α * H_B + β * H_SB
end

## Free evolution kernel
function Uf(H_tot, ρS, ρB)
  U = exp(-im * H_tot)
  ρtot = U * (ρS ⊗ ρB) * U'
  #= ptrace(ρtot, 2) =# # Switch to QuantumInformation.jl
  ptrace(ρtot, [2,2], 2) #Remove if diamond norm is used 
end

## DD sequence map
mutable struct DDchannel
  ham
  P
  Q
  f
  opL
  opR
  function DDchannel(ham, P, Q)
    f = exp(-im * ham)
    opL = reduce( *, map(x -> Q * x * f, P))
    opR = reduce( *, map(x -> f' * x' * Q', reverse(P)))
    return new(ham, P, Q, f, opL, opR)
  end 
end

function (channel::DDchannel)(ρS, ρB)  
  opL, opR = channel.opL, channel.opR
  ρtot = opL * (ρS ⊗ ρB) * opR
  #= ptrace(ρtot, 2) =# # Switch to QuantumInformation.jl
  ptrace(ρtot, [2,2], 2) #Remove if diamond norm is used 
end

# Switch to QuantumInformation.jl
#= function fidelity(ρ1,ρ2) =#
#=   ρ1_sqrt = √ρ1 =#
#=   F = tr(√(ρ1_sqrt * ρ2 * ρ1_sqrt))^2 =#
#=   imag(F) < 1e-12 && return real(F) =#
#= end =#

#= function fidelity_pure(ρ1,ρ2) =#
#=   F = tr(ρ1 * ρ2) =#
#=   imag(F) < 1e-12 && return real(F) =#
#= end =#

# Fidelity wrapper
function fidelity(channel::Function, ρ0)
  #= F = fidelity_pure(channel(ρ0), ρ0) =#
  F = fidelity(channel(ρ0), ρ0)
  imag(F) < 1e-12 && return real(F)
end

function tr_distance(ρ1, ρ2)
  Δ = ρ1 - ρ2
  0.5 * tr(√(Δ^2)) 
end

function tr_distance(channel::Function, ρ0)
  dist = tr_distance(channel(ρ0),ρ0)
  imag(dist) < 1e-12 && return real(dist)
end

function max_tr_distance(channel::Function, ρ)
  dist_ = 0
  for i = 1:100
    ρ = rand_pure_state(2)
    dist = trace_distance(channel(ρ), ρ) #Switch to QuantumInformation.jl for trace distance
    dist > dist_ && (dist_ = dist)
  end
  return dist_
end

# Check the warning from Convex.jl
function diamond_norm(channel::Function, ρ)
  su, sp, sm, sd = [1 0; 0 0], [0 1; 0 0], [0 0;1 0], [0 0; 0 1]
  basis1 = (su, sp, sm, sd)
  basis2 = (su, sm, sp, sd)
  sup_op = map2sup(channel, basis1, basis2) 
  #= println(sup_op) =#
  kraus_ops = to_kraus(sup_op)
  #= println(kraus_ops) =#
  c1 = KrausOperators([kraus_ops[1], kraus_ops[2], kraus_ops[3], kraus_ops[4]])
  c2 = KrausOperators([[1.0 0im;0im 1.0]])
  #= println(c1) =#
  #= error("Stop") =# 
  diamond_norm = diamond_distance(c1, c2)
end

function unitary_noise(η)
  σ = get_pauli()
  #= σ = [Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1)] =#
  K1 = sum(rand(Float64,4) .* σ)
  normalize!(K1, norm)
  Q = exp(-1im * η * K1) ⊗ σ[1] 
end

function F_condition(F_free, F_dd) 
  max_free = maximum(F_free)
  max_dd = maximum(F_dd)
  min_dd = minimum(F_dd)
  dF = min_dd - max_free
  return dF, max_free, min_dd
end

function T_condition(F_free, F_dd) 
  max_free = maximum(F_free)
  min_free = minimum(F_free)
  max_dd = maximum(F_dd)
  min_dd = minimum(F_dd)
  dF = min_free - max_dd 
  return dF, max_free, max_dd
end

function maxmax_condition(F_free, F_dd) 
  max_free = maximum(F_free)
  max_dd = maximum(F_dd)
  dF = max_free - max_dd 
  return dF, max_free, max_dd
end

function A_condition(F_free, F_dd)
  max_free = maximum(F_free)
  min_free = minimum(F_free)
  max_dd = maximum(F_dd)
  min_dd = minimum(F_dd)
  if all(F_dd .>= F_free)
    dF = 1
  else 
    dF = 0 
  end
  return dF, max_free, max_dd
end

function A_condition_map(F_free, F_dd)
  max_free = maximum(F_free)
  min_free = minimum(F_free)
  max_dd = maximum(F_dd)
  min_dd = minimum(F_dd)
  if all( sqrt.(1 .- F_dd .+ 1e-7) .>= sqrt.(1 .- F_free .+ 1e-7) )
    dF = 1
  else 
    dF = 0 
  end
  return dF, max_free, max_dd
end


function A_condition_tr(F_free, F_dd)
  max_free = maximum(F_free)
  min_free = minimum(F_free)
  max_dd = maximum(F_dd)
  min_dd = minimum(F_dd)
  if all(F_dd .<= F_free)
    dF = 1
  else 
    dF = 0 
  end
  return dF, max_free, max_dd
end


function model_fixed(α, β ,η, measure, condition, rand_states, n_samples = 10000)
  ## Define the Hamiltonian
  ham = model(α, β, 527)

  F_free = Vector{Float64}(undef,n_samples)
  F_dd = Vector{Float64}(undef,n_samples)
  # Bath state
  ρB= Matrix{Float64}(I,2,2)/2.0 

  for i = 1:n_samples
    σ = Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1) 
    P = map(x-> x ⊗ σ[1], [σ[4],σ[2],σ[4],σ[2]])
    Q = unitary_noise(η) 
    XZXZ = DDchannel(ham, P, Q)

    # System state
    ρS = rand_states(2)

    Ufree(ρS) = Uf(ham, ρS, ρB)
    F_free[i]= measure(Ufree, ρS) 
    Udd(ρS) = XZXZ(ρS, ρB)
    F_dd[i]= measure(Udd, ρS) 
  end

  return condition(F_free, F_dd)
end




function model_random(α, β ,η, measure, condition, rand_states, n_samples = 10000 )

  F_free = Vector{Float64}(undef,n_samples)
  F_dd = Vector{Float64}(undef,n_samples)
  # Bath state
  #= ρB= rand_dm(2,10) =#
  ρB= Matrix{Float64}(I,2,2)/2.0 
  for i = 1:n_samples
    ham = model(α, β)
    σ = get_pauli()
    #= σ = Matrix{ComplexF64}(I,2,2), σx(1,1), σy(1,1), σz(1,1) =# 
    P = map(x-> x ⊗ σ[1], [σ[4],σ[2],σ[4],σ[2]]) #  Z<-X<-Z<-X 
    #= Q = I =# 
    Q = unitary_noise(η) 
    XZXZ = DDchannel(ham, P, Q)

    # System state
    ρS = rand_states(2)


    Ufree(ρS) = Uf(ham, ρS, ρB)
    F_free[i]= measure(Ufree, ρS) 
    Udd(ρS) = XZXZ(ρS, ρB)
    F_dd[i]= measure(Udd, ρS) 
  end

  return condition(F_free, F_dd)
end


function run(example, η, measure, condition, rand_states, n_samples=1000)
  α = collect(0:0.01:0.6)
  β = collect(0.01:0.01:0.6)
  dF = zeros(length(β),length(α)) 
  F_free = zeros(length(β),length(α)) 
  F_dd = zeros(length(β),length(α)) 
  @showprogress for i = 1:length(β)
    for j = 1:length(α)
      dF[i,j],F_free[i,j],F_dd[i,j] = example(α[j],β[i],η, measure, condition, rand_states, n_samples)
    end
  end
  return dF
end


