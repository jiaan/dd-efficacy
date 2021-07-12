include("fidelity.jl")
using Plot


#= run(example, η, measure, condition, rand_states, n_samples=1000) =#

# example: model_random, model_fixed
# η: noise in pulses
# measure: fidelity, tr_distance, max_tr_distance
# condition: F_condition, A_condition, T_condition,  A_condtion_tr 
# rand_states: rand_pure_state, rand_mixed_state


Δ = run(model_random, 0, fidelity, F_condition, rand_pure_state, 1)
Δ = reverse(Δ,dims=1) # flip y-axis up down 

imshow(Δ)
clim([0,1e-15])
show()
