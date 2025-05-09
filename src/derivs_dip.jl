using QuantumPropagators.Generators: Operator, evaluate
using QuantumPropagators.Generators_dip: Generator_dip, _make_generator_dip, evaluate
using QuantumPropagators.Amplitudes: LockedAmplitude, ShapedAmplitude



@doc raw"""
Get the derivative of the generator ``G`` w.r.t. the control ``ϵ(t)``.

```julia
μ  = get_control_deriv(generator, control)
```

returns `nothing` if the `generator` (Hamiltonian or Liouvillian) does not
depend on `control`, or a generator

```math
μ = \frac{∂G}{∂ϵ(t)}
```

otherwise. For linear control terms, `μ` will be a static operator, e.g. an
`AbstractMatrix` or an [`Operator`](@ref). For non-linear controls, `μ` will be
time-dependent, e.g. a [`Generator`](@ref). In either case,
[`evaluate`](@ref) should be used to evaluate `μ` into a constant operator
for particular values of the controls and a particular point in time.

For constant generators, e.g. an [`Operator`](@ref), the result is always
`nothing`.
"""

"""
function get_control_deriv(generator::Generator, control)
    terms = []
    drift_offset = length(generator.ops) - length(generator.amplitudes)
    for (i, ampl) in enumerate(generator.amplitudes)
        ∂a╱∂ϵ = get_control_deriv(ampl, control)
        if ∂a╱∂ϵ == 0.0
            continue
        elseif ∂a╱∂ϵ == 1.0
            mu_op = generator.ops[i+drift_offset]
            push!(terms, mu_op)
        else
            mu_op = generator.ops[i+drift_offset]
            push!(terms, (mu_op, ∂a╱∂ϵ))
        end
    end
    if length(terms) == 0
        return nothing
    else
        return _make_generator(terms...)
    end
end
"""

#TODO build this function 
function get_control_deriv(generator::Generator_dip, control)
    terms = []
    drift_offset = length(generator.ops) - length(generator.dresses)
    amplitudes = generator.amplitudes
    dresses = generator.dresses


    for (j, dress) in enumerate(generator.dresses)

        #derivfunction = [] #???

        for (i, ampl) in enumerate(generator.amplitudes)

            ∂a╱∂ϵ = get_control_deriv(ampl, control)
            if ∂a╱∂ϵ == 0.0
                continue
            else
                ∂dj╱∂ai =  get_dress_deriv(dress, ampl)

                #derivfuction += ∂a╱∂ϵ * ∂dj╱∂ai
            end


        end
    end


    
    for (i, ampl) in enumerate(generator.amplitudes)

        for (j, dress) in enumerate(generator.dresses)

            ∂dj╱∂ai = get_dress_deriv(dress, ampl)

        end

        ∂a╱∂ϵ = get_control_deriv(ampl, control)
        if ∂a╱∂ϵ == 0.0
            continue
        elseif ∂a╱∂ϵ == 1.0
            mu_op = generator.ops[i+drift_offset]
            push!(terms, mu_op)
        else
            mu_op = generator.ops[i+drift_offset]
            push!(terms, (mu_op, ∂a╱∂ϵ))
        end
    end
    if length(terms) == 0
        return nothing
    else
        return _make_generator_dip(terms...; ampl_vec=amplitudes)
    end
end


function get_dress_deriv(dress, ampl)
    return 0.0
end