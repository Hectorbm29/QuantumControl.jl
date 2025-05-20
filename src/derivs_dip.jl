using QuantumPropagators.Generators: Operator, evaluate
using QuantumPropagators.Generators_dip: Generator_dip, _make_generator_dip, evaluate


"""
Get the derivative of the generator ``G_dip`` w.r.t. the control ``ϵ(t)``.

```julia
μ  = get_control_deriv(generator_dip, control)
```

returns `nothing` if the `generator` (Hamiltonian or Liouvillian) does not
depend on `control`, or a generator

```math
μ = \frac{∂G_dip}{∂ϵ(t)}
```

otherwise. The generator is a `Generator_dip` object, which is a generalization of
the `Generator` object to include the dressing operators. To achieve this, the analytical
derivatives of the dressing operators with respect to the amplitudes are required. 
"""

function get_control_deriv(generator::Generator_dip, control)
    terms = []
    drift_offset = length(generator.ops) - length(generator.dresses)
    amplitudes = generator.amplitudes

    if generator.dresses_derivatives == nothing
        return nothing
    end

    for (i, dress) in enumerate(generator.dresses) # iteration over the dresses (different operators)
        constant_terms = 0.0
        function_terms = Any[]

        for (j, ampl) in enumerate(generator.amplitudes) # Each dress depends on all the amplitudes

            ∂d╱∂a = get_dress_deriv(generator, dress, ampl)
            if ∂d╱∂a == 0.0
                continue
            end

            ∂a╱∂ϵ = get_control_deriv(ampl, control)
            if ∂a╱∂ϵ == 0.0
                continue
            end

            # If both are numbers we can just multiply them and add them to the constant terms
            if typeof(∂d╱∂a) <: Number && typeof(∂a/∂ϵ) <: Number
                constant_terms += ∂d╱∂a * ∂a╱∂ϵ
                continue
            elseif typeof(∂d╱∂a) <: Number # if not we need to create a callable constant
                ∂d╱∂a = create_callable_constant(∂d╱∂a)
            elseif typeof(∂a╱∂ϵ) <: Number
                ∂a╱∂ϵ = create_callable_constant(∂a╱∂ϵ)
            end

            chan_rule_deriv = mulfunctions(∂d╱∂a, ∂a╱∂ϵ)

            push!(function_terms, chan_rule_deriv)

        end

        if constant_terms != 0.0
            push!(terms, constant_terms * generator.ops[i+drift_offset])
        end

        if length(function_terms) > 0
            dress_op = generator.ops[i+drift_offset]
            if length(function_terms) == 1
                dress_terms = function_terms[1]

            else
                dress_terms = addfunctions(function_terms)
            end
            op = generator.ops[i+drift_offset]
            push!(terms, (op , dress_terms))
        end

    end


    if length(terms) == 0
        return nothing
    else
        return _make_generator_dip(terms...; ampl_vec=amplitudes)
    end
end

function mulfunctions(a, b)
    c = (args...) -> a(args...) * b(args...)
    return c 
end

function __addfunctions(a, b)
    c = (args...) -> a(args...) + b(args...)
    return c
end

function addfunctions(vector_of_functions)

    c = __addfunctions(vector_of_functions[1], vector_of_functions[2])
    if length(vector_of_functions) == 2
        return c
    end
    for i in 3:length(vector_of_functions)
        c = __addfunctions(c, vector_of_functions[i])
    end
    return c
end

function create_callable_constant(a)
    c = (args...) -> a
    return c
end

function get_dress_deriv(generator::Generator_dip, dress, ampl)

    dresses = generator.dresses
    ampls = generator.amplitudes
    d_index = findfirst(x -> x === dress, dresses)
    a_index = findfirst(x -> x === ampl, ampls)

    d_derivs = generator.dresses_derivatives 

    if isnothing(d_derivs)
        error("The generator does not have analytical derivatives for the dresses\n")
    else
        return d_derivs[d_index, a_index]
    end

end