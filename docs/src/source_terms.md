# User-Provided Source Terms

HallThruster allows users to provide additional source terms if they want. We make use of this function internally when we perform order verification studies (see [Verification](@ref), as well as our order verification tests in /tests/order_verification.jl). This may also be useful when implementing additional physics. Users may provide seperate source terms for each of the solved equations. The corresponding fields in the `Config` struct are

- `source_neutrals`
- `source_potential`
- `source_electron_energy`
- `source_ion_continuity`
- `source_ion_momentum`

The first three act similarly. HallThruster expects a function which takes `(U, params, i)` as an argument, where `U` is the state matrix, `params` is the NamedTuple of paramters, and `i` is the cell index. The source term should then return a scalar, which is added to the right-hand side of the corresponding equation.

As we allow for multiple ion charge states, `source_ion_continuity` and `source_ion_momentum` should be of type `Tuple` or `Vector` and contain one function for each charge state, with the same expected signature and return type as above.

!!! warning "Concrete vs abstract types"
    For performance, it is important that Julia is able to infer the return type of your source term. Therefore, you must ensure that the values passed to `source_ion_continuity` and `source_ion_momentum` are concretely-typed. For example, this would not be ideal:

    ```julia
    source_ion_momentum = [
        (U, params, i) -> 2.0,
        (U, params, i) -> 3.0,
        (U, params, i) -> -1
    ]
    ```

    because the type of this term is not concrete:

    ```julia-repl
    julia> typeof(source_ion_momentum)
    Vector{Function} (alias for Array{Function, 1})

    julia> isconcretetype(ans)
    false
    ```

    Better would be to use FunctionWrappers or to implement a callable singleton type, as shown below:

    ```julia
    struct MySource
        number::Float64
    end

    (s::MySource)(U, params, i) = s.number

    source_ion_momentum_concrete = [
        MySource(2.0),
        MySource(3.0),
        MySource(-1),
    ]
    ```

    This has identical behavior to the first example, but is a concrete type:

    ```julia-repl
    julia> typeof(source_ion_momentum_concrete)
    Vector{MySource} (alias for Array{MySource, 1})

    julia> isconcretetype(ans)
    true
    ```