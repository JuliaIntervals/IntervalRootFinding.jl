
import IntervalArithmetic: diam, isinterior
import Base: iterate, eltype, IteratorSize, copy

export branch_and_prune, Bisection, Newton, RootSearch, SearchStrategy
export start, next, done, copy, step!, eltype, iteratorsize

diam(x::Root) = diam(x.interval)

"""
    SearchStrategy(constructor, store!, retrieve!)

Type describing the chosen strategy determining the order in which the
intervals are processed. Given a type, the `constructor` function must return
a container for elements of that type on which `store!` and `retrieve!` act.
The function `store!(container, element)` must add `element` to `container` and
`retrieve!(container)` must return the next element to be processed and delete
it from `set`.

The container type returned by `constructor` must support a `deepcopy` function
to allow to copy the `RootSearchState` generated.
"""
struct SearchStrategy
    constructor::Function
    store!::Function
    retrieve!::Function
end

"""
    function SearchStrategy(CONTAINER::Type, push!::Function, retrieve!::Function)

Constructor for a `SearchStrategy` constructing a container for elements of type
`EL` with the sythax `CONTAINER{EL}()`.
"""
function SearchStrategy(CONTAINER::Type, push!::Function, retrieve!::Function)
    SearchStrategy(EL -> CONTAINER{EL}(), push!, retrieve!)
end

SearchStrategy() = SearchStrategy(Vector, push!, pop!)

"""
    RootSearch{R <: Union{Interval,IntervalBox}, S <: Contractor, T <: Real}

Type implementing the `Base.Iterator` interface to the branch and prune routine.
Returns the `RootSearchState` at each iteration. Note: Each iteration mutates
the `RootSearchState`. Use `copy(state::RootSearchState)` to create an
independent instance if necessary.
"""
struct RootSearch{R <: Union{Interval,IntervalBox}, C <: Contractor, S <: SearchStrategy, T <: Real}
    region::R
    contractor::C
    strategy::S
    tolerance::T
end

function RootSearch(region::R, contractor::C, tol::T) where {R <: Union{Interval,IntervalBox}, C <: Contractor, T <: Real}
    RootSearch(region, contractor, SearchStrategy(), tol)
end

eltype(::Type{RS}) where {RS <: RootSearch} = RootSearchState # Warning: Not a concrete type
iteratorsize(::Type{RS}) where {RS <: RootSearch} = Base.SizeUnknown()


struct RootSearchState{V, VR}
    working::V  # Should be a container of the form  CONT{T}
    outputs::VR  # Should ba a container of root of the form CONT{Root{T}}
end

function RootSearchState(rs::RootSearch)
    return RootSearchState(rs.region, rs.strategy)
end

function RootSearchState(region::R, strat::S) where {R <: Union{Interval,IntervalBox}, S <: SearchStrategy}
    working = strat.constructor(R)
    outputs = strat.constructor(Root{R})
    strat.store!(working, region)
    return RootSearchState(working, outputs)
end

function RootSearchState(region::T) where {T<:Union{Interval,IntervalBox}}
    working = [region]
    outputs = Root{T}[]

    sizehint!(working, 1000)
    sizehint!(outputs, 100)

    RootSearchState(working, outputs)
end

copy(state::RootSearchState) =
    RootSearchState(deepcopy(state.working), deepcopy(state.outputs))


function start(iter::RootSearch)
    state = RootSearchState(iter)
    return state
end


"""
    step!(state::RootSearchState, contractor, tolerance)

Progress `state` by treating one of its `working` regions. Note: `state.working`
is always modified. If a root is found, it is added to `state.outputs`.
"""
function step!(state::RootSearchState, contractor, searchstrat, tolerance)
    X = searchstrat.retrieve!(state.working)
    status, output = contractor(X, tolerance)
    if status == :empty
        return nothing
    elseif status == :unique
        searchstrat.store!(state.outputs, Root(output, :unique))
    elseif diam(output) < tolerance
        searchstrat.store!(state.outputs, Root(output, :unknown))
    else # branch
        X1, X2 = bisect(X)
        searchstrat.store!(state.working, X1, X2)
    end
    return nothing
end

function next(iter::RootSearch, state::RootSearchState)
    step!(state, iter.contractor, iter.strategy, iter.tolerance)
    return state, state
end


# done(iter::RootSearch, state::RootSearchState) = isempty(state.working)

"""
    branch_and_prune(X, contractor, tol=1e-3)

Generic branch and prune routine for finding isolated roots using the `contract`
function as the contractor.

Inputs:
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that determines the status of a given box `X`. It
    returns the new box and a symbol indicating the status. Current possible
    values are of type `Bisection` or `Newton`

"""
function branch_and_prune(X, contractor, tol=1e-3)
    iter = RootSearch(X, contractor, tol)
    local output
    # complete iteration
    for state in iter
        output = state.outputs
    end
    return output
end

export recursively_branch_and_prune

function recursively_branch_and_prune(h, X, contractor=BisectionContractor, final_tol=1e-14)
    tol = 2
    roots = branch_and_prune(h, X, IntervalRootFinding.BisectionContractor, tol)

    while tol > 1e-14
       tol /= 2
       roots = branch_and_prune(h, roots, IntervalRootFinding.BisectionContractor, tol)
    end

    return roots
end


const IntervalLike{T} = Union{Interval{T}, IntervalBox{T}}
const NewtonLike = Union{Type{Newton}, Type{Krawczyk}}


"""
    roots(f, X, contractor, tol=1e-3)
    roots(f, deriv, X, contractor, tol=1e-3)

Uses a generic branch and prune routine to find in principle all isolated roots of a function
`f:R^n → R^n` in a box `X`, or a vector of boxes.

Inputs:
- `f`: function whose roots will be found
- `X`: `Interval` or `IntervalBox`
- `contractor`: function that, when applied to the function `f`, determines
    the status of a given box `X`. It returns the new box and a symbol indicating
    the status. Current possible values are `Bisection`, `Newton` and `Krawczyk`
- `deriv` ; explicit derivative of `f` for `Newton` and `Krawczyk`
"""
function roots(f, X::IntervalLike{T}, ::Type{Bisection}, tol::Float64=1e-3) where {T}
    branch_and_prune(X, Bisection(f), tol)
end

function roots(f, X::Interval{T}, C::NewtonLike, tol::Float64=1e-3) where {T}

    deriv = x -> ForwardDiff.derivative(f, x)

    roots(f, deriv, X, C, tol)
end

function roots(f, deriv, X::Interval{T}, C::NewtonLike, tol::Float64=1e-3) where {T}
    branch_and_prune(X, C(f, deriv), tol)
end

function roots(f, X::IntervalBox{T}, C::NewtonLike, tol::Float64=1e-3) where {T}

    deriv = x -> ForwardDiff.jacobian(f, x)

    roots(f, deriv, X, C, tol)
end

function roots(f, deriv, X::IntervalBox{T}, C::NewtonLike, tol::Float64=1e-3) where {T}
    branch_and_prune(X, C(f, deriv), tol)
end

roots(f, r::Root, contractor::Type{C}, tol::Float64=1e-3) where {C<:Contractor} =
    roots(f, r.interval, contractor, tol)
roots(f, deriv, r::Root, contractor::Type{C}, tol::Float64=1e-3) where {C<:Contractor} =
    roots(f, deriv, r.interval, contractor, tol)


# Acting on a Vector:

# TODO: Use previous status information about roots:
roots(f, V::Vector{Root{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor} =
    vcat(roots.(f, V, contractor, tol)...)
roots(f, deriv, V::Vector{Root{T}}, contractor::Type{C}, tol::Float64=1e-3) where {T, C<:Contractor} =
    vcat(roots.(f, deriv, V, contractor, tol)...)



# Complex:

function roots(f, Xc::Complex{Interval{T}}, contractor::Type{C},
        tol::Float64=1e-3) where {T, C<:Contractor}

    g = realify(f)
    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, Y, contractor, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function roots(f, Xc::Complex{Interval{T}}, C::NewtonLike, tol::Float64=1e-3) where {T}

    g = realify(f)

    g_prime = x -> ForwardDiff.jacobian(g, x)

    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, g_prime, Y, C, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

function roots(f, deriv, Xc::Complex{Interval{T}}, C::NewtonLike, tol::Float64=1e-3) where {T}

    g = realify(f)

    g_prime = realify_derivative(deriv)

    Y = IntervalBox(reim(Xc)...)
    rts = roots(g, g_prime, Y, C, tol)

    return [Root(Complex(root.interval...), root.status) for root in rts]
end

# Default
roots(f::F, deriv::F, X, tol::Float64=1e-15) where {F<:Function} = roots(f, deriv, X, Newton, tol)
roots(f::F, X, tol::Float64=1e-15) where {F<:Function} = roots(f, X, Newton, tol)
