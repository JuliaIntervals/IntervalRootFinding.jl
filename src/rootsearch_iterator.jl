import Base: start, next, done, copy, eltype, iteratorsize

export RootSearch, SearchStrategy
export start, next, done, copy, step!, eltype, iteratorsize

"""
    SearchStrategy{CON}(store!, retrieve!)

Type describing the strategy followed to chose the order in which intervals are processed during a `RootSearch`, using a container of type `CON` to store
intervals.

The container type `CON` must support the syntax `CON{EL}()` to create an empty
container for elements of type `EL`.

# Fields:
    - `store!(container, element)`: add `element` to `container`
    - `retrieve!(container)`: return the next element to be processed and remove
        it from `container`
"""
struct SearchStrategy{CON}
    store!::Function
    retrieve!::Function
end

SearchStrategy(CON::Type, store!::Function, retrieve!::Function) = SearchStrategy{CON}(store!, retrieve!)

SearchStrategy() = SearchStrategy(Vector, push!, pop!)

"""
    RootSearch{R <: Union{Interval,IntervalBox}, C <: Contractor, CON, T <: Real}

Type implementing the `Base.Iterator` interface to the branch and prune routine.
Returns the `RootSearchState` at each iteration. Note: Each iteration mutates
the `RootSearchState`. Use `copy(state::RootSearchState)` to create an
independent instance if necessary.

# Fields:
    - `region`: Region in which zeros are searched
    - `contractor`: Contractor used to determine the status of a region
    - `strategy`: Strategy determining the order in which regions are processed
        and the type of container used to store them
    - `tolerance`: Absolute tolerance of the search
"""
struct RootSearch{R <: Union{Interval,IntervalBox}, C <: Contractor, CON, T <: Real}
    region::R
    contractor::C
    strategy::SearchStrategy{CON}
    tolerance::T
end

function RootSearch(region::R, contractor::C, tol::T) where {R <: Union{Interval,IntervalBox}, C <: Contractor, T <: Real}
    RootSearch(region, contractor, SearchStrategy(), tol)
end

eltype(::Type{RS}) where {R, C, T, CON, RS <: RootSearch{R, C, CON, T}} = RootSearchState{CON{R}, CON{Root{R}}}
iteratorsize(::Type{RS}) where {RS <: RootSearch} = Base.SizeUnknown()

"""
    RootSearchState{V, VR}

State type for a `RootSearch` iterator.

The type of the containers used is determinde by the `SearchStrategy` used by
the `RootSearch`.
"""
struct RootSearchState{V, VR}
    working::V  # Should be a container of the form  CON{T}
    outputs::VR  # Should ba a container of root of the form CON{Root{T}}
end

function RootSearchState(rs::RootSearch)
    return RootSearchState(rs.region, rs.strategy)
end

function RootSearchState(region::R, strat::SearchStrategy{CON}) where {R <: Union{Interval,IntervalBox}, CON}
    working = CON{R}()
    outputs = CON{Root{R}}()
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


"""
    copy(state::RootSearchState)

Return an independant copy of `state`. The underlying container type used to
store regions must support the `deepcopy` function.
"""
copy(state::RootSearchState) =
    RootSearchState(deepcopy(state.working), deepcopy(state.outputs))


function start(iter::RootSearch{R, C, CON, T}) where {R, C, CON, T}
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

done(iter::RootSearch, state::RootSearchState) = isempty(state.working)
