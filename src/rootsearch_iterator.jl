import Base: start, next, done, copy, eltype, iteratorsize
import Base: getindex, setindex!, delete!

export BBSearch, SearchStrategy
export BreadthFirstSearch, DepthFirstSearch
export start, next, done, copy, step!, eltype, iteratorsize

abstract type AbstractWorkingNode end

struct WorkingNode <: AbstractWorkingNode
    parent::Int
    children::Vector{Int}
end

struct WorkingLeaf{DATA} <: AbstractWorkingNode
    data::DATA
    parent::Int
    status::Symbol
end

function WorkingLeaf(data::DATA, parent::Int) where {DATA}
    WorkingLeaf{DATA}(data, parent, :working)
end

function WorkingNode(leaf::WorkingLeaf, child1::Int, child2::Int)
    WorkingNode(parent_id(leaf), Int[child1, child2])
end

struct WorkingTree{DATA}
    nodes::Dict{Int, Union{WorkingNode, WorkingLeaf{DATA}}}
    working_leafs::Vector{Int}
end

function WorkingTree(rootdata::DATA) where {DATA}
    rootleaf = WorkingLeaf(rootdata, 0)
    WorkingTree{DATA}(Dict{Int, Union{WorkingNode, WorkingLeaf{DATA}}}(1 => rootleaf), Int[1])
end

parent_id(node::AbstractWorkingNode) = node.parent
parent_id(wt::WorkingTree, node::AbstractWorkingNode) = parent_id(node)
parent_id(wt::WorkingTree, id::Int) = parent_id(wt[id])
parent(wt::WorkingTree, node::AbstractWorkingNode) = wt[parent_id(node)]
parent(wt::WorkingTree, id::Int) = wt[parent_id(wt, id)]
child_ids(node::WorkingNode) = node.children
data(leaf::WorkingLeaf) = leaf.data

nnodes(wt::WorkingTree) = length(wt.nodes)
newid(wt::WorkingTree) = maximum(keys(wt.nodes)) + 1
data(wt::WorkingTree) = [data(val) for val in values(wt.nodes) if isa(val, WorkingLeaf)]

# Index operations
getindex(wt::WorkingTree, id) = wt.nodes[id]
setindex!(wt::WorkingTree, id, val) = setindex!(wt.nodes, id, val)
delete!(wt::WorkingTree, id) = delete!(wt.nodes, id)

function discard_leaf!(wt::WorkingTree, id::Int)
    leaf = wt[id]
    recursively_delete_child!(wt, parent_id(leaf), id)
end

function recursively_delete_child!(wt, id_parent, id_child)
    parent = wt[id_parent]
    cc = child_ids(parent)
    filter!(id -> id == id_child, cc)
    if isempty(cc) && parent_id(parent) != 0
        recursively_delete_child!(wt, parent_id(parent), id_parent)
    end
    delete!(wt, id_child)
end


"""
    SearchStrategy{KEY}

Abstract type for the strategies followed to chose the order in which elements
are processed during a `BBSearch`.
"""
abstract type SearchStrategy{KEY} end

struct BreadthFirstSearch <: SearchStrategy{Void} end
struct DepthFirstSearch <: SearchStrategy{Void} end

struct KeySearch{KEY <: Function} <: SearchStrategy{KEY}
    keyfunc::KEY
end

get_leaf_id!(strat::BreadthFirstSearch, wt::WorkingTree) = shift!(wt.working_leafs)
get_leaf_id!(strat::DepthFirstSearch, wt::WorkingTree) = pop!(wt.working_leafs)
get_leaf_id!(strat::KeySearch, wt::WorkingTree) = shift!(wt.working_leafs)

function insert_leaf!(strat::Union{BreadthFirstSearch, DepthFirstSearch},
                      wt::WorkingTree{DATA}, leaf::WorkingLeaf{DATA}) where {DATA}
    id = newid(wt)
    wt[id] = leaf
    push!(wt.working_leafs, id)
    return id
end

function insert_leaf!(strat::KeySearch, wt::WorkingTree{DATA}, leaf::WorkingLeaf{DATA}) where {DATA}
    id = newid(wt)
    wt[id] = leaf
    keys = strat.keyfunc.(wt.working_leafs)
    current = strat.keyfunc(leaf)

    # Keep the working_leafs sorted
    insert!(wt.working_leafs, searchsortedfirst(keys, current), id)
    return id
end

"""
    BBSearch{DATA, PFUNC <: Function, BFUNC <: Function, KEY}

Type implementing the `Base.Iterator` interface to the branch and prune routine.
Returns the `WorkingTree` at each iteration. Note: Each iteration mutates
the `WorkingTree`.

# Fields:
    - `inital_element`
    - `process`: Function deciding the status of an element and optionnally
        refining it. Must return an element and the action to be performed on it.
        Valid action are `:bisect`, `:discard` and `:store`.
    - `bisect`: Function bisecting an element.
    - `strategy`: `SearchStrategy` determining the order in which the elements
        are processed.
"""
struct BBSearch{DATA, PFUNC <: Function, BFUNC <: Function, KEY}
    initial_element::DATA
    process::PFUNC
    bisect::BFUNC
    strategy::SearchStrategy{KEY}
end

function BBSearch(init, process::Function, bisect::Function, S::Type{STRAT}) where {STRAT <: Union{BreadthFirstSearch, DepthFirstSearch}}
    return BBSearch(init, process, bisect, S())
end

eltype(::Type{BBS}) where {DATA, PFUNC, BFUNC, KEY, BBS <: BBSearch{DATA, PFUNC, BFUNC, KEY}} = WorkingTree{DATA}
iteratorsize(::Type{BBS}) where {BBS <: BBSearch} = Base.SizeUnknown()


function start(iter::BBSearch{DATA, PFUNC, BFUNC, KEY}) where {DATA, PFUNC, BFUNC, KEY}
    return WorkingTree(iter.initial_element)
end

"""
    step!(state::WorkingTree, contractor, tolerance)

Progress `state` by treating one of its `working` working_leafs. Note: the
working tree `wt` is always modified.
"""
function step!(wt::WorkingTree, search)
    strat = search.strategy
    id = get_leaf_id!(strat, wt)
    X = wt[id]
    action, newdata = search.process(data(X))
    if action == :store
        wt[id] = WorkingLeaf(newdata, parent_id(X), :final)
    elseif action == :bisect
        parent = wt[id]
        child1, child2 = search.bisect(newdata)
        leaf1 = WorkingLeaf(child1, id, :working)
        leaf2 = WorkingLeaf(child2, id, :working)
        id1 = insert_leaf!(strat, wt, leaf1)
        id2 = insert_leaf!(strat, wt, leaf2)
        wt[id] = WorkingNode(parent, id1, id2)
    elseif action == :discard
        discard_leaf!(wt, id)
    else
        warn("Branch and bound: process function of the search object return unkown action: $action, element $X is ignored. Valid actions are :store, :bisect and :discard.")
    end
end

function next(iter::BBSearch, state::WorkingTree)
    step!(state, iter)
    return state, state
end

done(iter::BBSearch, state::WorkingTree) = isempty(state.working_leafs)
