import Base: copy, eltype, iterate, IteratorSize
import Base: getindex, setindex!, delete!

export BBSearch, SearchStrategy
export BreadthFirstBBSearch, DepthFirstBBSearch
export copy, eltype, iterate, IteratorSize

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
    BBSearch{DATA}

Branch and bound search interface in element of type DATA.

This interface provide an iterable that perform the search.

# Methods to implement:
    - `root_element(::BBSearch)`: return the element with which the searc is started
    - `process(::BBSearch, elem::DATA)`: return a symbol representing the action
        to perform with the element `elem` and a `elem` itself (possibly refined)
    - `bisect(::BBSearch, elem::DATA)`: return two elements build by bisecting `elem`

# Actions returned by the process function
    - `:store`: the element is considered as final and is stored, it will not be
        further processed
    - `:bisect`: the element is bisected and each of the two resulting part will
        be further processed
    - `:discard`: the element is discarded from the tree, allowing to free memory
"""
abstract type BBSearch{DATA} end

abstract type BreadthFirstBBSearch{DATA} <: BBSearch{DATA} end
abstract type DepthFirstBBSearch{DATA} <: BBSearch{DATA} end
abstract type KeyBBSearch{DATA} <: BBSearch{DATA} end

get_leaf_id!(::BreadthFirstBBSearch, wt::WorkingTree) = shift!(wt.working_leafs)
get_leaf_id!(::DepthFirstBBSearch, wt::WorkingTree) = pop!(wt.working_leafs)
get_leaf_id!(::KeyBBSearch, wt::WorkingTree) = shift!(wt.working_leafs)

function insert_leaf!(::Union{BreadthFirstBBSearch{DATA}, DepthFirstBBSearch{DATA}},
                      wt::WorkingTree{DATA}, leaf::WorkingLeaf{DATA}) where {DATA}
    id = newid(wt)
    wt[id] = leaf
    push!(wt.working_leafs, id)
    return id
end

function insert_leaf!(::KS, wt::WorkingTree{DATA}, leaf::WorkingLeaf{DATA}) where {DATA, KS <: KeyBBSearch{DATA}}
    id = newid(wt)
    wt[id] = leaf
    keys = keyfunc.(KS, wt.working_leafs)
    current = keyfunc(KS, leaf)

    # Keep the working_leafs sorted
    insert!(wt.working_leafs, searchsortedfirst(keys, current), id)
    return id
end

eltype(::Type{BBS}) where {DATA, BBS <: BBSearch{DATA}} = WorkingTree{DATA}
IteratorSize(::Type{BBS}) where {BBS <: BBSearch} = Base.SizeUnknown()

function iterate(search::BBSearch{DATA},
                 wt::WorkingTree=WorkingTree(root_element(search))) where {DATA}
    isempty(wt.working_leafs) && return nothing

    id = get_leaf_id!(search, wt)
    X = wt[id]
    action, newdata = process(search, data(X))
    if action == :store
        wt[id] = WorkingLeaf(newdata, parent_id(X), :final)
    elseif action == :bisect
        parent = wt[id]
        child1, child2 = bisect(search, newdata)
        leaf1 = WorkingLeaf(child1, id, :working)
        leaf2 = WorkingLeaf(child2, id, :working)
        id1 = insert_leaf!(search, wt, leaf1)
        id2 = insert_leaf!(search, wt, leaf2)
        wt[id] = WorkingNode(parent, id1, id2)
    elseif action == :discard
        discard_leaf!(wt, id)
    else
        warn("Branch and bound: process function of the search object return unkown action: $action, element $X is ignored. Valid actions are :store, :bisect and :discard.")
    end
    return wt, wt
end
