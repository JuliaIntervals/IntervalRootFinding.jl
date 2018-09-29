import Base: copy, eltype, iterate, IteratorSize
import Base: getindex, setindex!, delete!

export BBSearch, SearchStrategy
export BreadthFirstBBSearch, DepthFirstBBSearch
export copy, eltype, iterate, IteratorSize

abstract type AbstractBBNode end

"""
    BBNode <: AbstractBBNode

Intermediate node of a `BBTree`. Does not contain any data by itself,
only redirect toward its children.
"""
struct BBNode <: AbstractBBNode
    parent::Int
    children::Vector{Int}
end

"""
    BBLeaf{DATA} <: AbstractBBLeaf

Leaf node of a `BBTree`. Contains both data and a status indicating if it
will be further processed by the branch and bound search.
"""
struct BBLeaf{DATA} <: AbstractBBNode
    data::DATA
    parent::Int
    status::Symbol
end

function BBLeaf(data::DATA, parent::Int) where {DATA}
    BBLeaf{DATA}(data, parent, :working)
end

function BBNode(leaf::BBLeaf, child1::Int, child2::Int)
    BBNode(parent_id(leaf), Int[child1, child2])
end

"""
    BBTree{DATA}

Tree storing the data used and produced by a branch and bound search in a
structured way.

Node can be accessed using their index using the bracket syntax `wt[node_id]`.
Support the iterator interface. The element yielded by the iteration are
tuples `(node_id, lvl)` where `lvl` is the depth of the node in the tree.
"""
struct BBTree{DATA}
    nodes::Dict{Int, Union{BBNode, BBLeaf{DATA}}}
    working_leafs::Vector{Int}
end

function BBTree(rootdata::DATA) where {DATA}
    rootleaf = BBLeaf(rootdata, 0)
    BBTree{DATA}(Dict{Int, Union{BBNode, BBLeaf{DATA}}}(1 => rootleaf), Int[1])
end

show(io::IO, wn::BBNode) = print(io, "BBNode with children $(wn.children)")
function show(io::IO, wl::BBLeaf)
    print(io, "BBLeaf (:$(wl.status)) with data ($(wl.data))")
end

function show(io::IO, wt::BBTree{DATA}) where {DATA}
    println(io, "Working tree with $(nnodes(wt)) elements of type $DATA")
    println(io, "Indices: ", keys(wt.nodes) |> collect |> sort)
    println(io, "Structure:")
    for (id, lvl) in wt
        println(io, "  "^lvl * "[$id] $(wt[id])")
    end
end

parent_id(node::AbstractBBNode) = node.parent
parent_id(wt::BBTree, node::AbstractBBNode) = parent_id(node)
parent_id(wt::BBTree, id::Int) = parent_id(wt[id])
parent(wt::BBTree, node::AbstractBBNode) = wt[parent_id(node)]
parent(wt::BBTree, id::Int) = wt[parent_id(wt, id)]
children_ids(node::BBNode) = node.children
children_ids(wt::BBTree, id::Int) = children_ids(wt[id])
data(leaf::BBLeaf) = leaf.data
root(wt::BBTree) = wt[1]

nnodes(wt::BBTree) = length(wt.nodes)
newid(wt::BBTree) = maximum(keys(wt.nodes)) + 1
data(wt::BBTree) = [data(val) for val in values(wt.nodes) if isa(val, BBLeaf)]

# Index operations
getindex(wt::BBTree, id) = wt.nodes[id]
setindex!(wt::BBTree, id, val) = setindex!(wt.nodes, id, val)
delete!(wt::BBTree, id) = delete!(wt.nodes, id)

"""
    discard_leaf!(wt::BBTree, id::Int)

Delete the `BBLeaf` with index `id` and all its ancestors to which it is
the last descendant.
"""
function discard_leaf!(wt::BBTree, id::Int)
    leaf = wt[id]
    recursively_delete_child!(wt, parent_id(leaf), id)
end

function recursively_delete_child!(wt, id_parent, id_child)
    parent = wt[id_parent]
    cc = children_ids(parent)
    filter!(id -> id != id_child, cc)
    if isempty(cc) && parent_id(parent) != 0
        recursively_delete_child!(wt, parent_id(parent), id_parent)
    end
    delete!(wt, id_child)
end

function iterate(wt::BBTree, (id, lvl)=(0, 0))
    id, lvl = next_id(wt, id, lvl)
    lvl == 0 && return nothing
    return (id, lvl), (id, lvl)
end

function next_id(wt::BBTree, id, lvl)
    lvl == 0 && return (1, 1)
    node = wt[id]
    isa(node, BBNode) && return (children_ids(node)[1], lvl + 1)
    return next_sibling(wt, id, lvl)
end

function next_sibling(wt::BBTree, sibling, lvl)
    parent = parent_id(wt, sibling)
    parent == 0 && return (0, 0)
    children = children_ids(wt, parent)
    maximum(children) == sibling && return next_sibling(wt, parent, lvl - 1)
    id = minimum(filter(x -> x > sibling, children))
    return (id, lvl)
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

"""
    KeyBBSearch{DATA} <: BBSearch{DATA}

Interface to a branch and bound search that use a key function to decide which
element to process first. The search process first the element with the largest
key as computed by `keyfunc(ks::KeyBBSearch, elem)`.

!WARNING: Untested.
"""
abstract type KeyBBSearch{DATA} <: BBSearch{DATA} end

"""
    get_leaf_id!(::BBSearch, wt::BBTree)

Return the id of the next leaf that will be processed and remove it from the
list of working leafs.

Must be define for custom searches that are direct subtype of `BBSearch`.
"""
get_leaf_id!(::BreadthFirstBBSearch, wt::BBTree) = popfirst!(wt.working_leafs)
get_leaf_id!(::DepthFirstBBSearch, wt::BBTree) = pop!(wt.working_leafs)
get_leaf_id!(::KeyBBSearch, wt::BBTree) = popfirst!(wt.working_leafs)

"""
    insert_leaf!(::BBSearch, wt::BBTree, leaf::BBLeaf)

Insert the id of a new leaf that has been produced by bisecting an older leaf
into the list of working leafs.

Must be define for custom searches that are direct subtype of `BBSearch`.
"""
function insert_leaf!(::Union{BreadthFirstBBSearch{DATA}, DepthFirstBBSearch{DATA}},
                      wt::BBTree{DATA}, leaf::BBLeaf{DATA}) where {DATA}
    id = newid(wt)
    wt[id] = leaf
    push!(wt.working_leafs, id)
    return id
end

function insert_leaf!(::KS, wt::BBTree{DATA}, leaf::BBLeaf{DATA}) where {DATA, KS <: KeyBBSearch{DATA}}
    id = newid(wt)
    wt[id] = leaf
    keys = keyfunc.(KS, wt.working_leafs)
    current = keyfunc(KS, leaf)

    # Keep the working_leafs sorted
    insert!(wt.working_leafs, searchsortedfirst(keys, current), id)
    return id
end

eltype(::Type{BBS}) where {DATA, BBS <: BBSearch{DATA}} = BBTree{DATA}
IteratorSize(::Type{BBS}) where {BBS <: BBSearch} = Base.SizeUnknown()

function iterate(search::BBSearch{DATA},
                 wt::BBTree=BBTree(root_element(search))) where {DATA}
    isempty(wt.working_leafs) && return nothing

    id = get_leaf_id!(search, wt)
    X = wt[id]
    action, newdata = process(search, data(X))
    if action == :store
        wt[id] = BBLeaf(newdata, parent_id(X), :final)
    elseif action == :bisect
        parent = wt[id]
        child1, child2 = bisect(search, newdata)
        leaf1 = BBLeaf(child1, id, :working)
        leaf2 = BBLeaf(child2, id, :working)
        id1 = insert_leaf!(search, wt, leaf1)
        id2 = insert_leaf!(search, wt, leaf2)
        wt[id] = BBNode(parent, id1, id2)
    elseif action == :discard
        discard_leaf!(wt, id)
    else
        warn("Branch and bound: process function of the search object return unkown action: $action, element $X is ignored. Valid actions are :store, :bisect and :discard.")
    end
    return wt, wt
end
