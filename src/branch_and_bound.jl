import Base: copy, eltype, iterate, IteratorSize
import Base: getindex, setindex!, delete!

export BBSearch, BreadthFirstBBSearch, DepthFirstBBSearch
export data
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

Leaf node of a `BBTree` that contains some data. Its status is either
    - `:working`: the leaf will be further processed.
    - `:final`: the leaf won't be touched anymore.
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
    BBNode(leaf.parent, Int[child1, child2])
end

"""
    BBTree{DATA}

Tree storing the data used and produced by a branch and bound search in a
structured way.

Nodes and leaves can be accessed using their index using the bracket syntax
`wt[node_id]`. However this is slow, as nodes and leaves are stored separately.

Support the iterator interface. The element yielded by the iteration are
tuples `(node_id, lvl)` where `lvl` is the depth of the node in the tree.
"""
struct BBTree{DATA}
    nodes::Dict{Int, BBNode}
    leaves::Dict{Int, BBLeaf{DATA}}
    working_leaves::Vector{Int}
end

function BBTree(rootdata::DATA) where {DATA}
    rootleaf = BBLeaf(rootdata, 0)
    BBTree{DATA}(Dict{Int, BBNode}(), Dict(1 => rootleaf), Int[1])
end

show(io::IO, wn::BBNode) = print(io, "Node with children $(wn.children)")

function show(io::IO, wl::BBLeaf)
    print(io, "Leaf (:$(wl.status)) with data $(wl.data)")
end

function show(io::IO, wt::BBTree{DATA}) where {DATA}
    println(io, "Working tree with $(nnodes(wt)) elements of type $DATA")

    if nnodes(wt) > 0
        println(io, "Indices: ", vcat(collect(keys(wt.nodes)), collect(keys(wt.leaves))) |> sort)
        println(io, "Structure:")
        for (id, lvl) in wt
            println(io, "  "^lvl * "[$id] $(wt[id])")
        end
    end
end

# Root node has id 1 and parent id 0
root(wt::BBTree) = wt[1]
is_root(wt::BBTree, id::Int) = (id == 1)

"""
    nnodes(wt::BBTree)

Number of nodes (including leaves) in a `BBTree`.
"""
nnodes(wt::BBTree) = length(wt.nodes) + length(wt.leaves)

"""
    data(leaf::BBLeaf)

Return the data stored in the leaf.
"""
data(leaf::BBLeaf) = leaf.data

"""
    data(wt::BBTree)

Return all the data stored in a `BBTree` as a list. The ordering of the elements
is arbitrary.
"""
data(wt::BBTree) = data.(values(wt.leaves))

function newid(wt::BBTree)
    k1 = keys(wt.nodes)
    k2 = keys(wt.leaves)

    if length(k1) > 0
        m1 = maximum(k1)
    else
        m1 = 0
    end

    if length(k2) > 0
        m2 = maximum(k2)
    else
        m2 = 0
    end

    return max(m1, m2) + 1
end

# Index operations (slower than manipulating the node directly in the correct
# dictionary)
function getindex(wt::BBTree, id)
    haskey(wt.nodes, id) && return wt.nodes[id]
    haskey(wt.leaves, id) && return wt.leaves[id]
    error("getindex failed: no index $id")  # TODO: make better error
end

setindex!(wt::BBTree, val::BBNode, id) = setindex!(wt.nodes, val, id)
setindex!(wt::BBTree, val::BBLeaf, id) = setindex!(wt.leaves, val, id)

function delete!(wt::BBTree, id)
    if haskey(wt.nodes, id)
        delete!(wt.nodes, id)
    elseif haskey(wt.leaves, id)
        delete!(wt.leaves, id)
    else
        error("delete! failed: no index $id")  # TODO: make better error
    end
end

"""
    discard_leaf!(wt::BBTree, id::Int)

Delete the `BBLeaf` with index `id` and all its ancestors to which it is
the last descendant.
"""
function discard_leaf!(wt::BBTree, id::Int)
    leaf = wt.leaves[id]
    delete!(wt.leaves, id)
    recursively_delete_parent!(wt, leaf.parent, id)
end

function recursively_delete_parent!(wt, id_parent, id_child)
    if !is_root(wt, id_child)
        parent = wt.nodes[id_parent]
        siblings = parent.children
        if length(parent.children) == 1  # The child has no siblings, so delete the parent
            delete!(wt.nodes, id_parent)
            recursively_delete_parent!(wt, parent.parent, id_parent)
        else  # The child has siblings so remove it from the children list
            deleteat!(parent.children, searchsortedfirst(parent.children, id_child))
        end
    end
end

function iterate(wt::BBTree, (id, lvl)=(0, 0))
    id, lvl = next_id(wt, id, lvl)
    lvl == 0 && return nothing
    return (id, lvl), (id, lvl)
end

function next_id(wt::BBTree, id, lvl)
    lvl == 0 && return (1, 1)
    node = wt[id]
    isa(node, BBNode) && return (node.children[1], lvl + 1)
    return next_sibling(wt, id, lvl)
end

function next_sibling(wt::BBTree, sibling, lvl)
    parent = wt[sibling].parent
    parent == 0 && return (0, 0)
    children = wt[parent].children
    maximum(children) == sibling && return next_sibling(wt, parent, lvl - 1)
    id = minimum(filter(x -> x > sibling, children))
    return (id, lvl)
end


"""
    BBSearch{DATA}

Branch and bound search interface in element of type DATA.

This interface provide an iterable that perform the search.

There is currently three types of search supported `BreadFirstBBSearch`,
`DepthFirstBBSearch` and `KeyBBSearch`, each one processing the element of the
tree in a different order. When subtyping one of these, the following methods
must be implemented:
    - `root_element(::BBSearch)`: return the element with which the search is started
    - `process(::BBSearch, elem::DATA)`: return a symbol representing the action
        to perform with the element `elem` and an object of type `DATA` reprensenting
        the state of the element after processing (may return `elem` unchanged).
    - `bisect(::BBSearch, elem::DATA)`: return two elements of type `DATA` build
        by bisecting `elem`

Subtyping `BBSearch` directly allows to have control over the order in which
the elements are process. To do this the following methods must be implemented:
    - `root_element(::BBSearch)`: return the first element to be processed. Use
        to build the initial tree.
    - `get_leaf_id!(::BBSearch, wt::BBTree)`: return the id of the next leaf that
        will be processed and remove it from the list of working leaves of `wt`.
    - `insert_leaf!(::BBSearch, wt::BBTree, leaf::BBLeaf)`: insert a leaf in the
        list of working leaves.

# Valid symbols returned by the process function
    - `:store`: the element is considered as final and is stored, it will not be
        further processed
    - `:bisect`: the element is bisected and each of the two resulting part will
        be processed
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

!!WARNING: Untested.
"""
abstract type KeyBBSearch{DATA} <: BBSearch{DATA} end

"""
    root_element(search::BBSearch)

Return the initial element of the search. The `BBTree` will be build around it.

Can be define for custom searches that are direct subtype of `BBSearch`, default
behavior is to fetch the field `initial` of the search.
"""
root_element(search::BBSearch) = search.initial

"""
    get_leaf_id!(::BBSearch, wt::BBTree)

Return the id of the next leaf that will be processed and remove it from the
list of working leaves.

Must be define for custom searches that are direct subtype of `BBSearch`.
"""
get_leaf_id!(::BreadthFirstBBSearch, wt::BBTree) = popfirst!(wt.working_leaves)
get_leaf_id!(::DepthFirstBBSearch, wt::BBTree) = pop!(wt.working_leaves)
get_leaf_id!(::KeyBBSearch, wt::BBTree) = popfirst!(wt.working_leaves)

"""
    insert_leaf!(::BBSearch, wt::BBTree, leaf::BBLeaf)

Insert the id of a new leaf that has been produced by bisecting an older leaf
into the list of working leaves.

Must be define for custom searches that are direct subtype of `BBSearch`.
"""
function insert_leaf!(::Union{BreadthFirstBBSearch{DATA}, DepthFirstBBSearch{DATA}},
                      wt::BBTree{DATA}, leaf::BBLeaf{DATA}) where {DATA}
    id = newid(wt)
    wt.leaves[id] = leaf
    push!(wt.working_leaves, id)
    return id
end

function insert_leaf!(::KS, wt::BBTree{DATA}, leaf::BBLeaf{DATA}) where {DATA, KS <: KeyBBSearch{DATA}}
    id = newid(wt)
    wt.leaves[id] = leaf
    keys = keyfunc.(KS, wt.working_leaves)
    current = keyfunc(KS, leaf)

    # Keep the working_leaves sorted
    insert!(wt.working_leaves, searchsortedfirst(keys, current), id)
    return id
end

eltype(::Type{BBS}) where {DATA, BBS <: BBSearch{DATA}} = BBTree{DATA}
IteratorSize(::Type{BBS}) where {BBS <: BBSearch} = Base.SizeUnknown()

function iterate(search::BBSearch{DATA},
                 wt::BBTree=BBTree(root_element(search))) where {DATA}

    isempty(wt.working_leaves) && return nothing

    id = get_leaf_id!(search, wt)
    X = wt.leaves[id]
    action, newdata = process(search, data(X))
    if action == :store
        wt.leaves[id] = BBLeaf(newdata, X.parent, :final)
    elseif action == :bisect
        child1, child2 = bisect(search, newdata)
        leaf1 = BBLeaf(child1, id, :working)
        leaf2 = BBLeaf(child2, id, :working)
        id1 = insert_leaf!(search, wt, leaf1)
        id2 = insert_leaf!(search, wt, leaf2)
        wt.nodes[id] = BBNode(X, id1, id2)
        delete!(wt.leaves, id)
    elseif action == :discard
        discard_leaf!(wt, id)
    else
        error("Branch and bound: process function of the search object return " *
              "unknown action: $action for element $X. Valid actions are " *
              ":store, :bisect and :discard.")
    end
    return wt, wt
end
