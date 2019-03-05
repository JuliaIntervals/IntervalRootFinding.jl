import IntervalRootFinding: BBLeaf, BBNode, BBTree
import IntervalRootFinding: root, nnodes, data, newid, discard_leaf!

import IntervalRootFinding: BreadthFirstBBSearch
import IntervalRootFinding: process, bisect

@testset "Branch and bound tree" begin
    #= Manually create a simple tree with the following structure

        (1)
         |
        (2)
    =#
    leaf = BBLeaf("I'm a leaf", 1, :final)
    node = BBNode(0, [2])
    tree = BBTree(Dict(1 => node),
                  Dict(2 => leaf),
                  [2])

    # Check that printing does not error
    io = IOBuffer()
    println(io, leaf)
    println(io, node)
    println(io, tree)

    @test root(tree) == node
    @test nnodes(tree) == 2

    d = data(tree)
    @test length(d) == 1
    @test d[1] == "I'm a leaf"

    k = newid(tree)
    @test k ∉ 0:2  # The new ID must be new and not 0
    @test isa(k, Integer)

    discard_leaf!(tree, 2)
    # The last leaf was deleted so the tree should be empty
    @test nnodes(tree) == 0

    #= Tests on a slighty more complicated tree

          (1)
         /   \
       (2)   (8)
      /  \     \
    (3)  (5)   (9)
     |   | \     \
    (4) (6) (7)  (10)
    =#

    n1 = BBNode(0, [2, 8])
    n2 = BBNode(1, [3, 5])
    n3 = BBNode(2, [4])
    n4 = BBLeaf("Leaf 4", 3, :working)
    n5 = BBNode(2, [6, 7])
    n6 = BBLeaf("Leaf 6", 5, :final)
    n7 = BBLeaf("Leaf 7", 5, :working)
    n8 = BBNode(1, [9])
    n9 = BBNode(8, [10])
    n10 = BBLeaf("Leaf 10", 9, :working)

    tree = BBTree(Dict(1 => n1, 2 => n2, 3 => n3, 5 => n5, 8 => n8, 9 => n9),
                  Dict(4 => n4, 6 => n6, 7 => n7, 10 => n10),
                  [3, 4, 8])

    # Check that printing does not error
    io = IOBuffer()
    println(io, tree)

    @test newid(tree) ∉ 0:10
    @test length(data(tree)) == 4

    @test nnodes(tree) == 10

    discard_leaf!(tree, 6)
    @test nnodes(tree) == 9

    discard_leaf!(tree, 4)
    @test nnodes(tree) == 7

    discard_leaf!(tree, 10)
    @test nnodes(tree) == 4

    @test length(data(tree)) == 1
end


# Implement the interface for breadth first in a dummy way
struct DummySearch <: BreadthFirstBBSearch{Symbol} end

process(::DummySearch, s::Symbol) = s, s
bisect(::DummySearch, s::Symbol) = :store, :store

@testset "Interface functions" begin
    #= Build the following tree for testing
                    (1)
                   /   \
                (2)     (5: to_bisect)
               /   \
    (3: to_store)  (4: to_discard)
    =#
    rt = BBNode(0, [2, 5])
    node = BBNode(1, [3, 4])
    to_store = BBLeaf(:store, 2, :working)
    to_discard = BBLeaf(:discard, 2, :working)
    to_bisect = BBLeaf(:bisect, 1, :working)

    tree = BBTree(Dict(1 => rt, 2 => node),
                  Dict(3 => to_store, 4 => to_discard, 5 => to_bisect),
                  [3, 4, 5])

    search = DummySearch()

    # First iteration processes the to_store leaf
    tree, _ = iterate(search, tree)

    @test nnodes(tree) == 5
    @test tree[3].status == :final
    @test length(tree.working_leaves) == 2

    # Second iteration processes the to_discard leaf
    tree, _ = iterate(search, tree)

    @test nnodes(tree) == 4
    @test length(tree.working_leaves) == 1

    # Second iteration processes the to_bisect leaf
    tree, _ = iterate(search, tree)

    @test nnodes(tree) == 6
    @test length(tree.working_leaves) == 2
    @test isa(tree[5], BBNode)
end
