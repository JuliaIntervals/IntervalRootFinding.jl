include("../examples/smiley_examples.jl")

using .SmileyExample22, .SmileyExample52, .SmileyExample54, .SmileyExample55

const abstol = 1e-6

@testset "$(SmileyExample22.title)" begin
    for contractor in (Newton, Krawczyk)
        roots_found = roots(SmileyExample22.f, SmileyExample22.region ; contractor, abstol)
        @test length(roots_found) == 8
        @test all(isunique, roots_found)
        # no reference data for roots given
    end
end

@testset "$(example.title)" for example in (SmileyExample52, SmileyExample54)
    for contractor in (Newton, Krawczyk)
        roots_found = roots(example.f, example.region ; contractor, abstol)
        @test length(roots_found) == length(example.known_roots)
        @test all(isunique, roots_found)
        for rf in roots_found
            # check there is exactly one known root for each found root
            @test 1 == sum(example.known_roots) do rk
                !isempty_region(intersect_region(rk, root_region(rf)))
            end
        end
    end
end
