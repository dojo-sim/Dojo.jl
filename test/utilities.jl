@testset "Utilities" begin 
    @test length(Dojo.module_dir()) > 0
    @test Dojo.scn(1.0) == " 1.0e+0"
    @test Dojo.scn(0.0) == " 0.0e+0"
    @test Dojo.scn(Inf) == " Inf"
    @test Dojo.scn(123.1, 
        digits=0) == " 1e+2"
end


