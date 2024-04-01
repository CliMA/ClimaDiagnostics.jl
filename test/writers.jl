using Test

import ClimaDiagnostics.Writers

@testset "DictWriter" begin
    writer = Writers.DictWriter()

    # Test with some strings and floats instead of actual Fields and ScheduledDiagnostics
    Writers.write_field!(writer, 10.0, "mytest", nothing, nothing, 0.0)
    @test writer.dict["mytest"][0.0] == 10.0
    Writers.write_field!(writer, 20.0, "mytest", nothing, nothing, 2.0)
    @test writer.dict["mytest"][2.0] == 20.0
    Writers.write_field!(writer, 50.0, "mytest2", nothing, nothing, 8.0)
    @test writer.dict["mytest2"][8.0] == 50.0
end
