# test_wavelets.jl
using BasisFunctions
using Base.Test
using Wavelets.DWT: wavelet, scaling

function bf_wavelets_implementation_test()
    @testset begin
        b1 = DaubechiesWaveletBasis(3,2)
        b2 = CDFWaveletBasis(3,1,5)
        b = CDFWaveletBasis(1,1,3)

        eval_element(b, 1, 0.5)
        println("Timings should be in the order of 0.001 seconds, 22.01k allocations, 625.094 KB")
        @time begin
            for x in linspace(0,1,1000)
                BasisFunctions.unsafe_eval_element(b, 1, x)
            end
        end
        supports = ((0,1),(0,1),(0.0,0.5),(0.5,1.0),(0.0,0.25),(0.25,0.5),(0.5,0.75),(0.75,1.0));
        for i in ordering(b)
            @test left(b,i) == supports[value(i)][1]
            @test right(b,i) == supports[value(i)][2]
        end
        for i in ordering(b1)
            @test support(b1,i) == (0.,1.)
        end

        @test BasisFunctions.subdict(b1,1:1) == BasisFunctions.DaubechiesWaveletBasis(3,0)
        @test BasisFunctions.subdict(b2,1:3) == BasisFunctions.LargeSubdict(b2,1:3)
        @test b2[1:5] == BasisFunctions.LargeSubdict(b2,1:5)
        @test b2[1:4] == BasisFunctions.CDFWaveletBasis(3,1,2)
        @test BasisFunctions.dyadic_length(b1) == 2
        @test BasisFunctions.dyadic_length(b2) == 5
        @test length(b1) == 4
        @test length(b2) == 32
        @test promote_eltype(b1,Complex128) == DaubechiesWaveletBasis(3,2, Complex128)
        @test promote_eltype(b2,Complex128) == CDFWaveletBasis(3,1,5, Complex128)
        @test resize(b1,8) == BasisFunctions.DaubechiesWaveletBasis(3,3)
        @test BasisFunctions.name(b1) == "Basis of db3 wavelets"
        @test BasisFunctions.name(b2) == "Basis of cdf31 wavelets"
        @test BasisFunctions.native_index(b,1) == (scaling, 0, 0)
        @test BasisFunctions.native_index(b,2) == (wavelet, 0, 0)
        @test BasisFunctions.native_index(b,3) == (wavelet, 1, 0)
        @test BasisFunctions.native_index(b,4) == (wavelet, 1, 1)
        @test BasisFunctions.native_index(b,5) == (wavelet, 2, 0)
        @test BasisFunctions.native_index(b,6) == (wavelet, 2, 1)
        @test BasisFunctions.native_index(b,7) == (wavelet, 2, 2)
        @test BasisFunctions.native_index(b,8) == (wavelet, 2, 3)

        for i in 1:length(b)
            @test(linear_index(b,BasisFunctions.native_index(b,i))==i )
        end
        @test grid(b1) == PeriodicEquispacedGrid(4,0,1)
        @test grid(b2) == PeriodicEquispacedGrid(32,0,1)
        @test BasisFunctions.period(b1)==1.

        # test grid eval functions
        g = plotgrid(b,200)
        for i in ordering(b)
            tic(); e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g); t1 = toq();
            tic(); e2 = BasisFunctions._unsafe_eval_element_in_dyadic_grid(b, i, g); t2 = toq();
        end
        for i in ordering(b)
            tic(); e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g); t1 = toq();
            tic(); e2 = BasisFunctions._unsafe_eval_element_in_dyadic_grid(b, i, g); t2 = toq();
            @test e1 ≈ e2
            @test t2 < t1
        end
    end
end
bf_wavelets_implementation_test()

using Plots
b = CDFWaveletBasis(2,2,3)
gr()
plot(b)
