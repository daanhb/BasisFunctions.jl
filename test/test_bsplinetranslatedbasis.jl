function test_generic_periodicbsplinebasis(T)

  for B in (BSplineTranslatesBasis, SymBSplineTranslatesBasis,)
    tol = sqrt(eps(real(T)))
    n = 5
    b = B(n,3, T)
    @test left(b) == 0
    @test right(b)==1

    @test length(b)==5
    @test degree(b)==3
    @test is_basis(b)
    @test is_biorthogonal(b)
    @test !is_orthogonal(b)
    @test !is_orthonormal(b)
    @test !has_unitary_transform(b)

    @test BasisFunctions.left_of_compact_function(b) <= 0 <= BasisFunctions.right_of_compact_function(b)
    @test 0 < BasisFunctions.right_of_compact_function(b) - BasisFunctions.left_of_compact_function(b) < BasisFunctions.period(b)

    @test instantiate(B, 4, Float16)==B(4,3,Float16)
    @test promote_domaintype(b, Float16)==B(n,degree(b),Float16)
    @test promote_domaintype(b, complex(Float64))==B(n,degree(b),Complex128)
    @test resize(b, 20)==B(20,degree(b),T)
    @test BasisFunctions.grid(b)==PeriodicEquispacedGrid(n,0,1)
    @test BasisFunctions.period(b)==T(1)
    @test BasisFunctions.stepsize(b)==T(1//5)

    n = 3
    b=B(n,1,T)
    @test abs(sum(BasisFunctions.grammatrix(span(b)) - [2//3 1//6 1//6; 1//6 2//3 1//6;1//6 1//6 2//3]//n)) < tol
    @test abs(sum(BasisFunctions.dualgrammatrix(span(b)) - [5/3 -1/3 -1/3; -1/3 5/3 -1/3; -1/3 -1/3 5/3]*n)) < tol

    n = 8
    b=B(n,0,T)
    t = linspace(T(0),T(1))
    fp = map(x->BasisFunctions.eval_element(b,1,x),t)
    fd = 1/n*map(x->BasisFunctions.eval_dualelement(b,1,x),t)
    @test fp≈fd
  end
end

function test_translatedbsplines(T)
  tol = sqrt(eps(real(T)))
  n = 5
  bb = BSplineTranslatesBasis(n,1, T;scaled=true)
  b = BSplineTranslatesBasis(n,1, T)
  e = rand(n)
  @test norm(Gram(span(b))*e-Gram(span(bb))*e/n) < tol


  b = BSplineTranslatesBasis(n,3, T)

  @test BasisFunctions.name(b) == "Set of translates of a function (B spline of degree 3)"

  @test left(b,1)≈ 0
  @test left(b,2)≈1//5
  @test left(b,3)≈2//5
  @test left(b,4)≈3//5
  @test left(b,5)≈4//5
  @test right(b,1)≈4//5
  @test right(b,2)≈5//5
  @test right(b,3)≈6//5
  @test right(b,4)≈7//5
  @test right(b,5)≈8//5

  t = .001
  @test in_support(b,1,.5)
  @test !in_support(b,1,.8+t)
  @test !in_support(b,1,1.-t)
  @test in_support(b,3,.2-t)
  @test in_support(b,3,.4+t)
  @test !in_support(b,3,.2+t)
  @test !in_support(b,3,.4-t)

  @test BasisFunctions.compatible_grid(b, grid(b))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n+1,0,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0.1,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1.1))

  grid(BSplineTranslatesBasis(n,2, T)) == MidpointEquispacedGrid(n,0,1)
  @test degree(BSplineTranslatesBasis(5,2, T)) == 2
  b = BSplineTranslatesBasis(n,2,T)
  @test BasisFunctions.compatible_grid(b, grid(b))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n+1,0,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0.1,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1.1))

  # Test extension_operator and invertability of restriction_operator w.r.t. extension_operator.
  n = 8
  for degree in 0:3
    b = span(BSplineTranslatesBasis(n, degree, T))
    basis_ext = extend(b)
    r = restriction_operator(basis_ext, b)
    e = extension_operator(b, basis_ext)
    @test abs(sum(eye(n)-matrix(r*e))) < tol

    grid_ext = grid(basis_ext)
    L = evaluation_operator(b, grid_ext)
    e = random_expansion(b)
    z = L*e
    L2 = evaluation_operator(basis_ext, grid_ext) * extension_operator(b, basis_ext)
    z2 = L2*e
    @test maximum(abs.(z-z2)) < tol
  end


  if T == Float64
    for K in 0:3
      for s2 in 5:6
        s1 = s2<<1
        b1 = span(BSplineTranslatesBasis(s1,K,T))
        b2 = span(BSplineTranslatesBasis(s2,K,T))

        e1 = random_expansion(b1)
        e2 = random_expansion(b2)

        @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b2))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b2), grid(b2))*e2)
        @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b1))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b1), grid(b1))*e2)
        @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b1))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b1), grid(b1))*e1)
        @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b2))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b2), grid(b2))*e1)

        mr = matrix(restriction_operator(b1, b2))
        me = matrix(extension_operator(b2, b1))
        pinvme = pinv(me)
        r = rand(size(pinvme,2))
        @test pinvme*r ≈ mr*r
      end
    end

    for K in 0:3
      for s2 in 5:6
        s1 = s2<<1
        b1 = span(BSplineTranslatesBasis(s1,K,T; scaled=true))
        b2 = span(BSplineTranslatesBasis(s2,K,T; scaled=true))

        e1 = random_expansion(b1)
        e2 = random_expansion(b2)

        @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b2))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b2), grid(b2))*e2)
        @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b1))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b1), grid(b1))*e2)
        @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b1))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b1), grid(b1))*e1)
        @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b2))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b2), grid(b2))*e1)

        mr = matrix(restriction_operator(b1, b2))
        me = matrix(extension_operator(b2, b1))
        pinvme = pinv(me)
        r = rand(size(pinvme,2))
        @test pinvme*r ≈ mr*r
      end
    end
  end

  @test_throws AssertionError restriction_operator(span(BSplineTranslatesBasis(4,0,T)), span(BSplineTranslatesBasis(3,0,T)))
  @test_throws AssertionError extension_operator(span(BSplineTranslatesBasis(4,0,T)), span(BSplineTranslatesBasis(6,0,T)))
end


function test_translatedsymmetricbsplines(T)
  tol = sqrt(eps(real(T)))
  n = 5

  b = SymBSplineTranslatesBasis(n,1,T)
  bb = BSplineTranslatesBasis(n,1,T)
  @test norm((Gram(span(b))-Gram(span(bb)))*rand(n)) < tol

  b = SymBSplineTranslatesBasis(n,3, T)

  @test BasisFunctions.name(b) == "Set of translates of a function (symmetric B spline of degree 3)"


  @test left(b,1)≈ -4//10
  @test left(b,2)≈ -2//10
  @test left(b,3)≈ 0//10
  @test left(b,4)≈ 2//10
  @test left(b,5)≈ 4//10
  @test right(b,1)≈4//10
  @test right(b,2)≈6//10
  @test right(b,3)≈8//10
  @test right(b,4)≈10//10
  @test right(b,5)≈12//10

  t = .001
  @test in_support(b,1,.0)
  @test !in_support(b,1,.4+t)
  @test !in_support(b,1,.6-t)
  @test in_support(b,3,.8-t)
  @test in_support(b,3,.0+t)
  @test !in_support(b,3,.8+t)
  @test !in_support(b,3,.0-t)

  grid(SymBSplineTranslatesBasis(n,2, T)) == EquispacedGrid(n,0,1)
  @test BasisFunctions.compatible_grid(b, grid(b))
  @test BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n+1,0,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0.1,1))
  @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1.1))

  @test BasisFunctions.compatible_grid(b, grid(b))
  @test BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n+1,0,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0.1,1))
  @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1.1))

  # Test extension_operator and invertability of restriction_operator w.r.t. extension_operator.
  n = 8
  for degree in 1:2:3
    b = span(SymBSplineTranslatesBasis(n, degree, T))
    basis_ext = extend(b)
    r = restriction_operator(basis_ext, b)
    e = extension_operator(b, basis_ext)
    @test abs(sum(eye(n)-matrix(r*e))) < tol

    grid_ext = grid(basis_ext)
    L = evaluation_operator(b, grid_ext)
    e = random_expansion(b)
    z = L*e
    L2 = evaluation_operator(basis_ext, grid_ext) * extension_operator(b, basis_ext)
    z2 = L2*e
    @test maximum(abs.(z-z2)) < tol
  end



  for K in 1:2:3
    for s2 in 5:6
      s1 = s2<<1
      b1 = span(SymBSplineTranslatesBasis(s1,K))
      b2 = span(SymBSplineTranslatesBasis(s2,K))

      e1 = random_expansion(b1)
      e2 = random_expansion(b2)

      @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b2))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b2), grid(b2))*e2)
      @test coefficients(BasisFunctions.default_evaluation_operator(b2, gridspace(b1))*e2) ≈ coefficients(grid_evaluation_operator(b2, gridspace(b1), grid(b1))*e2)
      @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b1))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b1), grid(b1))*e1)
      @test coefficients(BasisFunctions.default_evaluation_operator(b1, gridspace(b2))*e1) ≈ coefficients(grid_evaluation_operator(b1, gridspace(b2), grid(b2))*e1)

      mr = matrix(restriction_operator(b1, b2))
      me = matrix(extension_operator(b2, b1))
      pinvme = pinv(me)
      r = rand(size(pinvme,2))
      @test pinvme*r ≈ mr*r
    end
  end

  @test_throws AssertionError restriction_operator(span(SymBSplineTranslatesBasis(4,1)), span(SymBSplineTranslatesBasis(3,1)))
  @test_throws AssertionError extension_operator(span(SymBSplineTranslatesBasis(4,1)), span(SymBSplineTranslatesBasis(6,1)))
end

function test_orthonormalsplinebasis(T)
  b = OrthonormalSplineBasis(5,2,Float64)
  b = OrthonormalSplineBasis(5,2,T)
  @test name(b) == "Set of translates of a function (B spline of degree 2) (orthonormalized)"
  @test instantiate(OrthonormalSplineBasis,5)==OrthonormalSplineBasis(5,3)

  G = sqrt(DualGram(span(b.superset)))
  e = zeros(eltype(G),size(G,1))
  e[1] = 1
  @test BasisFunctions.coefficients(b) ≈ G*e

  d = BasisFunctions.primalgramcolumn(span(b); atol=1e-3)
  @test d ≈ e
  @test typeof(Gram(span(b))) <: IdentityOperator

  n = 8
  for degree in 0:3
    b = span(OrthonormalSplineBasis(n, degree, T))
    basis_ext = extend(b)
    r = restriction_operator(basis_ext, b)
    e = extension_operator(b, basis_ext)
    @test eye(n) ≈ matrix(r*e)

    grid_ext = grid(basis_ext)
    L = evaluation_operator(b, grid_ext)
    e = random_expansion(b)
    z = L*e
    L2 = evaluation_operator(basis_ext, grid_ext) * extension_operator(b, basis_ext)
    z2 = L2*e
    @test 1+maximum(abs.(z-z2)) ≈ T(1)
  end
end

function test_discrete_orthonormalsplinebasis(T)
  b = DiscreteOrthonormalSplineBasis(5,2,Float64)
  b = DiscreteOrthonormalSplineBasis(5,2,T)
  @test name(b) == "Set of translates of a function (B spline of degree 2) (orthonormalized, discrete)"
  @test instantiate(DiscreteOrthonormalSplineBasis,5)==DiscreteOrthonormalSplineBasis(5,3)

  E = evaluation_operator(span(b))
  e = zeros(eltype(E),size(E,1))
  e[1] = 1
  @test (E'*E*e) ≈ 5*e
  @test typeof(E) <: CirculantOperator

  n = 8
  for degree in 0:3
    b = span(DiscreteOrthonormalSplineBasis(n, degree, T))
    basis_ext = extend(b)
    r = restriction_operator(basis_ext, b)
    e = extension_operator(b, basis_ext)
    @test eye(n) ≈ matrix(r*e)

    grid_ext = grid(basis_ext)
    L = evaluation_operator(b, grid_ext)
    e = random_expansion(b)
    z = L*e
    L2 = evaluation_operator(basis_ext, grid_ext) * extension_operator(b, basis_ext)
    z2 = L2*e
    @test 1+maximum(abs.(z-z2)) ≈ T(1)
  end
  for dgr in 0:4, oversampling in 1:4, n in 10:11
    b = DiscreteOrthonormalSplineBasis(n,dgr,T; oversampling=oversampling)

    e = zeros(T,n); e[1] = 1
    @test sqrt(DiscreteDualGram(span(BSplineTranslatesBasis(n,dgr,T)); oversampling=oversampling))*e≈BasisFunctions.coefficients(b)
    @test DiscreteGram(span(b); oversampling=oversampling)*e≈e
  end
end
using QuadGK
function test_dualsplinebasis(T)
  n = 10; degree = 2;
  tol = max(sqrt(eps(real(T))), 1e-16)
  b = BSplineTranslatesBasis(n,degree)
  bb = BasisFunctions.dual(b; rtol=tol, atol=tol)
  @test dual(bb) == b
  e = coefficients(random_expansion(span(b)))

  @test Gram(span(b); atol=tol, rtol=tol)*e ≈ DualGram(span(bb); atol=tol, rtol=tol)*e
  @test Gram(span(bb); atol=tol, rtol=tol)*e ≈ DualGram(span(b); atol=tol, rtol=tol)*e
  @test BasisFunctions.dualgramcolumn(span(b); rtol=tol, atol=tol) ≈ BasisFunctions.coefficients(bb)
  @test QuadGK.quadgk(x->b[1](x)*bb[1](x),left(b), right(b); rtol=tol, atol=tol)[1] - T(1) < sqrt(tol)
  @test QuadGK.quadgk(x->b[1](x)*bb[2](x),left(b), right(b); rtol=tol, atol=tol)[1] - T(1) < sqrt(tol)
end

function test_discrete_dualsplinebasis(T)
  for degree in 0:4, n in 10:11, oversampling=1:3
    b = BSplineTranslatesBasis(n,degree, T)
    d = discrete_dual(b; oversampling=oversampling)

    e = zeros(T,n); e[1] = 1
    @test DiscreteDualGram(span(b); oversampling=oversampling)*e≈BasisFunctions.coefficients(d)
    E1 = matrix(discrete_dual_evaluation_operator(span(b); oversampling=oversampling))
    oss = []
    # For even degrees points on a coarse grid do not overlap with those on e fine grid.
    if isodd(degree)
      oss = [1, oversampling]
    else
      oss = [oversampling,]
    end
    for os in oss
      E2 = matrix(evaluation_operator(span(d); oversampling = os))
      E2test = E1[1:Int(oversampling/os):end,:]
      @test E2test*e ≈ E2*e
    end
  end
end

# exit()
# using Base.Test
# using BasisFunctions
# @testset begin test_discrete_dualsplinebasis(Float64) end
#
# @testset begin test_dualsplinebasis(Float64) end
# @testset begin test_discrete_orthonormalsplinebasis(Float64) end
# @testset begin test_orthonormalsplinebasis(Float64) end
# @testset begin test_translatedbsplines(Float64) end
# @testset begin test_translatedsymmetricbsplines(Float64) end
# @testset begin test_generic_periodicbsplinebasis(Float64) end

# using Plots
# using BasisFunctions
# n = 7
# k = 0
# b = BSplineTranslatesBasis(n,k)
# t = linspace(-0,2,200)
# f = map(x->eval_element(b,1,x), t)
# using Plots
# plot!(t,f)
# @which eval_element(b,1,.25)
# BasisFunctions.fun(b)
# f = map(BasisFunctions.fun(b), t)
# plot!(t,f)
# f = map(x->BasisFunctions.eval_expansion(b,ones(n),x),t)
# # f = map(x->BasisFunctions.eval_expansion(b,[1,0,0,0,0,0,0,0,0,0],x),t)
# plot!(t,f,ylims=[-n,2n])
#
