
#####################
# Generic transforms
#####################

# A function set can have several associated transforms. The default transform is
# associated with the grid of the set, e.g. the FFT and the DCTII for Chebyshev expansions
# which convert between coefficient space and value space. In this case, the
# transform maps coefficients to or from a GridBasis.
#
# transform_operator takes two arguments, a source and destination set, in order
# to allow for different transforms.
#
# We assume that the transform itself is unitary (ASK up to a normalization). In order to compute an approximation
# to a function from function values, the transform is typically preceded and followed
# by an additional computation (e.g. the first Chebyshev coefficient is halved after the DCT).
# These additional operations are achieved by the operator returned by the
# transform_operator_pre and transform_operator_post functions: pre acts on the
# coefficients of the source space, post on the coefficients of the dest space.
#
# An important case is where the destination is a GridBasis. In that case,
# the routines transform_from_grid, transform_from_grid_pre and transform_from_grid_post
# are invoked, with the grid as a third argument. When the source is a GridBasis,
# similar routines with "to_grid" in their names are invoked.
# Dictionary's can intercept these functions and define their coefficient-to-value
# and value-to-coefficient transforms.


# The default transform space is the space associated with the grid of the set
transform_dict(s::Dictionary; options...) = gridbasis(s)

for op in (:transform_operator, :transform_operator_pre, :transform_operator_post)
    # With only one argument, use the default transform space
    @eval $op(src::Dictionary; options...) =
        $op(src, transform_dict(src; options...); options...)

    # Convert complexified spaces to uncomplexified ones
    @eval $(op)(src::ComplexifiedDict, dest::ComplexifiedDict; T = op_eltype(src, dest), options...) =
        $(op)(superdict(src), superdict(dest); T = T, options...)
end

# If the destination is a GridBasis, invoke the "to_grid" routines
for op in ( (:transform_operator, :transform_to_grid),
            (:transform_operator_pre, :transform_to_grid_pre),
            (:transform_operator_post, :transform_to_grid_post) )
    @eval $(op[1])(src::Dictionary, dest::GridBasis; options...) =
        $(op[2])(src, dest, grid(dest); options...)
    # Convenience function: convert a grid to a grid space
    @eval $(op[1])(src::Dictionary, grid::AbstractGrid; options...) =
        $(op[1])(src, gridbasis(grid, coefficienttype(src)); options...)
end

# If the source is a GridBasis, invoke the "from_grid" routines
for op in ( (:transform_operator, :transform_from_grid),
            (:transform_operator_pre, :transform_from_grid_pre),
            (:transform_operator_post, :transform_from_grid_post) )
    @eval $(op[1])(src::GridBasis, dest::Dictionary; options...) =
        $(op[2])(src, dest, grid(src); options...)
    # Convenience function: convert a grid to a grid space
    @eval $(op[1])(src::AbstractGrid, dest::Dictionary; options...) =
        $(op[1])(gridbasis(src, codomaintype(dest)), dest; options...)
end

# Pre and post operations are identity by default.
for op in (:transform_from_grid_pre, :transform_to_grid_pre)
    @eval $op(src, dest, grid; T = op_eltype(src,dest), options...) =
        IdentityOperator(T, src)
end

for op in (:transform_from_grid_post, :transform_to_grid_post)
    @eval $op(src, dest, grid; T = op_eltype(src,dest), options...) =
        IdentityOperator(T, dest)
end

# Return all three of them in a tuple
transform_operators(sets::Dictionary...; options...) =
    (transform_operator_pre(sets...; options...),
     transform_operator(sets...; options...),
     transform_operator_post(sets...; options...))

# Return the full operation: Post * Trans * Pre
function full_transform_operator(sets::Dictionary...; options...)
    Pre,T,Post = transform_operators(sets...; options...)
    Post * T * Pre
end
