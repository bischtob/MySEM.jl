abstract type AbstractEmbedding end

struct Embedding{𝒯} <: AbstractEmbedding
    manifold::𝒯
end

function CartesianEmbedding(; manifold)
    return Embedding(manifold)
end

function (embedding::CartesianEmbedding{Interval})(x)
    a = embedding.manifold.a
    b = embedding.manifold.b
    @assert a ≤ x ≤ b

    return [x]
end

function (embedding::CartesianEmbedding{Rectangle})(x, y)
    a = embedding.manifold.a
    b = embedding.manifold.b
    c = embedding.manifold.c
    d = embedding.manifold.d
    @assert a ≤ x ≤ b
    @assert c ≤ y ≤ d

    return [x, y]
end

function (embedding::SphericalEmbedding{Sphere})(λ, ϕ)
    @assert -π ≤ λ ≤ π
    @assert -π/2 ≤ ϕ ≤ π/2
    a = embedding.manifold.radius

    return [a*cos(λ)*sin(ϕ), ]
end