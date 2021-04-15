abstract type Manifold end
abstract type FlatManifold end
abstract type CurvedManifold end

"""
  Interval
"""
struct Interval{𝒯} <: FlatManifold
  a::𝒯
  b::𝒯
  periodicity::Bool
end

function Interval(; a, b, periodicity)
  @assert a < b
  return Interval(a, b, periodicity)
end

"""
  Circle
"""
struct Circle{𝒯} <: CurvedManifold
  radius::𝒯
end


function Circle(; radius)
  return Circle(radius = radius)
end

"""
  isperiodic
"""
isperiodic(::Manifold) = error("This manifold type is not implemented.")
isperiodic(manifold::Interval) = manifold.periodicity
isperiodic(::Circle) = true

"""
  isflat
"""
isflat(::Manifold) = error("This manifold type is not implemented.")
isflat(::FlatManifold) = true
isflat(::CurvedManifold) = false
iscurved(manifold::Manifold) = !isflat(manifold)

"""
  dim
"""
dim(::Manifold) = error("This manifold type is not implemented.")
dim(::Interval) = 1
dim(::Circle) = 1

"""
  embedding_dim
"""
embedding_dim(::Manifold) = error("This manifold type is not implemented.")
embedding_dim(manifold::FlatManifold) = dim(manifold)
embedding_dim(manifold::Circle) = 2