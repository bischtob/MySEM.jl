abstract type Domain end

"""
"""
struct Interval{𝒯} <: Domain
  min::𝒯
  max::𝒯
  periodicity::Bool
end

"""
"""
struct Circle{𝒯} <: Domain
  radius::𝒯
  periodicity::true
end

"""
"""
struct ProductDomain{𝒯} <: Domain
  subdomains::𝒯
end

"""
"""
function Interval(; min, max, periodic)
  @assert min < max
  return Interval(min, max, periodic)
end

"""
"""
function Circle(; radius)
  @assert radius > 0
  return Circle(radius)
end

"""
"""
function Rectangle(; a, b, c, d)
  interval₁ = Interval(min = a, max = b)
  interval₂ = Interval(min = c, max = d)
  return ProductDomain((interval₁, interval₂)) 
end

"""
"""
function Annulus(; radius, height)
  circle = Circle(radius = radius)
  interval = Interval(min = radius, max = radius + height)
  return ProductDomain((annulus, interval))
end

"""
"""
function isperiodic(domain::Domain)
  return domain.periodicity
end

function isperiodic(domain::ProductDomain)
  return [(typeof(item), isperiodic(item)) for item in domain.subdomains]
end

"""
"""
function ndims(::Domain)
  return 1
end

function ndims(domain::ProductDomain)
  return length(domain.subdomains)
end
