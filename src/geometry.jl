"""
"""
struct Position{𝒯}
  domain::𝒯
end

function Position(; domain)
  return Position(domain)
end

function (position::Position{Interval})(x)
  
  return [x]
end

function (position::Position{Circle})(λ)
  a = position.domain.radius
  return [a*cos(λ), a*sin(λ)]
end

function (position::Position{Rectangle})(x, y)
  return [x, y]
end

function (position::Position{Annulus})(λ, r)
  rmax = position.domain.subdomains.interval.max
  rmin = position.domain.subdomains.interval.min
  @assert rmin ≤ r ≤ rmax

  return [r*cos(λ), r*sin(λ)]
end

"""
"""
struct Metric{𝒯}
  domain::𝒯
end

function Metric(; domain)
  return Metric(domain)
end

function (g::Metric{Interval})(x)
  return [1]
end

function (g::Metric{Circle})(λ)
  a = g.domain.radius
  return [a^2]
end

function (g::Metric{Rectangle})(x, y)
  return [[1, 0] [0, 1]]
end

function (g::Metric{Annulus})(λ, r)
  rmax = g.domain.subdomains.interval.max
  rmin = g.domain.subdomains.interval.min
  @assert rmin ≤ r ≤ rmax

  return [[1 0] [0 r^2]]
end

"""
"""
struct VolumeElement{𝒯}
  domain::𝒯
end

function VolumeElement(; domain)
  return VolumeElement(domain)
end

function (dV::VolumeElement{Interval})(x)
  return 1
end

function (dV::VolumeElement{Circle})(λ)
  a = dV.domain.radius
  return a^2
end

function (dV::VolumeElement{Rectangle})(x, y)
  return 1
end

function (dV::VolumeElement{Annulus})(λ, r)
  rmax = dV.domain.subdomains.interval.max
  rmin = dV.domain.subdomains.interval.min
  @assert rmin ≤ r ≤ rmax

  return r^2
end
