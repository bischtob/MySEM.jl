"""
"""
struct PositionVector{𝒯}
  domain::𝒯
end

function PositionVector(; domain)
  return PositionVector(domain)
end

function (position::PositionVector{Interval})(x)
  return [x]
end

function (position::PositionVector{Circle})(λ)
  a = position.domain.radius
  return [a*cos(λ), a*sin(λ)]
end

function (position::PositionVector{Rectangle})(x, y)
  return [x, y]
end

function (position::PositionVector{Annulus})(λ, r)
  rmax = position.domain.subdomains.interval.max
  rmin = position.domain.subdomains.interval.min
  @assert rmin ≤ r ≤ rmax

  return [r*cos(λ), r*sin(λ)]
end

"""
"""
struct MetricTensor{𝒯}
  domain::𝒯
end

function MetricTensor(; domain)
  return MetricTensor(domain)
end

function (g::MetricTensor{Interval})(x)
  return [1]
end

function (g::MetricTensor{Circle})(λ)
  a = g.domain.radius
  return [a^2]
end

function (g::MetricTensor{Rectangle})(x, y)
  return [[1, 0] [0, 1]]
end

function (g::MetricTensor{Annulus})(λ, r)
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

function (dv::VolumeElement{Interval})(x)
  return 1
end

function (dv::VolumeElement{Circle})(λ)
  a = dV.domain.radius
  return a^2
end

function (dv::VolumeElement{Rectangle})(x, y)
  return 1
end

function (dv::VolumeElement{Annulus})(λ, r)
  rmax = dv.domain.subdomains.interval.max
  rmin = dv.domain.subdomains.interval.min
  @assert rmin ≤ r ≤ rmax

  return r^2
end