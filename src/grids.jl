function create_grid(Np, Npoin_cg, Ne, x, DFloat)

    coord_cg=zeros(DFloat,Npoin_cg)
    # intma_cg=zeros(Int64,Np,Ne)
    # periodicity_cg=zeros(Int64,Npoin_cg)

    #Set some constants
    #Always assume that the domain is divided in equispaced intervals
    xmin=-1
    xmax=+1
    dx=(xmax-xmin)/Ne

    #Generate Grid Points
    # I=1
    # coord_cg[1]=xmin
    # for e=1:Ne
    #     x0=xmin + (e-1)*dx
    #     for i=2:Np
    #         I+=1
    #         coord_cg[I]= x0 + 0.5 * ( x[i] + 1 ) * dx
    #     end
    # end

    grid_points = copy(coord_cg)
    nelements = Ne
    npoints = Np
    ξ = copy(x)
    x_min = xmin

    grid_points[1] = x_min
    for i in 1:nelements
      # calculate left element edge location
      # and index offset
      i_left = (i - 1) * (npoints - 1)
      x_left = (i - 1) * dx + x_min
      
      # loop of number of interpolation points
      # within each element starting at the second
      # point to avoid duplication
      for j in 2:npoints
        println((i, i_left, j))
        # calculate interpolation point location
        # relative to left element edge
        x_relative = 0.5 * (ξ[j] + 1) * dx
        grid_points[i_left + j] = x_left + x_relative
      end
    end

    coord_cg = copy(grid_points)

    # #Generate INTMA
    # I=1
    # for e=1:Ne
    #     intma_cg[1,e]=I
    #     for i=2:Np
    #         I+=1
    #         intma_cg[i,e]=I
    #     end
    # end

    # #Generate Periodicity_CG array
    # for I=1:Npoin_cg
    #     periodicity_cg[I]=I
    # end
    # periodicity_cg[Npoin_cg]=1

    return coord_cg#, intma_cg, periodicity_cg)
end

# Np = poly_order + 1 (dof)
# Npoin_cg = num_elements * poly_order + 1
# ne = num_elements
# x = interpolation points on reference element (-1,1) there should be P+1

# |. . . , . . . , . . .| = Ne * (P + 1) - (Ne-1) = Ne * P + Ne - Ne + 1 = Ne * P + 1

function lobatto_gauss(P::Integer,DFloat)

    xgl=zeros(DFloat,P)
    wgl=zeros(DFloat,P)

    p=P-1; #Order of the Polynomials
    ph=floor(typeof(P), (p+1)/2 )

    for i=1:ph
        x=cos( (2*i-1)*pi/(2*p+1) )
        for k=1:20
            (L0,L0_1,L0_2)=legendre_poly(p,x); #Compute Nth order Derivatives of Legendre Polys
            #Get new Newton Iteration
            dx=-(1-x^2)*L0_1/(-2*x*L0_1 + (1-x^2)*L0_2)
            x=x+dx
            if (abs(dx) < 1.0e-20)
                break
            end
        end
        xgl[p+2-i]=x
        wgl[p+2-i]=2/(p*(p+1)*L0^2)
    end

    #Check for Zero Root
    if (p+1 != 2*ph)
        x=0
        (L0,L0_1,L0_2)=legendre_poly(p,x)
        xgl[ph+1]=x
        wgl[ph+1]=2/(p*(p+1)*L0^2)
    end

    #Find remainder of roots via symmetry
    for i=1:ph
        xgl[i]=-xgl[p+2-i]
        wgl[i]=+wgl[p+2-i]
    end
    return xgl,wgl
end #function

function legendre_poly(p,x)

    L1, L1_1, L1_2 = 0,0,0
    L0, L0_1, L0_2 = 1,0,0

    for i=1:p
        L2=L1;L2_1=L1_1;L2_2=L1_2
        L1=L0;L1_1=L0_1;L1_2=L0_2
        a=(2*i-1)/i
        b=(i-1)/i
        L0=a*x*L1 - b*L2
        L0_1=a*(L1 + x*L1_1) - b*L2_1
        L0_2=a*(2*L1_1 + x*L1_2) - b*L2_2
    end
    return (L0,L0_1,L0_2)
end

"""
  Domain
"""
abstract type AbstractDomain end

struct IntervalDomain{𝒯} <: AbstractDomain
    min::𝒯
    max::𝒯
end

function IntervalDomain(; min, max)
    return IntervalDomain(min, max)
end

"""
  Topology
  "Informally, a topology tells how elements of a set of points relate spatially to each other." (wikipedia.org)
"""
abstract type AbstractElementTopology end

struct IntervalTopology{𝒯} <: AbstractElementTopology
    interface_locations::𝒯
end

function RegularIntervalTopology(; domain::IntervalDomain, nelements::Int)
    min = domain.min
    max = domain.max
    interface_location = Array(range(min, max, length = nelments+1))

    return IntervalTopology(interface_locations)
end

"""
  Grid
"""
abstract type AbstractGrid end

struct SpectralElementGrid{𝒯,𝒰,𝒱} <: AbstractGrid
    points::𝒯
    #periodicity
    #storage_map
    quadrature_style::𝒰
    continuity_style::𝒱
end

function SpectralElementGrid(; topology, quadrature_style, continuity_style)
    points = _create_grid_points(continuity_style, topology, quadrature_style)
    return SpectralElementGrid(points, quadrature_style, continuity_style)
end

function _create_grid_points(::ContinuousStyle, topology::RegularInterval, quadrature_style::𝒯) where {𝒯}
    # Unpack relevant parameters
    x_faces   = topology.interface_locations
    nelements = length(topology.interface_locations) - 1
    npoints   = quadrature_style.polynomial_order + 1

    # Create interpolation points
    ξ = create_interpolation_points(quadrature_style)

    grid_points[1] = x_min
    for i in 1:nelements
      # calculate left element edge location
      # and index offset
      dx     = x_faces[i+1] - x_faces[i]
      i_left = (i - 1) * (npoints - 1)
      x_left = (i - 1) * dx + x_min
      
      # loop of number of interpolation points
      # within each element starting at the second
      # point to avoid duplication
      for j in 2:npoints
        println((i, i_left, j))
        # calculate interpolation point location
        # relative to left element edge
        x_relative = 0.5 * (ξ[j] + 1) * dx
        grid_points[i_left + j] = x_left + x_relative
      end
    end

    return grid_points
end