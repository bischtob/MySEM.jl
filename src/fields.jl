abstract type Field end
abstract type DiscreteField end

struct ScalarField <: Field end
struct VectorField <: Field end
struct MatrixField <: Field end

struct DiscretizedField <: Field
  
end

