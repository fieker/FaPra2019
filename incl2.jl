#############################
using Hecke,Primes,Markdown
mutable struct QuadForm
  a::fmpz
  b::fmpz
  c::fmpz
  function QuadForm(a, b, c)
    return new(fmpz(a), fmpz(b), fmpz(c))
  end
end

function Hecke.discriminant(f::QuadForm)
  return f.b^2 - 4*f.a*f.c
end
function Base.show(io::IO, f::QuadForm)
  print(io, "my first form: <$(f.a), $(f.b), $(f.c)>")
end

include("findforms.jl")
include("formalg.jl")
