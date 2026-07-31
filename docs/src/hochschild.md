```@meta
DocTestSetup = :(using AbstractAlgebra, HodgeDiamonds)
```

# Hochschild homology

```@docs
HochschildHomology
from_list
from_positive
from_polynomial(::AbstractAlgebra.LaurentPolyRingElem)
polynomial(::HochschildHomology)
dimension(::HochschildHomology)
euler(::HochschildHomology)
Base.getindex(::HochschildHomology, ::Integer)
symmetric_power(::HochschildHomology, ::Integer)
sym
```
