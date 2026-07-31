```@meta
DocTestSetup = :(using AbstractAlgebra, HodgeDiamonds)
```

# Hochschild homology

```@docs
HochschildHomology
from_positive
polynomial(::HochschildHomology)
dimension(::HochschildHomology)
euler(::HochschildHomology)
Base.getindex(::HochschildHomology, ::Integer)
symmetric_power(::HochschildHomology, ::Integer)
sym
```
