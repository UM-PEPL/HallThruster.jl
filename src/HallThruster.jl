module HallThruster

export func

"""
    func(x)

Returns double the number `x` plus `1`.

```jldoctest;setup = :(using HallThruster)
func(1)

# output

3
```
"""
func(x) = 2x + 1

end # module
