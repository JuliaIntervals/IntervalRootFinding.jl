# IA v1.0 comments

## Disparitions

- `where_bisect`
- `Region`
- `contains_zero`

## Weird errors

```
julia> intersect(X, X)
ERROR: ArgumentError: `union!` is purposely not supported for intervals. See instead `hull`
```

```
julia> 0 in X
ERROR: ArgumentError: `==` is purposely not supported for intervals. See instead `isequal_interval` `hull`
```

```
julia> isinf(X)
ERROR: ArgumentError: `isnan` is purposely not supported for intervals. See instead `isnai`
```