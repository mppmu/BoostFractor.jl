# Implementation notes on dancer() function

The following outlines the indexing of the arrays.


```
Boundaries:         [ 1,       2,      3,      4,     ...                ]
               MIRROR |       |       |       |      
               MIRROR |       |       |       |        ...
               MIRROR |       |       |       |      
Regions:              [  1    ,   2   ,   3   ,   4   ,  ...                ]

```

Arrays refering to **boundaries** are:
```
# in general:
SetupBoundaries.r
```

Arrays refering to **regions** are:
```
# in general:
SetupBoundaries.distance
               .eps
               .relative_tilt_x
               .relative_tilt_y
               .relative_surfaces

# in dancer() function:
fields
```

Before the call of the propagator the `rightmoving` fields are those leaving a boundary to the right, yet upropagated. The `leftmoving` fields are those leaving a boundary to the left, yet upropagated. Therefore in every gap in each step the `rightmoving` fields go from left to right (left leg) and vice versa.

In each outer iteration (while loop, "measure of music") we do all propagations and reflection / transmissions for the whole system once.

In each sub-iteration (for loop, "legs") we iterate over _boundaries_. At each boundary we "pull" the fields to the boundary using the propagator. Afterwards we calculate the reflection and transmission.

To save memory the index of `rightmoving` and `leftmoving` is swapped after each iteration, such that the reflected fields stay in the same place in memory and only one copy for one transmitted field has to be made. This is done at the end of each iteration of the outer loop.
