# Assignment model

The model features an assignment model that can be used to build assignment
schemes and can be extended. As of writing, plans for the assignment model is
described in the next section.

## Overview

Main logic is done by the action of an `AbstractCriterion` on a set of arguments
to yield an `AbstractAssignment`. While `AbstractCriterion` objects are callable,
they are not subtypes of `Base.Callable`. So far, the only defined criteria
subtypes are `DiscreteCriterion` and `CompoundCriterion`.

The expectated behavior is that `AbstractAssignment` subtypes behave _similarly_
as [enum] "flags" and that the usual basic operations on integers: `+`, `-`,
`&`, and `|` can be applied.
Criterion types are parametrized by the assignment type they yield, and collections
of criteria with the same yielded assignment type can be "squashed" to create
compound criteria.

## Assignments and the assignment macro

Assignments can be treated as flags (assignment flags) and can be combined to
yield complex assignments. Assignment types are subtypes of `AbstractAssignment`,
and are expected to have property `val`. For example, the following is a minimal
valid definition:

```julia
struct DSSPAlphaType <: AbstractAssignment
    val::Int32
end
const NoAssign   = DSSPAlphaType(0)
const ShortTurn = DSSPAlphaType(1)
const AlphaTurn = DSSPAlphaType(2)
const PiTurn    = DSSPAlphaType(4)
```

The above code defined `DSSPAlphaType` as an assignment with integer field `val`.
Top-level constants are defined so that they can be used within the program/module.

!!! tip "Top-level assignments"
    
    It is recommended to declare assignments as top-level elements for clarity.

The above seven line declaration can be shortened to a more compact form using the
`ProtStructs.@Assignment` macro. The behavior of the macro is not yet final, and
hence is not exported:

```julia
ProtStructs.@Assignment DSSPAlphaType NoAssign=0 ShortTurn=1 AlphaTurn=2 PiTurn=4
```

## Discrete criteria

A very common way of assigning something to an object is that it is given when the
object follows a certain condition. Otherwise, a default assignment is given. Since
this criterion acts as a yes-no switch, it is called a _discrete criterion_.

As examples, suppose you want to know when a specific ``n``-turn is present.
To learn more about the ``n``-turn variants, check the docs for 
[nturn](../API.md#ProtStructs.nturn).
`DiscreteCriterion` objects are constructed as follows:

```julia
criterion1 = DiscreteCriterion(alphaturn, AlphaTurn, NoAssign)
criterion2 = DiscreteCriterion(shortturn, ShortTurn, NoAssign)
criterion3 = DiscreteCriterion(piturn, PiTurn, NoAssign)
```

## Squashing and compound criteria

The collection of criteria above can be "squashed" into a `CompoundCriterion` as
follows:

```julia
alphasquash = squash(criterion1, criterion2, criterion3)
```

This should work when the input criteria can accept the same numbers and types of
arguments.

## Future plans

The above described behavior is experimental, and will change depending on future
needs. Other features to be included are the following:
+ Support for blocks in assignment declaration using macro
+ Other forms of criteria aside from discrete and compound
+ Assignment blocks and layers

