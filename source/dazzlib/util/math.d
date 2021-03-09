/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.math;

import std.algorithm;
import std.traits :
    isIntegral,
    isNumeric;


///
enum RoundingMode : byte
{
    ///
    floor,
    ///
    round,
    ///
    ceil,
}

/**
    Round x upward according to base, ie. returns the next integer larger or
    equal to x which is divisible by base.

    Returns: x rounded upward according to base.
*/
Integer ceil(Integer)(in Integer x, in Integer base) pure nothrow if (isIntegral!Integer)
{
    return x % base == 0
        ? x
        : (x / base + 1) * base;
}

///
unittest
{
    assert(ceil(8, 10) == 10);
    assert(ceil(32, 16) == 32);
    assert(ceil(101, 100) == 200);
}

/**
    Round x downward according to base, ie. returns the next integer smaller or
    equal to x which is divisible by base.

    Returns: x rounded downward according to base.
*/
Integer floor(Integer)(in Integer x, in Integer base) pure nothrow if (isIntegral!Integer)
{
    return (x / base) * base;
}

///
unittest
{
    assert(floor(8, 10) == 0);
    assert(floor(32, 16) == 32);
    assert(floor(101, 100) == 100);
}

/// Returns the absolute difference between two numbers.
Num absdiff(Num)(in Num a, in Num b) pure nothrow if (isNumeric!Num)
{
    return a > b
        ? a - b
        : b - a;
}

///
unittest
{
    assert(absdiff(2UL, 3UL) == 1UL);
    assert(absdiff(-42, 13) == 55);
    assert(absdiff(2.5, 5) == 2.5);
}


/// Returns the result of `ceil(a / b)` but uses integer arithmetic only.
Integer ceildiv(Integer)(in Integer a, in Integer b) pure nothrow if (isIntegral!Integer)
{
    Integer resultSign = (a < 0) ^ (b < 0) ? -1 : 1;

    return resultSign < 0 || a % b == 0
        ? a / b
        : a / b + resultSign;
}

///
unittest
{
    assert(ceildiv(0, 3) == 0);
    assert(ceildiv(1UL, 3UL) == 1UL);
    assert(ceildiv(2L, 3L) == 1L);
    assert(ceildiv(3U, 3U) == 1U);
    assert(ceildiv(4, 3) == 2);
    assert(ceildiv(-4, 4) == -1);
    assert(ceildiv(-4, 3) == -1);
}
