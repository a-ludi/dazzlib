/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.math;

import std.algorithm;
import std.traits;


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


/// Convert to given type without errors by bounding values to target type
/// limits.
IntTo boundedConvert(IntTo, IntFrom)(IntFrom value) pure nothrow @safe @nogc
    if (isIntegral!IntTo && isIntegral!IntFrom)
{
    static if (isSigned!IntFrom == isSigned!IntTo)
    {
        if (IntTo.min <= value && value <= IntTo.max)
            return cast(IntTo) value;
        else if (IntTo.min > value)
            return IntTo.min;
        else
            return IntTo.max;
    }
    else static if (isSigned!IntFrom)
    {
        static assert(isUnsigned!IntTo);

        if (value < 0)
            return IntTo.min;
        else if (cast(Unsigned!IntFrom) value < IntTo.max)
            return cast(IntTo) value;
        else
            return IntTo.max;
    }
    else
    {
        static assert(isUnsigned!IntFrom && isSigned!IntTo);

        if (value < cast(Unsigned!IntTo) IntTo.max)
            return cast(IntTo) value;
        else
            return IntTo.max;
    }
}

///
unittest
{
    assert((0).boundedConvert!uint == 0u);
    assert((42).boundedConvert!uint == 42u);
    assert((int.max).boundedConvert!uint == cast(uint) int.max);
    assert((-1).boundedConvert!uint == 0u);
}

unittest
{
    import std.meta;

    alias IntTypes = AliasSeq!(
        byte,
        ubyte,
        short,
        ushort,
        int,
        uint,
        long,
        ulong,
        //cent,
        //ucent,
    );

    static foreach (alias IntFrom; IntTypes)
        static foreach (alias IntTo; IntTypes)
        {
            assert(IntFrom(0).boundedConvert!IntTo == IntTo(0));
            assert(IntFrom(42).boundedConvert!IntTo == IntTo(42));
        }

    // IntFrom = byte
    assert(boundedConvert!byte      (byte.max) == byte.max);
    assert(boundedConvert!ubyte     (byte.max) == byte.max);
    assert(boundedConvert!short     (byte.max) == byte.max);
    assert(boundedConvert!ushort    (byte.max) == byte.max);
    assert(boundedConvert!int       (byte.max) == byte.max);
    assert(boundedConvert!uint      (byte.max) == byte.max);
    assert(boundedConvert!long      (byte.max) == byte.max);
    assert(boundedConvert!ulong     (byte.max) == byte.max);
    assert(boundedConvert!byte      (byte.min) == byte.min);
    assert(boundedConvert!ubyte     (byte.min) == ubyte.min);
    assert(boundedConvert!short     (byte.min) == byte.min);
    assert(boundedConvert!ushort    (byte.min) == ushort.min);
    assert(boundedConvert!int       (byte.min) == byte.min);
    assert(boundedConvert!uint      (byte.min) == uint.min);
    assert(boundedConvert!long      (byte.min) == byte.min);
    assert(boundedConvert!ulong     (byte.min) == ulong.min);

    // IntFrom = ubyte
    assert(boundedConvert!byte      (ubyte.max) == byte.max);
    assert(boundedConvert!ubyte     (ubyte.max) == ubyte.max);
    assert(boundedConvert!short     (ubyte.max) == ubyte.max);
    assert(boundedConvert!ushort    (ubyte.max) == ubyte.max);
    assert(boundedConvert!int       (ubyte.max) == ubyte.max);
    assert(boundedConvert!uint      (ubyte.max) == ubyte.max);
    assert(boundedConvert!long      (ubyte.max) == ubyte.max);
    assert(boundedConvert!ulong     (ubyte.max) == ubyte.max);

    // IntFrom = int
    assert(boundedConvert!byte      (int.max) == byte.max);
    assert(boundedConvert!ubyte     (int.max) == ubyte.max);
    assert(boundedConvert!short     (int.max) == short.max);
    assert(boundedConvert!ushort    (int.max) == ushort.max);
    assert(boundedConvert!int       (int.max) == int.max);
    assert(boundedConvert!uint      (int.max) == int.max);
    assert(boundedConvert!long      (int.max) == int.max);
    assert(boundedConvert!ulong     (int.max) == int.max);
    assert(boundedConvert!byte      (int.min) == byte.min);
    assert(boundedConvert!ubyte     (int.min) == ubyte.min);
    assert(boundedConvert!short     (int.min) == short.min);
    assert(boundedConvert!ushort    (int.min) == ushort.min);
    assert(boundedConvert!int       (int.min) == int.min);
    assert(boundedConvert!uint      (int.min) == uint.min);
    assert(boundedConvert!long      (int.min) == int.min);
    assert(boundedConvert!ulong     (int.min) == uint.min);

    // IntFrom = uint
    assert(boundedConvert!byte      (uint.max) == byte.max);
    assert(boundedConvert!ubyte     (uint.max) == ubyte.max);
    assert(boundedConvert!short     (uint.max) == short.max);
    assert(boundedConvert!ushort    (uint.max) == ushort.max);
    assert(boundedConvert!int       (uint.max) == int.max);
    assert(boundedConvert!uint      (uint.max) == uint.max);
    assert(boundedConvert!long      (uint.max) == uint.max);
    assert(boundedConvert!ulong     (uint.max) == uint.max);
}
