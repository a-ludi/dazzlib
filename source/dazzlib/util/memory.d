/**
    Provides convenient function for working without the GC.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.memory;

import core.exception;
static import core.stdc.stdlib;
import core.stdc.string;


/// Requests an aligned block of *unmanaged* memory from the OS. If allocation
/// fails, this function will call `onOutOfMemory` which is expected to throw
/// an `OutOfMemoryError`.
T[] malloc(T)(size_t numElements) @nogc
{
    auto memoryPtr = core.stdc.stdlib.malloc(numElements * T.sizeof);
    if (memoryPtr is null)
        onOutOfMemoryError(memoryPtr);

    return (cast(T*) memoryPtr)[0 .. numElements];
}


/// Requests an aligned block of *unmanaged* memory from the OS. If allocation
/// fails, this function will call `onOutOfMemory` which is expected to throw
/// an `OutOfMemoryError`.
T* mallocObject(T, Args...)(Args args) @nogc
{
    auto object = malloc!T(1).ptr;
    *object = T(args);

    return object;
}


/// Returns a zero-terminated copy of `s` in *unmanaged* memory. If allocation
/// fails, this function will call `onOutOfMemory` which is expected to throw
/// an `OutOfMemoryError`.
char* stringzCopy(const(char)[] s) @nogc
{
    auto zstr = malloc!char(s.length + 1);
    memcpy(zstr.ptr, s.ptr, s.length);
    zstr[$ - 1] = 0;

    return zstr.ptr;
}
