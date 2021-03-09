/**
    Safe-guarded IO operations.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.safeio;

import std.exception;
import std.stdio;


/// Like File.rawRead but enforces that buffer is completely filled.
void rawReadAll(T)(File file, T[] buffer, lazy string msg = "invalid read: premature end of file")
{
    enforce(file.rawRead(buffer).length == buffer.length, msg);
}


/// Like File.rawRead but reads just a scalar value and enforces that this
/// value is actually read.
void rawReadScalar(T)(File file, ref T scalar, lazy string msg = "invalid read: premature end of file")
{
    rawReadAll(file, (&scalar)[0 .. 1], msg);
}
