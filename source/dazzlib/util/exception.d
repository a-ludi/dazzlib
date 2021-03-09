/**
    Defines execeptions for use in dazzlib.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.exception;

import dazzlib.core.c.DB;
import std.exception;


///
class DazzlibException : Exception
{
    ///
    mixin basicExceptionCtors;
}


///
alias dazzlibEnforce = enforce!DazzlibException;


///
void catchErrorMessage(alias work, Args...)(Args args)
{
    clearCurrentError();
    work(args);
    dazzlibEnforce(!hasError, currentError);
}
