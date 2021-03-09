/**
    High-level access to DB/DAM data.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.db;

import dazzlib.core.c.DB;
import std.algorithm;
import std.format;
import std.path;

/// File suffixes of essential auxiliary .db files.
enum auxiliaryDbFileSuffixes = [".bps", ".idx"];

/// File suffixes of essential auxiliary .dam files.
enum auxiliaryDamFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .db file extension.
enum dbFileExtension = ".db";

/// Constant holding the .dam file extension.
enum damFileExtension = ".dam";

/// The Dazzler tools require sequence of a least minSequenceLength base pairs.
enum minSequenceLength = 14;


/// Structure that holds the names of essential files associated with a DB/DAM.
struct EssentialDbFiles
{
    string stub;
    string basePairs;
    string index;
    string headers;

    this(string stub) pure nothrow @safe
    {
        this.stub = stub;
        assert(isDb || isDam, "must use with Dazzler DB");

        this.basePairs = auxiliaryFile(".bps");
        this.index = auxiliaryFile(".idx");
        this.headers = isDam?  auxiliaryFile(".hdr") : null;
    }


    ///
    @property bool isDb() const pure nothrow @safe
    {
        try
        {
            return stub.endsWith(dbFileExtension);
        }
        catch (Exception e)
        {
            assert(0, "unexpected exception: " ~ e.msg);

            return false;
        }
    }

    ///
    unittest
    {
        assert(EssentialDbFiles("/path/to/test.db").isDb);
        assert(!EssentialDbFiles("/path/to/test.dam").isDb);
    }


    ///
    @property bool isDam() const pure nothrow @safe
    {
        try
        {
            return stub.endsWith(damFileExtension);
        }
        catch (Exception e)
        {
            assert(0, "unexpected exception: " ~ e.msg);

            return false;
        }
    }

    ///
    unittest
    {
        assert(!EssentialDbFiles("/path/to/test.db").isDam);
        assert(EssentialDbFiles("/path/to/test.dam").isDam);
    }


    /// Return directory part of DB stub.
    @property string dbdir() const pure nothrow @safe
    {
        return stub.dirName;
    }

    ///
    unittest
    {
        assert(EssentialDbFiles("/path/to/test.db").dbdir == "/path/to");
    }


    /// Return base name part of DB stub without extension.
    @property string dbname() const pure nothrow @safe
    {
        return stub.baseName.stripExtension;
    }

    ///
    unittest
    {
        assert(EssentialDbFiles("/path/to/test.db").dbname == "test");
    }


    /// Return auxiliary file with given suffix.
    string auxiliaryFile(string suffix) const pure nothrow @safe
    {
        return buildPath(dbdir, auxiliaryDbFilePrefix ~ dbname ~ suffix);
    }

    ///
    unittest
    {
        assert(
            EssentialDbFiles("/path/to/test.db").auxiliaryFile(".tan.anno") ==
            "/path/to/.test.tan.anno"
        );
    }


    /// Return an array of all files in alphabetical order (stub is always
    /// first).
    string[] list() const pure nothrow @safe
    {
        if (isDb)
            return [stub, basePairs, index];
        else
            return [stub, basePairs, headers, index];
    }
}

///
unittest
{
    auto dbFiles = EssentialDbFiles("/path/to/test.db");

    assert(dbFiles.stub == "/path/to/test.db");
    assert(dbFiles.basePairs == "/path/to/.test.bps");
    assert(dbFiles.index == "/path/to/.test.idx");
}

///
unittest
{
    auto dbFiles = EssentialDbFiles("/path/to/test.dam");

    assert(dbFiles.stub == "/path/to/test.dam");
    assert(dbFiles.basePairs == "/path/to/.test.bps");
    assert(dbFiles.headers == "/path/to/.test.hdr");
    assert(dbFiles.index == "/path/to/.test.idx");
}
