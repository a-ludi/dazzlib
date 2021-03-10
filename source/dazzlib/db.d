/**
    High-level access to DB/DAM data.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.db;

import dazzlib.basictypes;
import dazzlib.core.c.DB;
import dazzlib.util.exception;
import dazzlib.util.math;
import std.algorithm;
import std.format;
import std.path;
import std.stdio;
import std.string;
import std.typecons;

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


///
alias TrimDb = Flag!"trimDb";


class DazzDb
{
    public alias Flag = DAZZ_DB.Flag;

    enum DbType : byte
    {
        undefined = -1,
        db = 0,
        dam = 1,
    }

    private DAZZ_DB dazzDb;
    private DbType _dbType;


    /// Construct from `dbFile` by opening and optionally trimming the DB.
    /// Errors in the underlying C routines will be promoted to
    /// `DazzlibException`s.
    ///
    /// Throws: DazzlibException on errors
    this(string dbFile, TrimDb trimDb = Yes.trimDb)
    {
        auto result = Open_DB(dbFile.toStringz, &this.dazzDb);
        dazzlibEnforce(result >= 0, currentError.idup);

        this._dbType = cast(DbType) result;

        if (trimDb)
            catchErrorMessage!Trim_DB(&dazzDb);
    }


    /// Construct from existing `DAZZ_DB` object and invalidate the passed
    /// object. This will leave `dbType` `undefined`.
    this(ref DAZZ_DB dazzDb)
    {
        this.dazzDb = dazzDb;
        // invalidate original DAZZ_DB object
        dazzDb = DAZZ_DB();
    }


    ~this()
    {
        Close_DB(&dazzDb);
    }


    /// Set to the type of DB dertermined by Open_DB.
    @property DbType dbType() const pure nothrow @safe
    {
        return _dbType;
    }


    /// Total number of reads in untrimmed DB
    @property id_t numReadsUntrimmed() const pure @safe
    {
        return dazzDb.ureads.boundedConvert!(typeof(return));
    }

    /// Total number of reads in trimmed DB
    @property id_t numReadsTrimmed() const pure @safe
    {
        return dazzDb.treads.boundedConvert!(typeof(return));
    }

    /// Minimum read length in block (-1 if not yet set)
    @property arithmetic_t numReadsUntrimmed() const pure nothrow @safe
    {
        return dazzDb.cutoff;
    }

    /// DB_ALL | DB_ARROW
    @property BitFlags!Flag flags() const pure nothrow @safe
    {
        return typeof(return)(dazzDb.allarr);
    }

    /// frequency of A, C, G, T, respectively
    @property float[4] baseFrequency() const pure nothrow @safe
    {
        return dazzDb.freq;
    }

    /// length of maximum read (initially over all DB)
    @property coord_t maxReadLength() const pure nothrow @safe
    {
        return dazzDb.maxlen.boundedConvert!(typeof(return));
    }

    /// total # of bases (initially over all DB)
    @property size_t totalBps() const pure nothrow @safe
    {
        return dazzDb.totlen.boundedConvert!(typeof(return));
    }

    /// # of reads in actively loaded portion of DB
    @property id_t numReads() const pure nothrow @safe
    {
        return dazzDb.nreads.boundedConvert!(typeof(return));
    }

    /// DB has been trimmed by cutoff/all
    @property bool isTrimmed() const pure nothrow @safe
    {
        return dazzDb.trimmed != 0;
    }

    /// DB block (if > 0), total DB (if == 0)
    @property id_t block() const pure nothrow @safe
    {
        return dazzDb.part.boundedConvert!(typeof(return));
    }

    /// Index of first read in block (without trimming)
    @property id_t firstReadUntrimmedPtr() const pure nothrow @safe
    {
        return dazzDb.ufirst.boundedConvert!(typeof(return));
    }

    /// Index of first read in block (with trimming)
    @property id_t firstReadTrimmedPtr() const pure nothrow @safe
    {
        return dazzDb.tfirst.boundedConvert!(typeof(return));
    }

    /// Root name of DB for .bps, .qvs, and tracks
    @property const(char)[] dbName() const pure
    {
        return dazzDb.path.fromStringz;
    }

    /// Are reads loaded in memory?
    @property bool areReadsLoaded() const pure nothrow @safe
    {
        return dazzDb.loaded != 0;
    }

    /// Array of DAZZ_READ
    @property DAZZ_READ[] reads() pure nothrow @trusted
    {
        return dazzDb.reads[0 .. numReads];
    }

    /// Linked list of loaded tracks
    @property DAZZ_TRACK* tracksPtr() pure nothrow @safe
    {
        return dazzDb.tracks;
    }
}


/// Validate DB by opening the DB once.
///
/// Returns: `null` if DB is valid; otherwise error message.
string validateDb(string dbFile, Flag!"allowBlock" allowBlock)
{
    try
    {
        auto dazzDb = new DazzDb(dbFile, No.trimDb);

        if (!allowBlock && dazzDb.block != 0)
            return "operation not allowed on a block";

        if (dazzDb.dbType == DazzDb.DbType.dam)
            // check for presence of headers file
            cast(void) File(EssentialDbFiles(dbFile).headers, "r");

        destroy(dazzDb);

        return null;
    }
    catch (Exception e)
    {
        return e.msg;
    }
}
