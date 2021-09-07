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
import dazzlib.util.safeio;
import std.algorithm;
import std.ascii;
import std.conv;
import std.format;
import std.mmfile;
import std.path;
import std.range;
import std.stdio;
import std.string;
import std.typecons;


///
public import dazzlib.core.c.DB : TrackKind;

///
public import dazzlib.core.c.DB : SequenceFormat;


/// Collection of file extensions for DB-related files.
struct DbExtension
{
    enum string db = ".db";
    enum string dam = ".dam";
    enum string basePairs = ".bps";
    enum string index = ".idx";
    enum string headers = ".hdr";
}

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

        this.basePairs = auxiliaryFile(DbExtension.basePairs);
        this.index = auxiliaryFile(DbExtension.index);
        this.headers = isDam?  auxiliaryFile(DbExtension.headers) : null;
    }


    ///
    @property bool isDb() const pure nothrow @safe
    {
        try
        {
            return stub.endsWith(DbExtension.db);
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
            return stub.endsWith(DbExtension.dam);
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


    ///
    string trackAnnotationFile(string trackName, id_t block = 0)
    {
        if (block == 0)
            return auxiliaryFile(format!".%s.anno"(trackName));
        else
            return auxiliaryFile(format!".%d.%s.anno"(block, trackName));
    }


    ///
    string trackDataFile(string trackName, id_t block = 0)
    {
        if (block == 0)
            return auxiliaryFile(format!".%s.data"(trackName));
        else
            return auxiliaryFile(format!".%d.%s.data"(block, trackName));
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


///
class DazzDb
{
    ///
    public alias Flag = DAZZ_DB.Flag;

    ///
    enum DbType : byte
    {
        undefined = -1,
        db = 0,
        dam = 1,
    }

    private DAZZ_DB dazzDb;
    private string _dbFile;
    private DbType _dbType;
    private DazzTrack[string] _trackIndex;


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

        if (dbFile.endsWith(DbExtension.dam, DbExtension.db))
            this._dbFile = dbFile;
        else
            this._dbFile = dbType == DbType.dam?  DbExtension.dam : DbExtension.db;

        if (trimDb)
            catchErrorMessage!Trim_DB(&dazzDb);
    }


    /// Construct from existing `DAZZ_DB` object and invalidate the passed
    /// object. This will leave `dbType` `undefined`.
    protected this(ref DAZZ_DB dazzDb)
    {
        this.dazzDb = dazzDb;
        // invalidate original DAZZ_DB object
        dazzDb = DAZZ_DB();
    }


    ~this()
    {
        Close_DB(&dazzDb);
    }


    ///
    @property EssentialDbFiles files() const pure nothrow @safe
    {
        return EssentialDbFiles(_dbFile);
    }


    /// Set to the type of DB dertermined by Open_DB.
    @property DbType dbType() const pure nothrow @safe @nogc
    {
        return _dbType;
    }


    /// Total number of reads in untrimmed DB
    @property id_t numReadsUntrimmed() const pure @safe @nogc
    {
        return dazzDb.ureads.boundedConvert!(typeof(return));
    }

    /// Total number of reads in trimmed DB
    @property id_t numReadsTrimmed() const pure @safe @nogc
    {
        return dazzDb.treads.boundedConvert!(typeof(return));
    }

    /// Minimum read length in block (-1 if not yet set)
    @property arithmetic_t cutoff() const pure nothrow @safe @nogc
    {
        return dazzDb.cutoff;
    }

    /// DB_ALL | DB_ARROW
    @property BitFlags!Flag flags() const pure nothrow @safe @nogc
    {
        return typeof(return)(dazzDb.allarr);
    }

    /// frequency of A, C, G, T, respectively
    @property float[4] baseFrequency() const pure nothrow @safe @nogc
    {
        return dazzDb.freq;
    }

    /// length of maximum read (initially over all DB)
    @property coord_t maxReadLength() const pure nothrow @safe @nogc
    {
        return dazzDb.maxlen.boundedConvert!(typeof(return));
    }

    /// total # of bases (initially over all DB)
    @property size_t totalBps() const pure nothrow @safe @nogc
    {
        return dazzDb.totlen.boundedConvert!(typeof(return));
    }

    /// # of reads in actively loaded portion of DB
    @property id_t numReads() const pure nothrow @safe @nogc
    {
        return dazzDb.nreads.boundedConvert!(typeof(return));
    }

    /// DB has been trimmed by cutoff/all
    @property bool isTrimmed() const pure nothrow @safe @nogc
    {
        return dazzDb.trimmed != 0;
    }

    /// DB block (if > 0), total DB (if == 0)
    @property id_t block() const pure nothrow @safe @nogc
    {
        return dazzDb.part.boundedConvert!(typeof(return));
    }

    /// Index of first read in block (without trimming)
    @property id_t firstReadUntrimmedPtr() const pure nothrow @safe @nogc
    {
        return dazzDb.ufirst.boundedConvert!(typeof(return));
    }

    /// Index of first read in block (with trimming)
    @property id_t firstReadTrimmedPtr() const pure nothrow @safe @nogc
    {
        return dazzDb.tfirst.boundedConvert!(typeof(return));
    }

    /// Root name of DB for .bps, .qvs, and tracks
    @property const(char)[] dbName() const pure
    {
        return dazzDb.path.fromStringz;
    }

    /// Are reads loaded in memory?
    @property bool areReadsLoaded() const pure nothrow @safe @nogc
    {
        return dazzDb.loaded != 0;
    }


    /// Allocate and return a buffer big enough for the largest read in the
    /// DB, leaving room for an initial delimiter character.
    auto makeReadBuffer() const pure nothrow @safe
    {
        return makeReadBuffer(maxReadLength);
    }


    /// Allocate and return a buffer for the given number of base pairs,
    /// leaving room for an initial delimiter character.
    static auto makeReadBuffer(size_t numBaseBairs) pure nothrow @safe
    {
        auto fullBuffer = new char[numBaseBairs + 4UL];

        // make value at index -1 accessible (required by DAZZLER routines)
        return fullBuffer[1 .. $];
    }

    /// Load the i'th (sub)read in this DB. A new buffer will be allocated if
    /// not given.
    ///
    /// Note, the byte at `buffer.ptr - 1` will be set to a delimiter
    /// character and must be allocated! Use `makeReadBuffer` for that.
    char[] loadRead(
        id_t readId,
        SequenceFormat seqFormat = SequenceFormat.init,
        char[] buffer = [],
    ) const
    {
        const readLength = reads[readId - 1].rlen.boundedConvert!coord_t;
        if (buffer.length == 0)
            buffer = makeReadBuffer(readLength);

        const result = Load_Read(&dazzDb, readId - 1, buffer.ptr, seqFormat);
        dazzlibEnforce(result == 0, currentError.idup);

        return buffer[0 .. readLength];
    }

    /// ditto
    char[] loadRead(
        id_t readId,
        coord_t begin,
        coord_t end,
        SequenceFormat seqFormat = SequenceFormat.init,
        char[] buffer = [],
    ) const
    in (begin < end, format!"Selected empty read slice: %d .. %d"(begin, end))
    in (end <= reads[readId - 1].rlen, format!"end is out of bounds: %d > %d"(end, reads[readId - 1].rlen))
    {
        const subreadLength = end - begin;
        if (buffer.length == 0)
            buffer = makeReadBuffer(subreadLength);

        const result = Load_Subread(&dazzDb, readId - 1, begin, end, buffer.ptr, seqFormat);
        dazzlibEnforce(result !is null, currentError.idup);
        assert(
            buffer.ptr <= result && result < buffer.ptr + 4,
            "sequence pointer out of bounds",
        );
        const offset = result - buffer.ptr;

        return buffer[offset .. subreadLength + offset];
    }

    unittest
    {
        import dazzlib.util.tempfile;
        import dazzlib.util.testdata;
        import std.exception;
        import std.file;
        import std.algorithm;

        auto tmpDir = mkdtemp("./.unittest-XXXXXX");
        scope (exit)
            rmdirRecurse(tmpDir);

        auto dbFile = buildPath(tmpDir, "test.db");

        writeTestDb(dbFile);

        auto dazzDb = new DazzDb(dbFile);

        foreach (i, expSequence; testSequences)
            assert(dazzDb.loadRead(cast(id_t) i + 1, SequenceFormat.asciiLower) == expSequence);

        foreach (i, expSequence; testSequences)
        {
            const begin = cast(coord_t) (i % 2 == 0 ? i                  : 0);
            const end   = cast(coord_t) (i % 2 == 0 ? expSequence.length : i);

            assert(dazzDb.loadRead(
                cast(id_t) i + 1,
                begin,
                end,
                SequenceFormat.asciiLower,
            ) == expSequence[begin .. end]);
        }
    }

    /// Array of DAZZ_READ
    @property inout(DAZZ_READ)[] reads() inout pure nothrow @trusted @nogc
    {
        return dazzDb.reads[0 .. numReads];
    }

    ///
    auto tracks() pure nothrow @safe @nogc
    {
        return _trackIndex.byValue();
    }


    /// Get track with trackName wrapped in a DazzTrack. This will open the
    /// track if it was not opened. The track data is not loaded automatically.
    DazzTrack getTrack(DataType = int)(string trackName)
    out(dazzTrack; dazzTrack !is null, "track is not of expected type")
    {
        if (trackName in _trackIndex)
        {
            auto dazzTrack = _trackIndex[trackName];

            final switch (dazzTrack.type)
            {
                case DazzTrack.Type.annotation:
                    return cast(AnnotationDazzTrack!DataType) dazzTrack;
                case DazzTrack.Type.mask:
                    return cast(MaskTrack) dazzTrack;
                case DazzTrack.Type.intData:
                    return cast(DataDazzTrack!(int32, DataType)) dazzTrack;
                case DazzTrack.Type.longData:
                    return cast(DataDazzTrack!(int64, DataType)) dazzTrack;
                case DazzTrack.Type.unsupported:
                    return null;
            }
        }

        auto trackPtr = Open_Track(&dazzDb, trackName.toStringz);
        dazzlibEnforce(trackPtr !is null, currentError.idup);

        auto rawSize = readTrackHead(trackName).recordSize;
        auto dazzTrackType = DazzTrack.type(trackPtr, rawSize);
        DazzTrack dazzTrack = dazzTrackType.predSwitch(
            DazzTrack.Type.annotation, cast(DazzTrack) new AnnotationDazzTrack!DataType(trackPtr, rawSize),
            DazzTrack.Type.mask, cast(DazzTrack) new MaskTrack(trackPtr, rawSize),
            DazzTrack.Type.intData, cast(DazzTrack) new DataDazzTrack!(int32, DataType)(trackPtr, rawSize),
            DazzTrack.Type.longData, cast(DazzTrack) new DataDazzTrack!(int64, DataType)(trackPtr, rawSize),
            null
        );

        _trackIndex[trackName] = dazzTrack;

        return dazzTrack;
    }

    unittest
    {
        import dazzlib.util.tempfile;
        import dazzlib.util.testdata;
        import std.exception;
        import std.file;
        import std.algorithm;

        auto tmpDir = mkdtemp("./.unittest-XXXXXX");
        scope (exit)
            rmdirRecurse(tmpDir);

        auto dbFile = buildPath(tmpDir, "test.db");
        enum maskName = "test-mask";

        writeTestDb(dbFile);
        writeTestMask(dbFile, maskName);

        auto dazzDb = new DazzDb(dbFile);
        auto maskTrack = dazzDb.getTrack(maskName);

        assert(maskTrack.name == maskName);
        assert(!maskTrack.isDataLoaded);
        assert(maskTrack.type == DazzTrack.Type.mask);
        assert(maskTrack.maskTrack !is null);
        assert(maskTrack.dataTrack!(int64, int32) !is null);
        assert(maskTrack.maskTrack.dataSegments[] == testMaskData);
        assert(equal!equal(
            maskTrack.maskTrack.intervals(),
            testMaskData.map!(readData => readData
                .chunks(2)
                .map!(dataPair => MaskTrack.Interval(dataPair[0], dataPair[1]))
            ),
        ));
    }


    protected auto readTrackHead(string trackName)
    {
        int[2] trackHead;
        File(files.trackAnnotationFile(trackName), "r").rawReadAll(trackHead[]);

        return tuple!("trackLength", "recordSize")(trackHead[0], trackHead[1]);
    }


    /// Validate track with trackName.
    void validateTrack(string trackName, TrackKind kind, TrackFor trackFor = TrackFor.trimmed) const
    {
        string corrupted(string explanation, Args...)(Args args)
        {
            return format!("track files for %s are corrupted: " ~ explanation)(trackName, args);
        }

        auto annoFile = File(files.trackAnnotationFile(trackName), "r");

        int trackLength;
        int recordSize;

        annoFile.rawReadScalar(
            trackLength,
            corrupted!"reached EOF while reading trackLength",
        );
        annoFile.rawReadScalar(
            recordSize,
            corrupted!"reached EOF while reading recordSize",
        );

        dazzlibEnforce(recordSize >= 0, corrupted!"negative recordSize");
        if (recordSize == 0)
            dazzlibEnforce(kind == TrackKind.mask, "track kind does not match expectation");
        else
            dazzlibEnforce(kind == TrackKind.custom, "track kind does not match expectation");

        int numReadsTrimmed;
        int numReadsUntrimmed;

        if (block > 0)
        {
            numReadsUntrimmed = *(cast(int*) dazzDb.reads -1);
            numReadsTrimmed = *(cast(int*) dazzDb.reads -2);
        }
        else
        {
            numReadsUntrimmed = dazzDb.ureads;
            numReadsTrimmed = dazzDb.treads;
        }

        dazzlibEnforce(
            trackLength.among(numReadsUntrimmed, numReadsTrimmed),
            corrupted!"trackLength matches neither trimmed nor untrimmed size",
        );
        auto observedTrackFor = numReadsTrimmed == numReadsUntrimmed
            ? TrackFor.any
            : trackLength == numReadsUntrimmed
                ? TrackFor.untrimmed
                : TrackFor.trimmed;

        dazzlibEnforce(
            trackFor & observedTrackFor,
            format!"track is for %s DB but expected %s"(observedTrackFor, trackFor),
        );
    }

    unittest
    {
        import dazzlib.util.tempfile;
        import dazzlib.util.testdata;
        import std.exception;
        import std.file;
        import std.algorithm;

        auto tmpDir = mkdtemp("./.unittest-XXXXXX");
        scope (exit)
            rmdirRecurse(tmpDir);

        auto dbFile = buildPath(tmpDir, "test.db");
        enum maskName = "test-mask";

        writeTestDb(dbFile);
        writeTestMask(dbFile, maskName);

        auto dazzDb = new DazzDb(dbFile);
        assertNotThrown!DazzlibException(dazzDb.validateTrack(
            maskName,
            TrackKind.mask,
            TrackFor.any,
        ));

        assertThrown!DazzlibException(dazzDb.validateTrack(
            maskName,
            TrackKind.custom,
            TrackFor.any,
        ));
    }


    /// Get the FASTA header of `reads[readIdx]`
    string getContigHeader(size_t readIdx)
    {
        dazzlibEnforce(dbType == DbType.dam, "may only be called on a DAM");

        static MmFile headerFile;

        if (headerFile is null)
            headerFile = new MmFile(cast(string) (dbName ~ DbExtension.headers));

        // read MAX_NAME bytes at reads[readIdx].coff
        auto headerOffset = reads[readIdx].coff;
        auto headerBytes = cast(char[]) headerFile[headerOffset .. min(headerOffset + MAX_NAME, $)];

        // chop off the leading '>'
        if (headerBytes.length > 0 && headerBytes[0] == '>')
            headerBytes = headerBytes[1 .. $];

        return headerBytes.until(newline).to!string;
    }


    /// Get number of reads in (un)trimmed dbFile without opening the entire DB.
    static id_t numReads(string what = "trimmed")(string dbFile) if (what.among("trimmed", "untrimmed"))
    {
        auto dazzStub = Read_DB_Stub(dbFile.toStringz, StubPart.blocks);
        dazzlibEnforce(dazzStub !is null, "DB stub cannot be read; is the DB split?");

        scope (exit)
            Free_DB_Stub(dazzStub);

        static if (what == "trimmed")
            auto indexPtr = dazzStub.tblocks;
        else static if (what == "untrimmed")
            auto indexPtr = dazzStub.ublocks;
        else
            static assert(0, "`what` must be \"trimmed\" or \"untrimmed\"");

        auto blockIndex = indexPtr[0 .. dazzStub.nblocks];
        auto numReads = blockIndex[$ - 1] - blockIndex[0];

        return numReads;
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


/// Wrapper around a `DAZZ_TRACK*`.
class DazzTrack
{
    ///
    enum Type : ubyte
    {
        /// Unsupported configuration.
        unsupported,
        /// There are numReads records of recordSize bytes each.
        annotation,
        /// There is variable length data vector of pairs of 32bit integers
        /// indexed by 64bit integers.
        mask,
        /// There is variable length data indexed by 32bit integers.
        intData,
        /// There is variable length data indexed by 64bit integers.
        longData,
    }


    // Provide names for special values of the `DAZZ_TRACK.size` field.
    private enum SizeOf : int
    {
        mask = 0,
        intData = 4,
        longData = 8,
    }


    private DAZZ_TRACK* _dazzTrack;
    private int _rawSize;


    protected this(DAZZ_TRACK* dazzTrack, int rawSize) nothrow @nogc
    {
        this._dazzTrack = dazzTrack;
        this._rawSize = rawSize;
    }


    this(DazzTrack baseDazzTrack) nothrow @nogc
    {
        this(baseDazzTrack.dazzTrack, baseDazzTrack._rawSize);
    }


    protected @property ref inout(DAZZ_TRACK*) dazzTrack() inout pure nothrow @nogc @safe
    {
        assert(_dazzTrack !is null);

        return _dazzTrack;
    }


    /// Symbolic name of track
    @property const(char)[] name() const pure nothrow
    {
        return dazzTrack.name.fromStringz();
    }

    /// Size of track records
    @property size_t recordSize() const pure nothrow @safe @nogc
    {
        return dazzTrack.size.boundedConvert!(typeof(return));
    }

    // Value the should be writte to the track header.
    protected @property int rawSize() const pure nothrow @safe @nogc
    {
        return _rawSize;
    }

    /// Number of reads in track
    @property id_t numReads() const pure nothrow @safe @nogc
    {
        return dazzTrack.nreads.boundedConvert!(typeof(return));
    }


    /// Get type of track
    Type type() const pure nothrow @nogc
    {
        return type(dazzTrack, _rawSize);
    }


    protected static Type type(const DAZZ_TRACK* dazzTrack, int rawSize) pure nothrow @safe @nogc
    {
        if (dazzTrack is null)
            return Type.unsupported;

        if (dazzTrack.data is null)
            return Type.annotation;
        else if (rawSize == SizeOf.mask)
            return Type.mask;
        else if (rawSize == SizeOf.intData)
            return Type.intData;
        else if (rawSize == SizeOf.longData)
            return Type.longData;
        else
            return Type.unsupported;
    }

    /// Is track data loaded in memory?
    @property bool isDataLoaded() const pure nothrow @safe @nogc
    {
        return dazzTrack.loaded != 0;
    }

    /// Largest read data segment in bytes
    @property size_t maxDataSegmentBytes() const pure nothrow @safe @nogc
    {
        return dazzTrack.dmax.boundedConvert!(typeof(return));
    }


    AnnotationDazzTrack!T annotationTrack(T)() nothrow @nogc
    {
        if (type == Type.annotation)
            return cast(typeof(return)) this;
        else
            return null;
    }


    DataDazzTrack!(A, D) dataTrack(A, D)() nothrow @nogc
    {
        if (type.among(Type.mask, Type.intData, Type.longData))
            return cast(typeof(return)) this;
        else
            return null;
    }


    MaskTrack maskTrack() nothrow @nogc
    {
        if (type == Type.mask)
            return cast(typeof(return)) this;
        else
            return null;
    }
}



///
class AnnotationDazzTrack(T) : DazzTrack
{
    protected this(DAZZ_TRACK* dazzTrack, int rawSize) nothrow @nogc
    {
        super(dazzTrack, rawSize);

        assert(
            type == Type.annotation,
            "cannot derive " ~ typeof(this).stringof ~ ": type must be annotation",
        );
        assert(recordSize == T.sizeof, "record size does not match the the size of " ~ T.stringof);
    }


    inout(T)[] annotations() inout nothrow @nogc
    {
        return (cast(inout(T)*) dazzTrack.anno)[0 .. numReads];
    }
}

// trigger compilation
private shared AnnotationDazzTrack!int __annotationDazzTrackInt;


///
class DataDazzTrack(T, Data=byte) if (is(T == int32) || is(T == int64)) : DazzTrack
{
    protected this(DAZZ_TRACK* dazzTrack, int rawSize) nothrow @nogc
    {
        super(dazzTrack, rawSize);

        assert(type.among(Type.mask, Type.intData, Type.longData),
            "cannot derive " ~ typeof(this).stringof ~ ": type must be intData or longData",
        );
        assert(recordSize == T.sizeof, "record size does not match the size of " ~ T.stringof);
    }


    /// Load the all data into a single block of memory.
    void loadAllData()
    {
        auto result = Load_All_Track_Data(dazzTrack);
        dazzlibEnforce(result == 0, currentError.idup);
    }


    protected inout(T)[] annotations() inout nothrow @nogc
    {
        return (cast(inout(T)*) dazzTrack.anno)[0 .. numReads + 1];
    }


    /// Return an array-like object with typed and bounds-checked access to
    /// data.
    ///
    /// Supported operations:
    /// ```
    /// /// Return a two-dimensional array of the data.
    /// Data[][] opIndex()
    ///
    /// /// Return the data for given read index.
    /// Data[] opIndex(size_t readIdx)
    ///
    /// /// Copy the data for given read index into buffer and return the
    /// /// remaining unused part of buffer.
    /// Data[] copy(size_t readIdx, Data[] buffer) @nogc
    /// ```
    auto dataSegments() nothrow @nogc
    {
        return DataSegments!(typeof(this))(this);
    }

    /// ditto
    auto dataSegments() const nothrow @nogc
    {
        return DataSegments!(typeof(this))(this);
    }


    static struct DataSegments(This)
    {
        enum constAccess = is(This == const(This));
        enum mutableAccess = !constAccess;

        This dazzTrack;

        static if (mutableAccess)
            Data[][] opIndex()
            {
                dazzTrack.loadAllData();

                return getDataField();
            }

        static if (constAccess)
            const(Data)[][] opIndex() const
            {
                enforceDataLoaded();

                return getDataField();
            }


        @property size_t opDollar() const pure nothrow @safe @nogc
        {
            return dazzTrack.numReads;
        }

        alias length = opDollar;


        static if (mutableAccess)
            Data[] opIndex(size_t readIdx)
            {
                if (dazzTrack.isDataLoaded)
                    return getDataSliceFor(readIdx);
                else
                {
                    auto dataBuffer = new Data[segmentSize(readIdx)];
                    auto bufferRest = copyDataFor(readIdx, dataBuffer);
                    assert(bufferRest.length == 0);

                    return dataBuffer;
                }
            }


        const(Data)[] opIndex(size_t readIdx) const
        {
            enforceDataLoaded();

            return getDataSliceFor(readIdx);
        }


        static if (mutableAccess)
            Data[] copy(size_t readIdx, Data[] buffer) @nogc
            {
                if (dazzTrack.isDataLoaded)
                    return .copy(getDataSliceFor(readIdx), buffer);
                else
                    return copyDataFor(readIdx, buffer);
            }


        static if (constAccess)
            Data[] copy(size_t readIdx, Data[] buffer) const
            {
                enforceDataLoaded();

                return .copy(getDataSliceFor(readIdx), buffer);
            }


        protected size_t segmentSize(size_t readIdx) const nothrow @nogc
        {
            assert(
                segmentPtr(readIdx) <= segmentPtr(readIdx + 1),
                "segment begins after it ends",
            );

            return boundedConvert!size_t(segmentPtr(readIdx + 1) - segmentPtr(readIdx));
        }


        protected size_t segmentPtr(size_t readIdx) const nothrow @nogc
        {
            assert(readIdx <= length);
            auto beginPtr = dazzTrack.annotations[readIdx];
            assert(
                beginPtr % Data.sizeof == 0,
                "illegal data pointer: not aligned to multiple of " ~ Data.sizeof.stringof,
            );

            return beginPtr / Data.sizeof;
        }


        protected bool validateSegment(size_t readIdx) const nothrow @nogc
        {
            // just invoke assertions in segmentSize
            assert(segmentSize(readIdx) >= 0);

            return true;
        }



        protected inout(Data)[] getDataSliceFor(size_t readIdx) inout nothrow
        {
            assert(validateSegment(readIdx));
            assert(dazzTrack.isDataLoaded);
            auto dataBasePtr = cast(inout(Data)*) dazzTrack.dazzTrack.data;

            return dataBasePtr[segmentPtr(readIdx) .. segmentPtr(readIdx + 1)];
        }


        protected inout(Data)[][] getDataField() inout nothrow
        {
            assert(dazzTrack.isDataLoaded);

            auto dataField = new inout(Data)[][dazzTrack.numReads];

            foreach (readIdx, ref data; dataField)
                data = getDataSliceFor(readIdx);

            return dataField;
        }


        static if (mutableAccess)
            protected Data[] copyDataFor(size_t readIdx, Data[] buffer)
            {
                assert(validateSegment(readIdx));

                auto numBytesRead = Load_Track_Data(dazzTrack.dazzTrack, readIdx.boundedConvert!int, buffer.ptr);
                assert(numBytesRead == segmentSize(readIdx));

                return buffer[numBytesRead .. $];
            }


        protected void enforceDataLoaded() const
        {
            dazzlibEnforce(dazzTrack.isDataLoaded, "data must be loaded before const access");
        }
    }
}


// trigger compilation
private shared DataDazzTrack!int32 __dataDazzTrackInt32;
// trigger compilation
private shared DataDazzTrack!int64 __dataDazzTrackInt64;

///
class MaskTrack : DataDazzTrack!(int64, int32)
{
    static struct Interval
    {
        coord_t begin;
        coord_t end;

        invariant
        {
            assert(begin <= end, "interval end < begin");
        }

        @property coord_t length() const pure nothrow @safe @nogc
        {
            return end - begin;
        }
    }

    protected this(DAZZ_TRACK* dazzTrack, int rawSize) nothrow @nogc
    {
        super(dazzTrack, rawSize);

        assert(type == Type.mask, "cannot derive " ~ typeof(this).stringof ~ ": type must be mask");
    }


    /// Returns a lazy random-access range of ranges where each elements is
    /// the result of calling `intervals(id_t)`.
    ///
    /// See_also: toIntervals
    auto intervals() inout
    {
        return iota(numReads).map!(readIdx => intervals(readIdx));
    }


    /// Returns a lazy random-access range of intervals for readIdx.
    ///
    /// See_also: toIntervals
    auto intervals(id_t readIdx) inout
    {
        return toIntervals(dataSegments[readIdx]);
    }


    /// Returns a lazy random-access range of `Interval`s for dataSegment.
    static auto toIntervals(const int[] dataSegment)
    {
        return dataSegment
            .map!(boundedConvert!(coord_t, int))
            .chunks(2)
            .map!((chunk) {
                assert(chunk.length == 2, "corrupted mask data");
                auto interval = Interval(chunk[0], chunk[1]);
                assert(&interval);

                return interval;
            });
    }
}


enum TrackFor : int
{
    untrimmed = 0b01,
    trimmed = 0b10,
    any = untrimmed | trimmed,
}


/// Validate DB by opening the DB once.
///
/// Returns: `null` if DB is valid; otherwise error message.
//string validateDbTrack(string dbFile, TrackFor trackFor = TrackFor.any)
//{
//    new DazzDb(dbFile).validateTrack()
//}

/// Convert sequence to numeric/ascii representation.
alias toNumeric = toAlphabet!"\0\1\2\3\4";

/// ditto
alias toLowerAscii = toAlphabet!"acgt\0";

/// ditto
alias toUpperAscii = toAlphabet!"ACGT\0";


/// Convert sequence to given alphabet. Note, works only for numeric- or
/// ascii-encoded input sequences.
void toAlphabet(char[5] alphabet)(char[] sequence) pure nothrow @nogc
in (
    (*(sequence.ptr - 1)).among('\0', '\4', alphabet[4]) &&
    *(sequence.ptr - 1) == *(sequence.ptr + sequence.length),
    "sequence must be terminated"
)
out (;
    *(sequence.ptr - 1) == alphabet[4] &&
    *(sequence.ptr - 1) == *(sequence.ptr + sequence.length),
    "transformed sequence must be terminated"
)
{
    // set terminators
    *(sequence.ptr - 1) = alphabet[4];
    *(sequence.ptr + sequence.length) = alphabet[4];

    foreach (ref c; sequence)
        switch (c)
        {
            case 0:
            case 'a':
            case 'A':
                c = alphabet[0];
                break;
            case 1:
            case 'c':
            case 'C':
                c = alphabet[1];
                break;
            case 2:
            case 'g':
            case 'G':
                c = alphabet[2];
                break;
            case 3:
            case 't':
            case 'T':
                c = alphabet[3];
                break;
            default:
                c = alphabet[0];
                assert(0, "illegal char");
        }
}

unittest
{
    enum asciiLowerSeq = "\0tctcacctcgatccccctagcactctaatattacagcgttca\0";
    enum asciiUpperSeq = "\0TCTCACCTCGATCCCCCTAGCACTCTAATATTACAGCGTTCA\0";
    enum numericSeq = "43131011312031111130210131300303301021233104".map!"cast(char) (a - '0')".array;
    enum numericAsciiSeq = "Z313101131203111113021013130030330102123310Z";

    char[] seq;
    static foreach (iseq; [asciiLowerSeq, asciiUpperSeq, numericSeq])
    {
        seq = iseq.dup[1 .. $ - 1];
        seq.toLowerAscii();
        assert(seq == asciiLowerSeq[1 .. $ - 1]);

        seq = iseq.dup[1 .. $ - 1];
        seq.toUpperAscii();
        assert(seq == asciiUpperSeq[1 .. $ - 1]);

        seq = iseq.dup[1 .. $ - 1];
        seq.toNumeric();
        assert(seq == numericSeq[1 .. $ - 1]);

        seq = iseq.dup[1 .. $ - 1];
        seq.toAlphabet!"0123Z"();
        assert(seq == numericAsciiSeq[1 .. $ - 1]);
    }
}


/// Convert sequence from given alphabet to numeric.
void fromAlphabet(char[5] alphabet)(char[] sequence) pure nothrow @nogc
in (
    *(sequence.ptr - 1) == alphabet[4] &&
    *(sequence.ptr - 1) == *(sequence.ptr + sequence.length),
    "sequence must be terminated"
)
out (;
    *(sequence.ptr - 1) == 4 &&
    *(sequence.ptr - 1) == *(sequence.ptr + sequence.length),
    "transformed sequence must be terminated"
)
{
    // set terminators
    *(sequence.ptr - 1) = 4;
    *(sequence.ptr + sequence.length) = 4;

    foreach (ref c; sequence)
        switch (c)
        {
            case alphabet[0]:
                c = 0;
                break;
            case alphabet[1]:
                c = 1;
                break;
            case alphabet[2]:
                c = 2;
                break;
            case alphabet[3]:
                c = 3;
                break;
            default:
                c = 0;
                assert(0, "illegal char");
        }
}

unittest
{
    enum numericAsciiSeq = "Z313101131203111113021013130030330102123310Z";
    enum numericSeq = "43131011312031111130210131300303301021233104".map!"cast(char) (a - '0')".array;

    auto seq = numericAsciiSeq.dup[1 .. $ - 1];
    seq.fromAlphabet!"0123Z"();
    assert(seq == numericSeq[1 .. $ - 1]);
}
