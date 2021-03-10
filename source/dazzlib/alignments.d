/**
    High-level access to alignment data.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.alignments;

import dazzlib.basictypes;
import dazzlib.core.c.DB;
import dazzlib.core.c.align_;
import dazzlib.util.exception;
import dazzlib.util.safeio;
import dazzlib.util.math;
import std.algorithm;
import std.array;
import std.conv;
import std.format;
import std.range;
import std.stdio;
import std.string;
import std.traits;
import std.typecons;


///
public import dazzlib.core.c.align_ : LocalAlignmentFlag;


/// Typesafe structure for storing `LocalAlignmentFlag`s.
alias LocalAlignmentFlags = BitFlags!LocalAlignmentFlag;


/// Contig characteristics
struct Contig
{
    /// Contig id (1-based indexing)
    id_t id;
    /// Contig length; zero means undefined.
    coord_t length;
}


/// A locus is defined as the contig ID and the right-open, zero-based
/// coordinate interval and may contain the contig length.
struct Locus
{
    ///
    Contig contig;
    ///
    coord_t begin;
    ///
    coord_t end;


    @property coord_t length() const pure nothrow @safe
    {
        return end - begin;
    }


    @property void boundedBegin(coord_t begin) pure nothrow @safe
    {
        this.begin = min(begin, contig.length);
    }


    @property void boundedEnd(coord_t end) pure nothrow @safe
    {
        this.end = min(end, contig.length);
    }


    bool beginsWithin(coord_t allowance) const pure nothrow @safe
    {
        return begin <= allowance;
    }


    bool endsWithin(coord_t allowance) const pure nothrow @safe
    {
        return end + allowance >= contig.length;
    }


    bool isFullyContained(coord_t allowance) const pure nothrow @safe
    {
        return beginsWithin(allowance) && endsWithin(allowance);
    }
}


/// A single trace point (expanded to 16bit encoding).
struct TracePoint
{
    ///
    trace_point_t numDiffs;
    ///
    trace_point_t numBasePairs;
}


/// Stores coordinates on contigA and contigB that correspond to each other
/// according to a given trace path.
struct TranslatedTracePoint
{
    ///
    coord_t contigA;
    ///
    coord_t contigB;
}


///
struct Trace
{
    ///
    Locus contigA;
    ///
    Locus contigB;
    ///
    trace_point_t tracePointSpacing;
    ///
    const(TracePoint)[] tracePoints;


    ///
    TranslatedTracePoint translateTracePoint(string contig)(
        coord_t contigPos,
        RoundingMode roundingMode,
    ) const pure if (contig.among("contigA", "contigB"))
    {
        assert(mixin(contig ~ `.begin <= contigPos && contigPos <= ` ~ contig ~ `.end`));

        auto tracePointIndex = tracePointsUpTo!contig(contigPos, roundingMode);
        auto contigBPos = contigB.begin + tracePoints[0 .. tracePointIndex]
                .map!"a.numBasePairs"
                .sum;
        auto contigAPos = tracePointIndex == 0
            ? contigA.begin
            : tracePointIndex < tracePoints.length
                ? floor(contigA.begin, tracePointSpacing) + cast(coord_t) (tracePointIndex * tracePointSpacing)
                : contigA.end;

        return TranslatedTracePoint(contigAPos, contigBPos);
    }


    ///
    auto tracePointsUpTo(string contig)(
        coord_t contigAPos,
        RoundingMode roundingMode,
    ) const pure nothrow if (contig == "contigA")
    {
        assert(contigA.begin <= contigAPos && contigAPos <= contigA.end);

        auto firstTracePointRefPos = contigA.begin;
        auto secondTracePointRefPos = floor(firstTracePointRefPos, tracePointSpacing) + tracePointSpacing;
        auto secondFromLastTracePointRefPos = floor(contigA.end - 1, tracePointSpacing);

        final switch (roundingMode)
        {
        case RoundingMode.floor:
            if (contigAPos < secondTracePointRefPos)
                return 0;
            if (contigAPos < contigA.end)
                return 1 + (contigAPos - secondTracePointRefPos) / tracePointSpacing;
            else
                return tracePoints.length;
        case RoundingMode.round:
            assert(0, "unimplemented");
        case RoundingMode.ceil:
            if (firstTracePointRefPos == contigAPos)
                return 0;
            if (contigAPos <= secondTracePointRefPos)
                return 1;
            else if (contigAPos <= secondFromLastTracePointRefPos)
                return 1 + ceildiv(contigAPos - secondTracePointRefPos, tracePointSpacing);
            else
                return tracePoints.length;
        }
    }

    unittest
    {
        auto trace = Trace(
            Locus(Contig(), 50, 2897),
            Locus(Contig(), 50, 2905),
            100,
            [TracePoint(1, 50), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(7, 105)],
        );
        coord_t contigPos = 79;
        auto roundingMode = RoundingMode.ceil;
        auto tpsUpTo = trace.tracePointsUpTo!"contigA"(contigPos, roundingMode);

        assert(0 <= tpsUpTo && tpsUpTo <= trace.tracePoints.length);
    }

    /// ditto
    auto tracePointsUpTo(string contig)(
        coord_t contigBPos,
        RoundingMode roundingMode,
    ) const pure nothrow if (contig == "contigB")
    {
        assert(contigB.begin <= contigBPos && contigBPos <= contigB.end);

        if (contigBPos == contigB.begin)
            return 0;
        else if (contigBPos == contigB.end)
            return tracePoints.length;

        auto tracePointPositions = chain(only(contigB.begin), tracePoints.map!"a.numBasePairs")
            .cumulativeFold!"a + b"
            .enumerate;

        final switch (roundingMode)
        {
        case RoundingMode.floor:
            return tracePointPositions
                .find!(pair => contigBPos < pair.value)
                .front
                .index - 1;
        case RoundingMode.round:
            assert(0, "unimplemented");
        case RoundingMode.ceil:
            return tracePointPositions
                .find!(pair => contigBPos <= pair.value)
                .front
                .index;
        }
    }
}


/// D-friendly version of `dazzlib.core.c.align_.Alignment`.
struct LocalAlignment
{
    ///
    Locus contigA;
    ///
    Locus contigB;
    ///
    LocalAlignmentFlags flags;
    ///
    trace_point_t tracePointSpacing;
    ///
    TracePoint[] tracePoints;


    /// If not already disabled set the disabled flag to `disable`.
    bool disableIf(lazy bool disable) pure
    {
        if (!flags.disabled)
            flags.disabled = disable;

        return flags.disabled;
    }


    ///
    @property Trace trace() const pure nothrow @safe
    {
        return Trace(
            contigA,
            contigB,
            tracePointSpacing,
            tracePoints,
        );
    }


    /**
        Generate a cartoon of this alignment relative to `contig`. Every
        character in the cartoon corresponds to one trace point distance.

        Returns: a cartoon of this alignment
    */
    static string cartoon(string contig)(in LocalAlignment[] localAlignments...)
    {
        if (localAlignments.length == 0)
            return "";

        alias _locus = la => mixin("la." ~ contig);
        alias _contig = la => _locus(la).contig;

        auto referenceContig = _contig(localAlignments[0]);
        auto bpsPerChar = localAlignments[0].tracePointSpacing;

        dazzlibEnforce(
            localAlignments.all!(la => _contig(la) == referenceContig),
            "all alignment chains must share the same reference contig",
        );

        auto cartoonLine(in LocalAlignment la)
        {
            auto skipBps = contig == "contigA" || !la.flags.complement
                ? _locus(la).begin
                : referenceContig.length - _locus(la).end;

            return chain(
                ' '.repeat(skipBps / bpsPerChar),
                '-'.repeat(ceildiv(_locus(la).length, bpsPerChar)),
            );
        }

        return cast(string) chain(
            '-'.repeat(ceildiv(referenceContig.length, bpsPerChar)),
            only('\n'),
            localAlignments.map!cartoonLine.joiner(only('\n')),
        ).array;
    }

    ///
    unittest
    {
        auto las = [
            LocalAlignment(
                Locus(Contig(1, 10), 0, 3),
                Locus(Contig(1, 10), 0, 3),
                LocalAlignmentFlags(),
                1,
            ),
            LocalAlignment(
                Locus(Contig(1, 10), 4, 5),
                Locus(Contig(1, 10), 4, 5),
                LocalAlignmentFlags(),
                1,
            ),
            LocalAlignment(
                Locus(Contig(1, 10), 5, 8),
                Locus(Contig(1, 10), 0, 3),
                LocalAlignmentFlags(LocalAlignmentFlag.complement),
                1,
            ),
            LocalAlignment(
                Locus(Contig(1, 10), 9, 10),
                Locus(Contig(1, 10), 4, 5),
                LocalAlignmentFlags(LocalAlignmentFlag.complement),
                1,
            ),
        ];

        assert(cartoon!"contigA"(las) == "----------\n" ~
                                         "---\n" ~
                                         "    -\n" ~
                                         "     ---\n" ~
                                         "         -");
        assert(cartoon!"contigB"(las) == "----------\n" ~
                                         "---\n" ~
                                         "    -\n" ~
                                         "       ---\n" ~
                                         "     -");
    }
}


/// Returns true if 16bits are required for encoding the trace at
/// tracePointSpacing.
bool isLargeTraceType(Int)(const Int tracePointSpacing) pure nothrow @safe if (isIntegral!Int)
{
    return tracePointSpacing > TRACE_XOVR;
}


/// Returns number of bytes used to encode a trace with tracePointSpacing.
uint tracePointBytes(Int)(const Int tracePointSpacing) pure nothrow @safe if (isIntegral!Int)
{
    return isLargeTraceType(tracePointSpacing)?  trace_point_t.sizeof : small_trace_point_t.sizeof;
}


/// Controls how to manage the buffer when reading data into memory.
enum BufferMode : ubyte
{
    /// Keep a single buffer and keep overwriting it with every new record.
    overwrite,
    /// Allocate a new buffer for every record. Use an `Appender` if the
    /// number of records is unknown.
    dynamic,
    /// Write all records to a continuous stretch of memory preallocated by
    /// the caller. The memory may be uninitialized since it will always be
    /// written to before any read occurs.
    preallocated,
    /// Do not use memory at all instead just skip over the records.
    skip,
}


/// Input range that reads `LocalAlignment`s from a LAS file.
class LocalAlignmentReader
{
    protected
    {
        File las;
        BufferMode bufferMode;
        coord_t[] aLengths;
        coord_t[] bLengths;
        bool isLargeTraceType;
        id_t _numLocalAlignments;
        trace_point_t _tracePointSpacing;
        id_t numLocalAlignmentsLeft;
        LocalAlignment currentLA;
        Overlap overlapHead;
        TracePoint[] tracePointBuffer;
        TracePoint[] fullTracePointBuffer;
    }


    this(const string lasFile, BufferMode bufferMode, TracePoint[] tracePointBuffer = [])
    {
        this(
            lasFile,
            cast(string) null,
            cast(string) null,
            bufferMode,
            tracePointBuffer
        );
    }


    this(
        const string lasFile,
        string dbA,
        string dbB,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer = [],
    )
    {
        this(
            lasFile,
            contigLengths(dbA),
            dbB !is null && dbA != dbB
                ? contigLengths(dbB)
                : contigLengths(dbA),
            bufferMode,
            tracePointBuffer
        );
    }


    private this(
        const string lasFile,
        coord_t[] aLengths,
        coord_t[] bLengths,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer,
    )
    {
        this.las = File(lasFile, "rb");
        this.bufferMode = bufferMode;
        this.tracePointBuffer = tracePointBuffer;
        if (bufferMode.among(BufferMode.dynamic, BufferMode.skip))
            assert(
                tracePointBuffer.length == 0,
                "preallocated buffer was supplied but buffer mode does not one: "~
                "please remove the buffer or adjust the mode according to your needs."
            );
        else
            assert(
                tracePointBuffer.length > 0,
                "buffer mode requires a preallocated buffer but none was supplied: "~
                "please supply a preallocated buffer or adjust the mode according to your needs."
            );
        this.fullTracePointBuffer = tracePointBuffer;
        this.aLengths = aLengths;
        this.bLengths = bLengths;

        initialize();
    }


    private coord_t[] contigLengths(string db)
    {
        if (db is null)
            return [];

        DAZZ_DB dazzDb;
        dazzlibEnforce(Open_DB(db.toStringz, &dazzDb) >= 0, currentError);
        scope (exit)
            Close_DB(&dazzDb);
        catchErrorMessage!Trim_DB(&dazzDb);

        return dazzDb
            .reads[0 .. dazzDb.nreads]
            .map!(read => read.rlen.to!coord_t)
            .array;
    }


    private void initialize()
    {
        readHeader();
        this.numLocalAlignmentsLeft = this.numLocalAlignments;

        if (!empty)
        {
            ++numLocalAlignmentsLeft;
            popFront();
        }
    }


    /// Reset the reader to the first local alignment.
    void reset()
    {
        las.rewind();
        tracePointBuffer = fullTracePointBuffer;
        initialize();
    }


    ///
    @property id_t numLocalAlignments() const pure nothrow @safe
    {
        return _numLocalAlignments;
    }


    ///
    @property trace_point_t tracePointSpacing() const pure nothrow @safe
    {
        return _tracePointSpacing;
    }


    ///
    @property bool empty() const pure nothrow @safe
    {
        return numLocalAlignmentsLeft == 0;
    }


    /// Number of local alignments left to read.
    @property size_t length() const pure nothrow @safe
    {
        return numLocalAlignmentsLeft;
    }


    ///
    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty " ~ typeof(this).stringof);

        --numLocalAlignmentsLeft;
        if (!empty)
            readLocalAlignment();
    }

    ///
    @property LocalAlignment front() pure nothrow @safe
    {
        assert(!empty, "Attempting to fetch the front of an empty " ~ typeof(this).stringof);

        return currentLA;
    }


    /// True if trace points are skipped rather than read into memory.
    @property bool skipTracePoints() const pure nothrow @safe
    {
        return bufferMode == BufferMode.skip;
    }


    /// Number of trace points in `front`. This is useful if `skipTracePoints`
    /// is true.
    @property size_t currentNumTracePoints() const pure @safe
    {
        return (overlapHead.path.tlen / 2).to!size_t;
    }

protected:

    void readHeader()
    {
        readNumLocalAlignments();
        readTracePointSpacing();
    }


    void readNumLocalAlignments()
    {
        long rawNumLocalAlignments;

        las.rawReadScalar(rawNumLocalAlignments, unexpectedEOF!"numLocalAlignments");

        this._numLocalAlignments = rawNumLocalAlignments.to!id_t;
    }


    void readTracePointSpacing()
    {
        int rawTracePointSpacing;

        las.rawReadScalar(rawTracePointSpacing, unexpectedEOF!"tracePointSpacing");

        this._tracePointSpacing = rawTracePointSpacing.to!trace_point_t;
        this.isLargeTraceType = .isLargeTraceType(tracePointSpacing);
    }


    void readLocalAlignment()
    {
        readOverlapHead();
        fillInOverlapHead();
        if (aLengths.length > 0)
            fillInContigLengths();
        auto traceLength = getTraceVectorLength();

        if (skipTracePoints)
            skipTraceVector(traceLength);
        else
            readTraceVector(traceLength);

        if (!skipTracePoints)
        {
            if (!isLargeTraceType)
                rewriteTracePointBufferToLargeTracePointType();

            auto numTracePoints = overlapHead.path.tlen / 2;
            currentLA.tracePoints = tracePointBuffer[0 .. numTracePoints];
            if (bufferMode == BufferMode.preallocated)
                tracePointBuffer = tracePointBuffer[numTracePoints .. $];
        }
    }


    void readOverlapHead()
    {
        auto overlapBytes = (cast(void*) &overlapHead)[
            typeof(Overlap.path.trace).sizeof ..
            Overlap.sizeof
        ];
        las.rawReadAll(overlapBytes, unexpectedEOF!"overlapHead");
    }


    void fillInOverlapHead() pure nothrow @safe
    {
        currentLA.contigA.contig.id = overlapHead.aread + 1;
        currentLA.contigA.begin = overlapHead.path.abpos;
        currentLA.contigA.end = overlapHead.path.aepos;
        currentLA.contigB.contig.id = overlapHead.bread + 1;
        currentLA.contigB.begin = overlapHead.path.bbpos;
        currentLA.contigB.end = overlapHead.path.bepos;
        currentLA.flags = LocalAlignmentFlags(cast(LocalAlignmentFlag) overlapHead.flags);
        currentLA.flags.unchained = !(currentLA.flags.start || currentLA.flags.next);
        currentLA.tracePointSpacing = tracePointSpacing;
    }


    void fillInContigLengths() pure @safe
    {
        if (aLengths.length > 0)
        {
            dazzlibEnforce(
                currentLA.contigA.contig.id - 1 < aLengths.length,
                format!"contigA.contig.id out of bounds: %d >= %d"(currentLA.contigA.contig.id - 1, aLengths.length),
            );
            currentLA.contigA.contig.length = aLengths[currentLA.contigA.contig.id - 1];
        }

        if (bLengths.length)
        {
            dazzlibEnforce(
                currentLA.contigB.contig.id - 1 < bLengths.length,
                format!"contigB.contig.id out of bounds: %d >= %d"(currentLA.contigB.contig.id - 1, bLengths.length),
            );
            currentLA.contigB.contig.length = bLengths[currentLA.contigB.contig.id - 1];
        }
    }


    size_t getTraceVectorLength() const pure @safe
    {
        auto traceVectorLength = overlapHead.path.tlen * tracePointBytes(tracePointSpacing);
        dazzlibEnforce(
            traceVectorLength % 2 == 0,
            "illegal value for tlen: must be multiple of 2",
        );

        return traceVectorLength;
    }


    void readTraceVector(size_t traceLength)
    {
        if (bufferMode == BufferMode.dynamic)
            tracePointBuffer = uninitializedArray!(TracePoint[])(overlapHead.path.tlen / 2);
        auto rawBuffer = getRawTracePointBuffer(traceLength);
        las.rawReadAll(rawBuffer, unexpectedEOF!"tracePoints");
    }


    void skipTraceVector(size_t traceLength)
    {
        import core.stdc.stdio : SEEK_CUR;

        las.seek(traceLength, SEEK_CUR);
    }


    void rewriteTracePointBufferToLargeTracePointType()
    {
        auto traceVectorLength = overlapHead.path.tlen;
        // make sure enough memory is allocated
        auto smallTraceVector = getRawTracePointBuffer(traceVectorLength);
        // select appropriate portion of tracePointBuffer
        auto largeTracePoints = tracePointBuffer.ptr[0 .. traceVectorLength / 2];
        // convert numbers in reverse order as to avoid overwriting
        foreach_reverse (i, ref largeTracePoint; largeTracePoints)
        {
            largeTracePoint = TracePoint(
                smallTraceVector[2 * i],
                smallTraceVector[2 * i + 1],
            );
        }
    }


    @property ubyte[] getRawTracePointBuffer(size_t traceLength) pure nothrow
    {
        // prepare bytes buffer
        auto bufferBytes = tracePointBuffer.length * TracePoint.sizeof;
        auto rawBuffer = (cast(ubyte*) tracePointBuffer.ptr)[0 .. bufferBytes];
        // reduce size (triggers bounds check)
        rawBuffer = rawBuffer[0 .. traceLength];

        return rawBuffer;
    }


    @property string unexpectedEOF(string what)() const pure @safe
    {
        enum msgFormat = "error reading LAS file `%s`: premature end of file, expected " ~ what;

        return format!msgFormat(las.name);
    }
}

static assert(isInputRange!LocalAlignmentReader);
static assert(hasLength!LocalAlignmentReader);
static assert(is(ElementType!LocalAlignmentReader == LocalAlignment));


unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;
    import std.algorithm;
    import std.path;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto lasFile = buildPath(tmpDir, "test.las");
    auto computedLocalAlignments = testLocalAlignments
        .map!((la) {
            la.contigA.contig.length = 0;
            la.contigB.contig.length = 0;

            return la;
        });

    writeTestLas(lasFile);

    auto recoveredLocalAlignments = new LocalAlignmentReader(
        lasFile,
        BufferMode.dynamic,
    );

    assert(equal(computedLocalAlignments, recoveredLocalAlignments));
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;
    import std.algorithm;
    import std.path;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto lasFile = buildPath(tmpDir, "test.las");
    auto dbFile = buildPath(tmpDir, "test.db");

    writeTestDb(dbFile);
    writeTestLas(lasFile);

    auto recoveredLocalAlignments = new LocalAlignmentReader(
        lasFile,
        dbFile,
        dbFile,
        BufferMode.dynamic,
    );

    assert(equal(testLocalAlignments, recoveredLocalAlignments));
}


/// Construct a range that lazily reads local alignments from lasFile. If
/// dbA and/or dbB are given, the contig lengths will be filled in. If no
/// tracePointBuffer is given, it will be constructed internally according
/// to bufferMode.
///
/// See_also: LocalAlignmentReader
LocalAlignmentReader localAlignmentReader(
    in string lasFile,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return localAlignmentReader(lasFile, null, null, bufferMode, tracePointBuffer);
}

/// ditto
LocalAlignmentReader localAlignmentReader(
    in string lasFile,
    in string dbA,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return localAlignmentReader(lasFile, dbA, null, bufferMode, tracePointBuffer);
}

/// ditto
LocalAlignmentReader localAlignmentReader(
    in string lasFile,
    in string dbA,
    in string dbB,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    if (
        bufferMode.among(BufferMode.overwrite, BufferMode.preallocated) &&
        tracePointBuffer.length == 0
    )
    {
        auto alignmentStats = AlignmentStats.inferFrom(lasFile);

        if (bufferMode == BufferMode.overwrite)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentStats.maxTracePoints,
            );
        else if (bufferMode == BufferMode.preallocated)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentStats.numTracePoints,
            );
    }

    return new LocalAlignmentReader(
        lasFile,
        dbA,
        dbB,
        bufferMode,
        tracePointBuffer,
    );
}


/// Returns true if lasFile is empty
bool isLasEmpty(in string lasFile)
{
    auto reader = localAlignmentReader(lasFile);

    return reader.numLocalAlignments == 0;
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;

    auto lasFile = mkstemp("./.unittest-XXXXXX", ".las").name;
    scope (exit)
        remove(lasFile);

    writeTestLas(lasFile);

    assert(!isLasEmpty(lasFile));
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;

    auto lasFile = mkstemp("./.unittest-XXXXXX", ".las").name;
    scope (exit)
        remove(lasFile);

    writeEmptyTestLas(lasFile);

    assert(isLasEmpty(lasFile));
}


/// Validate LAS by reading the header. By default, empty LAS files are
/// rejected because they are usually useless. This can be turned off using
/// `allowEmpty`.
///
/// Returns: `null` if DB is valid; otherwise error message.
string validateLas(in string lasFile, Flag!"allowEmpty" allowEmpty = No.allowEmpty)
{
    try
    {
        auto reader = localAlignmentReader(lasFile);

        if (!reader.empty || allowEmpty)
            return null;
        else
            return "empty LAS file";
    }
    catch (Exception e)
    {
        return e.msg;
    }
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;

    auto lasFile = mkstemp("./.unittest-XXXXXX", ".las").name;
    scope (exit)
        remove(lasFile);

    writeTestLas(lasFile);

    assert(validateLas(lasFile) is null);
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;

    auto lasFile = mkstemp("./.unittest-XXXXXX", ".las").name;
    scope (exit)
        remove(lasFile);

    writeEmptyTestLas(lasFile);

    assert(validateLas(lasFile) !is null);
    assert(validateLas(lasFile, Yes.allowEmpty) is null);
}


void writeAlignments(R)(const string lasFile, R localAlignments)
    if (isInputRange!R && is(const(ElementType!R) == const(LocalAlignment)))
{
    auto las = File(lasFile, "wb");

    long numLocalAlignments = 0; // will be overwritten at the end
    auto tracePointSpacing = cast(int) AlignmentStats.inferTracePointSpacingFrom(localAlignments);
    las.rawWrite([numLocalAlignments]);
    las.rawWrite([tracePointSpacing]);

    foreach (localAlignment; localAlignments)
    {
        las.writeLocalAlignment(localAlignment, tracePointSpacing);
        ++numLocalAlignments;
    }

    las.rewind();
    las.rawWrite([numLocalAlignments]);

    las.close();
}


unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.algorithm;
    import std.file;
    import std.path;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto dbFile = buildPath(tmpDir, "test.db");
    auto lasFile = buildPath(tmpDir, "test.las");

    writeTestDb(dbFile);
    lasFile.writeAlignments(testLocalAlignments);
    auto recoveredLocalAlignments = new LocalAlignmentReader(
        lasFile,
        dbFile,
        dbFile,
        BufferMode.dynamic,
    );

    assert(equal(testLocalAlignments, recoveredLocalAlignments));
}


private auto writeLocalAlignment(
    File las,
    const LocalAlignment localAlignment,
    const int tracePointSpacing,
)
{
    Overlap overlap;

    // set read IDs
    overlap.aread = localAlignment.contigA.contig.id - 1;
    overlap.bread = localAlignment.contigB.contig.id - 1;
    // set flags
    overlap.flags = cast(LocalAlignmentFlag) cast(int) localAlignment.flags;
    // remove non-standard unchained flag. This is derived from the other
    // flags when reading the LAS.
    overlap.flags &= ~(LocalAlignmentFlag.unchained);

    // set coordinates
    overlap.path.abpos = localAlignment.contigA.begin;
    overlap.path.aepos = localAlignment.contigA.end;
    overlap.path.bbpos = localAlignment.contigB.begin;
    overlap.path.bepos = localAlignment.contigB.end;

    writeOverlap(las, overlap, localAlignment.tracePoints, tracePointSpacing);
}


private void writeOverlap(
    File las,
    ref Overlap overlap,
    const TracePoint[] tracePoints,
    const int tracePointSpacing,
)
{
    // set trace vector length
    overlap.path.tlen = to!int(2*tracePoints.length);
    // set diffs
    overlap.path.diffs = tracePoints.map!"a.numDiffs".sum;

    // write overlap "header"
    auto overlapBytes = (cast(void*) &overlap)[
        typeof(overlap.path.trace).sizeof ..
        Overlap.sizeof
    ];
    las.rawWrite(overlapBytes);

    // write trace vector
    if (tracePoints.length > 0)
    {
        if (isLargeTraceType(tracePointSpacing))
            las.rawWrite(tracePoints);
        else
            // rewrite trace to use `ubyte`s
            las.rawWrite(
                tracePoints
                    .map!(tp => [
                        tp.numDiffs.to!small_trace_point_t,
                        tp.numBasePairs.to!small_trace_point_t,
                    ])
                    .joiner
                    .takeExactly(2*tracePoints.length)
                    .array
            );
    }
}


struct AlignmentStats
{
    /// Total number of local alignments (disregarding chaining)
    size_t numLocalAlignments;

    /// Total number of unchained alignments
    size_t numUnchainedAlignments;

    /// Total number of alignment chains
    size_t numAlignmentChains;

    /// Maximum number of local alignments per chain
    size_t maxLocalAlignments;

    /// Maximum total number of local alignments per contig
    size_t maxLocalAlignmentsPerContig;

    /// Total number of trace points
    size_t numTracePoints;

    /// Maximum number of trace points per local alignment
    size_t maxTracePoints;

    /// Maximum number of trace points per contig
    size_t maxTracePointsPerContig;

    /// Trace point spacing
    trace_point_t tracePointSpacing;


    /// Infer stats from given range by traversing it once.
    static AlignmentStats inferFrom(R)(R localAlignments)
        if (isInputRange!R && is(const(ElementType!R) == const(LocalAlignment)))
    {
        AlignmentStats headerData;

        headerData.tracePointSpacing = inferTracePointSpacingFrom(localAlignments);

        id_t lastContig;
        id_t numLocalAlignmentsSinceLastContig;
        id_t numTracePointsSinceLastContig;
        id_t numLocalAlignmentsInChain = 1;
        foreach (localAlignment; localAlignments)
        {
            if (lastContig != localAlignment.contigA.contig.id)
            {
                // udpate maxLocalAlignmentsPerContig
                headerData.maxLocalAlignmentsPerContig = max(
                    headerData.maxLocalAlignmentsPerContig,
                    numLocalAlignmentsSinceLastContig,
                );
                numLocalAlignmentsSinceLastContig = 0;

                // udpate maxTracePointsPerContig
                headerData.maxTracePointsPerContig = max(
                    headerData.maxTracePointsPerContig,
                    numTracePointsSinceLastContig,
                );
                numTracePointsSinceLastContig = 0;
            }

            ++headerData.numLocalAlignments;
            if (localAlignment.flags.start)
                ++headerData.numAlignmentChains;
            else if (localAlignment.flags.unchained)
                ++headerData.numUnchainedAlignments;
            headerData.maxLocalAlignments = max(
                headerData.maxLocalAlignments,
                numLocalAlignmentsInChain,
            );

            if (localAlignment.flags.start || localAlignment.flags.unchained)
                numLocalAlignmentsInChain = 0;
            ++numLocalAlignmentsInChain;

            size_t currentNumTracePoints = localAlignment.tracePoints.length;

            static if (is(typeof(localAlignments.currentNumTracePoints)))
                if (currentNumTracePoints == 0)
                    currentNumTracePoints = localAlignments.currentNumTracePoints;

            headerData.numTracePoints += currentNumTracePoints;

            headerData.maxTracePoints = max(
                headerData.maxTracePoints,
                currentNumTracePoints,
            );

            ++numLocalAlignmentsSinceLastContig;
            numTracePointsSinceLastContig += currentNumTracePoints;
            lastContig = localAlignment.contigA.contig.id;
        }

        return headerData;
    }


    /// Infer stats from lasFile range by traversing it once.
    static AlignmentStats inferFrom(string lasFile)
    {
        auto lasScanner = new LocalAlignmentReader(lasFile, BufferMode.skip);

        return inferFrom(lasScanner);
    }


    /// Infer trace point spacing from the range by peeking at the first
    /// element if possible; otherwise returns 100, the default spacing.
    static trace_point_t inferTracePointSpacingFrom(R)(R alignments)
        if (isInputRange!R && is(const(ElementType!R) == const(LocalAlignment)))
    {
        if (alignments.empty)
            return 100;

        return alignments.front.tracePointSpacing;
    }
}

unittest
{
    import dazzlib.util.testdata;

    assert(AlignmentStats.inferFrom(testLocalAlignments) == AlignmentStats(
        4,      // numLocalAlignments
        0,      // numUnchainedAlignments
        2,      // numAlignmentChains
        2,      // maxLocalAlignments
        2,      // maxLocalAlignmentsPerContig
        4,      // numTracePoints
        1,      // maxTracePoints
        2,      // maxTracePointsPerContig
        100,    // tracePointSpacing
    ));
}

unittest
{
    import dazzlib.util.tempfile;
    import dazzlib.util.testdata;
    import std.file;

    auto lasFile = mkstemp("./.unittest-XXXXXX", ".las");
    scope (exit)
        remove(lasFile.name);
    lasFile.file.close();

    writeTestLas(lasFile.name);

    assert(AlignmentStats.inferFrom(lasFile.name) == AlignmentStats(
        4,      // numLocalAlignments
        0,      // numUnchainedAlignments
        2,      // numAlignmentChains
        2,      // maxLocalAlignments
        2,      // maxLocalAlignmentsPerContig
        4,      // numTracePoints
        1,      // maxTracePoints
        2,      // maxTracePointsPerContig
        100,    // tracePointSpacing
    ));
}
