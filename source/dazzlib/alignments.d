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
import dazzlib.db;
import dazzlib.util.exception;
import dazzlib.util.safeio;
import dazzlib.util.math;
import std.algorithm;
import std.array;
import std.ascii : toLower;
import std.conv;
import std.format;
import std.range;
import std.stdio;
import std.string;
import std.traits;
import std.typecons;


///
public import dazzlib.core.c.align_ : LocalAlignmentFlag;
public import dazzlib.util.math : RoundingMode;


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


    invariant (isValidInterval());
    invariant (isIntervalInBounds());


    @property coord_t length() const pure nothrow @safe @nogc
    {
        return end - begin;
    }


    @property void boundedBegin(coord_t begin) pure nothrow @safe @nogc
    {
        this.begin = min(begin, contig.length);
    }


    @property void boundedEnd(coord_t end) pure nothrow @safe @nogc
    {
        this.end = min(end, contig.length);
    }


    bool beginsWithin(coord_t allowance) const pure nothrow @safe @nogc
    {
        return begin <= allowance;
    }


    bool endsWithin(coord_t allowance) const pure nothrow @safe @nogc
    {
        return end + allowance >= contig.length;
    }


    bool isFullyContained(coord_t allowance) const pure nothrow @safe @nogc
    {
        return beginsWithin(allowance) && endsWithin(allowance);
    }


    private bool _isValid() const pure nothrow @safe @nogc
    {
        return isValidInterval() && isIntervalInBounds;
    }


    private bool isValidInterval() const pure nothrow @safe @nogc
    {
        return begin <= end;
    }


    private bool isIntervalInBounds() const pure nothrow @safe @nogc
    {
        return contig.length == 0 || end <= contig.length;
    }
}


/// Validates the locus without triggering the invariant.
bool isValid(const Locus self) pure nothrow @safe @nogc
{
    return self._isValid();
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
    coord_t translateCoordLinear(string contig)(coord_t contigPos) const pure nothrow @safe @nogc
        if (contig.among("contigA"))
    {
        const translatedTracePoint = translateTracePoint!contig(contigPos, RoundingMode.floor);

        // cap the linear interpolation at the next trace point and at the
        // end coordinate of the alignment
        return min(
            translatedTracePoint.contigB + (contigPos - translatedTracePoint.contigA),
            translatedTracePoint.contigB + tracePointSpacing - 1,
            contigB.end - 1,
        );
    }


    ///
    TranslatedTracePoint translateTracePoint(string contig)(
        coord_t contigPos,
        RoundingMode roundingMode,
    ) const pure nothrow @safe @nogc if (contig.among("contigA", "contigB"))
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
    ) const pure nothrow @safe @nogc if (contig == "contigA")
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


    ///
    @property diff_t numDiffs() const pure nothrow @safe @nogc
    {
        return tracePoints.map!"a.numDiffs".sum;
    }


    ///
    @property size_t numTracePoints() const pure nothrow @safe @nogc
    in (tracePointSpacing > 0, "tracePointSpacing requried to compute tracePointCount")
    out (n; tracePoints.length == 0 || tracePoints.length == n,
        "computed tracePointCount does not match")
    {
        const firstTracePoint = floor(contigA.begin, tracePointSpacing);
        const lastTracePoint = ceil(contigA.end, tracePointSpacing);

        return max(1, (lastTracePoint - firstTracePoint) / tracePointSpacing);
    }


    int opCmp(ref const LocalAlignment other) const pure nothrow @safe @nogc
    {
        long cmp;

        enum compare(string field) = q{
            cmp = cast(long) this.}~field~q{ -
                  cast(long) other.}~field~q{;

            if (cmp != 0)
                return cmp.boundedConvert!int;
        };

        mixin(compare!"contigA.contig.id");
        mixin(compare!"contigB.contig.id");
        mixin(compare!"flags.complement");
        mixin(compare!"contigA.begin");
        mixin(compare!"contigA.end");
        mixin(compare!"contigB.begin");
        mixin(compare!"contigB.end");
        mixin(compare!"numDiffs");

        return 0;
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


class AlignmentSpec
{
    private Align_Spec* alignSpec;
    private alias alignSpec this;
    private float[4] _baseFrequency;

    /**
        Construct a new alignment specification.

        Params:
            averageCorrelation =
                The average correlation (1 - 2*error_rate) for the sought
                alignments. For Pacbio data we set this to .70 assuming an
                average of 15% error in each read.
            tracePointSpacing =
                The spacing interval for keeping trace points and segment
                differences.
            baseFrequency =
                A 4-element vector where afreq[0] = frequency of A, f(A),
                freq[1] = f(C), freq[2] = f(G), and freq[3] = f(T). This
                vector is part of the header of every DAZZ database.
            reach =
                If set extend the alignment to the boundary when reasonable,
                otherwise the terminate only at suffix-positive points.
    */
    this(
        double averageCorrelation,
        trace_point_t tracePointSpacing,
        float[4] baseFrequency,
        Flag!"reach" reach,
    )
    in (0.5 <= averageCorrelation && averageCorrelation <= 1.0, "averageCorrelation must be within [0.5, 1.0]")
    in (0 < tracePointSpacing, "tracePointSpacing must be positive")
    in (baseFrequency[].all!"0.0 <= a && a <= 1.0", "base frequencies must be within [0, 1]")
    in (absdiff(baseFrequency[].sum, 1.0) < 1e-6, "base frequencies sum up to 1")
    {
        this._baseFrequency = baseFrequency;
        this.alignSpec = New_Align_Spec(
            averageCorrelation,
            tracePointSpacing.boundedConvert!int,
            this._baseFrequency.ptr,
            cast(bool) reach,
        );
    }


    ~this()
    {
        if (alignSpec !is null)
            Free_Align_Spec(alignSpec);
    }


    @property double averageCorrelation() const nothrow @trusted @nogc
    {
        if (alignSpec is null)
            return typeof(return).init;

        return Average_Correlation(alignSpec);
    }


    @property trace_point_t tracePointSpacing() const nothrow @trusted @nogc
    {
        if (alignSpec is null)
            return typeof(return).init;

        return Trace_Spacing(alignSpec).boundedConvert!trace_point_t;
    }


    @property float[4] baseFrequency() const nothrow @trusted @nogc
    {
        return _baseFrequency;
    }


    @property Flag!"reach" reach() const nothrow @trusted @nogc
    {
        if (alignSpec is null)
            return typeof(return).init;

        return cast(typeof(return)) Overlap_If_Possible(alignSpec);
    }
}


/// Many routines like computeLocalAlignment need working storage that is more
/// efficiently reused with each call, rather than being allocated anew with
/// each call. Each *thread* can create a AlignmentWorkData object and this
/// object holds and retains the working storage for routines of this module
/// between calls to the routines.
class AlignmentWorkData
{
    private Work_Data* workData;
    private alias workData this;

    /// Construct opaque work data for alignment operations.
    this()
    {
        this.workData = New_Work_Data();
    }


    ~this()
    {
        if (workData !is null)
            Free_Work_Data(workData);
    }
}


///
LocalAlignment[2] computeLocalAlignment(
    const(char)[] aSequence,
    const(char)[] bSequence,
    coord_t antiBase,
    arithmetic_t antiLow,
    arithmetic_t antiHigh,
    const AlignmentSpec alignmentSpecification,
)
{
    return computeLocalAlignment(
        aSequence,
        bSequence,
        antiBase,
        antiLow,
        antiHigh,
        alignmentSpecification,
        new AlignmentWorkData(),
    );
}

/// ditto
LocalAlignment[2] computeLocalAlignment(
    const(char)[] aSequence,
    const(char)[] bSequence,
    coord_t antiBase,
    arithmetic_t antiLow,
    arithmetic_t antiHigh,
    const AlignmentSpec alignmentSpecification,
    AlignmentWorkData workData,
)
in (aSequence.length > 0 && aSequence[0].among(0, 1, 2, 3), "aSequence format must be numeric")
in (*(aSequence.ptr - 1) == 4 && *(aSequence.ptr + aSequence.length) == 4, "aSequence must be terminated")
in (bSequence.length > 0 && bSequence[0].among(0, 1, 2, 3), "bSequence format must be numeric")
in (*(bSequence.ptr - 1) == 4 && *(bSequence.ptr + bSequence.length) == 4, "bSequence must be terminated")
in (antiLow <= antiHigh, "invalid anti-diagonal: antiLow > antiHigh")
in (
    0 <= (cast(arithmetic_t) antiBase + antiLow)/2 &&
    (cast(arithmetic_t) antiBase + antiHigh)/2 <= aSequence.length,
    "anti-diagonal is out of bounds of read A"
)
in (
    0 <= (cast(arithmetic_t) antiBase - antiHigh)/2 &&
    (cast(arithmetic_t) antiBase - antiLow)/2 <= bSequence.length,
    "anti-diagonal is out of bounds of read B"
)
out (las; las[0].flags.disabled == las[1].flags.disabled)
out (las; las[0].contigA == las[1].contigB && las[0].contigB == las[1].contigA)
{
    // TODO handle complement correctly
    Path aPath;
    Alignment alignment;
    alignment.path = &aPath;
    alignment.aseq = cast(char*) aSequence.ptr;
    alignment.alen = aSequence.length.boundedConvert!int;
    alignment.bseq = cast(char*) bSequence.ptr;
    alignment.blen = bSequence.length.boundedConvert!int;

    auto bPathPtr = Local_Alignment(
        &alignment,
        workData,
        alignmentSpecification,
        antiLow,
        antiHigh,
        antiBase,
        -1,
        -1,
    );
    dazzlibEnforce(bPathPtr !is null, currentError.idup);
    auto bPath = *bPathPtr;

    enum Direction { forward, backward }
    const shouldDisable = (ref const Path path) =>
        path.tlen == 0 ||
        (path.tlen == 2 && (cast(const(trace_point_t)*) path.trace)[0 .. 2] == [0, 0]) ||
        (path.abpos == path.aepos && path.bbpos == path.bepos);
    const makeLocalAlignment = (ref const Path path, Direction direction) => LocalAlignment(
        Locus(
            Contig(0, (direction == Direction.forward
                ? alignment.alen
                : alignment.blen).boundedConvert!coord_t),
            path.abpos.boundedConvert!coord_t,
            path.aepos.boundedConvert!coord_t,
        ),
        Locus(
            Contig(0, (direction == Direction.forward
                ? alignment.blen
                : alignment.alen).boundedConvert!coord_t),
            path.bbpos.boundedConvert!coord_t,
            path.bepos.boundedConvert!coord_t,
        ),
        LocalAlignmentFlags(alignment.flags) | (shouldDisable(path)
            ? LocalAlignmentFlag.disabled
            : LocalAlignmentFlag.init
        ),
        alignmentSpecification.tracePointSpacing,
        path.tlen > 0
            ? (cast(const(TracePoint)*) path.trace)[0 .. path.tlen / 2].dup
            : [],
    );

    typeof(return) localAlignments = [
        makeLocalAlignment(aPath, Direction.forward),
        makeLocalAlignment(bPath, Direction.backward),
    ];

    return localAlignments;
}

unittest
{
    enum tracePointSpacing = 3;
    char[] aSequence = [4, 0, 1, 2, 3, 2, 1, 0, 4];
    char[] bSequence = [4, 0, 1, 2, 0, 2, 1, 0, 4];
    auto localAlignments = computeLocalAlignment(
        aSequence[1 .. $ - 1],
        bSequence[1 .. $ - 1],
        0,
        0,
        0,
        new AlignmentSpec(0.9, tracePointSpacing, [0.25f, 0.25f, 0.25f, 0.25f], No.reach),
    );

    assert(localAlignments == [
        LocalAlignment(
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            LocalAlignmentFlags(),
            tracePointSpacing,
            [
                TracePoint(0, 3),
                TracePoint(1, 3),
                TracePoint(0, 1),
            ]),
        LocalAlignment(
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            LocalAlignmentFlags(),
            tracePointSpacing,
            [
                TracePoint(0, 3),
                TracePoint(1, 3),
                TracePoint(0, 1),
            ]
        ),
    ]);
}

unittest
{
    enum tracePointSpacing = 3;
    char[] aSequence = [4, 0, 1, 2, 3, 2, 1, 0, 4];
    char[] bSequence = [4, 0, 0, 1, 2, 3, 2, 1, 0, 4];
    auto localAlignments = computeLocalAlignment(
        aSequence[1 .. $ - 1],
        bSequence[1 .. $ - 1],
        2,
        -2,
        2,
        new AlignmentSpec(0.9, tracePointSpacing, [0.25f, 0.25f, 0.25f, 0.25f], No.reach),
    );

    assert(localAlignments == [
        LocalAlignment(
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            Locus(
                Contig(0, 8),
                1,
                8,
            ),
            LocalAlignmentFlags(),
            tracePointSpacing,
            [
                TracePoint(0, 3),
                TracePoint(0, 3),
                TracePoint(0, 1),
            ]),
        LocalAlignment(
            Locus(
                Contig(0, 8),
                1,
                8,
            ),
            Locus(
                Contig(0, 7),
                0,
                7,
            ),
            LocalAlignmentFlags(),
            tracePointSpacing,
            [
                TracePoint(0, 2),
                TracePoint(0, 3),
                TracePoint(0, 2),
            ]
        ),
    ]);
}


/// Returns true if 16bits are required for encoding the trace at
/// tracePointSpacing.
bool isLargeTraceType(Int)(const Int tracePointSpacing) pure nothrow @safe @nogc if (isIntegral!Int)
{
    return tracePointSpacing > TRACE_XOVR;
}


/// Returns number of bytes used to encode a trace with tracePointSpacing.
uint tracePointBytes(Int)(const Int tracePointSpacing) pure nothrow @safe @nogc if (isIntegral!Int)
{
    return isLargeTraceType(tracePointSpacing)?  trace_point_t.sizeof : small_trace_point_t.sizeof;
}


/// Controls how to manage the buffer when reading data into memory.
enum BufferMode : ubyte
{
    /// Keep a single buffer and keep overwriting it with every new record.
    overwrite,
    /// Keep two buffers and keep overwriting them in turn. This is useful
    /// in combination with e.g. `std.algorithm.iteration.chunkBy`.
    doubleBuffer,
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
        if (bufferMode == BufferMode.doubleBuffer)
            assert(
                tracePointBuffer.length % 2 == 0,
                "buffer mode `doubleBuffer` requires double-sized buffer"
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

        auto dazzDb = new DazzDb(db, Yes.trimDb);

        return dazzDb
            .reads
            .map!(read => read.rlen.boundedConvert!coord_t)
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
    @property id_t numLocalAlignments() const pure nothrow @safe @nogc
    {
        return _numLocalAlignments;
    }


    ///
    @property trace_point_t tracePointSpacing() const pure nothrow @safe @nogc
    {
        return _tracePointSpacing;
    }


    ///
    @property bool empty() const pure nothrow @safe @nogc
    {
        return numLocalAlignmentsLeft == 0;
    }


    /// Number of local alignments left to read.
    @property size_t length() const pure nothrow @safe @nogc
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

protected:

    enum headerSize = long.sizeof + int.sizeof;

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
        if (bufferMode == BufferMode.doubleBuffer)
        {
            // keep switching buffers
            if (&tracePointBuffer[0] == &fullTracePointBuffer[0])
                tracePointBuffer = fullTracePointBuffer[$/2 .. $];
            else
                tracePointBuffer = fullTracePointBuffer[0 .. $/2];
        }

        readOverlapHead();
        fillInOverlapHead();
        if (aLengths.length > 0)
            fillInContigLength!'A'();
        if (bLengths.length > 0)
            fillInContigLength!'B'();

        if (skipTracePoints)
            skipTraceVector();
        else
            readTraceVector();

        if (!skipTracePoints)
        {
            if (!isLargeTraceType)
                rewriteTracePointBufferToLargeTracePointType();

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

        dazzlibEnforce(
            overlapHead.path.tlen % 2 == 0,
            "illegal value for tlen: must be multiple of 2",
        );
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


    void fillInContigLength(char read)() pure @safe if (read.among('A', 'B'))
    {
        auto contig = mixin("&currentLA.contig" ~ read ~ ".contig");
        auto contigLengths = mixin(read.toLower ~ "Lengths");

        if (contigLengths.length > 0)
        {
            dazzlibEnforce(
                contig.id - 1 < contigLengths.length,
                format!("contig"~read~".contig.id out of bounds: %d >= %d")(
                    contig.id - 1, contigLengths.length,
                ),
            );
            contig.length = contigLengths[contig.id - 1];
        }
    }


    void readTraceVector()
    {
        if (bufferMode == BufferMode.dynamic)
            tracePointBuffer = uninitializedArray!(TracePoint[])(numTracePoints);
        auto rawBuffer = getRawTracePointBuffer();
        las.rawReadAll(rawBuffer, unexpectedEOF!"tracePoints");
    }


    void skipTraceVector()
    {
        import core.stdc.stdio : SEEK_CUR;

        las.seek(traceVectorLength, SEEK_CUR);
    }


    void rewriteTracePointBufferToLargeTracePointType() pure nothrow
    {
        // make sure enough memory is allocated
        auto smallTraceVector = getRawTracePointBuffer();
        // select appropriate portion of tracePointBuffer
        auto largeTracePoints = tracePointBuffer.ptr[0 .. numTracePoints];
        // convert numbers in reverse order as to avoid overwriting
        foreach_reverse (i, ref largeTracePoint; largeTracePoints)
        {
            largeTracePoint = TracePoint(
                smallTraceVector[2 * i],
                smallTraceVector[2 * i + 1],
            );
        }
    }


    @property ubyte[] getRawTracePointBuffer() pure nothrow
    {
        // prepare bytes buffer
        auto bufferBytes = tracePointBuffer.length * TracePoint.sizeof;
        auto rawBuffer = (cast(ubyte*) tracePointBuffer.ptr)[0 .. bufferBytes];
        // reduce size (triggers bounds check)
        rawBuffer = rawBuffer[0 .. traceVectorLength];

        return rawBuffer;
    }


    @property size_t traceVectorLength() const pure nothrow @safe @nogc
    {
        return overlapHead.path.tlen * tracePointBytes(tracePointSpacing);
    }


    @property size_t numTracePoints() const pure nothrow @safe @nogc
    {
        assert(overlapHead.path.tlen % 2 == 0, "illegal value for tlen: must be multiple of 2");

        return overlapHead.path.tlen /2 ;
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


/// Enhaces `LocalAlignmentReader` by `popFrontExactly` and `reset(size_t)`.
class IndexedLocalAlignmentReader
{
    LocalAlignmentReader reader;
    alias reader this;
    const(LasIndex) index;


    this(
        const string lasFile,
        const LasIndex index,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer = [],
    )
    {
        this(
            lasFile,
            index,
            cast(string) null,
            cast(string) null,
            bufferMode,
            tracePointBuffer
        );
    }


    this(
        const string lasFile,
        const LasIndex index,
        string dbA,
        string dbB,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer = [],
    )
    {
        this.reader = new LocalAlignmentReader(lasFile, dbA, dbB, bufferMode, tracePointBuffer);
        this.index = index;
        assert(index.localAlignmentIndex.length - 1 == this.reader.numLocalAlignments);
    }


    /// Reset the reader to the i'th (0-based) local alignment.
    void reset(size_t i)
    {
        assert(i <= numLocalAlignments, "index out of bounds");

        reader.las.seek(index.localAlignmentIndex[i]);
        tracePointBuffer = fullTracePointBuffer;
        numLocalAlignmentsLeft = numLocalAlignments - cast(id_t) i;
        readLocalAlignment();
    }


    /// Advance the range by numSkip elements without reading them.
    void popFrontExactly(size_t numSkip)
    {
        assert(numSkip <= numLocalAlignmentsLeft, "not enough elements left to pop");

        const currentIndex = numLocalAlignments - numLocalAlignmentsLeft;
        assert(
            index.localAlignmentIndex[currentIndex + 1] == las.tell(),
            "index does not match LAS file",
        );

        reader.las.seek(index.localAlignmentIndex[currentIndex + numSkip]);
        numLocalAlignmentsLeft -= numSkip;
        readLocalAlignment();
    }

    alias popFrontN = popFrontExactly;
}

static assert(isInputRange!IndexedLocalAlignmentReader);
static assert(hasLength!IndexedLocalAlignmentReader);
static assert(is(ElementType!IndexedLocalAlignmentReader == LocalAlignment));


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

    auto recoveredLocalAlignments = new IndexedLocalAlignmentReader(
        lasFile,
        LasIndex.inferFrom(lasFile),
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

    auto recoveredLocalAlignments = new IndexedLocalAlignmentReader(
        lasFile,
        LasIndex.inferFrom(lasFile),
        dbFile,
        dbFile,
        BufferMode.dynamic,
    );

    assert(equal(testLocalAlignments, recoveredLocalAlignments));
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

    auto recoveredLocalAlignments = new IndexedLocalAlignmentReader(
        lasFile,
        LasIndex.inferFrom(lasFile),
        dbFile,
        dbFile,
        BufferMode.dynamic,
    );
    recoveredLocalAlignments.popFrontExactly(2);

    assert(equal(testLocalAlignments[2 .. $], recoveredLocalAlignments));

    recoveredLocalAlignments.reset(2);

    assert(equal(testLocalAlignments[2 .. $], recoveredLocalAlignments));
}


/// Speed up sparse accesses to LAS files.
struct LasIndex
{
    private
    {
        size_t[] localAlignmentIndex;
        size_t[] readIndex;
    }


    @property bool isInitialized() const pure nothrow @safe @nogc
    {
        return localAlignmentIndex.length > 0;
    }


    /// Construct index from LAS file or local alignments range by traversing it once.
    static LasIndex inferFrom(in string lasFile)
    {
        return inferFrom(localAlignmentReader(lasFile));
    }

    /// ditto
    static LasIndex inferFrom(LocalAlignmentReader lasReader)
    {
        auto lasIndex = LasIndex(uninitializedArray!(size_t[])(lasReader.numLocalAlignments + 1));

        size_t readIdx;
        id_t lastReadId;
        lasReader.reset();
        lasIndex.localAlignmentIndex[0] = LocalAlignmentReader.headerSize;
        foreach (i, localAlignment; lasReader.enumerate)
        {
            const filePos = lasReader.las.tell();
            lasIndex.localAlignmentIndex[i + 1] = filePos;

            if (localAlignment.contigA.contig.id != lastReadId)
            {
                lasIndex.storeReadIndex(readIdx++, i);
                lastReadId = localAlignment.contigA.contig.id;
            }
        }
        lasIndex.storeReadIndex(readIdx++, lasReader.numLocalAlignments);
        lasIndex.readIndex = lasIndex.readIndex[0 .. readIdx];

        return lasIndex;
    }


    private void storeReadIndex(size_t i, size_t value) pure nothrow @safe
    {
        if (readIndex.length == 0)
            readIndex = uninitializedArray!(size_t[])(10_000);

        while (i >= readIndex.length)
            readIndex.length += readIndex.length/3;

        readIndex[i] = value;
    }


    ///
    size_t numReads() const pure nothrow @safe @nogc
    {
        return readIndex.length - 1;
    }


    /// Get the index slice for the i'th (0-based) group of local alignments
    /// grouped by A-read.
    size_t[2] areadSlice(size_t i) const pure nothrow @safe @nogc
    {
        assert(i + 1 < readIndex.length);
        typeof(return) slice;
        slice[] = readIndex[i .. i + 2];

        return slice;
    }
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

    writeTestLas(lasFile);

    auto lasIndex = LasIndex.inferFrom(lasFile);

    assert(lasIndex.localAlignmentIndex == [12, 54, 96, 138, 180]);
    assert(lasIndex.readIndex == [0, 2, 4]);
    assert(lasIndex.areadSlice(0) == [0, 2]);
    assert(lasIndex.areadSlice(1) == [2, 4]);
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
    prepareBuffer(tracePointBuffer, bufferMode, AlignmentStats.inferFrom(lasFile));

    return new LocalAlignmentReader(
        lasFile,
        dbA,
        dbB,
        bufferMode,
        tracePointBuffer,
    );
}


/// Construct a pseudo-random access range that lazily reads local alignments
/// from lasFile. The range offers two additional operations `reset`, which
/// resets the range to a specific position, and `popFrontExactly`, which
/// skips a given number of records efficiently. If dbA and/or dbB are given,
/// the contig lengths will be filled in. If no tracePointBuffer is given, it
/// will be constructed  internally according to bufferMode.
///
/// See_also: IndexedLocalAlignmentReader, LasIndex
IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return indexedLocalAlignmentReader(lasFile, LasIndex.init, null, null, bufferMode, tracePointBuffer);
}

IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    LasIndex index,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return indexedLocalAlignmentReader(lasFile, index, null, null, bufferMode, tracePointBuffer);
}

/// ditto
IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    in string dbA,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return indexedLocalAlignmentReader(lasFile, LasIndex.init, dbA, null, bufferMode, tracePointBuffer);
}

/// ditto
IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    LasIndex index,
    in string dbA,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return indexedLocalAlignmentReader(lasFile, index, dbA, null, bufferMode, tracePointBuffer);
}

/// ditto
IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    in string dbA,
    in string dbB,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return indexedLocalAlignmentReader(lasFile, LasIndex.init, dbA, dbB, bufferMode, tracePointBuffer);
}

/// ditto
IndexedLocalAlignmentReader indexedLocalAlignmentReader(
    in string lasFile,
    LasIndex index,
    in string dbA,
    in string dbB,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    prepareBuffer(tracePointBuffer, bufferMode, AlignmentStats.inferFrom(lasFile));
    prepareIndex(index, lasFile);

    return new IndexedLocalAlignmentReader(
        lasFile,
        index,
        dbA,
        dbB,
        bufferMode,
        tracePointBuffer,
    );
}


private void prepareBuffer(
    ref TracePoint[] tracePointBuffer,
    BufferMode bufferMode,
    lazy AlignmentStats alignmentStats,
) pure @safe
{
    if (
        bufferMode.among(BufferMode.overwrite, BufferMode.doubleBuffer, BufferMode.preallocated) &&
        tracePointBuffer.length == 0
    )
    {
        if (bufferMode == BufferMode.overwrite)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentStats.maxTracePoints,
            );
        else if (bufferMode == BufferMode.doubleBuffer)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                2 * alignmentStats.maxTracePoints,
            );
        else if (bufferMode == BufferMode.preallocated)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentStats.numTracePoints,
            );
    }
}


private void prepareIndex(ref LasIndex index, string lasFile)
{
    if (!index.isInitialized)
        index = LasIndex.inferFrom(lasFile);
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
                        tp.numDiffs.boundedConvert!small_trace_point_t,
                        tp.numBasePairs.boundedConvert!small_trace_point_t,
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
    alias maxLocalAlignmentsPerRead = maxLocalAlignmentsPerContig;

    /// Total number of trace points
    size_t numTracePoints;

    /// Maximum number of trace points per local alignment
    size_t maxTracePoints;

    /// Maximum number of trace points per contig
    size_t maxTracePointsPerContig;
    alias maxTracePointsPerRead = maxTracePointsPerContig;

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

        alias updateMaxPerContig = () {
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
        };

        foreach (localAlignment; localAlignments)
        {
            if (lastContig != localAlignment.contigA.contig.id)
                updateMaxPerContig();

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

            headerData.numTracePoints += localAlignment.numTracePoints;
            headerData.maxTracePoints = max(
                headerData.maxTracePoints,
                localAlignment.numTracePoints,
            );

            ++numLocalAlignmentsSinceLastContig;
            numTracePointsSinceLastContig += localAlignment.numTracePoints;
            lastContig = localAlignment.contigA.contig.id;
        }

        updateMaxPerContig();

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
