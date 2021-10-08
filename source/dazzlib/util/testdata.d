/**
    Generate data for unit testing.

    Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dazzlib.util.testdata;

version (unittest):

import dazzlib.alignments;
import dazzlib.basictypes;
import dazzlib.db;
import dazzlib.core.c.DB;
import std.algorithm;
import std.array;
import std.ascii;
import std.base64;
import std.format;
import std.random;
import std.range;
import std.stdio;
import std.string;


private string removeWhitespace(string base64)
{
    return base64.tr(whitespace, "", "d");
}


@property string[] testSequencesUntrimmed() pure
{
    // get a reproducable pseudo-random sequence
    auto random = Random(0x0a46e94eU);
    auto fakeSeq = (coord_t length) => cast(string) iota(length)
        .map!(i => uniform(0, 4, random).predSwitch(
            0, 'a',
            1, 'c',
            2, 'g',
               't',
        ))
        .array;

    return [
        fakeSeq(5),
        fakeSeq(15),
        fakeSeq(16),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(31),
        fakeSeq(33),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
        fakeSeq(42),
    ];
}


@property string[] testSequencesTrimmed() pure
{
    return testSequencesUntrimmed[1 .. $];
}


@property string[] testSequencePrologsUntrimmed() pure
{
    return [
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "faked",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
        "test",
    ];
}


@property string[] testSequencePrologsTrimmed() pure
{
    return testSequencePrologsUntrimmed[1 .. $];
}


@property string[] testFastaRecords() pure
{
    auto fakeRead = (id_t id, string well, string sequence) => format!">%s/%d/0_%d RQ=0.85\n%s"(
        well,
        id,
        id,
        sequence,
    );

    return zip(testSequencePrologsUntrimmed, testSequencesUntrimmed)
        .enumerate(1)
        .map!(enumSeq => fakeRead(enumSeq.index, enumSeq.value[0], enumSeq.value[1]))
        .array;
}

unittest
{
    import std.stdio;

    //writefln!"%-(%s\n%)"(testFastaRecords());
}


void writeTestDb(string stubFilename)
{
    enum dbStubBase64 = `
        ZmlsZXMgPSAgICAgICAgIDIKICAgICAgICAgMTQgc3Rkb3V0IGZha2VkCiAgICAgICAgIDI3IHN0ZG91dCB0ZXN0CmJs
        b2NrcyA9ICAgICAgICAgMQpzaXplID0gICAyMDAwMDAwMDAgY3V0b2ZmID0gICAgICAgIDEwIGFsbCA9IDAKICAgICAg
        ICAgMCAgICAgICAgIDAKICAgICAgICAyNyAgICAgICAgMjYK
    `.removeWhitespace;
    enum dbBasePairsBase64 = `
        JsBqoH0c+KMrdF2NVckdwzxJvTawul2DMmpK2bUrSECP/yWjO9vBKQIIAOSfu9l7mqJfo0CAVj3hd4/02KClqtBAe58U
        KBkcW2eBsK8k/gv7SPlkwUUQa9t++vsel7uHI5D2rzeasP9TZt2+cPNmfQmDneePIwnwzDjVpAiCP0ugNdBATc2MQGjy
        TLh+0IKYRWLLKIXx0J3gRXS924Zac1CtK8BMQBXsnSrvFHLdgKih6JKd7EpaFwzgTbu2TwZD6RQWQeQEoiDT/0CxPqV6
        NIU3b9FaEE2XqIn50J0THjXwVjFV93khnBrW8hCOVUZ0ycGDveCiQGEqqW97bLUZcdMQ8PKP6wjJAzGl/AA=
    `.removeWhitespace;
    enum dbIndexBase64 = `
        GwAAABoAAAAKAAAAAAAAAAAAdz4AgIM+AAB+PgAAgj4qAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAP//
        ////////VQgAAAAAAAACAAAADwAAAAAAAAAAAAAAAgAAAAAAAAD//////////1UIAAAAAAAAAwAAABAAAAAAAAAAAAAA
        AAYAAAAAAAAA//////////9VCAAAAAAAAAQAAAAqAAAAAAAAAAAAAAAKAAAAAAAAAP//////////VQgAAAAAAAAFAAAA
        KgAAAAAAAAAAAAAAFQAAAAAAAAD//////////1UIAAAAAAAABgAAACoAAAAAAAAAAAAAACAAAAAAAAAA//////////9V
        CAAAAAAAAAcAAAAqAAAAAAAAAAAAAAArAAAAAAAAAP//////////VQgAAAAAAAAIAAAAKgAAAAAAAAAAAAAANgAAAAAA
        AAD//////////1UIAAAAAAAACQAAACoAAAAAAAAAAAAAAEEAAAAAAAAA//////////9VCAAAAAAAAAoAAAAqAAAAAAAA
        AAAAAABMAAAAAAAAAP//////////VQgAAAAAAAALAAAAKgAAAAAAAAAAAAAAVwAAAAAAAAD//////////1UIAAAAAAAA
        DAAAACoAAAAAAAAAAAAAAGIAAAAAAAAA//////////9VCAAAAAAAAA0AAAAqAAAAAAAAAAAAAABtAAAAAAAAAP//////
        ////VQgAAAAAAAAOAAAAKgAAAAAAAAAAAAAAeAAAAAAAAAD//////////1UIAAAAAAAADwAAACoAAAAAAAAAAAAAAIMA
        AAAAAAAA//////////9VCAAAAAAAABAAAAAqAAAAAAAAAAAAAACOAAAAAAAAAP//////////VQgAAAAAAAARAAAAKgAA
        AAAAAAAAAAAAmQAAAAAAAAD//////////1UIAAAAAAAAEgAAACoAAAAAAAAAAAAAAKQAAAAAAAAA//////////9VCAAA
        AAAAABMAAAAqAAAAAAAAAAAAAACvAAAAAAAAAP//////////VQgAAAAAAAAUAAAAHwAAAAAAAAAAAAAAugAAAAAAAAD/
        /////////1UIAAAAAAAAFQAAACEAAAAAAAAAAAAAAMIAAAAAAAAA//////////9VCAAAAAAAABYAAAAqAAAAAAAAAAAA
        AADLAAAAAAAAAP//////////VQgAAAAAAAAXAAAAKgAAAAAAAAAAAAAA1gAAAAAAAAD//////////1UIAAAAAAAAGAAA
        ACoAAAAAAAAAAAAAAOEAAAAAAAAA//////////9VCAAAAAAAABkAAAAqAAAAAAAAAAAAAADsAAAAAAAAAP//////////
        VQgAAAAAAAAaAAAAKgAAAAAAAAAAAAAA9wAAAAAAAAD//////////1UIAAAAAAAAGwAAACoAAAAAAAAAAAAAAAIBAAAA
        AAAA//////////9VCAAAAAAAAA==
`.removeWhitespace;

    auto dbFiles = EssentialDbFiles(stubFilename);

    File(dbFiles.stub, "wb").rawWrite(Base64.decode(dbStubBase64));
    File(dbFiles.basePairs, "wb").rawWrite(Base64.decode(dbBasePairsBase64));
    File(dbFiles.index, "wb").rawWrite(Base64.decode(dbIndexBase64));
}


@property coord_t[][] testMaskData() pure
{
    auto fakeIntervals = (id_t id, coord_t length) => id % 2 == 0
        ? [id, length]
        : [0, length/2, length/2 + 1, length];

    return [
        fakeIntervals(1, 15),
        fakeIntervals(2, 16),
        fakeIntervals(3, 42),
        fakeIntervals(4, 42),
        fakeIntervals(5, 42),
        fakeIntervals(6, 42),
        fakeIntervals(7, 42),
        fakeIntervals(8, 42),
        fakeIntervals(9, 42),
        fakeIntervals(10, 42),
        fakeIntervals(11, 42),
        fakeIntervals(12, 42),
        fakeIntervals(13, 42),
        fakeIntervals(14, 42),
        fakeIntervals(15, 42),
        fakeIntervals(16, 42),
        fakeIntervals(17, 42),
        fakeIntervals(18, 42),
        fakeIntervals(19, 31),
        fakeIntervals(20, 33),
        fakeIntervals(21, 42),
        fakeIntervals(22, 42),
        fakeIntervals(23, 42),
        fakeIntervals(24, 42),
        fakeIntervals(25, 42),
        fakeIntervals(26, 42),
    ];
}


void writeTestMask(string stubFilename, string maskName)
{
    auto dbFiles = EssentialDbFiles(stubFilename);
    auto annoFile = File(dbFiles.auxiliaryFile("."~maskName~".anno"), "wb");

    auto numReads = cast(int) testMaskData.length;
    auto size = 0;
    int currentContig = 1;

    annoFile.rawWrite([numReads, size]);
    annoFile.rawWrite([int64(0)]);
    annoFile.rawWrite(testMaskData
        .cumulativeFold!((dataPointer, intervalsData) =>
            dataPointer + cast(int64) (intervalsData.length * int.sizeof)
        )(int64(0))
        .array
    );

    auto dataFile = File(dbFiles.auxiliaryFile("."~maskName~".data"), "wb");
    dataFile.rawWrite(testMaskData.join());

    //enum annoBase64 = ``.removeWhitespace;
    //enum dataBase64 = ``.removeWhitespace;

    //auto dbFiles = EssentialDbFiles(stubFilename);

    //File(dbFiles.auxiliaryFile("."~maskName~".anno"), "wb").rawWrite(Base64.decode(annoBase64));
    //File(dbFiles.auxiliaryFile("."~maskName~".data"), "wb").rawWrite(Base64.decode(dataBase64));
}


@property LocalAlignment[] testLocalAlignments() pure nothrow @safe
{
    enum trace_point_t tracePointSpacing = 100;

    return [
        LocalAlignment(
            Locus(Contig(1, 15), 3, 4),
            Locus(Contig(2, 16), 5, 6),
            LocalAlignmentFlags(
                LocalAlignmentFlag.start,
                LocalAlignmentFlag.best,
            ),
            tracePointSpacing,
            [TracePoint(7, 1)],
        ),
        LocalAlignment(
            Locus(Contig(1, 15), 12, 13),
            Locus(Contig(2, 16), 14, 15),
            LocalAlignmentFlags(
                LocalAlignmentFlag.next,
            ),
            tracePointSpacing,
            [TracePoint(16, 1)],
        ),
        LocalAlignment(
            Locus(Contig(19, 31), 21, 22),
            Locus(Contig(20, 33), 23, 24),
            LocalAlignmentFlags(
                LocalAlignmentFlag.start,
                LocalAlignmentFlag.complement,
            ),
            tracePointSpacing,
            [TracePoint(25, 1)],
        ),
        LocalAlignment(
            Locus(Contig(19, 31), 30, 31),
            Locus(Contig(20, 33), 32, 33),
            LocalAlignmentFlags(
                LocalAlignmentFlag.next,
                LocalAlignmentFlag.complement,
            ),
            tracePointSpacing,
            [TracePoint(0, 1)],
        ),
    ];
}


void writeTestLas(string filename)
{
    enum lasBase64 = `
        BAAAAAAAAABkAAAAAgAAAAcAAAADAAAABQAAAAQAAAAGAAAAFAAAAAAAAAABAAAAAAAAAAcBAgAAABAAAAAMAAAADgAA
        AA0AAAAPAAAACAAAAAAAAAABAAAAAAAAABABAgAAABkAAAAVAAAAFwAAABYAAAAYAAAABQAAABIAAAATAAAAAAAAABkB
        AgAAAAAAAAAeAAAAIAAAAB8AAAAhAAAACQAAABIAAAATAAAAAAAAAAAB
    `.removeWhitespace;

    File(filename, "wb").rawWrite(Base64.decode(lasBase64));
}


void writeEmptyTestLas(string filename)
{
    enum emptyLasBase64 = `AAAAAAAAAABkAAAA`;

    File(filename, "wb").rawWrite(Base64.decode(emptyLasBase64));
}


