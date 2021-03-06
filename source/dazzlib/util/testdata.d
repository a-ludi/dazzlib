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
import std.range;
import std.stdio;
import std.string;


private string removeWhitespace(string base64)
{
    return base64.tr(whitespace, "", "d");
}


@property string[] testFastaRecords() pure
{
    auto fakeSeq = (coord_t length) => cast(string) iota(length).map!"'a'".array;
    auto fakeRead = (id_t id, coord_t length) => format!">faked/%d/0_%d RQ=0.85\n%s"(
        id,
        id,
        fakeSeq(length),
    );

    return [
        fakeRead(1, 15),
        fakeRead(2, 16),
        fakeRead(3, 42),
        fakeRead(4, 42),
        fakeRead(5, 42),
        fakeRead(6, 42),
        fakeRead(7, 42),
        fakeRead(8, 42),
        fakeRead(9, 42),
        fakeRead(10, 42),
        fakeRead(11, 42),
        fakeRead(12, 42),
        fakeRead(13, 42),
        fakeRead(14, 42),
        fakeRead(15, 42),
        fakeRead(16, 42),
        fakeRead(17, 42),
        fakeRead(18, 42),
        fakeRead(19, 31),
        fakeRead(20, 33),
        fakeRead(21, 42),
        fakeRead(22, 42),
        fakeRead(23, 42),
        fakeRead(24, 42),
        fakeRead(25, 42),
        fakeRead(26, 42),
    ];
}


void writeTestDb(string stubFilename)
{
    enum dbStubBase64 = `
        ZmlsZXMgPSAgICAgICAgIDEKICAgICAgICAgMjYgc3Rkb3V0IGZha2VkCmJsb2NrcyA9ICAgICAgICAgMQpzaXplID0g
        ICAyMDAwMDAwMDAgY3V0b2ZmID0gICAgICAgICAwIGFsbCA9IDAKICAgICAgICAgMCAgICAgICAgIDAKICAgICAgICAy
        NiAgICAgICAgMjYK
    `.removeWhitespace;
    enum dbBasePairsBase64 = `
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    `.removeWhitespace;
    enum dbIndexBase64 = `
        GgAAABoAAAAAAAAAAAAAAAAAgD8AAAAAAAAAAAAAAAAqAAAAAAAAAPsDAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAPAAAAAAAAAAAAAAAAAAAAAAAAAP//
        ////////VQgAAAAAAAACAAAAEAAAAAAAAAAAAAAABAAAAAAAAAD//////////1UIAAAAAAAAAwAAACoAAAAAAAAAAAAA
        AAgAAAAAAAAA//////////9VCAAAAAAAAAQAAAAqAAAAAAAAAAAAAAATAAAAAAAAAP//////////VQgAAAAAAAAFAAAA
        KgAAAAAAAAAAAAAAHgAAAAAAAAD//////////1UIAAAAAAAABgAAACoAAAAAAAAAAAAAACkAAAAAAAAA//////////9V
        CAAAAAAAAAcAAAAqAAAAAAAAAAAAAAA0AAAAAAAAAP//////////VQgAAAAAAAAIAAAAKgAAAAAAAAAAAAAAPwAAAAAA
        AAD//////////1UIAAAAAAAACQAAACoAAAAAAAAAAAAAAEoAAAAAAAAA//////////9VCAAAAAAAAAoAAAAqAAAAAAAA
        AAAAAABVAAAAAAAAAP//////////VQgAAAAAAAALAAAAKgAAAAAAAAAAAAAAYAAAAAAAAAD//////////1UIAAAAAAAA
        DAAAACoAAAAAAAAAAAAAAGsAAAAAAAAA//////////9VCAAAAAAAAA0AAAAqAAAAAAAAAAAAAAB2AAAAAAAAAP//////
        ////VQgAAAAAAAAOAAAAKgAAAAAAAAAAAAAAgQAAAAAAAAD//////////1UIAAAAAAAADwAAACoAAAAAAAAAAAAAAIwA
        AAAAAAAA//////////9VCAAAAAAAABAAAAAqAAAAAAAAAAAAAACXAAAAAAAAAP//////////VQgAAAAAAAARAAAAKgAA
        AAAAAAAAAAAAogAAAAAAAAD//////////1UIAAAAAAAAEgAAACoAAAAAAAAAAAAAAK0AAAAAAAAA//////////9VCAAA
        AAAAABMAAAAfAAAAAAAAAAAAAAC4AAAAAAAAAP//////////VQgAAAAAAAAUAAAAIQAAAAAAAAAAAAAAwAAAAAAAAAD/
        /////////1UIAAAAAAAAFQAAACoAAAAAAAAAAAAAAMkAAAAAAAAA//////////9VCAAAAAAAABYAAAAqAAAAAAAAAAAA
        AADUAAAAAAAAAP//////////VQgAAAAAAAAXAAAAKgAAAAAAAAAAAAAA3wAAAAAAAAD//////////1UIAAAAAAAAGAAA
        ACoAAAAAAAAAAAAAAOoAAAAAAAAA//////////9VCAAAAAAAABkAAAAqAAAAAAAAAAAAAAD1AAAAAAAAAP//////////
        VQgAAAAAAAAaAAAAKgAAAAAAAAAAAAAAAAEAAAAAAAD//////////1UIAAAAAAAA
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


