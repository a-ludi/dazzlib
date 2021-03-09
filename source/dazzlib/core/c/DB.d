/**
    Compressed data base module.  Auxiliary routines to open and manipulate a
    data base for which the sequence and read information are separated into
    two separate files, and the sequence is compressed into 2-bits for each
    base. Support for tracks of additional information, and trimming according
    to the current partition.  Eventually will also support compressed quality
    information.

    Copyright: © 2014 Dr. Eugene W. Myers <gene.myers@gmail.com>. All rights
               reserved.
               © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>,
             Dr. Eugene W. Myers <gene.myers@gmail.com>
*/
module dazzlib.core.c.DB;

import dazzlib.core.c.QV;
import core.stdc.stdio;
import std.string;
import std.typecons;

extern (C):


/// Name of calling program which will be reported in error messages.
extern __gshared immutable(char)* Prog_Name;

// set default Prog_Name
extern (D) static this() { Prog_Name = "dazzlib"; }


/// Error buffer for Dazzler routines
extern __gshared char[1000] Ebuffer;


@property char[] currentError() nothrow
{
    return fromStringz(Ebuffer.ptr);
}


/// By default, auxiliary DB files start with a `.` so they are "hidden".
/// Define version `ShowAuxiliaryDbFiles` if you don't want this.
version (ShowAuxiliaryDbFiles)
    enum auxiliaryDbFilePrefix = "";
else
    enum auxiliaryDbFilePrefix = ".";


/// 8-byte unsigned integer type.
alias uint8 = ubyte;
/// 16-byte unsigned integer type.
alias uint16 = ushort;
/// 32-byte unsigned integer type.
alias uint32 = uint;
/// 64-byte unsigned integer type.
alias uint64 = ulong;
/// 8-byte signed integer type.
alias int8 = byte;
/// 16-byte signed integer type.
alias int16 = short;
/// 32-byte signed integer type.
alias int32 = int;
/// 64-byte signed integer type.
alias int64 = long;
/// 32-byte float type.
alias float32 = float;
/// 64-byte float type.
alias float64 = double;
/// Type for C-style strings.
alias cstring = immutable(char)*;


/**
    Represent a single read or contig.

    Fields have different interpretations if a .db versus a .dam
*/
struct DAZZ_READ
{
    ///  Well # (DB), Contig # (DAM)
    int origin;
    /// Length of the sequence (Last pulse = fpulse + rlen)
    int rlen;
    /// First pulse (DB), left index of contig in scaffold (DAM)
    int fpulse;
    /// Offset (in bytes) of compressed read in 'bases' file, or offset of
    /// uncompressed bases in memory block
    int64 boff;
    /// Offset (in bytes) of compressed quiva streams in '.qvs' file (DB),
    /// Offset (in bytes) of scaffold header string in '.hdr' file (DAM)
    /// 4 compressed shorts containing snr info if an arrow DB.
    int64 coff;
    /// QV of read + flags above (DB only)
    int rawFlags;


    /// Mask for 3-digit quality value
    enum qvMask = 0x03ff;


    /// Read flags.
    enum Flag : int
    {
        /// This is the second or later of a group of subreads from a given insert
        ccs = 0x0400,
        /// This is the "best" subread of a given insert (may be the only 1)
        best = 0x0800,
    }


    @property int qv() const pure nothrow @safe
    {
        return rawFlags & qvMask;
    }


    @property void qv(int value) pure nothrow @safe
    {
        rawFlags &= (rawFlags & ~qvMask) | (value & qvMask);
    }


    @property BitFlags!Flag flags() const pure nothrow @safe
    {
        return BitFlags!Flag(cast(Flag) (rawFlags & ~qvMask));
    }


    @property void flags(BitFlags!Flag value) pure nothrow @safe
    {
        rawFlags &= ((cast(int) value) & ~qvMask) | (rawFlags & qvMask);
    }
}


///  A track can be of 3 types:
///
///  - `data == NULL`: there are nreads 'anno' records of size 'size'.
///  - `data != NULL && size == 4`: anno is an array of `nreads+1` int's and
///    data[anno[i]..anno[i+1]) contains the variable length data
///  - `data != NULL && size == 8`: anno is an array of `nreads+1` int64's and
///    data[anno[i]..anno[i+1]) contains the variable length data
struct DAZZ_TRACK
{
    /// Link to next track
    DAZZ_TRACK* next;
    /// Symbolic name of track
    char* name;
    /// Size in bytes of anno records
    int size;
    /// Number of reads in track
    int nreads;
    /// over [0,nreads]: read i annotation: int, int64, or 'size' records
    void* anno;
    /// length of track data for read i (if data != NULL)
    int* alen;
    /// data[anno[i] .. anno[i]+alen[i[) is data for read i (if data != NULL)
    void* data;
    /// Is track data loaded in memory?
    int loaded;
    /// Largest read data segment in bytes
    int64 dmax;
}

/// The tailing part of a .anno track file can contain meta-information
/// produced by the command that produced the track.  For example, the
/// coverage, or good/bad parameters for trimming, or even say a histogram of
/// QV values.  Each item is an array of 'nelem' 64-bit ints or floats
/// ('vtype' = DB_INT or DB_REAL), has a 'name' string that describes it, and
/// an indicator as to whether the values should be equal accross all block
/// tracks, or summed accross all block tracks (by Catrack).  'value' points
/// at the array of values
struct DAZZ_EXTRA
{
    /// Element type.
    enum Type : int
    {
        /// Element type is `int64`
        int_ = 0,
        /// Element type is `float64`
        real_ = 1,
    }

    /// Accumulator when combining tracks of several blocks.
    enum Accumulator : int
    {
        /// Data must be exactly the same for all blocks.
        exact = 0,
        /// Sum data of all blocks.
        sum = 0,
    }

    /// Vector element type
    Type vtype;
    /// Number of elements
    int nelem;
    /// Accumulator operation
    Accumulator accum;
    /// Extra name
    char* name;
    /// Vector of `nelem` elements of `vtype`.
    void* value;
}

/// The information for accessing QV streams is in a DAZZ_QV record that is a
/// "pseudo-track" named ".@qvs" and is always the first track record in the
/// list (if present). Since normal track names cannot begin with a .
/// (this is enforced), this pseudo-track is never confused with a normal
/// track.
struct DAZZ_QV
{
    /// Link to next QV track
    DAZZ_QV* next;
    /// Name of the track
    char* name;
    /// # of coding tables
    int ncodes;
    /// array [0..ncodes-1] of coding schemes (see QV.h)
    QVcoding* coding;
    /// for i in [0,db->nreads-1]: read i should be decompressed with scheme
    /// coding[table[i]]
    uint16* table;
    /// the open file pointer to the .qvs file
    FILE* quiva;
}

/// The information for accessing Arrow streams is in a DAZZ_ARW record that
/// is a "pseudo-track" named ".@arw" and is always the first track record in
/// the list (if present). Since normal track names cannot begin with a .
/// (this is enforced), this pseudo-track is never confused with a normal
/// track.
struct DAZZ_ARROW
{
    /// Link to next Arrow stram
    DAZZ_ARROW* next;
    /// Name of the track
    char* name;
    /// offset in file or memory of arrow vector for read i
    int64* aoff;
    /// FILE * to the .arw file if not loaded, memory block otherwise
    void* arrow;
    /// Are arrow vectors loaded in memory?
    int loaded;
}

/// Every DB is referred to by an ASCII stub file with extension .db or .dam.
/// This file contains the information about the SMRT cells in the DB and the
/// current division of the DB into blocks for HPC processing. This file can
/// be read into the following data structure:
struct DAZZ_STUB
{
    /// Number of files/SMRT cells in DB
    int nfiles;
    /// [0..nfiles) = # of reads from cell
    int* nreads;
    /// [0..nfiles) = file name of cell
    char** fname;
    /// [0..nfiles) = fasta header prolog for cell
    char** prolog;
    /// Keep only best read from each well?
    int all;
    /// Trim reads less than cutoff
    int cutoff;
    /// Target size for blocks
    int64 bsize;
    /// Number of blocks for DB
    int nblocks;
    /// [0..nblcoks] = index of 1st read in block in untrimmed DB
    int* ublocks;
    /// [0..nblcoks] = index of 1st read in block in trimmed DB
    int* tblocks;
}

/// Suppose DB is the name of an original database. Then there will be files
/// .DB.idx, .DB.bps, .DB.qvs, and files .DB.<track>.anno and DB.<track>.data
/// where <track> is a track name (not containing a . !).
///
/// A DAM is basically a DB except that:
///
/// 1. there are no QV's, instead .coff points to the '\0' terminated fasta
///    header of the read in an additional file: .DB.hdr
/// 2. .origin contains the contig # of the read within a fasta entry
///    (assembly sequences contain N-separated contigs), and .fpulse the first
///    base of the contig in the fasta entry
///
/// The DB record holds all information about the current state of an active
/// DB including an array of DAZZ_READS, one per read, and a linked list of
/// DAZZ_TRACKs the first of which is always a DAZZ_QV pseudo-track (if the
/// QVs have been loaded).
///
/// Fields maxlen, totlen and nreads are set with respect to "active" part of
/// DB (all vs block, untrimmed vs trimmed).
///
/// In order to avoid forcing users to have to rebuild all thier DBs to
/// accommodate the addition of fields for the size of the actively loaded
/// trimmed and untrimmed blocks, an additional read record is allocated in
/// "reads" when a DB is loaded into memory (reads[-1]) and the two desired
/// fields are crammed into the first two integer spaces of the record.
struct DAZZ_DB
{
    // DB flags.
    enum Flag : int
    {
        /// DB is an arrow DB
        arrow = 0x2,
        /// all wells are in the trimmed DB
        all = 0x1,
    }

    /// Total number of reads in untrimmed DB
    int ureads;
    /// Total number of reads in trimmed DB
    int treads;
    /// Minimum read length in block (-1 if not yet set)
    int cutoff;
    /// DB_ALL | DB_ARROW
    Flag allarr;
    /// frequency of A, C, G, T, respectively
    float[4] freq;

    /// length of maximum read (initially over all DB)
    int maxlen;
    /// total # of bases (initially over all DB)
    int64 totlen;
    /// # of reads in actively loaded portion of DB
    int nreads;
    /// DB has been trimmed by cutoff/all
    int trimmed;
    /// DB block (if > 0), total DB (if == 0)
    int part;
    /// Index of first read in block (without trimming)
    int ufirst;
    /// Index of first read in block (with trimming)
    int tfirst;

    /// Root name of DB for .bps, .qvs, and tracks
    char* path;
    /// Are reads loaded in memory?
    int loaded;
    /// file pointer for bases file (to fetch reads from), or memory pointer to uncompressed block of all sequences.
    void* bases;
    /// Array [-1..nreads] of DAZZ_READ
    DAZZ_READ* reads;
    /// Linked list of loaded tracks
    DAZZ_TRACK* tracks;
}


/// Longest file name or fasta header line
enum MAX_NAME = 10000;


/// DB STUB FILE FORMAT = NFILE FDATA^nfile NBLOCK PARAMS BDATA^nblock
enum DbStubFileFormat : string
{
    /// number of files
    nfile = "files = %9d\n",
    /// last read index + 1, fasta name, file prolog
    fdata = "  %9d %s %s\n",
    /// number of blocks
    nblock = "blocks = %9d\n",
    /// block size, len cutoff, all in well
    params = "size = %11d cutoff = %9d all = %1d\n",
    /// First read index (untrimmed), first read index (trimmed)
    bdata = " %9d %9d\n",
    /// combined format in D's format syntax
    combined = nfile ~ "%(%(" ~ fdata ~ "%)%)" ~ nblock ~ params ~ "%(%(" ~ bdata ~ "%)%)",
}


/// Used to tell Read_DB_Stub which part of the DB stub to read.
enum StubPart : int
{
    ///
    nreads = 0x1,
    ///
    files = 0x2,
    ///
    prologs = 0x4,
    ///
    blocks = 0x8,
    ///
    all = nreads | files | prologs | blocks,
}


/// Read the specified contents of the DB stub file at "path" and return it
/// encoded in a DAZZ_STUB structure. This is allocated by the routine. "path"
/// is assumed to be the complete name of the file. If all flags are off, then
/// just the scalar parts of the stub are returned (i.e. nfiles, all, cutoff,
/// bsize, nblocks). Returns NULL if an error occured in INTERACTIVE mode.
DAZZ_STUB* Read_DB_Stub (cstring path, StubPart what);

unittest
{
    import dazzlib.util.tempfile;
    import std.algorithm;
    import std.conv;
    import std.file;
    import std.typecons;
    import std.string;

    auto testStub = mkstemp(".unittest-XXXXXX", ".db");
    scope (exit)
        remove(testStub.name);

    enum nfiles = 3;
    enum fdata = [
        tuple(1, "file1", "prolog1"),
        tuple(2, "file2", "prolog2"),
        tuple(3, "file3", "prolog3"),
    ];
    enum nblocks = 5;
    enum bsize = 200000000;
    enum cutoff = 1000;
    enum all = 1;
    enum bdata = [
        tuple(0, 1),
        tuple(2, 3),
        tuple(4, 5),
        tuple(6, 7),
        tuple(8, 9),
        tuple(10, 11),
    ];
    testStub.file.writef!(cast(string) DbStubFileFormat.combined)(
        nfiles,
        fdata,
        nblocks,
        bsize,
        cutoff,
        all,
        bdata,
    );
    testStub.file.close();

    auto dbStub = Read_DB_Stub(toStringz(testStub.name), StubPart.all);
    scope (exit)
        Free_DB_Stub(dbStub);

    assert(dbStub !is null, currentError);
    assert(dbStub.nfiles == nfiles);
    assert(equal(
        dbStub.nreads[0 .. dbStub.nfiles],
        fdata.map!"a[0]",
    ));
    assert(equal(
        dbStub.fname[0 .. dbStub.nfiles].map!fromStringz,
        fdata.map!"a[1]",
    ));
    assert(equal(
        dbStub.prolog[0 .. dbStub.nfiles].map!fromStringz,
        fdata.map!"a[2]",
    ));
    assert(dbStub.all == all);
    assert(dbStub.cutoff == cutoff);
    assert(dbStub.bsize == bsize);
    assert(dbStub.nblocks == nblocks);
    assert(equal(
        dbStub.ublocks[0 .. dbStub.nblocks + 1],
        bdata.map!"a[0]",
    ));
    assert(equal(
        dbStub.tblocks[0 .. dbStub.nblocks + 1],
        bdata.map!"a[1]",
    ));
}


/// Read the DB stub file "path" and extract the read index range [*first,*last)
/// for block n, for the trimmed DB if trim is set, the untrimmed DB otherwise.
/// If n is out of range first and last will be set to -1.  Returns 0 unless
/// an error occurs in INTERACTIVE mode in which case it returns 1.
int Fetch_Block_Range (cstring path, int trim, int n, int* first, int* last);


/// Free a DAZZ_STUB data structure returned by Read_DB_Stub
void Free_DB_Stub (DAZZ_STUB* stub);


/// Open the given database or dam, "path", into the supplied DAZZ_DB record
/// "db". If the name has a part # in it then just the part is opened. The
/// index array is allocated (for all or just the part) and read in.
///
/// Returns:
///     - `-1`: The DB could not be opened for a reason reported by the routine
///       to EPLACE
///     - `0`: Open of DB proceeded without mishap
///     - `1`: Open of DAM proceeded without mishap
int Open_DB (cstring path, DAZZ_DB* db);


/// Trim the DB or part thereof and all loaded tracks according to the cutoff
/// and all settings of the current DB partition.  Reallocate smaller memory
/// blocks for the information kept for the retained reads.
void Trim_DB (DAZZ_DB* db);


/// Return the size in bytes of the given DB
int64 sizeof_DB (DAZZ_DB* db);


/// For the DB or DAM "path" = "prefix/root.[db|dam]", find all the files for
/// that DB, i.e. all those of the form "prefix/[.]root.part" and call actor
/// with the complete path to each file pointed at by path, and the suffix of
/// the path by extension.  The . proceeds the root name if the defined
/// constant HIDE_FILES is set.  Always the first call is with the path
/// "prefix/root.[db|dam]" and extension "db" or "dam". There will always be
/// calls for "prefix/[.]root.idx" and "prefix/[.]root.bps". All other calls
/// are for *tracks* and so this routine gives one a way to know all the
/// tracks associated with a given DB.
///
/// Returns:
///     - `-1`: if the path could not be found
///     - `1`: if an error (reported to EPLACE) occured and INTERACTIVE is
///       defined.
///     - `0`: Otherwise.
int List_DB_Files (cstring path, void function (cstring path, cstring extension) actor);


/// Shut down an open 'db' by freeing all associated space, including tracks
/// and QV structures, and any open file pointers.  The record pointed at by
/// db however remains (the user supplied it and so should free it).
void Close_DB (DAZZ_DB* db);


/// Allocate and return a buffer big enough for the largest read in 'db'.
/// **NB** free(x-1) if x is the value returned as *prefix* and suffix '\0'
/// (4)-byte are needed by the alignment algorithms. If cannot allocate memory
/// then return NULL if INTERACTIVE is defined, or print error to stderr and
/// exit otherwise.
char* New_Read_Buffer (DAZZ_DB* db);


/// Load into 'read' the i'th read in 'db'.  As a lower case ascii string if
/// ascii is 1, an upper case ascii string if ascii is 2, and a numeric string
/// over 0(A), 1(C), 2(G), and 3(T) otherwise. A '\0' (or 4) is prepended and
/// appended to the string so it has a delimeter for traversals in either
/// direction.  A non-zero value is returned if an error occured and
/// INTERACTIVE is defined.
int Load_Read (DAZZ_DB* db, int i, char* read, int ascii);


/// Load into 'read' the subread [beg,end] of the i'th read in 'db' and return
/// a pointer to the the start of the subinterval (not necessarily equal to
/// read !!! ). As a lower case ascii string if ascii is 1, an upper case
/// ascii string if ascii is 2, and a numeric string over 0(A), 1(C), 2(G),
/// and 3(T) otherwise. A '\0' (or 4) is prepended and appended to the string
/// holding the substring so it has a delimeter for traversals in either
/// direction. A NULL pointer is returned if an error occured and INTERACTIVE
/// is defined.
char* Load_Subread (DAZZ_DB* db, int i, int beg, int end, char* read, int ascii);


/// Allocate a block big enough for all the uncompressed read sequences and
/// read and uncompress the reads into it, reset the 'boff' in each read
/// record to be its in-memory offset, and set the bases pointer to point at
/// the block after closing the bases file.  Return with a zero, except when
/// an error occurs and INTERACTIVE is defined in which case return wtih 1.
int Load_All_Reads (DAZZ_DB* db, int ascii);


/// If the Arrow pseudo track is not already in db's track list, then load it
/// and set it up. The database reads must not have been loaded with
/// Load_All_Reads yet. -1 is returned if a .arw file is not present, and 1
/// is returned if an error (reported to EPLACE) occured and INTERACTIVE is
/// defined.  Otherwise a 0 is returned.
int Open_Arrow (DAZZ_DB* db);


/// Exactly the same as Load_Read, save the arrow information is loaded, not
/// the DNA sequence, and there is only a choice between numeric (0) or
/// ascii (1);
int Load_Arrow (DAZZ_DB* db, int i, char* read, int ascii);


/// Allocate a block big enough for all the uncompressed Arrow vectors, read
/// them into it, reset the 'off' in each arrow record to be its in-memory
/// offset, and set the arrow pointer to point at the block after closing the
/// arrow file. If ascii is non-zero then the arrows are converted to 0123
/// ascii, otherwise the arrows are left as numeric strings over [0-3].
int Load_All_Arrows (DAZZ_DB* db, int ascii);


/// Remove the Arrow pseudo track, all space associated with it, and close
/// the .arw file.
void Close_Arrow (DAZZ_DB*);


/// Type of track
enum TrackKind : int
{
    ///
    custom = 0,
    ///
    mask = 1,
}


/// Look up the file and header in the file of the indicated track. In
/// addition, if opened (0 or 1 returned), then kind points at an integer
/// indicating the type of track.
///
/// Returns:
///     - '1': Track is for trimmed DB
///     - '0': Track is for untrimmed DB
///     - '-1': Track is not the right size of DB either trimmed or untrimmed
///     - '-2': Could not find the track
int Check_Track (DAZZ_DB* db, cstring track, TrackKind* kind);


/// If track is not already in the db's track list, then allocate all the
/// storage for the anno index, read it in from the appropriate file, add it
/// to the track list, and return a pointer to the newly created DAZZ_TRACK
/// record.  If the track does not exist or cannot be opened for some
/// reason, then NULL is returned if INTERACTIVE is defined.  Otherwise the
/// routine prints an error message to stderr and exits if an error occurs,
/// and returns with NULL only if the track does not exist.
DAZZ_TRACK* Open_Track (DAZZ_DB* db, cstring track);


/// Allocate a data buffer large enough to hold the longest read data block
/// that will occur in the track.  If cannot allocate memory then return
/// NULL if INTERACTIVE is defined, or print error to stderr and exit
/// otherwise.
void* New_Track_Buffer (DAZZ_TRACK* track);

/// Load into 'data' the read data block for read i's "track" data.  Return
/// the length of the data in bytes, unless an error occurs and INTERACTIVE
/// is defined in which case return wtih -1.
int Load_Track_Data (DAZZ_TRACK* track, int i, void* data);


/// Allocate a block big enough for all the track data and read the data
/// into it, reset the 'off' in each anno pointer to be its in-memory
/// offset, and set the data pointer to point at the block after closing the
/// data file.  Return with a zero, except when an error occurs and
/// INTERACTIVE is defined in which case return wtih 1.
int Load_All_Track_Data (DAZZ_TRACK* track);


/// Assumming file pointer for afile is correctly positioned at the start of
/// an extra item, and aname is the name of the .anno file, decode the value
/// present and place it in extra if extra->nelem == 0, otherwise reduce the
/// value just read into extra according according to the directive given by
/// accum'.  Leave the read pointer at the next extra or end-of-file.
///
/// Returns:
///     - `1` if at the end of file,
///     - `0` if item was read and folded correctly,
///     - `-1` if there was a system IO or allocation error (if interactive), and
///     - `-2` if the new value could not be reduced into the current value of extra (interactive)
int Read_Extra (FILE* afile, cstring aname, DAZZ_EXTRA* extra);


/// Write extra record to end of file afile and advance write pointer If
/// interactive, then return non-zero on error, if batch, then print and halt
/// if an error
int Write_Extra (FILE* afile, DAZZ_EXTRA* extra);


/// If track is on the db's track list, then it is removed and all storage
/// associated with it is freed.
void Close_Track (DAZZ_DB* db, DAZZ_TRACK* track);


/// If QV pseudo track is not already in db's track list, then load it and
/// set it up. The database must not have been trimmed yet.  -1 is returned
/// if a .qvs file is not present, and 1 is returned if an error (reported to
/// EPLACE) occured and INTERACTIVE is defined.  Otherwise a 0 is returned.
int Open_QVs (DAZZ_DB* db);


/// QV line index.
enum QvLineIndex : size_t
{
    /// The deletion QVs are x[DEL_QV] if x is the buffer returned by
    /// New_QV_Buffer
    delQv = 0,
    /// The deleted characters
    delTag = 1,
    /// The insertion QVs
    insQv = 2,
    /// The substitution QVs
    subQv = 3,
    /// The merge QVs
    mrgQv = 4,
}


/// Allocate a set of 5 vectors large enough to hold the longest QV stream
/// that will occur in the database.  If cannot allocate memory then return
/// NULL if INTERACTIVE is defined, or print error to stderr and exit
/// otherwise.
char** New_QV_Buffer (DAZZ_DB* db);

/// Load into 'entry' the 5 QV vectors for i'th read in 'db'.  The deletion
/// tag or characters are converted to a numeric or upper/lower case ascii
/// string as per ascii.  Return with a zero, except when an error occurs and
/// INTERACTIVE is defined in which case return wtih 1.
int Load_QVentry (DAZZ_DB* db, int i, char** entry, int ascii);


/// Remove the QV pseudo track, all space associated with it, and close the
/// .qvs file.
void Close_QVs (DAZZ_DB* db);
