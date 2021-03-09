/**
    Local alignment module.  Routines for finding local alignments given a
    seed position, representing such an l.a. with its interval and a set of
    pass-thru points, so that a detailed alignment can be efficiently computed
    on demand.

    All routines work on a numeric representation of DNA sequences, i.e. 0 for
    A, 1 for C, 2 for G, and 3 for T.


    ### INTERACTIVE vs BATCH version

    The defined constant INTERACTIVE (set in DB.h) determines whether an
    interactive or batch version of the routines in this library are compiled.
    In batch mode, routines print an error message and exit.  In interactive
    mode, the routines place the error message in EPLACE (also defined in
    DB.h) and return an error value, typically NULL if the routine returns a
    pointer, and an unusual integer value if the routine returns an integer.

    Below when an error return is described, one should understand that this
    value is returned only if the routine was compiled in INTERACTIVE mode.


    Copyright: © 2013 Dr. Eugene W. Myers <gene.myers@gmail.com>. All rights
               reserved.
               © 2021 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>,
             Dr. Eugene W. Myers <gene.myers@gmail.com>
*/
module dazzlib.core.c.align_;

import dazzlib.core.c.DB;
import core.stdc.stdio;

extern (C):

/// If the trace spacing is not more than this value, then can and do compress
/// traces pts to 8-bit unsigned ints
enum TRACE_XOVR = 125;


/**
    Path abstraction.

    Coordinates are *between* characters where 0 is the tick just before the
    first char, 1 is the tick between the first and second character, and so
    on.  Our data structure is called a Path refering to its conceptualization
    in an edit graph.

    A local alignment is specified by the point '(abpos,bbpos)' at which its
    path in the underlying edit graph starts, and the point '(aepos,bepos)' at
    which it ends. In otherwords A[abpos+1..aepos] is aligned to
    B[bbpos+1..bepos] (assuming X[1] is the *first* character of X).

    There are 'diffs' differences in an optimal local alignment between the
    beginning and end points of the alignment (if computed by Compute_Trace),
    or nearly so (if computed by Local_Alignment).

    Optionally, a Path can have additional information about the exact nature
    of the aligned substrings if the field 'trace' is not NULL.  Trace points
    to either an array of integers (if computed by a Compute_Trace routine),
    or an array of unsigned short integers (if computed by Local_Alignment).

    If computed by Local_Alignment 'trace' points at a list of 'tlen' (always
    even) short values:

        d_0, b_0, d_1, b_1, ... d_n-1, b_n-1, d_n, b_n

    to be interpreted as follows.  The alignment from (abpos,bbpos) to
    (aepos,bepos) passes through the n trace points for i in [1,n]:

        (a_i,b_i) where a_i = floor(abpos/TS)*TS + i*TS
                    and b_i = bbpos + (b_0 + b_1 + b_i-1)

    where also let a_0,b_0 = abpos,bbpos and a_(n+1),b_(n+1) = aepos,bepos.
    That is, the interior (i.e. i != 0 and i != n+1) trace points pass through
    every TS'th position of the aread where TS is the "trace spacing" employed
    when finding the alignment (see New_Align_Spec).  Typically TS is 100.
    Then d_i is the number of differences in the portion of the alignment
    between (a_i,b_i) and (a_i+1,b_i+1).  These trace points allow the
    Compute_Trace routines to efficiently compute the exact alignment between
    the two reads by efficiently computing exact alignments between
    consecutive pairs of trace points. Moreover, the diff values give one an
    idea of the quality of the alignment along every segment of TS symbols of
    the aread.

    If computed by a Compute_Trace routine, 'trace' points at a list of 'tlen'
    integers < i1, i2, ... in > that encodes an exact alignment as follows.  A
    negative number j indicates that a dash should be placed before A[-j] and
    a positive number k indicates that a dash should be placed before B[k],
    where A and B are the two sequences of the overlap.  The indels occur in
    the trace in the order in which they occur along the alignment.  For a
    good example of how to "decode" a trace into an alignment, see the code
    for the routine Print_Alignment.
*/
struct Path
{
    /// Trace vector
    void* trace;
    /// Number of integers in trace vector
    int tlen;
    /// Total number of diffs in trace
    int diffs;
    /// A-read begin position
    int abpos;
    /// B-read begin position
    int bbpos;
    /// A-read end position
    int aepos;
    /// B-read end position
    int bepos;
}

enum LocalAlignmentFlag : int
{
    /// B-sequence is complemented
    complement = 0x1,
    /// A-sequence is complemented, not B !  Only Local_Alignment notices
    aComplement = 0x2,
    /// LA is the first of a chain of 1 or more la's
    start = 0x4,
    /// LA is the next segment of a chain.
    next = 0x8,
    /// This is the start of the best chain; must be combined with start.
    best = 0x10,
    /// This LA should be ignored
    disabled = 0x20,
    /// (dazzlib extension) This LA is not part of any chain.
    unchained = 0x40,
}


/**
    Alignment abstraction.

    An alignment is modeled by an Alignment record, which in addition to a
    *pointer* to a 'path', gives pointers to the A and B sequences, their
    lengths, and indicates whether the B-sequence needs to be complemented
    ('comp' non-zero if so).  The 'trace' pointer of the 'path' subrecord can
    be either NULL, a list of pass-through points, or an exact trace depending
    on what routines have been called on the record.

    One can (1) compute a trace, with Compute_Trace, either from scratch if
    'path.trace' = NULL, or using the sequence of pass-through points in
    trace, (2) print an ASCII representation of an alignment, or (3) reverse
    the roles of A and B, and (4) complement a sequence (which is a reversible
    process).

    If the alignment record shows the B sequence as complemented, *** THEN IT
    IS THE RESPONSIBILITY OF THE CALLER *** to make sure that bseq points at a
    complement of the sequence before calling Compute_Trace or
    Print_Alignment.  Complement_Seq complements the sequence a of length n.
    The operation does the complementation/reversal in place. Calling it a
    second time on a given fragment restores it to its original state.

    With the introduction of the DAMAPPER, we need to code chains of
    alignments between a pair of sequences.  The alignments of a chain are
    expected to be found in order either on a file or in memory, where the
    START_FLAG marks the first alignment and the NEXT_FLAG all subsequent
    alignmenst in a chain.  A chain of a single LA is marked with the
    START_FLAG. The BEST_FLAG marks one of the best chains for a pair of
    sequences.  The convention is that either every record has either a START-
    or NEXT-flag, or none of them do (e.g. as produced by daligner), so one
    can always check the flags of the first alignment to see whether or not
    the chain concept applies to a given collection or not.
*/
struct Alignment
{
    /// The alignment path/details.
    Path* path;
    /// Pipeline status and complementation flags
    LocalAlignmentFlag flags;
    /// Pointer to A sequence
    char* aseq;
    /// Pointer to B sequence
    char* bseq;
    /// Length of A sequence
    int alen;
    /// Length of B sequence
    int blen;
}

/// Complement the sequence in fragment aseq.  The operation does the
/// complementation/reversal in place.  Calling it a second time on a given
/// fragment restores it to its original state.
void Complement_Seq (char* a, int n);
/// ditto
extern (D) void Complement_Seq (ref char[] a)
{
    Complement_Seq(a.ptr, cast(int) a.length);
}

unittest
{
    enum sequence = "attacg";
    char[] seqStore = sequence.dup;

    Complement_Seq(seqStore.ptr, cast(int) seqStore.length);
    Complement_Seq(seqStore.ptr, cast(int) seqStore.length);

    assert(seqStore == sequence);
}


/// Many routines like Local_Alignment, Compute_Trace, and Print_Alignment
/// need working storage that is more efficiently reused with each call,
/// rather than being allocated anew with each call.  Each *thread* can create
/// a Work_Data object with New_Work_Data and this object holds and retains
/// the working storage for routines of this module between calls to the
/// routines.  If enough memory for a Work_Data is not available then NULL is
/// returned. Free_Work_Data frees a Work_Data object and all working storage
/// held by it.
alias Work_Data = void;

/// ditto
Work_Data* New_Work_Data ();

/// ditto
void Free_Work_Data (Work_Data* work);

/**
    Local_Alignment seeks local alignments of a quality determined by a number
    of parameters. These are coded in an Align_Spec object that can be created
    with New_Align_Spec and freed with Free_Align_Spec when no longer needed.

    Params:

        ave_corr =  the average correlatio n(1 - 2*error_rate) for the sought
                    alignments.  For Pacbio data we set this to .70 assuming
                    an average of 15% error in each read.
        trace_space = the spacing interval for keeping trace points and
                    segment differences (see description of 'trace' for Paths
                    above)
        freq =      a 4-element vector where afreq[0] = frequency of A, f(A),
                    freq[1] = f(C), freq[2] = f(G), and freq[3] = f(T). This
                    vector is part of the header of every DAZZ database
                    (see db.d).
        reach =     a boolean, if set alignment extend to the boundary when
                    reasonable, otherwise the terminate only at
                    suffix-positive points.

    If an alignment cannot reach the boundary of the d.p. matrix with this
    condition (i.e. overlap), then the last/first 30 columns of the alignment
    are guaranteed to be suffix/prefix positive at correlation ave_corr *
    g(freq) where g is an empirically measured function that increases from 1
    as the entropy of freq decreases.  If memory is unavailable or the freq
    distribution is too skewed then NULL is returned.

    You can get back the original parameters used to create an Align_Spec with
    the simple utility functions below.
*/
alias Align_Spec = void;

/// ditto
Align_Spec* New_Align_Spec (double ave_corr, int trace_space, float* freq, int reach);

/// ditto
void Free_Align_Spec (Align_Spec* spec);

/// ditto
int Trace_Spacing (Align_Spec* spec);

/// ditto
double Average_Correlation (Align_Spec* spec);

/// ditto
float* Base_Frequencies (Align_Spec* spec);

/// ditto
int Overlap_If_Possible (Align_Spec* spec);


/**
    Local_Alignment finds the longest significant local alignment between the
    sequences in 'align' subject to:

    a) the alignment criterion given by the Align_Spec 'spec',
    b) it passes through one of the points (anti+k)/2,(anti-k)/2 for k in
       [low,hgh] within the underlying dynamic programming matrix (i.e. the
       points on diagonals low to hgh on anti-diagonal anti or anti-1
       (depending on whether the diagonal is odd or even)),
    c) if lbord >= 0, then the alignment is always above diagonal low-lbord, and
    d) if hbord >= 0, then the alignment is always below diagonal hgh+hbord.

    The path record of 'align' has its 'trace' filled from the point of view
    of an overlap between the aread and the bread.  In addition a Path record
    from the point of view of the bread versus the aread is returned by the
    function, with this Path's 'trace' filled in appropriately.  The space for
    the returned path and the two 'trace's are in the working storage supplied
    by the Work_Data packet and this space is reused with each call, so if one
    wants to retain the bread-path and the two trace point sequences, then
    they must be copied to user-allocated storage before calling the routine
    again.  NULL is returned in the event of an error.

    Find_Extension is a variant of Local_Alignment that simply finds a local
    alignment that either ends (if prefix is non-zero) or begins (if prefix is
    zero) at the point (anti+diag)/2,(anti-diag)/2).  All other parameters are
    as before.  It returns a non-zero value only when INTERACTIVE is on and it
    cannot allocate the memory it needs. Only the path and trace with respect
    to the aread is returned.  This routine is experimental and may not
    persist in later versions of the code.
*/
Path* Local_Alignment (
    Alignment* align_,
    Work_Data* work,
    Align_Spec* spec,
    int low,
    int hgh,
    int anti,
    int lbord,
    int hbord);

/// ditto
int Find_Extension (
    Alignment* align_,
    Work_Data* work,
    Align_Spec* spec,
    int diag,
    int anti,
    int lbord,
    int hbord,
    int prefix);


/// Possible modes for "mode" parameter below)
enum ComputeTraceMode : int
{
    ///
    lowermost = -1,
    ///
    greediest = 0,
    ///
    uppermost = 1,
}


/**
    Given a legitimate Alignment object and associated trace point vector in
    'align->path.trace', Compute_Trace_X, computes an exact trace for the
    alignment and resets 'align->path.trace' to point at an integer array
    within the storage of the Work_Data packet encoding an exact optimal trace
    from the start to end points.  If the trace is needed beyond the next call
    to a routine that sets it, then it should be copied to an array allocated
    and managed by the caller.

    Compute_Trace_PTS computes a trace by computing the trace between
    successive trace points. It is much, much faster than Compute_Alignment
    below but at the tradeoff of not necessarily being optimal as pass-through
    points are not all perfect.  Compute_Trace_MID computes a trace by
    computing the trace between the mid-points of alignments between two
    adjacent pairs of trace points.  It is generally twice as slow as
    Compute_Trace_PTS, but it produces nearer optimal alignments.  Both these
    routines return 1 if an error occurred and 0 otherwise.
*/
int Compute_Trace_PTS (Alignment* align_, Work_Data* work, int trace_spacing, ComputeTraceMode mode);

/// ditto
int Compute_Trace_MID (Alignment* align_, Work_Data* work, int trace_spacing, ComputeTraceMode mode);


/**
    Compute_Trace_IRR (IRR for IRRegular) computes a trace for the given
    alignment where it assumes the spacing between trace points between both
    the A and B read varies, and futher assumes that the A-spacing is given in
    the short integers normally occupied by the differences in the alignment
    between the trace points.  This routine is experimental and may not
    persist in later versions of the code.
*/
int Compute_Trace_IRR (Alignment* align_, Work_Data* work, ComputeTraceMode mode);


/// Determine task of Compute_Alignment.
enum ComputeAlignmentTask : int
{
    /// 'path.trace' is set to point to an optimal trace of an optimal
    /// alignment. This must be preceeded by a call with diffOnly.
    plusAlign = 0,
    /// 'path.trace' is set to point at an integer array within the storage of
    /// the Work_Data packet encoding a trace point sequence for an optimal
    /// alignment. This must be preceeded by a call with diffOnly.
    plusTrace = 1,
    /// Only the difference of this alignment is computed and placed in the
    /// "diffs" field of align's path
    diffOnly = 2,
    /// Like plusAlign but implicitly computes diffOnly
    diffAlign = 3,
    /// Like plusTrace but implicitly computes diffOnly
    diffTrace = 4,
}


/**
    Compute Alignment determines the best alignment between the substrings
    specified by align. If the task is DIFF_ONLY, then only the difference of
    this alignment is computed and placed in the "diffs" field of align's
    path.  If the task is PLUS_TRACE or DIFF_TRACE, then 'path.trace' is set
    to point at an integer array within the storage of the Work_Data packet
    encoding a trace point sequence for an optimal alignment, whereas if the
    task is PLUS_ALIGN or DIFF_ALIGN, then it points to an optimal trace of an
    optimal alignment.  The PLUS tasks can only be called if the immmediately
    proceeding call was a DIFF_ONLY on the same alignment record and
    sequences, in which case a little efficiency is gained by avoiding the
    repetition of the top level search for an optimal mid-point.
*/
int Compute_Alignment (Alignment* align_, Work_Data* work, ComputeAlignmentTask task, int trace_spacing);


/**
    Alignment_Cartoon prints an ASCII representation of the overlap
    relationhip between the two reads of 'align' to the given 'file' indented
    by 'indent' space.  Coord controls the display width of numbers, it must
    be not less than the width of any number to be displayed.

    If the alignment trace is an exact trace, then one can ask Print_Alignment
    to print an ASCII representation of the alignment 'align' to the file
    'file'.  Indent the display by "indent" spaces and put "width" columns per
    line in the display.  Show "border" characters of sequence on each side of
    the aligned region.  If upper is non-zero then display bases in upper
    case.  If coord is greater than 0, then the positions of the first
    character in A and B in the given row is displayed with a field width
    given by coord's value.

    Print_Reference is like Print_Alignment but rather than printing exaclty
    "width" columns per segment, it prints "block" characters of the A
    sequence in each segment.  This results in segments of different lengths,
    but is convenient when looking at two alignments involving A as segments
    are guaranteed to cover the same interval of A in a segment.

    Both Print routines return 1 if an error occurred (not enough memory), and
    0 otherwise.

    Flip_Alignment modifies align so the roles of A and B are reversed.  If
    full is off then the trace is ignored, otherwise the trace must be to a
    full alignment trace and this trace is also appropriately inverted.
*/
void Alignment_Cartoon (FILE* file, Alignment* align_, int indent, int coord);

/// ditto
int Print_Alignment (
    FILE* file,
    Alignment* align_,
    Work_Data* work,
    int indent,
    int width,
    int border,
    int upper,
    int coord);

/// ditto
int Print_Reference (
    FILE* file,
    Alignment* align_,
    Work_Data* work,
    int indent,
    int block,
    int border,
    int upper,
    int coord);

/// ditto
void Flip_Alignment (Alignment* align_, int full);

/**
    Overlap abstraction.

    Externally, between modules an Alignment is modeled by an "Overlap"
    record, which

    a) replaces the pointers to the two sequences with their ID's in the DAZZ
       data bases,
    b) does not contain the length of the 2 sequences (must fetch from DB), and
    c) contains its path as a subrecord rather than as a pointer (indeed,
       typically the corresponding Alignment record points at the Overlap's
       path sub-record).

    The trace pointer is always to a sequence of trace points and can be
    either compressed (uint8) or uncompressed (uint16).  One can read and
    write binary records of an "Overlap".
*/
struct Overlap
{
    /// Path: begin- and end-point of alignment + diffs
    Path path;
    /// Pipeline status and complementation flags
    uint32 flags;
    /// Id # of A sequence (zero-based indexing)
    int aread;
    /// Id # of B sequence (zero-based indexing)
    int bread;
}


/// Reads the next Overlap record from stream 'input', not including the trace
/// (if any), and without modifying 'ovl's trace pointer.
int Read_Overlap (FILE* input, Overlap* ovl);

/// Reads the ensuing trace into the memory pointed at by the trace field of
/// 'ovl'.  It is assumed to be big enough to accommodate the trace where each
/// value take 'tbytes' bytes (1 if uint8 or 2 if uint16).
int Read_Trace (FILE* innput, Overlap* ovl, int tbytes);


/// Writes 'ovl' to stream 'output' followed by its trace vector (if any) that
/// occupies 'tbytes' bytes per value.  It returns non-zero if there was an
/// error writing.
int Write_Overlap (FILE* output, Overlap* ovl, int tbytes);


/// Prints an ASCII version of the contents of 'ovl' to stream 'output' where
/// the trace occupes 'tbytes' per value and the print out is indented from
/// the left margin by 'indent' spaces.
void Print_Overlap (FILE* output, Overlap* ovl, int tbytes, int indent);


/// Compress_TraceTo8 converts a trace fo 16-bit values to 8-bit values in
/// place, and Decompress_TraceTo16 does the reverse conversion.  If check is
/// set in a call to Compress then it checks whether the values fit in 8-bits,
/// and if not returns a non-zero result in interactive mode, or exits with an
/// error message in batch mode.
int Compress_TraceTo8 (Overlap* ovl, int check);

/// ditto
void Decompress_TraceTo16 (Overlap* ovl);


/// Checks that the number of trace points is correct and that the sum of the
/// b-read displacements equals the b-read alignment interval, assuming the
/// trace spacing is 'tspace'.  It reports an error message if there is a
/// problem and 'verbose' is non-zero.  The 'ovl' came from the file names
/// 'fname'.
int Check_Trace_Points (Overlap* ovl, int tspace, int verbose, char* fname);
