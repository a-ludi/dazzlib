/**
    Definitions about local alignments specific to the Marvel-suite.


    Copyright: © 2021 Martin Pippel <pippel@mpi-cbg.de>,
               Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>,
             Martin Pippel <pippel@mpi-cbg.de>
*/
module dazzlib.core.marvel.align_;

import dazzlib.core.c.DB;


enum MarvelLocalAlignmentFlag : uint32
{
    init = 0x0,

    /// B on reverse strand
    comp = 1 << 0,

    /// Overlap tagged for removal
    disabled = 1 << 1,

    /// repeat induced overlap
    repeat = 1 << 3,

    /// local alignment
    local = 1 << 4,

    /// too many differences
    diff = 1 << 5,

    /// stitched to another overlap
    stitch = 1 << 7,

    ///
    symdiscard = 1 << 8,

    /// overlap length
    olen = 1 << 9,

    /// read length
    rlen = 1 << 10,

    /// temporary flag, not written to disk
    temp = 1 << 11,

    /// containment
    cont = 1 << 15,

    /// gap
    gap = 1 << 16,

    /// trimmed
    trim = 1 << 17,

    /// overlap spans unique repeat modules junction
    module_ = 1 << 18,

    /// optional (risky) overlaps for touring
    optional = 1 << 19,
}
