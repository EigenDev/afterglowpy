#include <stdio.h>
#include <stdlib.h>
#include "include/interval.hpp"
/*
 * Here are 3 helper structs for performing global adaptive integration.
 * Each Mesh Is a priority queue implemented with as a simple
 * binary heap. Elements may be added, and the element with the worst error can
 * be retrieved.
 *
 * Each mesh is a container of Intervals. An Interval is a primitive struct
 * contains only its left and right bounds, the value of the integral and its
 * error over that interval.  The Interval3 and Interval5 are identical to
 * Interval but contain 3 (or 5) extra slots to contain function memorized
 * function values.
 *
 * The Mesh3 and Mesh5 are priority queues of Interval3s and Interval5s
 * respectively.
 */


namespace afterglowpy
{
    namespace mesh
    {
        void Mesh9::intervalWrite(Interval9 &i, FILE *stream)
        {
            fprintf(stream, "(%.3le, %.3le)  %.12le +/- %.3le   %d\n",
                    i.a, i.b, i.I, i.err, i.refinement);
            fprintf(stream, "   [%.3le %.3le %.3le %.3le %.3le %.3le"
                            " %.3le %.3le %.3le]\n", i.fa, i.fll, i.fl, i.flr,
                            i.fm, i.frl, i.fr, i.frr, i.fb);
        }
    } // namespace mesh
    
} // namespace afterglowpy
