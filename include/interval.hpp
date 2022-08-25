#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <vector>
#include <fstream>
/*
 * Here are 3 helper structs for performing global adaptive integration.
 * Each Mesh struct is a priority queue implemented with as a simple
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
        
        template <typename T>
        struct MeshBase 
        {
            MeshBase(int size = 4, int N = 0) : totalSize(size), N(N), heap(std::vector<T>(size)) {};

            // ~MeshBase<T>();
            void insert(T &interval);
            void extract(T &worst);
            double totalIntergral();
            double totalError();
            void write(char** buf);
            int meshCheck(); 
            void heapifyUp();
            void heapifyDown();
            void meshWrite(char **buf);
            size_t totalSize;
            size_t N;
            std::vector<T> heap;
        };

        struct Interval
        {
            double a;
            double b;
            double I;
            double err;
        };


        struct Mesh : MeshBase<Interval>
        {
            using MeshBase<Interval>::MeshBase;
        };


        struct Interval3
        {
            double a;
            double b;
            double I;
            double err;
            double fa;
            double fb;
            double fm;
        };


        struct Mesh3 : MeshBase<Interval3>
        {
            using MeshBase<Interval3>::MeshBase;
        };


        struct Interval5
        {
            double a;
            double b;
            double I;
            double err;
            double fa;
            double fb;
            double fl;
            double fm;
            double fr;
        };


        struct Mesh5 : MeshBase<Interval5>
        {
            using MeshBase<Interval5>::MeshBase;
        };


        struct Interval9
        {
            double a;
            double b;
            double I;
            double err;
            double fa;
            double fll;
            double fl;
            double flr;
            double fm;
            double frl;
            double fr;
            double frr;
            double fb;
            int refinement;
        };

        struct Mesh9 : MeshBase<Interval9>
        {
            using MeshBase<Interval9>::MeshBase;
            void intervalWrite(Interval9 &i, FILE *stream);
        };
    } // namespace mesh
    
} // namespace afterglowpy

#include "interval.tpp"
#endif
