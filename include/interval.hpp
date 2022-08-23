#ifndef AFTERGLOWPY_INTERVAL_HPP
#define AFTERGLOWPY_INTERVAL_HPP

#include <vector>
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
        
        struct Interval
        {
            double a;
            double b;
            double I;
            double err;
        };


        struct Mesh
        {
            Mesh();
            ~Mesh();
            
            void insert(Interval &interval);
            void extract(Interval &worst);
            double totalIntergral();
            double totalError();
            void write(char** buf);
            int meshCheck(); 
            size_t totalSize;
            size_t N;
            std::vector<Interval> heap;
            // struct Interval *heap;
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


        struct Mesh3
        {
            Mesh3::Mesh3();
            Mesh3::~Mesh3();
            
            void insert(Interval &interval);
            void extract(Interval &worst);
            double totalIntergral();
            double totalError();
            void write(char** buf);
            int meshCheck(); 
            size_t totalSize;
            size_t N;
            std::vector<Interval3> heap;
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


        struct Mesh5
        {
            Mesh5::Mesh5();
            Mesh5::~Mesh5();
            
            void insert(Interval &interval);
            void extract(Interval &worst);
            double totalIntergral();
            double totalError();
            void write(char** buf);
            int meshCheck(); 
            size_t totalSize;
            size_t N;
            struct Interval5 *heap;
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

        struct Mesh9
        {
            Mesh9::Mesh9();
            Mesh9::~Mesh9();
            
            void insert(Interval i&nterval);
            void extract(Interval &worst);
            double totalIntergral();
            double totalError();
            void write(char** buf);
            int meshCheck(); 
            size_t totalSize;
            size_t N;
            struct Interval9 *heap;
        };


        void meshInit(std::vector<Mesh> m);
        void meshFree(std::vector<Mesh> m);
        void meshInsert(std::vector<Mesh> m, std::vector<Interval> i);
        void meshExtract(std::vector<Mesh> m, std::vector<Interval> worst);
        double meshTotalIntegral(std::vector<Mesh> m);
        double meshTotalError(std::vector<Mesh> m);
        void meshHeapifyUp(std::vector<Mesh> m);
        void meshHeapifyDown(std::vector<Mesh> m);
        int meshCheck(std::vector<Mesh> m);
        void meshWrite(std::vector<Mesh> m, char **buf);

        void mesh3Init(std::vector<Mesh3> m);
        void mesh3Free(std::vector<Mesh3> m);
        void mesh3Insert(std::vector<Mesh3> m, std::vector<Interval3> i);
        void mesh3Extract(std::vector<Mesh3> m, std::vector<Interval3> worst);
        double mesh3TotalIntegral(std::vector<Mesh3> m);
        double mesh3TotalError(std::vector<Mesh3> m);
        void mesh3HeapifyUp(std::vector<Mesh3> m);
        void mesh3HeapifyDown(std::vector<Mesh3> m);
        int mesh3Check(std::vector<Mesh3> m);
        void mesh3Write(std::vector<Mesh3> m, char **buf);

        void mesh5Init(std::vector<Mesh5> m);
        void mesh5Free(std::vector<Mesh5> m);
        void mesh5Insert(std::vector<Mesh5> m, std::vector<Interval5> i);
        void mesh5Extract(std::vector<Mesh5> m, std::vector<Interval5> worst);
        double mesh5TotalIntegral(std::vector<Mesh5> m);
        double mesh5TotalError(std::vector<Mesh5> m);
        void mesh5HeapifyUp(std::vector<Mesh5> m);
        void mesh5HeapifyDown(std::vector<Mesh5> m);
        int mesh5Check(std::vector<Mesh5> m);
        void mesh5Write(std::vector<Mesh5> m, char **buf);

        void mesh9Init(std::vector<Mesh9> m);
        void mesh9Free(std::vector<Mesh9> m);
        void mesh9Insert(std::vector<Mesh9> m, std::vector<Interval9> i);
        void mesh9Extract(std::vector<Mesh9> m, std::vector<Interval9> worst);
        double mesh9TotalIntegral(std::vector<Mesh9> m);
        double mesh9TotalError(std::vector<Mesh9> m);
        void mesh9HeapifyUp(std::vector<Mesh9> m);
        void mesh9HeapifyDown(std::vector<Mesh9> m);
        int mesh9Check(std::vector<Mesh9> m);
        void mesh9Write(std::vector<Mesh9> m, char **buf);

        void interval9Write(std::vector<Interval9> i, FILE *stream);
    } // namespace mesh
    
} // namespace afterglowpy

#endif
