// #include "include/interval.hpp"
namespace afterglowpy
{
    namespace mesh
    {
        // template<typename T>
        // MeshBase<T>::MeshBase(int size = 4, int N = 0) : totalSize(size), N(N), heap(std::vector<T>(size)) {};

        template <typename T>
        void MeshBase<T>::insert(T &interval)
        {
            while(N >= totalSize) {
                T dummy;
                totalSize *= 2;
                heap.push_back(dummy);
            }
            N++;

            //Restore ordering
            heapifyUp();
        }

        template <typename T>
        void MeshBase<T>::extract(T &worst)
        {
            worst = heap[0];
            heap[0] = heap[N - 1];
            N--;

            heapifyDown();
        }

        template<typename T>
        double MeshBase<T>::totalIntergral()
        {
            double I = 0.0;
            size_t i = 0;
            for(size_t i = 0; i < N; i++)
                I += heap[i].I;
            return I;
        }

        template<typename T>
        double MeshBase<T>::totalError()
        {
            double err = 0.0;
            for(size_t i = 0; i < N; i++)
                err += heap[i].err;
            return err;
        }

        template<typename T>
        int MeshBase<T>::meshCheck()
        {
            for(size_t p = 0; p <= (N - 2) / 2; p++)
            {
                size_t c1 = 2 * p+1;
                size_t c2 = c1 + 1;

                if(c1 < N && heap[c1].err > heap[p].err)
                    return 0;
                
                if(c2 < N && heap[c2].err > heap[p].err)
                    return 0;
            }
            return 1;
        }

        template<typename T>
        void MeshBase<T>::heapifyUp()
        {
            size_t c =  N - 1;
            size_t p = (c - 1) / 2;

            while(c != 0 && heap[p].err < heap[c].err)
            {
                T tempP = heap[p];
                heap[p] = heap[c];
                heap[c] = tempP;
                c = p;
                p = (c - 1) / 2;
            }
        }

        template<typename T>
        void MeshBase<T>::heapifyDown()
        {
            size_t p  = 0;
            size_t c1 = 2 * p + 1;
            size_t c2 = c1 + 1;


            while(c1 < N)
            {
                //Find the child with largest error
                size_t c = c1;
                double e = heap[c1].err;

                if(c2 < N && heap[c2].err > e)
                {
                    c = c2;
                    e = heap[c2].err;
                }

                // If the child is already in order then we're done.
                if(e <= heap[p].err)
                    break;

                // Otherwise, swap with the child.
                T tempP = heap[p];
                heap[p] = heap[c];
                heap[c] = tempP;
                p = c;
                c1 = 2*p + 1;
                c2 = c1 + 1;
            }
        }

        template<typename T>
        void MeshBase<T>::meshWrite(char **buf)
        {
            *buf = (char *) malloc((N * 4 * 30 + 12) * sizeof(char));

            size_t i;
            int c = sprintf(*buf, "%lu", N);
            for(i=0; i < N; i++)
            {
                T *in = &(heap[i]);
                c += sprintf(*buf + c, " %.16e %.16e %.16e %.16e",
                            in->a, in->b, in->I, in->err);
            }
            *buf = (char *) realloc(*buf, (c+1) * sizeof(char));
        }
    } // namespace mesh
    
} // namespace afterglowpy
