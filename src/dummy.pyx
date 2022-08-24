
cdef extern from "dummy.cpp":
    void welp()

def pywelp():
    welp()