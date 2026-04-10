# distutils: language = c++

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "sh2dwave.h":
    pair[vector[vector[vector[double]]], vector[vector[double]]] wave_propagate(
        int nx, int nz, int nt, double dx, double dz, double dt,
        double f0, double t0, int xsrc, int zsrc,
        const vector[vector[double]]& vs,
        const vector[vector[double]]& rho,
        int use_absorb, int w, double a,
        int accuracy)

    pair[vector[double], vector[double]] SH_SEIS(
        int nx, int nz, int nt, double dx, double dz, double dt,
        double f0, double t0, int isrc, int jsrc, int ir, int jr,
        const vector[vector[double]]& vs,
        const vector[vector[double]]& rho,
        int use_absorb, int w, double a,
        int accuracy)


def py_wave_propagate(int nx, int nz, int nt, double dx, double dz, double dt,
                      double f0, double t0, int xsrc, int zsrc,
                      np.ndarray[np.double_t, ndim=2] vs,
                      np.ndarray[np.double_t, ndim=2] rho,
                      int use_absorb, int w, double a,
                      int accuracy):

    cdef vector[vector[double]] c_vs
    cdef vector[vector[double]] c_rho
    cdef vector[double] row_vs, row_rho
    cdef int i, j

    for i in range(nx):
        row_vs = vector[double]()
        row_rho = vector[double]()
        for j in range(nz):
            row_vs.push_back(vs[i, j])
            row_rho.push_back(rho[i, j])
        c_vs.push_back(row_vs)
        c_rho.push_back(row_rho)

    cdef pair[vector[vector[vector[double]]], vector[vector[double]]] result

    result = wave_propagate(
        nx, nz, nt, dx, dz, dt,
        f0, t0, xsrc, zsrc,
        c_vs, c_rho,
        use_absorb, w, a,
        accuracy
    )

    return np.array(result.first, dtype=np.float64), np.array(result.second, dtype=np.float64)


def py_sh_seis(int nx, int nz, int nt, double dx, double dz, double dt,
               double f0, double t0, int isrc, int jsrc, int ir, int jr,
               np.ndarray[np.double_t, ndim=2] vs,
               np.ndarray[np.double_t, ndim=2] rho,
               int use_absorb, int w, double a,
               int accuracy):

    cdef vector[vector[double]] c_vs
    cdef vector[vector[double]] c_rho
    cdef vector[double] row_vs, row_rho
    cdef int i, j

    for i in range(nx):
        row_vs = vector[double]()
        row_rho = vector[double]()
        for j in range(nz):
            row_vs.push_back(vs[i, j])
            row_rho.push_back(rho[i, j])
        c_vs.push_back(row_vs)
        c_rho.push_back(row_rho)

    cdef pair[vector[double], vector[double]] result

    result = SH_SEIS(
        nx, nz, nt, dx, dz, dt,
        f0, t0, isrc, jsrc, ir, jr,
        c_vs, c_rho,
        use_absorb, w, a,
        accuracy
    )

    return np.array(result.first, dtype=np.float64), np.array(result.second, dtype=np.float64)