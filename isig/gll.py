#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Create GLL interpolation points

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
import matplotlib.pyplot as plt
import numpy as np


def gauss_lobatto_legendre_quadruature_points_weights_fast(n=5, eps=1e-15):
    """
    compute Gauss-Lobatto-Legendre (GLL) quadrature weights and points [-1, 1]
    using Newton-Raphson optimization.

    :param n: number of integration points (order + 1)
    :type n: integer
    :param eps: accuracy of the nodes
    :type eps: float

    :returns: tuple of two numpy arrays of floats containing the points and
        weights
    """

    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess

    x = np.cos(np.pi * np.arange(n) / (n - 1))

    P = np.zeros((n, n))

    xold = 0

    while np.max(np.abs(x - xold)) > eps:

        xold = x

        # Compute P_(N) using the recursion relation
        P[:,0] = 1
        P[:,1] = x

        for k in np.arange(2, n):
            P[:, k] = ((2 * k - 1) * x * P[:,k-1] - (k - 1) * P[:, k-2]) / k


        # Compute its first and second derivatives and
        # update x using the Newton-Raphson method.
        x = xold - (x * P[:,n-1] - P[:,n-2]) / (n * P[:,n-1])

    w = 2. / ((n - 1) * n * P[:, n-1] ** 2)

    return x[::-1], w


#def get_gll_points(self, gll_points=None, glj_points=None):
#    """
#    gll_points and glj_points are the quadrature points in reference
#    coordinates in 1D
#    not optimized and only meant for debugging/plotting purposes
#    """
#    if gll_points is None or glj_points is None:
#        gll_points = self.gll_points.copy()
#        glj_points = self.glj_points.copy()
#
#    gll_xi, gll_eta = np.meshgrid(gll_points, gll_points)
#    glj_xi, glj_eta = np.meshgrid(glj_points, glj_points)
#    gll_xi_map = {False: gll_xi, True: glj_xi}
#
#    gll_local = np.zeros((2, self.nelem, self.npol + 1, self.npol + 1))
#
#    # can't vectorize the loop, because we rely on the fortran library for
#    # mapping
#    for ielem in np.arange(self.nelem):
#        nodes = self.points[self.connectivity[ielem]]
#        for jpol in np.arange(self.npol + 1):
#            for ipol in np.arange(self.npol + 1):
#                gll_local[:, ielem, jpol, ipol] = \
#                    finite_elem_mapping.mapping(
#                        gll_xi_map[self.is_axial[ielem]][jpol, ipol],
#                        gll_eta[jpol, ipol], nodes)
