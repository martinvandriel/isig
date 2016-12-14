#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Create GLL interpolation points

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
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
        P[:, 0] = 1
        P[:, 1] = x

        for k in np.arange(2, n):
            P[:, k] = ((2 * k - 1) * x * P[:, k-1] - (k - 1) * P[:, k-2]) / k

        # Compute its first and second derivatives and
        # update x using the Newton-Raphson method.
        x = xold - (x * P[:, n-1] - P[:, n-2]) / (n * P[:, n-1])

    w = 2. / ((n - 1) * n * P[:, n-1] ** 2)

    return x[::-1], w
