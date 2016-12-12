#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for quadrature points and weights.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2015
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lgpl.html)
"""
import numpy as np

from ..gll import gauss_lobatto_legendre_quadruature_points_weights_fast


def test_gauss_lobatto_legendre_quadruature_points_weights_fast():

    # reference data for n = 5
    ref_points = [-1., -21 ** 0.5 / 7., 0., 21 ** 0.5 / 7., 1.]
    ref_weights = [1./10., 49./90., 32./45., 49./90., 1./10.]

    points, weights = gauss_lobatto_legendre_quadruature_points_weights_fast(5)
    np.testing.assert_allclose(points, ref_points, atol=1e-15)
    np.testing.assert_allclose(weights, ref_weights, atol=1e-15)

    n = 100
    p, w = gauss_lobatto_legendre_quadruature_points_weights_fast(n)
    np.testing.assert_allclose(p, -p[::-1], atol=1e-15)
    np.testing.assert_allclose(w, w[::-1], atol=1e-15)
