#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Tests for spheroidal mapping.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
import numpy as np
from ..map_spheroid import map_spheroid
from pymesher import Skeleton


def test_map_spheroid():
    discontinuities = np.array([0., .9, 1.])
    nlayer = len(discontinuities) - 1
    hmax = np.ones(nlayer) * 0.5

    sk = Skeleton.create_spherical_nonconforming_mesh(
        discontinuities, hmax, refinement_factor=2, ndim=2, min_colat=0.,
        max_colat=45.)

    m = sk.get_unstructured_mesh()
    m.plot(show=False)

    gll_points = np.array([-1., 0., 1.])
    gll_x, gll_y = map_spheroid(gll_points, m.points[m.connectivity])

    px = np.array([0., 0., 0., 0.17558129, 0.18533581, 0.19509032, 0.34441509,
                   0.36354926, 0.38268343, 0.34441509, 0.36354926, 0.38268343,
                   0.50001321, 0.52779172, 0.55557023, 0.6363961, 0.67175144,
                   0.70710678])
    py = np.array([0.9, 0.95, 1., 0.88270675, 0.93174602, 0.98078528,
                   0.83149158, 0.87768556, 0.92387953, 0.83149158, 0.87768556,
                   0.92387953, 0.74832265, 0.78989613, 0.83146961, 0.6363961,
                   0.67175144, 0.70710678])

    np.testing.assert_allclose(px, gll_x.flatten(), atol=1e-15)
    np.testing.assert_allclose(py, gll_y.flatten(), atol=1e-15)
