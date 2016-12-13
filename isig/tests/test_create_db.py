#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
A new python script.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
import matplotlib.pyplot as plt
import numpy as np
from pymesher import Skeleton, models_1D
import os

from ..create_db import create_db
from ..gll import gauss_lobatto_legendre_quadruature_points_weights_fast as \
    get_gll
from ..map_spheroid import map_spheroid


def test_create_db():

    mod = models_1D.model.built_in('prem_ani')
    mod.plot()
    npol = 5
    discontinuities = mod.discontinuities #np.array([0., .9, 1.])
    discontinuities = discontinuities[[0, -3, -2, -1]]
    nlayer = len(discontinuities) - 1
    hmax = np.ones(nlayer) * 0.5

    sk = Skeleton.create_spherical_nonconforming_mesh(
        discontinuities, hmax, refinement_factor=2, ndim=2, min_colat=0.,
        max_colat=45.)

    m = sk.get_unstructured_mesh()

    gll = get_gll(npol)[0]
    gll_x, gll_y = map_spheroid(gll, m.points[m.connectivity])

    create_db('test.h5', mod, m.points, m.connectivity, gll_x, gll_y)
    os.remove('test.h5')
