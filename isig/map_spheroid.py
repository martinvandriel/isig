#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Spheroidal mapping.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
import matplotlib.pyplot as plt
import numpy as np


def map_spheroid(gll_points, nodes):
    """
    nodes(nelem, 4, 2)
    """
    r = np.sqrt((nodes ** 2).sum(axis=-1))
    theta = np.arctan2(nodes[...,0], nodes[...,1])

    nelem = nodes.shape[0]
    npol = gll_points.shape[0]
    points_x = np.zeros((nelem, npol, npol))
    points_y = np.zeros((nelem, npol, npol))

    eta, xi = np.meshgrid(gll_points, gll_points, indexing='ij')

    for i, (_r, _t) in enumerate(zip(r, theta)):
        points_x[i] = \
           (1 + eta) * _r[...,3] / 2 * np.sin(((1 - xi) * _t[...,3] + (1 + xi) * _t[...,2]) / 2)\
         + (1 - eta) * _r[...,0] / 2 * np.sin(((1 - xi) * _t[...,0] + (1 + xi) * _t[...,1]) / 2)

        points_y[i] = \
           (1 + eta) * _r[...,3] / 2 * np.cos(((1 - xi) * _t[...,3] + (1 + xi) * _t[...,2]) / 2) \
         + (1 - eta) * _r[...,0] / 2 * np.cos(((1 - xi) * _t[...,0] + (1 + xi) * _t[...,1]) / 2)

    return points_x, points_y
