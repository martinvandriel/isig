#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Generate a Instaseis input for Sismosphere

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
'''
import argparse
import numpy as np
from pymesher.models_1D import model
from pymesher.skeleton import Skeleton
from pymesher.global_numbering import get_global_lexi
import sys

from .map_spheroid import map_spheroid
from .gll import gauss_lobatto_legendre_quadruature_points_weights_fast as get_gll


if __name__ == "__main__":

    # handle negative digits in the arguments which could be misinterpreted as
    # options otherwise. See http://stackoverflow.com/a/21446783
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    parser = argparse.ArgumentParser(
        prog="python -m isig",
        description='Generate a Instaseis input for Sismosphere.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-m', '--model_file', type=str, required=True,
        help='path to 1D model in deck file format (or any other format '
             'compatible with salvus mesher)')

    parser.add_argument(
        '-p', '--period', type=float, default=50.,
        help='Shortest period to resolve.')

    parser.add_argument(
        '-d', '--max_depth', type=float, default=100.,
        help='Maximum source depth in km.')

    parser.add_argument(
        '--min_dist', type=float, default=0.,
        help='Minimum epicentral distance in degrees.')

    parser.add_argument(
        '--max_dist', type=float, default=180.,
        help='Maximum epicentral distance in degrees.')

    parser.add_argument(
        '-e', '--elements_per_wavelength', type=float, default=2.,
        help='Number of Elements per Wavelength.')

    parser.add_argument(
        '-n', '--npol', type=int, default=4,
        help='Polynomial order used for interpolation.')

    args = parser.parse_args()

    # input to SI and consistency checks
    max_depth = 1e3 * args.max_depth
    if args.min_dist < 0:
        raise ValueError('min_dist < 0')
    if args.max_dist < 0:
        raise ValueError('max_dist < 0')
    if args.max_depth< 0:
        raise ValueError('depth < 0')

    mod = model.read(args.model_file)

    discontinuities = mod.discontinuities
    hmax = mod.get_edgelengths(
        dominant_period=args.period,
        elements_per_wavelength=args.elements_per_wavelength)

    # adapt discontinuities and hmax for min_radius
    idx = discontinuities > 1. - max_depth / mod.scale
    ndisc = idx.sum() + 2

    discontinuities_new = np.zeros(ndisc)
    discontinuities_new[1] = 1. - max_depth / mod.scale
    discontinuities_new[-ndisc+2:] = discontinuities[idx]
    discontinuities = discontinuities_new

    hmax_new = np.ones(ndisc-1)
    hmax_new[-ndisc+2:] = hmax[-ndisc+2:]
    hmax = hmax_new

    sk = Skeleton.create_spherical_nonconforming_mesh(
        discontinuities, hmax, 2,
        refinement_factor=2,
        max_colat=args.max_dist,
        min_colat=args.min_dist,
        full_sphere=False)

    m = sk.get_unstructured_mesh()

    gll = get_gll(args.npol)[0]
    gll_x, gll_y = map_spheroid(gll, m.points[m.connectivity])

    r = np.sqrt(gll_x.ravel() ** 2 +  gll_y.ravel() ** 2)
    theta = np.arctan2(gll_x.ravel(), gll_y.ravel())

    theta_idx, ntheta = get_global_lexi(theta[np.newaxis,:])

    theta_unique = theta[np.unique(theta_idx, return_index=True)[1]]
    theta_unique

    for i in np.arange(ntheta):
        print np.rad2deg(theta_unique[i])
        print (1. - r[theta_idx==i]) * mod.scale / 1e3

    import matplotlib.pyplot as plt
    plt.scatter(gll_x, gll_y)
    plt.show()
