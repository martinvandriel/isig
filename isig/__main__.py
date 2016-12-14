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
import sys

from .create_db import create_db


DEFAULT_FILE_NAME = 'isig_<modelname>_<period>s_<depth>km'

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

    parser.add_argument('-dt', type=float, default=1., help='Sample rate.')

    parser.add_argument('-npts', type=int, default=3600,
                        help='Number of time samples.')

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
        '-n', '--npol', type=int, default=5,
        help='Polynomial order used for interpolation + 1.')

    parser.add_argument(
        '-o', '--output_filename', type=str, default=DEFAULT_FILE_NAME,
        help='Filename for the database will be <output_filename>.nc.')

    parser.add_argument(
        '--plot', dest='plot', action='store_true', default=False,
        help='Show plots of mesh and interpolation points.')

    args = parser.parse_args()

    # input to SI and consistency checks
    max_depth = 1e3 * args.max_depth
    if args.min_dist < 0:
        raise ValueError('min_dist < 0')
    if args.max_dist < 0:
        raise ValueError('max_dist < 0')
    if args.max_depthi < 0:
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

    if args.output_filename == DEFAULT_FILE_NAME:
        filename = 'isig_%s_%gs_%dkm.nc' % (mod.name, args.period,
                                            args.max_depth)
    else:
        filename = args.output_meshname + '.nc'

    gll_x, gll_y = create_db(filename, mod, m.points, m.connectivity,
                             npol=args.npol, dt=args.dt, npts=args.npts)
    # print mesh info
    info = [
        '=' * 78,
        'SUMMARY OF MESH PROPERTIES:',
        '',
        '  model name                 | %9s' % (mod.name,),
        '  period                     | %9.2f s' % (args.period,),
        '  elements per wavelength    | %9.2f' %
        (args.elements_per_wavelength,),
        '',
        '  time step dt               | %9.4f s' % (args.dt,),
        '  number of time samples     | %9d' % (args.npts,),
        '  number of elements         | %9d' % (m.nelem,),
        '  number of points           | %9d' % (gll_x.size,),
        '  estimated storage (uncomp) | %9.4f GB' % (
            gll_x.size * args.npts * 5 * 4. / 1024. ** 3,),
        '=' * 78]

    info_str = '\n'.join(info)
    print(info_str)

    if args.plot:
        import matplotlib.pyplot as plt
        m.plot(show=False)
        plt.scatter(gll_x, gll_y, color='r')
        plt.show()

    # create sorted and unique point set:
    # r = np.sqrt(gll_x.ravel() ** 2 +  gll_y.ravel() ** 2)
    # theta = np.arctan2(gll_x.ravel(), gll_y.ravel())

    # theta_idx, ntheta = get_global_lexi(theta[np.newaxis,:])

    # theta_unique = theta[np.unique(theta_idx, return_index=True)[1]]
    # theta_unique

    # for i in np.arange(ntheta):
    #     print np.rad2deg(theta_unique[i])
    #     print np.sort((1. - r[theta_idx==i]) * mod.scale / 1e3)
