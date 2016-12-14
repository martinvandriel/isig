#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
A new python script.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    None
'''
import getpass
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import h5netcdf
import socket
import sys
import os

from .gll import gauss_lobatto_legendre_quadruature_points_weights_fast as \
    get_gll
from .basis_polynomials import lagrange_basis_derivative_matrix
from .map_spheroid import map_spheroid


# group: Mesh {
#       dimensions:
#         elements = 9856 ;
#         control_points = 4 ;
#         npol = 5 ;
#       variables:
#         int midpoint_mesh(elements) ;              : not used in instaseis
#         int eltype(elements) ;                     : all 0
#         int axis(elements) ;                       : all 0
#         int fem_mesh(elements, control_points) ;   : OK
#         int sem_mesh(elements, npol, npol) ;       : OK
#         float mp_mesh_S(elements) ;                : OK
#         float mp_mesh_Z(elements) ;                : OK
#         double G0(npol) ;                          : not needed
#         double G1(npol, npol) ;                    : not needed
#         double G2(npol, npol) ;                    : OK
#         double gll(npol) ;                         : OK
#         double glj(npol) ;                         : not needed
#         float mesh_S(gllpoints_all) ;              : OK
#         float mesh_Z(gllpoints_all) ;              : OK
#         float mesh_vp(gllpoints_all) ;             : not needed
#         float mesh_vs(gllpoints_all) ;             : not needed
#         float mesh_rho(gllpoints_all) ;            : not needed
#         float mesh_lambda(gllpoints_all) ;         : not needed
#         float mesh_mu(gllpoints_all) ;             : OK
#         float mesh_xi(gllpoints_all) ;             : not needed
#         float mesh_phi(gllpoints_all) ;            : not needed
#         float mesh_eta(gllpoints_all) ;            : not needed
#         float mesh_Qmu(gllpoints_all) ;            : not needed
#         float mesh_Qka(gllpoints_all) ;            : not needed
#       }


def create_db(fname, model, points, connectivity, npol=5, dt=0.1,
              npts=1000):

    gll = get_gll(npol)[0]
    gll_x, gll_y = map_spheroid(gll, points[connectivity])

    nelem = connectivity.shape[0]
    nquad = 4

    with h5netcdf.File(fname, "w") as f:

        f.dimensions = {
            'gllpoints_all': gll_x.size,
            'snapshots': npts,
            'ipol': npol,
            'jpol': npol,
            'nvars': 5,
            'elements': nelem
        }

        # MERGED DATA VARIABLES
        f.create_variable('stf_dump', ('snapshots', ), 'float32')
        f.create_variable('stf_d_dump', ('snapshots', ), 'float32')
        f.create_variable('MergedSnapshots',
                          ('elements', 'nvars', 'jpol', 'ipol', 'snapshots'),
                          'float32')

        # MESH GROUP
        mesh_group = f.create_group("Mesh")
        mesh_group.dimensions = {
            # 'elements': nelem,
            'control_points': nquad,
            'npol': npol
        }

        sem_mesh = mesh_group.create_variable(
            'sem_mesh', ('elements', 'npol', 'npol'), 'int32')

        sem_mesh[...] = \
            np.arange(nelem * npol ** 2).reshape((nelem, npol, npol))

        fem_mesh = mesh_group.create_variable(
            'fem_mesh', ('elements', 'control_points'), 'int32')

        fem_mesh[:, 0] = sem_mesh[:, 0, 0]
        fem_mesh[:, 1] = sem_mesh[:, 0, -1]
        fem_mesh[:, 2] = sem_mesh[:, -1, -1]
        fem_mesh[:, 3] = sem_mesh[:, -1, 0]

        mp_mesh_S = mesh_group.create_variable(
            'mp_mesh_S', ('elements', ), 'float32')
        mp_mesh_Z = mesh_group.create_variable(
            'mp_mesh_Z', ('elements', ), 'float32')

        gll = np.zeros(1)
        mp_mesh_S[:], mp_mesh_Z[:] = [
            g.ravel() * model.scale
            for g in map_spheroid(gll, points[connectivity])]

        G2 = mesh_group.create_variable('G2', ('npol', 'npol'), float)
        gll = mesh_group.create_variable('gll', ('npol', ), float)

        gll[:] = get_gll(npol)[0]
        G2[...] = lagrange_basis_derivative_matrix(gll[:])

        mesh_S = mesh_group.create_variable(
            'mesh_S', ('gllpoints_all', ), 'float32')
        mesh_Z = mesh_group.create_variable(
            'mesh_Z', ('gllpoints_all', ), 'float32')

        mesh_S[:] = gll_x.ravel() * model.scale
        mesh_Z[:] = gll_y.ravel() * model.scale

        r = np.sqrt(gll_x.ravel() ** 2 + gll_y.ravel() ** 2)
        theta = np.arctan2(gll_x.ravel(), gll_y.ravel())
        mp = np.sqrt(mp_mesh_S[:] ** 2 + mp_mesh_Z[:] ** 2) / model.scale
        mp = mp.repeat(npol ** 2)

        var = mesh_group.create_variable('mesh_mu', ('gllpoints_all', ),
                                         'float32')
        var[:] = model.get_elastic_parameter('MU', r, mp)

        # GLOBAL ATTRIBUTES
        # @TODO: replace place holder and meaningless names
        f.attrs['dump type (displ_only, displ_velo, fullfields)'] = \
            'displ_only'
        # Don't make sense for Merged DBs
        # f.attrs['excitation_type'] =  ''
        # f.attrs['source type'] =  ''
        f.attrs['background model'] = model.name
        f.attrs['git commit hash'] = ''
        f.attrs['datetime'] = str(datetime.now())
        f.attrs['source time function'] = ''
        f.attrs['npol'] = npol
        f.attrs['file version'] = 0
        f.attrs['number of strain dumps'] = npts
        f.attrs['scalar source magnitude'] = 0.
        f.attrs['strain dump sampling rate in sec'] = dt
        f.attrs['source shift factor in sec'] = 0.
        f.attrs['source shift factor for deltat_coarse'] = 0
        f.attrs['npoints'] = gll_x.size
        f.attrs['attenuation'] = 1
        f.attrs['planet radius'] = model.scale
        f.attrs['dominant source period'] = 0.
        f.attrs['kernel wavefield rmin'] = r.min() * model.scale
        f.attrs['kernel wavefield rmax'] = r.max() * model.scale
        f.attrs['kernel wavefield colatmin'] = np.rad2deg(theta.min())
        f.attrs['kernel wavefield colatmax'] = np.rad2deg(theta.max())
        f.attrs['source depth in km'] = 0.
        f.attrs['nelem_kwf_global'] = 0

        f.attrs['compiler brand'] = ''
        f.attrs['compiler version'] = ''
        f.attrs['user name'] = getpass.getuser()
        f.attrs['host name'] = socket.gethostname()
        f.attrs['rundir'] = os.getcwdu()
        f.attrs['cmdl'] = 'python -m isig ' + ' '.join(sys.argv[1:])

    return gll_x, gll_y
