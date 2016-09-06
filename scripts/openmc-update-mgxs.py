#!/usr/bin/env python
"""Update OpenMC's deprecated multi-group cross section XML files to the latest
HDF5-based format.

Usage information can be obtained by running 'openmc-update-mgxs --help':

usage: openmc-update-mgxs [-h] in out

Update mgxs.xml files to the latest format. This will remove 'outside'
attributes/elements from lattices and replace them with 'outer' attributes. For
'cell' elements, any 'surfaces' attributes/elements will be renamed
'region'. Note that this script will not delete the given files; it will append
'.original' to the given files and write new ones.

positional arguments:
  in          Input mgxs xml file
  out         Output mgxs hdf5 file

optional arguments:
  -h, --help  show this help message and exit

"""

from __future__ import print_function
from shutil import move
import warnings
import xml.etree.ElementTree as ET

import argparse
import h5py
import numpy as np

import openmc.mgxs_library

description = """\
Update OpenMC's deprecated multi-group cross section XML files to the latest
HDF5-based format."""


def parse_args():
    """Read the input files from the commandline."""
    # Create argument parser
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='input XML file')
    parser.add_argument('-o', '--output', nargs='?', default='',
                        help='output file, in HDF5 format')
    parser.add_argument('-c', '--compression', type=int,
                        help='HDF5 Compression Level')
    args = vars(parser.parse_args())

    if args['output'] == '':
        filename = args['input'].name
        extension = filename[filename.rfind('.'):]
        if extension == '.xml':
            filename = filename[:filename.rfind('.')] + '.h5'
        args['output'] = filename

    # Parse and return commandline arguments.
    return args


if __name__ == '__main__':
    args = parse_args()

    # Parse the XML data.
    tree = ET.parse(args['input'])
    root = tree.getroot()

    if root.tag != 'library':
        raise ValueError("Invalid XML file type")

    # Get old metadata
    temp = tree.find('group_structure').text.strip()
    temp = np.array(temp.split())
    group_structure = temp.astype(np.float)
    energy_groups = openmc.mgxs.EnergyGroups(group_structure)
    temp = tree.find('inverse_velocities')
    if temp is not None:
        temp = temp.text.strip()
        temp = np.array(temp.split())
        inverse_velocities = temp.astype(np.float)
    else:
        inverse_velocities = None

    xsd = []
    names = []

    # Now move on to the cross section data itself
    for xsdata_elem in root.iter('xsdata'):
        name = xsdata_elem.find('name').text.strip()

        temperature = xsdata_elem.find('kT')
        if temperature is not None:
            temperature = \
                float(temperature.text.strip()) / openmc.data.K_BOLTZMANN
        else:
            temperature = 294.
        temperatures = [temperature]

        awr = xsdata_elem.find('awr')
        if awr is not None:
            awr = float(awr.text.strip())

        representation = xsdata_elem.find('representation')
        if representation is not None:
            representation = representation.text.strip()
        else:
            representation = 'isotropic'
        if representation == 'angle':
            n_azi = int(xsdata_elem.find('num_azimuthal').text.strip())
            n_pol = int(xsdata_elem.find('num_polar').text.strip())

        scatter_type = xsdata_elem.find('scatt_type')
        if scatter_type is not None:
            scatter_type = scatter_type.text.strip()
        else:
            scatter_type = 'legendre'

        order = int(xsdata_elem.find('order').text.strip())

        tab_leg = xsdata_elem.find('tabular_legendre')
        if tab_leg is not None:
            warnings.Warning('The tabular_legendre option has moved to the '
                             'settings.xml file and must be added manually')

        # Either add the data to a previously existing xsdata (if it is
        # for the same 'name' but a different temperature), or create a
        # new one.

        try:
            # It is in our list, so store that entry
            i = names.index(name)
        except:
            # It is not in our list, so add it
            i = -1
            xsd.append(openmc.XSdata(name, energy_groups,
                                     temperatures=temperatures,
                                     representation=representation))
            if awr is not None:
                xsd[-1].awr = awr
            if representation == 'angle':
                xsd[-1].num_azimuthal = n_azi
                xsd[-1].num_polar = n_pol
            xsd[-1].scatter_type = scatter_type
            xsd[-1].order = order
            names.append(name)

        if scatter_type == 'legendre':
            order_dim = order + 1
        else:
            order_dim = order

        if i != -1:
            xsd[i].add_temperature(temperature)

        temp = xsdata_elem.find('total')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            total = temp.astype(np.float)
            total = np.reshape(total, xsd[i].vector_shape)
            xsd[i].set_total(total, temperature)

        temp = xsdata_elem.find('absorption').text.strip()
        temp = np.array(temp.split())
        absorption = temp.astype(np.float)
        absorption = np.reshape(absorption, xsd[i].vector_shape)
        xsd[i].set_absorption(absorption, temperature)

        temp = xsdata_elem.find('scatter').text.strip()
        temp = np.array(temp.split())
        temp = temp.astype(np.float)
        scatter = np.reshape(temp, xsd[i].pn_matrix_shape)
        xsd[i].set_scatter_matrix(scatter, temperature)

        temp = xsdata_elem.find('multiplicity')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            temp = temp.astype(np.float)
            multiplicity = np.reshape(temp, xsd[i].matrix_shape)
            xsd[i].set_multiplicity_matrix(multiplicity, temperature)

        temp = xsdata_elem.find('fission')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            fission = temp.astype(np.float)
            fission = np.reshape(fission, xsd[i].vector_shape)
            xsd[i].set_fission(fission, temperature)

        temp = xsdata_elem.find('kappa_fission')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            kappa_fission = temp.astype(np.float)
            kappa_fission = np.reshape(kappa_fission, xsd[i].vector_shape)
            xsd[i].set_kappa_fission(kappa_fission, temperature)

        temp = xsdata_elem.find('chi')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            chi = temp.astype(np.float)
            chi = np.reshape(chi, xsd[i].vector_shape)
            xsd[i].set_chi(chi, temperature)
        else:
            chi = None

        temp = xsdata_elem.find('nu_fission')
        if temp is not None:
            temp = temp.text.strip()
            temp = np.array(temp.split())
            temp = temp.astype(np.float)
            if chi is not None:
                nu_fission = np.reshape(temp, xsd[i].vector_shape)
            else:
                nu_fission = np.reshape(temp, xsd[i].matrix_shape)
            xsd[i].set_nu_fission(nu_fission, temperature)

    # Build library as we go, but first we have enough to initialize it
    lib = openmc.MGXSLibrary(energy_groups)
    if inverse_velocities is not None:
        lib.inverse_velocities = inverse_velocities

    lib.add_xsdatas(xsd)

    lib.export_to_hdf5(args['output'], compression=args['compression'])
