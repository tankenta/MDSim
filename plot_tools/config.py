#!/usr/bin/env python
# -*- coding: utf-8 -*-

FIG_SIZE = (8, 6)

PLOT_CONF = {
        'total': {
            'file_name': 'total_energy.csv',
            'plot_pos': 321,
            'type': 'scatter',
            'label': 'total energy',
            'style': 'r',
            'size': 10,
        },
        'potential': {
            'file_name': 'potential_energy.csv',
            'plot_pos': 322,
            'type': 'scatter',
            'label': 'potential energy',
            'style': 'g',
            'size': 10,
        },
        'kinetic': {
            'file_name': 'kinetic_energy.csv',
            'plot_pos': 323,
            'type': 'scatter',
            'label': 'kinetic energy',
            'style': 'b',
            'size': 10,
        },
        'temp': {
            'file_name': 'temperature.csv',
            'plot_pos': 324,
            'type': 'scatter',
            'label': 'temperature',
            'style': 'm',
            'size': 10,
        },
        'rdf': {
            'file_name': 'RDF.csv',
            'plot_pos': 325,
            'type': 'plot',
            'label': 'RDF',
            'style': 'r',
        },
        'msd': {
            'file_name': 'MSD.csv',
            'plot_pos': 326,
            'type': 'scatter',
            'label': 'MSD',
            'style': 'c',
            'size': 10,
        }
}

