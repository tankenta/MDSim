#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import os
import sys

from matplotlib import pyplot as plt

from config import *

def csv2result_dict(csv_path):
    result_dict = dict()

    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        result_dict['label_x'] = header[0]
        result_dict['label_y'] = header[1]

        array_x = list()
        array_y = list()
        for row in reader:
            array_x.append(float(row[0]))
            array_y.append(float(row[1]))
        result_dict['array_x'] = array_x
        result_dict['array_y'] = array_y

    return result_dict

def plot_all_results(csv_dir):
    fig = plt.figure(figsize=FIG_SIZE)
    for key, conf in PLOT_CONF.items():
        csv_path = os.path.join(csv_dir, conf['file_name'])
        result_dict = csv2result_dict(csv_path)

        ax = fig.add_subplot(conf['plot_pos'])
        if conf['type'] == 'plot':
            ax.plot(
                    result_dict['array_x'], result_dict['array_y'],
                    conf['style'], label=conf['label'])
        else:
            ax.scatter(
                    result_dict['array_x'], result_dict['array_y'],
                    conf['size'], conf['style'], label=conf['label'])
        ax.set_xlabel(result_dict['label_x'])
        ax.set_ylabel(result_dict['label_y'])
        ax.xaxis.set_tick_params(
                which='both', direction='in', left=True, right=True)
        ax.yaxis.set_tick_params(
                which='both', direction='in', left=True, right=True)
        ax.legend()
    fig.tight_layout()
    plt.show()

def main():
    if len(sys.argv) < 2:
        print('Usage: ./plot_results.py csv_dir')
        sys.exit()
    csv_dir = sys.argv[1]
    plot_all_results(csv_dir)

if __name__ == '__main__':
    main()
