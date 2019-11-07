#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@email.arizona.edu>
Date   : 2019-11-07
Purpose: EBI Pipeline
"""

import argparse
import csv
import logging
import os
import re
import sys
import subprocess
import parallelprocs
from typing import Dict, List, TextIO, Optional, Sequence, Tuple, Any


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='EBI Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('experiment_id',
                        type=int,
                        nargs='+',
                        help='Experiment ID(s)')

    parser.add_argument('-d',
                        '--data_dir',
                        type=str,
                        metavar='DIR',
                        required=True,
                        help='Data directory')

    parser.add_argument('-m',
                        '--metadata',
                        type=argparse.FileType('r'),
                        metavar='FILE',
                        required=True,
                        help='Experiment metadata file')

    parser.add_argument('-D',
                        '--debug',
                        action='store_true',
                        help='Log at debug level')

    parser.add_argument('-l',
                        '--logfile',
                        metavar='str',
                        type=str,
                        default='.log',
                        help='Log filename')

    args = parser.parse_args()

    if not os.path.exists(args.data_dir):
        os.makedirs(args.data_dir)

    if not os.path.isdir(args.data_dir):
        parser.error(f'--data_dir "{args.data_dir}" not a directory')

    if not os.path.isabs(args.data_dir):
        args.data_dir = os.path.abspath(args.data_dir)

    return args


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    logging.basicConfig(filename=args.logfile,
                        filemode='w',
                        level=logging.DEBUG if args.debug else logging.INFO)

    for i, exp_id in enumerate(args.experiment_id, start=1):
        print(f'{i:3}: {exp_id}')
        process(exp_id, args.metadata, args.data_dir)

    print('Done.')


# --------------------------------------------------
def process(exp_id: int, fh: TextIO, data_dir: str) -> None:
    """Process experiment"""

    meta = find_meta(exp_id, fh)
    urls = list(map(lambda m: m['url'], meta))
    local = get_files(urls, data_dir)
    run_qc(local, data_dir)


# --------------------------------------------------
def find_meta(exp_id: int, fh: TextIO) -> Sequence[Dict[str, str]]:
    """Find metadata info"""

    reader = csv.DictReader(fh, delimiter=',')
    if not 'experiment_id' in reader.fieldnames:
        raise Exception(f'"{fh.name}" missing "experiment_id" field')

    meta = list(filter(lambda r: r['experiment_id'] == str(exp_id), reader))

    if not meta:
        raise Exception(f'experiment_id "{exp_id}" not found')

    return meta


# --------------------------------------------------
def get_files(files: List[str], data_dir: str) -> List[str]:
    """Put the files into the directory"""

    expected: List[str] = []
    todo = []
    raw_dir = os.path.join(data_dir, 'raw')
    if not os.path.isdir(raw_dir):
        os.makedirs(raw_dir)

    for file in files:
        path = os.path.join(raw_dir, os.path.basename(file))
        if path in expected:
            continue

        expected.append(path)

        if os.path.isfile(path):
            logging.debug('Already downloaded "%s"', file)
            continue

        logging.debug(f'Will get "{file}"')
        todo.append(f'iget -f {file} {path}')

    if todo:
        parallelprocs.run(todo, halt=1, num_procs=8)

    missing = list(filter(lambda f: not os.path.isfile(f), expected))
    if missing:
        raise Exception('Failed to download {}'.format(', '.join(missing)))

    return expected


# --------------------------------------------------
def run_qc(files: List[str], data_dir: str) -> List[str]:
    """Run QC on files"""

    qc_dir = os.path.join(data_dir, 'qc')
    if not os.path.isdir(qc_dir):
        os.makedirs(qc_dir)

    tmpl = ('{} | fastq_quality_filter -q 20 -p 80 -Q 33 '
            '| fastx_clipper -l 80 | fastx_collapser > {}')

    expected = []
    todo = []
    for file in files:
        basename = os.path.basename(file)
        gzipped = basename.endswith('.gz')
        match = re.match('(.+)\.f(?:ast)?q(?:\.gz)?$', basename)
        if match:
            basename = match.group(1)

        qc_file = os.path.join(qc_dir, basename + '.fa')
        expected.append(qc_file)

        if os.path.isfile(qc_file):
            logging.debug(f'QC file "{qc_file}" already exists')
            continue

        source = f'gunzip -c {file}' if gzipped else f'cat {file}'
        cmd = tmpl.format(source, qc_file)
        logging.debug(cmd)
        todo.append(cmd)

    if todo:
        parallelprocs.run(todo, halt=1, num_procs=8)

    missing = list(filter(lambda f: not os.path.isfile(f), expected))
    if missing:
        raise Exception('Failed to download {}'.format(', '.join(missing)))

    return expected


# --------------------------------------------------
def fst(t: Tuple) -> Any:
    """Return first of a tuple"""

    return t[0]


# --------------------------------------------------
def snd(t: Tuple) -> Any:
    """Return second of a tuple"""

    return t[1]


# --------------------------------------------------
if __name__ == '__main__':
    main()
