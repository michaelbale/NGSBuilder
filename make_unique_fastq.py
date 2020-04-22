#!/usr/bin/env python3

import sys
import os

def process_fq(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def read_fq(fn = None):
    if not os.path.exists(fn):
        raise SystemError("Error: File does not exist\n")
    records = []

    with open(fn, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4:
                records.append(process_fq(lines))
                lines = []
    return records

def reduce_fq(records=None):
    unique_recs = list({v['name']:v for v in records}.values())
    return( sorted(unique_recs, key = lambda i: i['name']) )

def write_fq(records=None, fn=None):
    if os.path.exists(fn):
        raise SystemError("Error: Supplied File Already Exists\n")
    
    with open(fn, 'w') as fh:
        for r in records:
            fh.write(r['name'] + '\n')
            fh.write(r['sequence'] + '\n')
            fh.write(r['optional'] + '\n')
            fh.write(r['quality'] + '\n')

def main(argv):

    in_file = argv[0]
    out_file = argv[1]
    fq = read_fq(fn = in_file)
    new_fq = reduce_fq(records = fq)
    write_fq(records = new_fq, fn = out_file)

if __name__ == "__main__":
    main(sys.argv[1:])
