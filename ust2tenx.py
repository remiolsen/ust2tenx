import sys
import re
import io
import os
import subprocess
import argparse
import json
import math
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import cycle

# 10X index kit to use as default
SIP02F8 = ["CATGAACA","TCACTCGC","AGCTGGAT","GTGACTTG"]
# "Random" 7 bp oligo
oligo = "TTGCGAG"
r1_nbc = ''.join(["N"] * 16)
r1_qual = ''.join(["A"] * 23)
i1_qual = ''.join(["A"] * 8)
gz_buf = 131072
try:
    subprocess.Popen(["pigz"], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
except OSError as e:
    if e.errno == os.errno.ENOENT:
        print("Could not find 'pigz' command in system environment")
        sys.exit(1)
    else:
        raise

ilmnbc = {}
ilmnbc_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "ilmn_indices.json")
try:
    with open(ilmnbc_file, "r") as f:
        ilmnbc = json.load(f)
except (IOError, ValueError):
    print("Error parsing the file 'ilmn_indices.json' will only use SI-P02-F8 barcode set")
    ilmnbc["SI-P02-F8"] = SIP02F8


### START DRY-principle violating code-block
def write_read1(ustmap, tenx, read, prefix, procs, tenx_kit):
    idx_loop = cycle(ilmnbc[tenx_kit])
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"R1_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    try:
                        ust_bc = tenx[ustmap[tarr[0]]]
                    except KeyError:
                        ust_bc = r1_nbc
                        #print("No barcode in file {} in read {}, inserting Ns".format(read, tarr[0]), file=sys.stderr)
                        pass
                    except IndexError:
                        ust_bc = r1_nbc
                        #print("No barcode in file {} in read {}, inserting Ns".format(read, tarr[0]), file=sys.stderr)
                        pass

                    outln = "@{} 1:N:0:{}\n".format(tarr[0], next(idx_loop))
                    outln += "{}{}{}\n".format(ust_bc, oligo, seq)
                    outln += "+\n"
                    outln += "{}{}\n".format(r1_qual, qual)
                    oz.stdin.write(outln.encode('utf-8'))

def write_read2(ustmap, tenx, read, prefix, procs, tenx_kit):
    idx_loop = cycle(ilmnbc[tenx_kit])
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"R2_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    outln = "@{} 2:N:0:{}\n".format(tarr[0], next(idx_loop))
                    outln += "{}\n".format(seq)
                    outln += "+\n"
                    outln += "{}\n".format(qual)
                    oz.stdin.write(outln.encode('utf-8'))

def write_i1(ustmap, tenx, read, prefix, procs, tenx_kit):
    idx_loop = cycle(ilmnbc[tenx_kit])
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"I1_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    idx_i = next(idx_loop)
                    outln = "@{} 1:N:0:{}\n".format(tarr[0], idx_i)
                    outln += "{}\n".format(idx_i)
                    outln += "+\n"
                    outln += "{}\n".format(i1_qual)
                    oz.stdin.write(outln.encode('utf-8'))
### END DRY-principle violating code-block


def main(tenxfile, r1, r2, i1, prefix, total_processes, minbc, tenx_kit, max_bc_split):
    idx = 0
    tenx_c = 0
    TENX_BC = []
    ustmap = {}
    ustmap_c = {}
    ustmap_r = {}
    ustmaps = [] # If using barcode splitting
    ustreads = [] # If using barcode splitting

    if tenx_kit not in ilmnbc.keys():
        print("Did not find 10X kit {}".format(tenx_kit))
        return

    with open(tenxfile, 'r') as f:
        for line in f:
            TENX_BC.append(line.strip())
            tenx_c += 1

    with subprocess.Popen(["pigz", "-d", "-c", "-p", str(total_processes), i1],
            stdout=subprocess.PIPE) as fz:
        with io.TextIOWrapper(fz.stdout, write_through=True) as f:
            for title, seq, _ in FastqGeneralIterator(f):
                rname = title.split()[0]
                if len(seq) == 18:
                    to_map = ustmap_r.get(seq, [])
                    to_map.append(rname)
                    ustmap_r[seq] = to_map

    for seq, reads in ustmap_r.items():
        count = len(reads)
        assert idx <= tenx_c, "Found more barcodes that available for 10X ({}), try using --min-bc argument".format(tenx_c)
        if count >= minbc:
            for read in reads:
                ustmap[read] = idx
            idx += 1

    print("Found:\t{} UST barcodes".format(len(ustmap_r.keys())))
    print("Made:\t{} 10X barcodes".format(idx))
    write_read1(ustmap, TENX_BC, r1, prefix, total_processes, tenx_kit)
    write_read2(ustmap, TENX_BC, r2, prefix, total_processes, tenx_kit)
    write_i1(ustmap, TENX_BC, r1, prefix, total_processes, tenx_kit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate UST (Tellseq) barcodes to 10X barcodes")
    parser.add_argument('--tenx-bc-file', '-t',type=str, default="{}/supernova/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt".format(os.path.dirname(os.path.realpath(__file__))), help="Path to text file with 10X barcodes (eg. 4M-with-alts-february-2016.txt)")
    parser.add_argument('--r1', '-1', type=str, required=True, help="Path to read 1")
    parser.add_argument('--r2', '-2', type=str, required=True, help="Path to read 2")
    parser.add_argument('--i1', '-i', type=str, required=True, help="Path to I1 read")
    parser.add_argument('--out-prefix', '-o', type=str, default="UST_OUT_S1_L001_", help="Prefix of the output fastq files (default: UST_OUT_S1_L001_)")
    parser.add_argument('--processes', '-p', type=str, default=2, help="Number of processes to spawn (default: 2)")
    parser.add_argument('--min-bc', '-m', type=int, default=1, help="Minumum barcode multiplicity to include it")
    parser.add_argument('--max-bc-split', '-s', type=int, default=0, help="Set a threshold of maximum 10X barcodes per file. If it exceeds this, additional files ('libraries') will be output. Leave this as 0 for only one set of output files.")
    parser.add_argument('--tenx-kit', '-k', default='SI-P02-F8', help="Which 10X barcode kit the converted files should use")
    args = parser.parse_args()
    sys.exit(main(args.tenx_bc_file, args.r1, args.r2, args.i1, args.out_prefix, args.processes, args.min_bc, args.tenx_kit, args.max_bc_split))
