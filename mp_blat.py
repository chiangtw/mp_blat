#! /usr/bin/env python

import argparse
import shlex
import shutil
import subprocess as sp
import os.path
import tempfile as tp
import concurrent.futures as cf
import logging
from functools import partial


logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


class Blat:
    def __init__(self, blat_bin='blat', blat_options=''):
        self.blat_bin = blat_bin
        self.blat_options = shlex.split(blat_options)

    def _run(self, ref_file, fa_file, out_file):
        logging.info(f'Starting the blat process on {fa_file}')
        blat_result = sp.run(
            [self.blat_bin, ref_file, fa_file, out_file] + self.blat_options,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            encoding='utf-8'
        )
        logging.info(f'The blat process completed. ({fa_file}) ')

        return blat_result

    def __call__(self, ref_file, fa_file, out_file):
        return self._run(ref_file, fa_file, out_file)


class FastaIndex:
    def __init__(self, name, length, offset, line_bases, line_width):
        self.name = name
        self.length = int(length)
        self.offset = int(offset)
        self.line_bases = int(line_bases)
        self.line_width = int(line_width)

    @property
    def end_pos(self):
        num_lines, r = divmod(self.length, self.line_bases)

        if r > 0:
            num_lines += 1

        end_pos = self.offset + self.length + num_lines

        return end_pos


class FastaFile:
    def __init__(self, fa_file):
        self.original_file = fa_file
        self.samtools_bin = samtools_bin
        self.index_file = None
        self.files = None

        self.original_filename = os.path.basename(self.original_file)

    def has_index(self):
        if os.path.exists(self.original_file + ".fai"):
            self.index_file = self.original_file + ".fai"
            return True
        else:
            return False

    def create_index(self, samtools_bin="samtools"):
        result = sp.run(
            [samtools_bin, 'faidx', self.original_file],
            stderr=sp.PIPE
        )

        if result.returncode == 0:
            self.index_file = self.original_file + ".fai"

        return result

    def split(self, num_files, work_dir='.', samtools_bin="samtools"):
        logging.info(f'Start to split the fasta file into {num_files} parts.')

        logging.info('Search for index file...')

        if not self.has_index():

            logging.info('Index file not found.')
            logging.info('Generating the index...')

            self.create_index(samtools_bin)

            logging.info('Done!')

        logging.info(f'Index file: {self.index_file}')

        logging.info('Calculating the size of each split block.')
        with open(self.index_file) as idx_in:
            fa_idx_data = []
            for line in idx_in:
                data = line.rstrip('\n').split('\t')
                idx_data = FastaIndex(*data)
                fa_idx_data.append(idx_data)

        # Get splicing pos
        num_for_each_file, remainder = divmod(len(fa_idx_data), num_files)
        num_list = [num_for_each_file] * num_files
        for i in range(remainder):
            num_list[i] += 1

        splice_pos = []
        idx = -1
        for num in num_list[:-1]:
            idx += num
            pos = fa_idx_data[idx].end_pos
            splice_pos.append(pos)

        # Get block size
        splice_pos_extend = [0] + splice_pos + [fa_idx_data[-1].end_pos]
        block_size = [
            p1 - p2
            for p1, p2 in zip(splice_pos_extend[1:], splice_pos_extend[:-1])
        ]

        # output data to files
        self.files = [
            os.path.join(work_dir, self.original_filename + f'.part_{i}')
            for i in range(1, num_files + 1)
        ]

        logging.info('Writing out the split fasta data to files.')
        with open(self.original_file) as fa_in:
            for file_, bsize in zip(self.files, block_size):
                with open(file_, 'w') as fa_out:
                    fa_out.write(fa_in.read(bsize))

        logging.info(f'The fasta file has been split into {num_files} parts.')


def mp_blat(reference_file,
            fasta_file,
            output_file,
            num_proc=1,
            tmp_path=".",
            blat_bin="blat",
            blat_options="",
            samtools_bin="samtools"):

    logging.info(
        f"Start mp_blat with {num_proc} process{'' if num_proc==1 else 'es'}"
    )

    _blat = Blat(blat_bin, blat_options)
    blat = partial(_blat, reference_file)

    if num_proc == 1:
        blat(fasta_file, output_file)
    else:
        # 1. split fa_file into n parts
        # 2. run blat with n processes
        # 3. merge blat results
        tmp_dir = tp.TemporaryDirectory(prefix='mp_blat_tmp.', dir='.')

        fa_file = FastaFile(fasta_file)
        fa_file.split(num_proc, tmp_dir.name, samtools_bin=samtools_bin)

        out_files = [file_ + '.psl' for file_ in fa_file.files]

        with cf.ProcessPoolExecutor(num_proc) as executor:
            executor.map(blat, fa_file.files, out_files)

        with open(output_file, 'wb') as out:
            # write headers
            with open(out_files[0], 'rb') as res_in:
                for _ in range(5):
                    out.write(res_in.readline())
                res_start_pos = res_in.tell()

            # write results
            for file_ in out_files:
                with open(file_, 'rb') as res_in:
                    res_in.seek(res_start_pos)
                    shutil.copyfileobj(res_in, out)

    logging.info(f"The process{'' if num_proc==1 else 'es'} completed!")


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("fasta")
    parser.add_argument("output")
    parser.add_argument("-p", "--num_proc", type=int, default=1)
    parser.add_argument("--tmp_path", type=str, default=".")
    parser.add_argument("--blat_bin", default="blat")
    parser.add_argument("--blat_options", type=str, default="")
    parser.add_argument("--samtools_bin", default="samtools")

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    mp_blat(
        args.reference,
        args.fasta,
        args.output,
        num_proc=args.num_proc,
        tmp_path=args.tmp_path,
        blat_bin=args.blat_bin,
        blat_options=args.blat_options,
        samtools_bin=args.samtools_bin
    )
