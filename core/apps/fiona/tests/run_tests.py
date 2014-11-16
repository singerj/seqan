#!/usr/bin/env python
"""Execute the tests for the fiona program.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import subprocess
import sys

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                    '..', '..', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests


class ResultChecker(object):
    """Check quality with compute_gain program."""

    def __init__(self, compute_gain, path_to_ref, path_to_sam,
                path_to_post, min_gain):
        self.compute_gain = compute_gain
        self.path_to_ref = path_to_ref
        self.path_to_sam = path_to_sam
        self.path_to_post = path_to_post
        self.min_gain = min_gain

    def __call__(self):
        cmd_line = [self.compute_gain,
                    '-g', self.path_to_ref,
                    '--pre', self.path_to_sam,
                    '--post', self.path_to_post]
        print ' '.join(cmd_line)
        process = subprocess.Popen(cmd_line, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        retcode = process.wait()
        if retcode != 0:
            raise app_tests.BadResultException('compute_gain did not return 0')
        report = process.stdout.read()
        gain = float(report.splitlines()[2].split()[0])
        if gain < self.min_gain:
            fmt = 'Gain too low. Expected >= %f, got %f'
            raise app_tests.BadResultException(fmt % (self.min_gain, gain))
        return True


def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for fiona'
    print '========================'
    print

    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/fiona/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_fiona = app_tests.autolocateBinary(
      binary_base, 'core/apps/fiona', 'fiona')
    path_to_fiona_illumina = app_tests.autolocateBinary(
      binary_base, 'core/apps/fiona', 'fiona_illumina')
    path_to_compute_gain = app_tests.autolocateBinary(
      binary_base, 'core/apps/fiona', 'compute_gain')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # We prepare a list of transforms to apply to the output files.  This is
    # used to strip the input/output paths from the programs' output to
    # make it more canonical and host independent.
    ph.outFile('-')  # To ensure that the out path is set.
    transforms = [
        app_tests.ReplaceTransform(
            os.path.join(ph.source_base_path,
                         'core/apps/fiona/tests') + os.sep,
            '', right=True),
        app_tests.ReplaceTransform(ph.temp_dir + os.sep, '', right=True),
        app_tests.NormalizeScientificExponentsTransform(),
        ]

    # ============================================================
    # Run on uniformly random DNA.
    # ============================================================

    # Note that instead of comparing the results with expected results, we
    # use a checker that computes the gain and compares it with a threshold.

    # Illumina Mode

    for i in [1, 2]:
        min_gain = {1: 40.0, 2: 50.0}
        conf = app_tests.TestConf(
            program=path_to_fiona_illumina,
            args=['-nt', '1',
                  '-i', str(i),
                  '-g', '10000',
                  ph.inFile('reads.illumina.fq'),
                  ph.outFile('reads.illumina.corrected.i%d.fa' % i)],
            redir_stdout=ph.outFile('reads.illumina.fq.i%d.stdout' % i),
            redir_stderr=ph.outFile('reads.illumina.fq.i%d.stderr' % i),
            check_callback=ResultChecker(
                path_to_compute_gain, ph.inFile('genome.10k.fa'),
                ph.inFile('reads.illumina.sam'),
                ph.outFile('reads.illumina.corrected.i%d.fa' % i),
                min_gain.get(i, 100.0)),
            to_diff=[(ph.inFile('reads.illumina.fq.i%d.stdout' % i),
                      ph.outFile('reads.illumina.fq.i%d.stdout' % i),
                      transforms),
                     (ph.inFile('reads.illumina.fq.i%d.stderr' % i),
                      ph.outFile('reads.illumina.fq.i%d.stderr' % i),
                      transforms),
                    ])
        conf_list.append(conf)

    # Indel Mode

    for i in [1, 2]:
        min_gain = {1: 70.0, 2: 85.0}
        conf = app_tests.TestConf(
            program=path_to_fiona,
            args=['-nt', '1',
                  '-i', str(i),
                  '-g', '10000',
                  ph.inFile('reads.454.fq'),
                  ph.outFile('reads.454.corrected.i%d.fa' % i)],
            redir_stdout=ph.outFile('reads.454.fq.i%d.stdout' % i),
            redir_stderr=ph.outFile('reads.454.fq.i%d.stderr' % i),
            check_callback=ResultChecker(
                path_to_compute_gain, ph.inFile('genome.10k.fa'),
                ph.inFile('reads.454.sam'),
                ph.outFile('reads.454.corrected.i%d.fa' % i),
                min_gain.get(i, 100.0)),
            to_diff=[(ph.inFile('reads.454.fq.i%d.stdout' % i),
                      ph.outFile('reads.454.fq.i%d.stdout' % i),
                      transforms),
                     (ph.inFile('reads.454.fq.i%d.stderr' % i),
                      ph.outFile('reads.454.fq.i%d.stderr' % i),
                      transforms),
                    ])
        conf_list.append(conf)

    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(conf.commandLineArgs())
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    # Cleanup.
    ph.deleteTempDir()

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))
