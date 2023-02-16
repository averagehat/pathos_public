'''
Usage:
    pathos_sheet.py <samplesheet> --config <config> --sampledir <sampledir> [--log <log>] [--outdir <outdir>] [-p] [--qsub]

Options:
    -c <config>, --config <config>
    --o <outdir>, --outdir <outdir>
    --log <log>
    --sampledir <sampledir>
    --qsub, -q
'''


from glob import glob
from docopt import docopt
import csv
import os
import pipeline
from path import Path
import yaml
import itertools
import sys
from functools import partial
import multiprocessing
# if there is more than one pair of read files need to handle that. I guess by concatenating.
# need to fix so that we join on conting ID. currently joining on something else.
# also ppl might do something with index


SEP=';'


def weave_files(sampledir, row):
  control_dirs = row['control_dirs'].split(SEP)
  sdir = partial(os.path.join, sampledir)
  read_dirs = row['read_dirs'].split(SEP)
  # control_dirs could be empty
  c1s, c2s = reads_from_dirs(map(sdir, control_dirs)) if control_dirs else ([], [])
  controls = list(itertools.chain(*zip(c1s, c2s)))
  # order of contams doesn't matter but we'll zip them anyway
  r1s, r2s = reads_from_dirs(map(sdir, read_dirs))
  # weave them into a sigle list of R1, R2^1, R1^2, R2^2 ....
  fastqs = list(itertools.chain(*zip(r1s, r2s)))
  return fastqs, controls

def reads_from_dirs(dirs):
  """Extract all R1/R2 files from a directory, ignores index (I1/I2) files"""
  # return list(itertools.chain(*map(lambda x: glob(os.path.join(x, '*_R[12]_*')), dirs)))
  r1s = list(itertools.chain(*map(lambda x: glob(os.path.join(x, '*_R1_*')), dirs)))
  r2s = list(itertools.chain(*map(lambda x: glob(os.path.join(x, '*_R2_*')), dirs)))
  return r1s, r2s

# run the main program


def main():
  args = docopt(__doc__, version='Version 1.0')

  # copy-pasted from pipeline.py :(
  cfg_f = args['--config']
  cfg_y = yaml.load(open(cfg_f))
  cfg = pipeline.Config(cfg_y)
  # it's probably better to have separate log files.
  if args['--log']:
    _log = Path(args['--log'])
    if _log.exists():
      print "Removing old log file %s" % _log
      _log.remove()
    log = open(args['--log'], 'a')
  else:
    log = sys.stdout

  sheet = open(args['<samplesheet>'])
  rows = csv.DictReader(sheet, fieldnames=['read_dirs', 'control_dirs'], delimiter='\t')
  rows=list(rows)
  base_out = args['--outdir'] or "." # os.getcwd()?
  sampledir = args['--sampledir']
  def p_run(rows):
      weave = partial(weave_files, sampledir)
      fqs_and_controls = list(map(weave, rows))
      run2_func = partial(pipeline.run2, cfg, log, base_out)
      pool = multiprocessing.Pool(len(rows))
      print "Launching %d processes.\n==========================\n\n" % len(rows)
      pool.map(run2_func, fqs_and_controls)
      pool.close()
      pool.join()
  if False: p_run(rows)
  if args['--qsub']:
    for i, row in enumerate(rows):
      #outdir = os.path.join(base_out, "sheet-sample-%d" % i)
      fastqs, controls = weave_files(sampledir, row)
      import tempfile
      import sh
      temp = tempfile.NamedTemporaryFile(prefix='pathos_sheet', suffix='qsub', delete=False)
      template = "{script} --fastq {fastqs} -c {cfg} -o {odir} --control {controls}"
      cmd = template.format(script='python /u/michael.panciera/CURRENT/pathos/pipeline.py',
                      fastqs=' '.join(fastqs),
                      controls=' '.join(controls),
                      cfg=args['--config'],
                      odir=base_out)
      temp.write("module load mpi\nmodule load bowtie\nmodule load blast\n")
      temp.write(cmd)
      temp.close()
      script = temp.name
      #print "qsub {script} -q batch -l nodes={node}:ppn={cores}".format(script=temp.name, node=amedpbswrair007.amed.ds.army.mil, cores=4)
      #print " -q batch -l nodes={node}:ppn={cores}".format(script=temp.name, node=amedpbswrair007.amed.ds.army.mil, cores=4)
      sample_num = row['read_dirs'].split(SEP)[0]
      sh.qsub(script,
              '-N',  "sheet-sample-%s" % sample_num,
              # "-M", "EMAIL HERE",
              # '-l', "nodes=1:ppn=8:mem=80514472881")
              '-l', "nodes=1:ppn=12",
              '-l', "mem=80514472881")
      print "Running %s" % script
  else:
    print "No --qsub flag, didn't run anything."
#  if p:
#
#  else:
#    for i, row in enumerate(rows):
#      cfg.outdir = os.path.join(base_out, "sheet-sample-%d" % i)
#      pipeline.run2(cfg, log, base_out, (fastqs, controls))


   # if controls:
   #     # TODO
   #   import sh
      #sh.python('pipeline.py', ' '.join(fastqs), '--control', ' '.join(controls), config=cfg_f, outdir=cfg.outdir)
   #   sh.python('pipeline.py', ' '.join(fastqs), '--config', cfg_f,  '--outdir', cfg.outdir, '--control', ' '.join(controls))
    #print 'python pipeline.py', ' '.join(fastqs), '--config', args['--config'], '-o', cfg.outdir, '--control', ' '.join(controls)



if __name__ == '__main__':
    main()
