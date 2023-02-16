'''
Usage:
    pipeline.py (--fastq <FASTQS>...) --config <config> [-o <outdir>] [--log <log>] [--control <CONTROLS>...]

Options:
    -f <FASTQS>, --fastq <FASTQS>
    -r <CONTROLS>, --control <CONTROLS>
    -c <config>, --config <config>
    --o <outdir>, --outdir <outdir>
    --log <log>, -l <log>
'''
import sh
from itertools import imap, izip, tee, chain, groupby, ifilter, starmap
from toolz.dicttoolz import keymap,valfilter,keyfilter,merge,assoc
from toolz.itertoolz import mapcat, partition
from funcy.seqs import _icut
from docopt import docopt
from path import Path
import types
import yaml
import sys
from functools import partial
import shutil
import os
from Bio import SeqIO
from operator import itemgetter as get
from collections import Counter
import csv
from ete2 import NCBITaxa
import plumbum
import pandas as pd
from schema import schema
from jsonschema import validate

#import numpy as np
#import matplotlib.pyplot as plt
# TODO: Log commands as doing them
# TODO: BLAST Contig results with mapped reads to get abundance:
# TODO: add contam optional arguments
#  - Does abyss record number of reads into contig? # no, it's just length and "kmer coverage"
#  - could duplicate contig blast entry for each read that maps to it and pass to krona

CONTAM_FLAG = 'CONTAM_MAPQ'
def lzw(sequence):
# https://github.com/betegonm/gen/blob/64aef21cfeefbf27b1e2bd6587c555d4df4f6913/gen.py#L294
  output = []
  table = dict(dict((chr(i), i) for i in range(256)))
  s = ''
  for ch in sequence:
    it = s + ch
    if it in table:
      s = it
    else:
      output.append(table[s])
      table[it] = len(table)
      s = ch
  output.append(table[s])
  return len(output)

#############
# Utilities #
#############
class Config:
    def __init__(self, entries):
        self.__dict__.update(entries)
        for k,v in self.__dict__.items():
            if type(v) == types.DictType:
                setattr(self, k, Config(v))
class Sh_(object):
    def __getattr__(self, attr):
        #cmd = getattr(sh, attr)
        def command(*args, **kwargs):
            #fixedargs = keymap("-{}".format, kwargs)
            bools = valfilter(lambda x: type(x) is bool, kwargs)
            vargs = keymap("-{}".format, keyfilter(lambda x: x not in ['_err', '_out'], valfilter(lambda x: not type(x) is bool, kwargs)))
            #bools.update(vargs)
            fixedargs = chain(vargs.items())
            getattr(sh, attr)(*(list(args) + list(fixedargs)), **bools)
        return command
sh_ = Sh_()

############
# Parts    #
############

def star(log, cfg, in1, in2):
  sh.STAR('--readFilesIn', in1, in2, #readFilesIn=unlist(in1, in2),
          genomeDir=cfg.star.starDB,
          outSAMtype="SAM",
          outReadsUnmapped="Fastx",
          outFileNamePrefix=os.path.join(cfg.outdir, STAR_PREFIX),
          _out=log, _err=log)

def pricefilter(log, cfg, in1, in2, o1, o2):
    cfg = cfg.pricefilter
    sh_.PriceSeqFilter('-fp', in1, in2,
                       '-op', o1, o2,
                       '-rqf', cfg.highQualPercent, cfg.highQualMin,
                       rnf=cfg.calledPercent)

def cdhitdup(log, cfg, r1, r2, o1, o2):
    sh_.cd_hit_dup(i=r1, i2=r2, o=o1, o2=o2, e=cfg.cdhitdup.minDifference, _err=log, _out=log)

# LZW!
def bowtie_sensitive(log, cfg, r1, r2, o1):
    args = {'1' : r1, '2' : r2,
                'very_sensitive_local' : True,
                'un_conc' : Path(o1).splitext()[0],
                'x' : cfg.bowtie2.bowtieDB,
            # '_err' : log, '_out' : log}
            '_err' : log}
    sh.bowtie2(**args)

def rapsearch(log, cfg, fq, out):
    out = out.splitext()[0] # rapsearch adds an m8 extension
    sh.rapsearch(o=out, d=cfg.rapsearch.rapsearchDB, q=fq, _err=log, _out=log)

def blastn(log, cfg, fq, out):
    print "attempting blast with %s %s" % (fq, out)
    #sh_.blastn(outfmt=6, db=cfg.ncbi.ntDB, query=fq, _err=log, _out=out)
    outformat = " ".join(["6"] + blast_columns)
    sh.blastn('-max_target_seqs', cfg.blastn.max_target_seqs, outfmt=outformat, db=cfg.ncbi.ntDB, query=fq, _err=log, _out=out, _long_prefix='-')

def krona(log, cfg, blast, out):
    sh.ktImportBLAST(blast, '-tax', cfg.ncbi.ktTaxonomy, o=out, _err=log, _out=log) # probably need config for kronadb!


def blastx(log, cfg, fq, out):
    outformat = " ".join(["6"] + blast_columns)
    sh.blastx('-max_target_seqs', '1', outfmt=outformat, db=cfg.ncbi.nrDB, query=fq, _err=log, _out=out, _long_prefix='-')

def abyss(log, cfg, r1, r2, out):
    out = Path(out)
    dir = out.dirname()
    f1 = dir.relpathto(r1)
    f2 = dir.relpathto(r2)
    prefix=out.basename().split('-')[0]
    #print f1, f2, 'name=%s' % prefix, 'k=%s' % 25
    #cmd("in='%s %s'" % (f1, f2),'name=%s' % prefix, 'k=%s' % 25, C=dir, _err=log, _out=log)
    # cmd.wait()
    abyss_log = os.path.join(dir, "abyss.log")
    kmer = cfg.assembly.kmer
    cmd_template = r"abyss-pe in='{f1} {f2}' name={prefix} k={kmer} -C {dir} > {abyss_log} 2>&1" 
    cmd = cmd_template.format(f1=f1,f2=f2, prefix=prefix, kmer=kmer, dir=dir, abyss_log=abyss_log)
    log.write('\nRunning:\n{command}\nwriting log to {abyss_log}'.format(command=cmd, abyss_log=abyss_log))
    os.system(cmd)
   

SEQID='qseqid'
#def lzw_filter(log, cfg, r1, r2, out1, out2):
def lzw_filter_single(min_complexity, x):
    un_comp_len = len(str(x.seq))
    comp_len = sum(imap(len, sh.gzip(f=True, _in=str(x.seq))))
    # sh.sh('lzw.sh', x) . . .
    complexity =  comp_len / float(un_comp_len)
    #print complexity
    return complexity >= min_complexity

def unzip(seq):
    t1, t2 = tee(seq)
    return imap(get(0), t1), imap(get(1), t2)

def filter_pair(func, r1, r2, o1, o2, format):
    CHUNK_SIZE = 100000
    fwd = SeqIO.parse(r1, format)
    rev = SeqIO.parse(r2, format)
    filtered = ((x, y) for (x, y) in izip(fwd, rev)
                if func(x) and func(y))
    chunks = _icut(False, CHUNK_SIZE/2, filtered)
    for chunk in chunks:
        with open(o1, 'a') as f1: # should this be 'a'?
            with open(o2, 'a') as f2: # 'a'? see above
                c1, c2 = unzip(chunk)
                SeqIO.write(c1, f1, format)
                SeqIO.write(c2, f2, format)
#    res1, res2 = unzip(filtered)
#    with open(o1, 'w') as f1:
#        with open(o2, 'w') as f2:
#            SeqIO.write(res1, f1, format)
#            SeqIO.write(res2, f2, format)

def lzw_filter_fastq(log, cfg, r1, r2, out1, out2):
    lzw_func = partial(lzw_filter_single, cfg.lzwfilter.minCompressionScore)
    filter_pair(lzw_func, r1, r2, out1, out2, 'fastq')

def sum_sam_by_ref(log, cfg, sam):
    res = sh.samtools.view(sam, S=True, F=260)
    count_d = {}
    contam_d  = {}
    for line in  res:
        fields = line.split('\t')
        ref, qname = fields[2], fields[0]
        count_d[ref] = count_d.get(ref, 0) + 1
        # contam_d[ref] = contam_d.get(ref, 0) + int(CONTAM_FLAG in qname)
        #contam_d[ref] = int(bool((qname.split(CONTAM_FLAG + '=')[-1])))
        contam_d[ref] = int((CONTAM_FLAG + '=') in qname)
    return count_d, contam_d

#    refs = imap(lambda x: x.split('\t')[2], res)
#    reduce(lambda z, x:  assoc(z, x[0],  \
#                               (z.get(x[0], (0, 0))[0] + 1, z.get(x[0], (0, 0))[1] \
#                               + (1 if CONTAM_FLAG in x[1] else 0))), xs, {})
#
#    return Counter(refs)

def readcounts_to_tsv(sam, out):
    counter, contam_counter = sum_sam_by_ref(None, None, sam)
    xs = []
    for ref, count in counter.iteritems():
        xs.append({SEQID : ref,
                  'read_count' : count,
                  'contam_percentage' : contam_counter[ref] / float(count) })
    with open(out, 'w') as o:
       writer = csv.DictWriter(o, [SEQID, 'read_count', 'contam_percentage'], delimiter='\t')
       writer.writeheader()
       for row in xs:
         writer.writerow(row)
       # writer = csv.writer(o, delimiter='\t')
       # writer.writerow([SEQID,'read_count'])
       # writer.writerows(counter.items())


#TODO:
# * use existing contig index
# * make sure there are contam reads as input
# * map with bowtie, same way as read_counts is generated
# * extract read_count via idxstats into a csv, using the SEQID label and "contam_mapped"
# * use the existing join function to join on SEQID after joining read_count
# * ! currently need to load bowtie and blast as modules


def fasta_add_metadata(fasta, tsv, fasta_out):
    meta = csv.DictReader(open(tsv), delimiter='\t')
    fa = list(SeqIO.parse(fasta, 'fasta'))
    ids = map(lambda x: x.id, fa)
    meta.sort(key=lambda x: ids.index(x[SEQID]))
    def annotate(d):
        del d[SEQID]
        return ' ; ' + ';'.join(dictmap("{0}={1}".format, d))
    for (rec, met) in izip(fa, meta):
        rec.id += annotate(meta)
    with open(fasta_out, 'w') as out:
            SeqIO.write(fa, out, 'fasta')
from StringIO import StringIO
# def bam_to_counter(bam):
def get_idxstats(bam):
    # assumes bam is sorted, and being indexed will speed up significantly
    sh.samtools.index(bam)
    idxstats_str = sh.samtools.idxstats(bam).stdout
    idxstats = csv.reader(StringIO(idxstats_str), delimiter='\t')
    return idxstats

def bam_to_mapped_tsv(bam, mapped_field, out):
    # needs to be sorted and indexed
    #idxstats = '\n'.join(idxstats)
    #idxstats_dict = csv.reader(StringIO(idxstats), delimiter='\t')
    # sequence name, sequence length, mapped, unmapped
    idxstats = get_idxstats(bam)
    name_and_mapped = [ { SEQID : x[0], mapped_field : x[2]}  for x in idxstats]
    with open(out, 'w') as o:
      writer = csv.DictWriter(o, [SEQID, mapped_field], delimiter='\t')
      writer.writeheader()
      writer.writerows(name_and_mapped)



#    ref_counts = dict([(ref, int(mapped)) for (ref, _, mapped, _) in idxstats_dict if ref != '*'])
#    return Counter(ref_counts)

def dup_blast(log, counter, blst, out):
    # counter, _ = sum_sam_by_ref(None, None, sam)
    log.write("Skipped Contigs:\n======\n")
    with open(out, 'w') as f:
        with open(blst, 'r') as blast:
            for k,_v in groupby(blast, lambda x: x.split('\t')[blast_columns.index('qseqid')]):
                if k not in counter:
                    #contig_length = next(v).split('\t')[3]
                    #log.write("Contig %s of length %s had no mapped reads.\n" % (k, contig_length))
                    log.write("Contig %s has no mapped reads.\n" % k)
                else:
                    #import ipdb; ipdb.set_trace()
                    def add_id(i, line):
                        # k is the contig ID
                        return '\t'.join(['%s-%d' % (k, i)] + line.split('\t')[1:])
                    v = list(_v)
                    repeat_count = (counter[k] / len(v)) or 1
                    # total number of duplicates will equal the read_count
                    repeated = v * repeat_count
                    with_unique_ids = list(starmap(add_id, enumerate(repeated)))

                    f.writelines(with_unique_ids) # with_unique_id = starmap(add_id, enumerate(v))
                    #f.writelines(list(v) * counter[k])
    log.write("======\n")



#                # because SeqIO is slow for writing single records
#        for (x, y) in izip(fwd, rev):
#            if func(x) and func(y):
#                '@
#                f1.write(x.format(form))
#                f2.write(y.format(form))



    #comp_len = sh.wc(sh.gzip(f=True, _in=str(x.seq)), c=True)

    #sh.run_abyss(r1, r2, 'name=%s' % prefix, 'k=%s' % 25, _err=log, _out=log)
#    args = ["in=%s %s" % (r1, r2), 'name=%s' % prefix, 'k=%s' % 25]
#    subprocess.call('abyss-pe ' + ' '.join(args))
    #subprocess.call('abyss-pe', ' '.join(args))
#    print ['abyss-pe'] + args
    #subprocess.call(['abyss-pe'] + args) #, stdout=log, stderr=log, shell=True)
#    ex = "abyss-pe in=%s %s" % (r1, r2) +  ' name=%s ' % prefix + ' k=%s' % 25
#    subprocess.call(ex, stdout=log, stderr=log, shell=True)

#    sh.abyss_pe('in=\'%s %s\'' % (r1, r2), 'name=%s' % prefix, 'k=%s' % 25, # '-n'  dryrun
#                _err=log, _out=log)
#    sh.abyss_pe("in='%s %s'" % (r1, r2), name=prefix, k=25, # '-n'  dryrun
#                _err=log, _out=log, _long_prefix='', _short_prefix='')

#    abyss-pe k=25 name=test     in='test-data/reads1.fastq test-data/reads2.fastq'


# taxid = 1056490

def blastdbcmd(**opts):
    cmd_opts = keymap('-{}'.format, opts).items()
    process = plumbum.local['blastdbcmd'][cmd_opts]
    print process
    for line in process.popen().iter_lines(retcode=None):
        yield line[0]

def get_taxid(db, seqids): # (Path, str) -> dict[str,str]
   #res = sh.blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % seqid, _long_prefix='-')
   max_ids = 1000/80
   if len(seqids) > max_ids:
      xs = [seqids[i*max_ids:(i+1)*max_ids] \
         for i in range(len(seqids) / max_ids)]
      xs.extend([seqids[sum(map(len, xs))-1:]])
   else: xs = [seqids]
   res = mapcat(lambda x: blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % ','.join(x)), xs)
   #res = blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % ','.join(xs))
   res = ifilter(bool, res)
   res = imap(lambda s: s.strip("'"), res)
   return dict(imap(unicode.split, res))

def taxonomy(ncbi, taxid):
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    def make_d(lineage, ranks, names):
        for lin in lineage:
            if ranks[lin] == 'no rank':
                continue
            yield (ranks[lin], names[lin])
    return dict(make_d(lineage, ranks, names))
def dictmap(f, d): return starmap(f, d.items())
# ['superkingdom','kingdom','superfamily','genus','family','order','class','phylum','species']
ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'superfamily', 'family', 'genus', 'species']
#@ blast_columns = "sseqid qlen pident length mismatch evalue qseqid sscinames scomnames sblastnames stitle staxids".split()
blast_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle sscinames scomnames sblastnames staxids".split()
# blast_columns = "sseqid qlen pident length mismatch evalue qseqid sscinames scomnames sblastnames stitle staxids".split()
csv_fields = blast_columns + ranks
# csv_fields = ['superkingdom','kingdom','superfamily','genus','qend','bitscore','family','evalue','gapopen','pid','send','order','class','phylum','alnlen','species','sseqid','qseqid','qstart','sstart']
csv_fields = blast_columns + ranks
def blast2summary_dict(db, blastpath, ete2_db): # (Path, Path) -> list[dict]

  """Reading in a blast output file, lookup all seqids to get taxids with a single blastdbcmd.
  Then, lookup the taxonomy using ETE2 via the taxid, and add that info to the blast info."""
  # rows = csv.DictReader(open(blastpath), delimiter='\t',fieldnames=[SEQID, 'sseqid','pid', 'alnlen','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
  rows = csv.DictReader(open(blastpath), delimiter='\t',fieldnames=blast_columns)
  rows = list(rows)
  seqids = map(get('sseqid'), rows)
  taxids = get_taxid(db, seqids)
  def get_gi(s):
    fields = s.split('|')
    if len(fields) > 1:
        return fields[1]
    else:
        raise ValueError("Seq ID %s is missing GI fields and '|'" % s)
  gis = imap(get_gi, seqids)
  #TODO: change matches to use something unique--not the TAXID! actually, why is it a dict
  # in the first place? it should be a list of dictionaries, and then map over
  # the dictionaries to merge them with the taxonomy info
  # this will replace the lines:
  # matches = . . .
  # items = . . .
  #matches = dict((taxids[gi], row) for gi, row in zip(gis,rows) if gi in taxids)
  ncbi = NCBITaxa(ete2_db) # downloads database and creates SQLite database if needed
 # items = dictmap(lambda tid,row: merge(row, taxonomy(ncbi, tid)), matches)
  matches = [assoc(row, 'taxid', taxids[gi]) for gi, row in zip(gis, rows) if gi in taxids]
  items = [merge(row1, taxonomy(ncbi, row1['taxid'])) for row1 in matches]
  res =  imap(partial(keyfilter, csv_fields.__contains__), items)
  return res
  #return {k:v for k, v in items if k in fields}

def blast2summary(db, blastpath, outpath, cfg): # (Path,Path,Path) -> None
    with_taxonomies = list(blast2summary_dict(db, blastpath, cfg.ete2_db))
    with open(outpath, 'w') as out:
       #writer = csv.DictWriter(out, head.keys(), delimiter='\t')
       writer = csv.DictWriter(out, csv_fields, delimiter='\t')
       #writer = csv.DictWriter(out, fieldnames=head.keys(), delimiter='\t')
       writer.writeheader()
       for row in with_taxonomies:
         writer.writerow(row)

def join_csv(fa, fb, on, out):
  # need to make sure that if it's not found in metadata, it still gets into
  # output, with 0 as the result, this is right-biased join or something
  # Right now, this works with how='outer'
  a = pd.read_csv(fa, sep='\t')
  b = pd.read_csv(fb, sep='\t')
  merged = a.merge(b, on=on, how='outer')
  merged.dropna(thresh=14, inplace=True)
  merged.to_csv(out, sep='\t', index=False)



from plumbum.cmd import samtools, gmap_build, gsnapl

def bam2fq(sam, outpath, pairflag):
  proc = samtools['view']['-f', pairflag, sam]
  #import ipdb; ipdb.set_trace()
  # TODO: fix
  # https://support.bioconductor.org/p/68040/
  # https://www.biostars.org/p/323175/
  # it would be good to test this too lol in a unit test
  # the contam number is just mapq
  with open(outpath, 'w') as outfile:
    #TODO: the following doesn't catch failures but errors with
    # a type error because line will be a list of [stdout, stderr] in the case
    # of error
    for line in proc.popen().iter_lines(retcode=None):
      if isinstance(line, list):
        line = line[0]
        if not line: continue
      if line.startswith('@'):
        continue
      x = line.split('\t')
      qname, ref, mapq, seq, qual = x[0], x[2], x[4], x[9], x[10]
      if ref != '*':
        qname += ":%s=%s" % (CONTAM_FLAG, mapq)
      entry = "@{}\n{}\n+\n{}\n".format(qname, seq, qual)
      outfile.write(entry)

# build index
# do mapping
#def map_contam(idxFolder, contams, sampR1, sampR2, sam, outR1, outR2):
#  idxName = "contam-idx"
#  #gmap_build['-d', "%s/%s" % (idxFolder, idxName), contR1, contR2]
#  gmap_build['-d', "%s/%s" % (idxFolder, idxName), *contams)]
#  gsnapl['-d', idxFolder, '-d', idxName, sampR1, sampR2, "-A", "sam", "-o", sam]
#  bam2fq(sam, outR1, '64')
#  bam2fq(sam, outR2, '128')

def map_contam(log, idxFolder, contams, sampR1, sampR2, sam, outR1, outR2):
#  idxName = "contam-idx"
#  sh.gmap_build('-d', "%s/%s" % (idxFolder, idxName), *contams)
#  sh.gsnapl('-d', idxFolder, '-d', idxName, sampR1, sampR2, "-A", "sam", "-o", sam)
#  bam2fq(sam, outR1, '64')
#  bam2fq(sam, outR2, '128')

   sh.bowtie2_build(','.join(contams), idxFolder)
   sh.bowtie2(**{'1' : sampR1, '2' : sampR2, 'x' : idxFolder, '_err' : log, '_out' : sam})
   bam2fq(sam, outR1, '64')
   bam2fq(sam, outR2, '128')

def seq_reaches_length(length, seq):
  return len(seq.seq) >= length

def filter_contigs(min_length, inc, outc):
  raw_seqs = SeqIO.parse(inc, 'fasta')
  filtered_seqs = ifilter(partial(seq_reaches_length, min_length), raw_seqs)
  with open(outc, 'w') as out:
    #SeqIO.write(filtered_seqs, out, 'fasta')
    for seq in filtered_seqs:
      if not seq.id.startswith('c'):
        seq.id = 'c' + seq.id
      SeqIO.write(seq, out, 'fasta')

def merge_fastqs(fastqs, r1, r2):
  if len(fastqs) > 2:
    x1s, x2s = fastqs[0::2], fastqs[1::2]
    with open(r1, 'wb') as o1:
      with open(r2, 'wb') as o2:
        for x1, x2 in zip(x1s, x2s):
          with open(x1, 'rb') as i1:
            shutil.copyfileobj(i1, o1) #, 1024*1024*10)
          with open(x2, 'rb') as i2:
            shutil.copyfileobj(i2, o2)
  else:
    r1, r2 = fastqs[0], fastqs[1]
  return r1, r2

############
# Pipeline #
############
STAR_PREFIX = 'STAR_'
def skip(input1, input2, o1 ,o2):
    shutil.copy(input1, o1)
    shutil.copy(input2, o2)

def run(cfg, input1, input2, contams, log=None):
  p = partial(os.path.join, cfg.outdir)
  _star1 = p("%sUnmapped.out.mate1" % STAR_PREFIX)
  _star2 = p("%sUnmapped.out.mate2" % STAR_PREFIX)
  star1 = p("star.r1.fq")
  star2 = p("star.r2.fq")
  psf1 =  p( "psf.r1.fq"             )
  psf2 =  p( "psf.r2.fq"             )

  cd1 =       p( "cd.r1.fq" )
  cd2 =       p( "cd.r2.fq" )

  lzw1 = p("lzw.r1")
  lzw2 = p("lzw.r2")

  _bowtie1 =   p( "bowtie.1.1" )
  _bowtie2 =   p( "bowtie.2.1" )

  unfiltered_contigs = p("abyss-contigs.fa")
  contigs = p("filtered-contigs.fa")
  contigs_sam = p('contigs.sam')

  ray_bam =  p('RAYOUT/bowtie2_mapping/out.bam')

  contig_nr = p('contigs.nr.blast')
  contig_nt = p('contigs.nt.blast')

  contig_nt_flagged = p('contigs.nt.flagged.tsv')

  dup_nt = p('contigs.nt.blast.dup')
  dup_nr = p('contigs.nr.blast.dup')
  contig_kronaNT = p('contigs.nt.html')
  contig_kronaNR = p('contigs.nr.html')
  contig_kronaNT_dup = p('contigs.nt.DUP.html')
  contig_kronaNR_dup = p('contigs.nr.DUP.html')

  contig_nt_tsv = p("contigs.nt.tsv")
  contig_nr_tsv = p("contigs.nr.tsv")
  nt_tsv = p('nt.tsv')
  nr_tsv = p('nr.tsv')
#  bowtie1 =   p( "bowtie.r1.fa" )
#  bowtie2 =   p( "bowtie.r2.fa" )
#  nr1     =   p( "rapsearch.r1.blast.m8" ) # rapsearch automatically adds .m8 extension
#  nr2     =   p( "rapsearch.r2.blast.m8" ) # rapsearch automatically adds .m8 extension
#
#  nt1 =       p( "r1.blast" )
#  nt2 =       p( "r2.blast" )
#  kronaNT1  = p( "r1.NT.html" )
#  kronaNT2  = p( "r2.NT.html" )
#  kronaNR1  = p( "r1.NR.html" )
#  kronaNR2  = p( "r2.NR.html" )
  with_tax_tsv = p('contigs.blast.tax.tsv')
  contigs_meta = p('contigs.meta.tsv')

  if not log:
    log = sys.stdout

  def logtime(s):
      from time import gmtime, strftime
      time = strftime("%Y-%m-%d %H:%M:%S", gmtime())
      message = "{}  {} started\n\n".format(time, s)
      log.write(message)

  logtime('star')
  if not cfg.star.skip:
    if need(_star1):
      star(log, cfg, input1, input2)
  
    if need(star1):
      shutil.copy(_star1, star1)
      shutil.copy(_star2, star2)
  else:  # building the STAR database requires a lot of memory. so does running it!
    if need(star1):
      shutil.copy(input1, star1)
      shutil.copy(input2, star2)

  logtime('priceseqfilter')
  if need(psf1):
    pricefilter(log, cfg, star1, star2, psf1, psf2)



  logtime('cd-hit-dup')
  if need(cd1):
    cdhitdup(log, cfg, psf1, psf2, cd1, cd2)

  logtime('lzw_filter_fastq')
  if need(lzw1):
    if not cfg.lzwfilter.minCompressionScore:
      skip(cd1, cd2, lzw1, lzw2)
    else:
      lzw_filter_fastq(log, cfg, cd1, cd2, lzw1, lzw2)

  logtime('bowtie_sensitive')
  if need(_bowtie1):
    bowtie_sensitive(log, cfg, lzw1, lzw2, _bowtie1)

  marked1 =   p( "contam_marked.R1.fastq" )
  marked2 =   p( "contam_marked.R2.fastq" )
  contam_sam =  p("contam.sam")
  #if contR1 and contR2:
  if need(marked1):
    if contams:
      idxFolder = p("contamIdx")
      #if not os.path.exists(idxFolder): os.mkdir(idxFolder)
      fasta_contams = []
      for c in contams:
        fasta_name = str( p( Path(c).basename().stripext() + '.fasta' )  )
        # gmap__build requires fasta, not fastq to build index from
        SeqIO.convert(c, 'fastq', fasta_name, 'fasta')
        fasta_contams.append(fasta_name)
      #map_contam(idxFolder, contR1, contR2, _bowtie1, _bowtie2, contam_sam, marked1, marked2)
      map_contam(log, idxFolder, fasta_contams,
                 _bowtie1, _bowtie2, contam_sam, marked1, marked2)

    else: # this isn't actually necessary; could use the _bowtie1 since they get absorbed in assembly anyway
      shutil.copy(_bowtie1, marked1)
      shutil.copy(_bowtie2, marked2)


  logtime('assembly')
  if need(contigs):
    if cfg.assembly.assembler == 'ray2':
      import subprocess
      subprocess.check_call(' '.join(['ray_script', str(marked1), str(marked2), str(unfiltered_contigs), str(contigs_sam), cfg.param_file]), shell=True)
      #sh.ray_script(marked1, marked2, contigs, contigs_sam)
    elif cfg.assembly.assembler == 'abyss':
      abyss(log, cfg, marked1, marked2, unfiltered_contigs)
    else:
      raise ValueError("Config Assembler %s not supported" % cfg.assembly.assembler)
    filter_contigs(cfg.assembly.minimum_contig_length, unfiltered_contigs, contigs)
  contigs_index = p('contigs-b2')
  if need(contigs_sam): # shouldn't happen w/ ray because ray_script copies over a bam
    # NOTE: used to be a mapping onto unfiltered contigs
    # sh.bowtie2_build(unfiltered_contigs, contigs_index)
    sh.bowtie2_build(contigs, contigs_index)
    sh.bowtie2(**{'1' : marked1, '2' : marked2, 'x' : contigs_index,
                  '_err' : log, '_out' : contigs_sam})
    # TODO: refactor below so that it works for R1 and NT. probably just drop
    # R2. what does pathdiscov do?

  logtime('blastn')
  if need(contig_nt):
    blastn(log, cfg, contigs, contig_nt)
  logtime('dup_nt')
  if need(dup_nt):
    counter, _ = sum_sam_by_ref(None, None, contigs_sam)
    dup_blast(log, counter, contig_nt, dup_nt)


#  logtime('blastx')
#  if need(contig_nr):
#    blastx(log, cfg, contigs, contig_nr)
#    dup_blast(log, contigs_sam, contig_nr, dup_nr)

  logtime('krona')
  if need(contig_kronaNT):
    krona(log, cfg, contig_nt, contig_kronaNT)
  if need(contig_kronaNT_dup):
    krona(log, cfg, dup_nt, contig_kronaNT_dup)
#  if need(contig_kronaNR):
#    krona(log, cfg, contig_nr, contig_kronaNR)

  logtime('blast2summary')
  if need(with_tax_tsv):
    blast2summary(cfg.ncbi.ntDB, contig_nt, with_tax_tsv, cfg) # this is shared by read files
    # the below is for contigs only.

  if need(contigs_meta):
      # can check for CONTAM flag
    readcounts_to_tsv(contigs_sam, contigs_meta)

  contig_nt_almost = p("contigs.nt.almost")
  if need(contig_nt_almost):
      # joins on SEQID, which collapses data.
    import diff_ranks
    diff_ranks.flag_annotated_blast(with_tax_tsv, contig_nt_flagged)
    join_csv(contigs_meta, contig_nt_flagged, SEQID, contig_nt_almost)

  if contams:
    contam_ctg_bam = p("control_to_contig.bam")
    if need(contam_ctg_bam):
      contam_ctg_sam = p("control_to_contig.sam")
      cr1, cr2 = (p("control_merged_R1.fq"), p("control_merged_R2.fq"))
      cr1, cr2 = merge_fastqs(contams, cr1, cr2)
        # actually need to merge contams too
      #TODO: remove the line below building the contigs_index again
      sh.bowtie2_build(contigs, contigs_index)
      sh.bowtie2(**{'1' : cr1, '2' : cr2, 'x' : contigs_index,
                    '_err' : log, '_out' : contam_ctg_sam})
      sh.samtools.sort(contam_ctg_sam, o=contam_ctg_bam)
      sh.samtools.index(contam_ctg_bam)
    contam_mapped_tsv = p("contam_mapped.tsv")
    bam_to_mapped_tsv(contam_ctg_bam, 'control_mapped', contam_mapped_tsv)
    join_csv(contig_nt_almost, contam_mapped_tsv, SEQID, contig_nt_tsv)
  else:
    shutil.copy(contig_nt_almost, contig_nt_tsv)

  logtime('finished!')
#  def mksummary(db, blast, s):
#    with_tax_tsv = '%s.blast.tax' % s
#    contigs_meta = '%s.meta' % s
#    blast2summary(db, blast, with_tax_tsv)
#    readcounts_to_tsv(contigs_sam, contigs_meta)
#    join_csv(contigs_meta, with_tax_tsv, SEQID, contig_nt_tsv)

    # TODO: fix below
    #fasta_add_metadata(contigs, contigs_meta, contigs_with_meta)

#  if need(bowtie1):
#    SeqIO.convert(_bowtie1, 'fastq', bowtie1, 'fasta')
#    SeqIO.convert(_bowtie2, 'fastq', bowtie2, 'fasta')
#
#  if need(nt1):
#    blastn(log, cfg, bowtie1, nt1)
#    blastn(log, cfg, bowtie2, nt2)
#
#  if need(nr1):
#    #rapsearch(log, cfg, bowtie1, nr1)
#    #rapsearch(log, cfg, bowtie2, nr2)
#    blastx(log, cfg, bowtie1, nr1)
#    blastx(log, cfg, bowtie2, nr2)
#
#  if need(kronaNT1):
#    krona(log, cfg, nt1, kronaNT1)
#    krona(log, cfg, nt2, kronaNT2)
#    krona(log, cfg, nr1, kronaNR1)
#    krona(log, cfg, nr2, kronaNR2)

need = lambda p: not os.path.exists(p)

def run2(cfg, log, base_out, fastqs_and_controls):
  '''If there is more than one pair of input files, merge them into a single pair (two files). r1 will have all the R1 files, r2 will have all the R2 files, in the order that they were passed to the commandline.
  NOTE: this requires that they were passed in the correct way!
  fastqs_and_controls is a tuple of the flie-lists'''
  fastqs, controls = fastqs_and_controls
  controls = controls if controls else []
  fastqs, controls = list(fastqs), list(controls)
  print 'FASTQS : %s' % fastqs, 'CONTROLS : %s' % controls
  assert len(fastqs) > 1 and (len(fastqs) % 2) == 0, "Must provide an even number of input files."
  vdbpm_id = Path(fastqs[0]).splitall()[-2]
  #print fastqs, control
  outdir = os.path.join(base_out, "sheet-%s" % vdbpm_id)
  cfg2 = Config(cfg.__dict__) # I'm not sure if this obj is shared with other processes
  cfg2.outdir = outdir
  p = partial(os.path.join, outdir)
  if not os.path.exists(outdir): os.makedirs(outdir)
  r1 = p("input_merged.r1.fq")
  r2 = p("input_merged.r2.fq")
  if need(r1) or need(r2):
    if len(fastqs) > 2:
      #if not os.path.exists(cfg.outdir): os.mkdir(cfg.outdir)
      x1s, x2s = fastqs[0::2], fastqs[1::2]
      with open(r1, 'wb') as o1:
        with open(r2, 'wb') as o2:
          for x1, x2 in zip(x1s, x2s):
            with open(x1, 'rb') as i1:
              shutil.copyfileobj(i1, o1) #, 1024*1024*10)
            with open(x2, 'rb') as i2:
              shutil.copyfileobj(i2, o2)
    else:
      f1 = os.path.abspath(fastqs[0])
      f2 = os.path.abspath(fastqs[1])
      os.symlink(f1, r1)
      os.symlink(f2, r2) 
      # r1, r2 = fastqs[0], fastqs[1]
  run(cfg2, r1, r2, controls, log)



def arg_parser():

  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('--fastq', nargs='+', required=True)
  parser.add_argument('-c', '--config', type=str, required=True)
  parser.add_argument('-o', '--outdir', type=str, required=True)
  parser.add_argument('--control', nargs='*', required=False)
  parser.add_argument('-l', '--log', type=str)
  return parser

def main():
  parser = arg_parser()
  args = parser.parse_args()
  #args = docopt(__doc__, version='Version 1.0')

  cfg = args.config
  cfg = yaml.load(open(cfg))
  validate(cfg, schema)
  cfg = Config(cfg)
  cfg.outdir = args.outdir or "."
  if args.log:
    _log = Path(args.log)
    if _log.exists():
      print "Removing old log file %s" % _log
      _log.remove()
    log = open(args.log, 'a')
  else:
    log = sys.stdout

  fastqs = args.fastq
  #c1, c2 = args['<c1>'], args['<c2>']
  #run(cfg, args['<r1>'], args['<r2>'], c1, c2, log)
  #TODO: fix . . .
  controls = args.control
  print args

  run2(cfg, log, cfg.outdir, (fastqs, controls))
#def run2(cfg, log, base_out, fastqs_and_controls):
  #run(cfg, r1, r2, controls, log)

#  try:
#    run(cfg, args['<r1>'], args['<r2>'], log)
#  except Exception as e:
#    log.write(str(e))
  if args.log: log.close()
  sys.exit(0)

if __name__ == '__main__': main()




