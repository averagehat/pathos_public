'''
Usage:
    summary.py <indir> <outdir>
'''
from docopt import docopt
import pandas as pd
from operator import add
# from itertools import accumulate # python 3 :(
import os
from StringIO import StringIO
from functools import partial
from Bio import SeqIO
from matplotlib import pyplot as plt
import sys
import numpy as np

LCD_RANK_FIELD_NAME = 'Lowest Common Rank'

def accumulate(iterable, func):
    it = iter(iterable)
    try:
        total = next(it)
    except StopIteration:
        return
    yield total
    for element in it:
        total = func(total, element)
        yield total

def n50(lens):
    half = sum(lens) / 2.0
    lns = sorted(lens, reverse=True)
    accumulated = accumulate(lns, add)
    whereOverHalf = [lns[i] for i,acc in enumerate(accumulated) if acc > half]
    return whereOverHalf[0]

DUMB_INDEX = [0]
# TODO: size doesn't do what I think it does!
def summ_contigs(tops):
  return pd.DataFrame(index = DUMB_INDEX, data={ 'N50' : n50(tops.qlen),
                       'contig_count' : len(tops),
                       'assembly_length' :  tops.qlen.sum(),
                       'species_count' : tops.species.unique().size })


def count(gen):
  return sum(1 for _ in gen)

def readCount(fp, format='fastq'):
    gen = SeqIO.parse(open(fp), format)
    return sum(1 for _ in gen)

#TODO: this needs to know about filepaths;
#
def summarize_nontaxon(files, indir):
  p = partial(os.path.join, indir)
  unfiltered_contigs = p(files.unfiltered_contigs)
  input1 = p(files.input1)
  star1 = p(files.star1)
  _bowtie1 = p(files._bowtie1)
  lzw1 = p(files.lzw1)
  contigs = p(files.contigs)
  contig_top_summary = p(files.contig_nt_tsv)
  # really all this needs is the filenames.
  # could return tops and summ_df to avoid any i/o here
  summ_df = pd.DataFrame(index = DUMB_INDEX)
  # actuallyt hese are totaly independent and should be treated as such
  summ_df['unfiltered_contig_count'] = readCount(unfiltered_contigs, 'fasta')
  summ_df['sample_input'] = readCount(input1)
  star_host_count = summ_df['sample_input'] - readCount(star1)
  bowtie_host_count = readCount(lzw1) - readCount(_bowtie1)
  summ_df['host_reads'] = star_host_count + bowtie_host_count
  # below is filtered value.
  summ_df['unblasted_contigs'] = readCount(contigs, 'fasta') - count(contig_top_summary) - 1 # for header lin
  return summ_df
  # return None

# type: rank
descendingRanksWithNoMatch = ['NoMatch', 'superkingdom', 'kingdom', 'phylum','class', 'order', \
                   'superfamily', 'family',   'genus', 'species']
descendingRanks = descendingRanksWithNoMatch[1:]
    # won't work b/c maximum not have "nomatch" as a return value
def max_rank(ranks): # -> rank
    min_idx =  max(map(descendingRanksWithNoMatch.index, ranks))
    return descendingRanksWithNoMatch[min_idx]

def matchedRank(rows, rank): # -> rank
    hasMatch = rows[rank].nunique() == 1
    return rank if hasMatch else "NoMatch"

def lcdRank(rows): # -> rank
    return max_rank(map(partial(matchedRank, rows), descendingRanks))

#TODO:  how does this handle e.g. NaN???
#TODO: take the current "summary file" as input and overwrite the different_ranks field, and give:
 # lowest_common_rank
 # top BLAST result based on maximum PID
# print(lcdRank(pd.DataFrame([dict(superkingdom='foo'), dict(superkingdom='bar')], columns=descendingRanks)))
# print(lcdRank(pd.DataFrame([dict(superkingdom='bar', species='foo'), dict(superkingdom='bar', species='baz')], columns=descendingRanks)))
# print(lcdRank(pd.DataFrame([dict(species='foo'), dict(species='foo')], columns=descendingRanks)))
# print(lcdRank(pd.DataFrame([dict(species='bar', genus='baz'), dict(species='foo', genus='baz')], columns=descendingRanks)))
#
def to_tops(df):
  gbs = df.groupby('qseqid')
  lcd_ranks = gbs.apply(lcdRank)
  tops = df.sort_values('pident', ascending=False).groupby('qseqid').first()
  tops[LCD_RANK_FIELD_NAME] = lcd_ranks
  tops.reset_index(inplace=True)
  return tops


#NOTE: Text overlaps in e.g. kingdom, when most only fit into one.
#NOTE: This doesn't give info when the parent rank is under-represented,
# i.e. it won't give us the top 5
def plot_rank(group, rank, rankpath):
    read_sorted =  group.read_count.sort_values(ascending=False)
    tooBig = read_sorted.size > 7
    if tooBig:
        others = pd.Series(read_sorted[5:].sum() , index=['others'])
        to_plot = pd.concat([read_sorted[:5], others])
    else:
        to_plot = read_sorted
    # to_plot.plot.pie(y='read_count')
    to_plot.plot.pie()
    #plt.ylabel(rank)
    plt.ylabel('') # remove to make more space for labels
    plt.savefig(rankpath + '.png')
    plt.close() # otherwise plot will look weird with "others" apparently not working and lots of labels


def make_rank_summaries(df, rankdir):
  if not os.path.exists(rankdir):
      os.mkdir(rankdir)
  rank_dfs = []
  ranks_str = "SuperKingdom | Kingdom | Phylum | Class | Order | SuperFamily | Family | Genus | Species"
  ranks = list(map(str.lower, map(str.strip, ranks_str.split('|'))))
  # TODO: do the below also for the whole of the top contigs.
  for rank in ranks:
    gb = df.groupby(rank)
    aggregated = gb.agg({'qseqid' : lambda x: ','.join(set(x)), # 'qlen' : sum,
                      'staxids' : lambda x: ','.join(np.vectorize(str)(x.unique())),
                       'read_count' : sum })
    # the value is a function on the series of just that thing.
    N50s = gb.qlen.apply(n50)
    aggregated['N50'] = N50s
    aggregated['contig_count'] = gb.size() # no this isn't correct, is it?
    aggregated['assembly_length'] = gb.qlen.sum()
    aggregated['species_count'] = gb.species.unique().agg(len)
    rankpath = os.path.join(rankdir, rank)
    columns = ['N50', 'contig_count', 'assembly_length', 'read_count', 'species_count',  'staxids', 'qseqid']
    aggregated.to_csv(rankpath + '.tsv', sep='\t', columns=columns)
    plot_rank(aggregated, rank, rankpath)
    aggregated['taxon'] = rank
    rank_dfs.append(aggregated)
  big_df = pd.concat(rank_dfs)
  big_summary_fp = os.path.join(rankdir, 'ranksummary.tsv')
  big_df.to_csv(big_summary_fp, sep='\t', columns= (columns + ['taxon']))
  return None

in_fp = 'contigs.15383.nt.tsv'
# def summarize(in_summary_fp, outdir):
import filenames
def summarize(indir, outdir):
  if not os.path.exists(outdir):
     os.mkdir(outdir)
  in_summary_fp = os.path.join(indir, filenames.contig_nt_tsv)
  in_df = pd.read_csv(open(in_summary_fp), sep='\t')
  rankdir = os.path.join(outdir, 'ranks/')
  make_rank_summaries(in_df, rankdir)
  tops = to_tops(in_df)
  p = partial(os.path.join, outdir)
  tops_fp = p('top_summary.tsv')
  tops.to_csv(tops_fp, sep='\t') # do I need to specify columns?
# summ_df can't be done without knowledge of other files.

  summ_df_ctg = summ_contigs(tops)
  summ_df_reads = summarize_nontaxon(filenames, indir)
  summ_df = summ_df_ctg.join(summ_df_reads)
  summ_fp = p('run_summary.tsv')
  #TODO: below writes the empty index in the csv. nah, fixed with index=False
  summ_df.to_csv(summ_fp, sep='\t', index=False)

def main():
  args = docopt(__doc__, version='Version 1.0')
  indir = args['<indir>']
  outdir = args['<outdir>']
  summarize(indir, outdir)

if __name__ == '__main__':
  main()
  # summarize(sys.argv[1], sys.argv[2])

# need to rename qlen->assemblyLength, add:
 # mean conting length
 # unblasted contig count
 #
 # possibly unfiltered stats




# don't really want sample->contaminant mapping count. We'll stop using it prob.
# STAR produces a sam file
# calculate number of host reads by subtracting readcounts from each other, or via samtools idxstats.
# then adding the results from the two host-mapping stages
# assumes its sorted/indexed. must be a BAM.
# get unmapped reads from 'unmapped' column of last (ref='*') row.

def get_idx_stats(indexed_bam):
  import sh
  idxstats_cmd=sh.samtools.idxstats(indexed_bam)
  idxstats_str = '\n'.join(idxstats_cmd)
  fields=['ref', 'ref_length', 'mapped', 'unmapped']
  idx_stats = pd.read_csv(StringIO(idxstats_str), sep='\t', header=None, names=fields)
  return idx_stats



