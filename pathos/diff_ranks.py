from operator import itemgetter as get
from collections import Counter
#import pandas as pd
import csv
from io import StringIO

# need to fix so that we join on conting ID. currently joining on something else.
test_string = u'''qseqid	read_count	contam_percentage	superkingdom	kingdom	superfamily	genus	qend	bitscore	family	evalue	gapopen	pid	send	order	class	phylum	alnlen	species	sseqid	qstart	sstart
15368	2.0	0.0	Bacteria			Bacillus	1.0	5e-08	Bacillaceae	4605131.0	3.0	93.18	4605174.0	Bacillales	Bacilli	Firmicutes	44.0	Bacillus cereus	gi|300373910|gb|CP001746.1|	0.0	44.0
15368	2.0	0.0	Bacteria			Bacillus	1.0	5e-08	Bacillaceae	4605131.0	3.0	93.18	4605174.0	Bacillales	Bacilli	Firmicutes	44.0	Bacillus cereus	gi|300373910|gb|CP001746.1|	0.0	44.0
4142	2.0	0.0	Bacteria			Bacillus	3.0	1e-14	Bacillaceae	3013062.0	0.0	100.0	3013016.0	Bacillales	Bacilli	Firmicutes	47.0	Bacillus kochii	gi|1237941648|gb|CP022983.1|	0.0	49.0
16708	6.0	0.0	Eukaryota			Skeletonema	5.0	2e-54	Skeletonemataceae	2857.0	6.0	95.65	2720.0	Thalassiosirales	Coscinodiscophyceae	Bacillariophyta	138.0	Skeletonema costatum	gi|1179910890|dbj|LC258372.1|	0.0	142.0
13020			Bacteria			Enterobacter	1.0	1e-07	Enterobacteriaceae	4603241.0	0.0	100.0	4603273.0	Enterobacteriales	Gammaproteobacteria	Proteobacteria	33.0	Enterobacter cloacae	gi|1395322938|gb|CP021167.1|	0.0	33.0
1476			Viruses			Flavivirus	1.0	7e-12	Flaviviridae	212.0	2.0	95.83	165.0				48.0	Kokobera virus	gi|1098496923|gb|KU059115.1|	0.0	48.0
'''

def test_differing_rank():
  test_csv = StringIO(test_string)
  #rows = pd.read_csv(test_csv, sep='\t')
  rows = list(csv.DictReader(test_csv, delimiter='\t'))
  match_rows = rows[:2]
  assert differing_rank(match_rows) == 'MATCH'

  species_rows = rows[:3]
  assert differing_rank(species_rows) == 'species'

  kingdom_rows = rows[0], rows[1], rows[2], rows[4]
  # some columns are empty: in this case kingdom. so with that unknown, we go down to phylum.
  assert differing_rank(kingdom_rows) == 'phylum', differing_rank(kingdom_rows) # map(get('kingdom'), kingdom_rows)

  superkingdom_rows = rows
  assert differing_rank(superkingdom_rows) == 'superkingdom'
  return "\ntest_differing_rank SUCCESS\n"


'''
Find the smallest rank where the blast results DO NOT match.
'''
def differing_rank(rows):
  # no domain
  ranks = reversed(['superkingdom',	'kingdom',	 'phylum',\
   'class', 'order', 'superfamily', 'family',	'genus', 'species'])
  differ_at = 'MATCH'
  for rank in ranks:
    row_ranks = map(get(rank), rows)
    uniq = set(row_ranks)
    if len(uniq) > 1:
      differ_at = rank
      # none of the below are used
      popular_rank, popular_count = Counter(row_ranks).most_common(1)[0]
      percent_match = popular_count / float(len(rows))

  return differ_at

# TODO
# groupby blast input by field 1 (contig name=c1, c2 .  . . .
# pick the top blast hit (should be first entry
# keep the info relevant to the top hit (pid, tax, etc.)
# collapse diff_ranks
# add diff_ranks as a field
#
# todo some point: normalize contig names
# meta tsv wasn't created ?
# readcount is a float somehow and it shouldn't be
# readcount is coming out blank sometimes in contigs.nt.tsv
# RAYOUT/contig_numreads.txt is not being used -- can't use it
# because it forgets what reads contribute to a contig
# and we need to calculate the contam-percentage
# should we just map the contaminants back to the contig?
# I'm getting contam=0 probably where I should get no flag at all.
#      if ref != '*':
#        qname += ":%s=%s" % (CONTAM_FLAG, mapq)
#  contam.sam is missing

from toolz.dicttoolz import merge
import itertools
from itertools import groupby
def flag(group_obj):
  group = list(group_obj[1])
  rank = differing_rank(group)
  add_rank = lambda x: merge(x, {'different_rank' : rank})
  return map(add_rank, group)

def flag_annotated_blast(input_fn, output_fn):
  with open(input_fn) as f_in, open(output_fn, 'w') as f_out:
    rows = list(csv.DictReader(f_in, delimiter='\t'))
    groups = groupby(rows, lambda x: x['qseqid'])
    flagged_rows = list(itertools.chain(*map(flag, groups)))
    writer = csv.DictWriter(f_out, flagged_rows[0].keys(), delimiter='\t')
    writer.writeheader()
    for row in flagged_rows:
      writer.writerow(row)



if __name__ == '__main__':
  print test_differing_rank()
