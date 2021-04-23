famd = fastaread('rcsb_pdb_1BBT.fasta')

seq1 = getfield(famd(1,1),'Sequence')
seq2 = getfield(famd(2,1),'Sequence')
seq3 = getfield(famd(3,1),'Sequence')
seq4 = getfield(famd(4,1),'Sequence')

data(1).Sequence = aa2nt(seq1);
data(1).Header = getfield(famd(1,1),'Header');
data(2).Sequence = aa2nt(seq2);
data(2).Header = getfield(famd(2,1),'Header');
data(3).Sequence = aa2nt(seq3);
data(3).Header = getfield(famd(3,1),'Header');
data(4).Sequence = aa2nt(seq4);
data(4).Header = getfield(famd(4,1),'Header');

fastawrite('nuc_seq.fasta', data)


