from Bio import SeqIO
from differences import diff
from sys import argv

prot_name: str = argv[1]
out = open(prot_name + ".tsv", "w")
out.write("index\tstart\tend\tref\talt\tlocation\tdate\tproteome_id\n")

records = SeqIO.parse(f"{prot_name}.fasta", "fasta")
ref = next(records)
ref_seq: str = ref.seq

gen_index = 1
mut_index = 1
for record in records:
    annotations = record.name.split("|")
    meta = ()
    try:
        meta = (annotations[1].split("/")[1], annotations[2], f"{prot_name}{gen_index}")
        gen_index += 1
    except Exception:
        continue

    tups = diff(ref_seq, record.seq)
    for tup in tups:
        out.write("\t".join([str(mut_index), str(tup[0]), str(tup[1]), tup[2], tup[3], *meta]) + "\n")
        mut_index += 1

out.close()
