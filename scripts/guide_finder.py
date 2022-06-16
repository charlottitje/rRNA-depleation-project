import argparse

parser = argparse.ArgumentParser(description='combines barrnap rRNA fasta files from multiple species in a single file to streamline ribocutter input and define a name for the combination of genomes for futre refference')

parser.add_argument("-i", type=str, help='path to rRNA fasta file')
parser.add_argument("-o",type=str, help='path to output fasta file')

args = parser.parse_args()

all_guides=[]

def rev_c(seq):
    """
    simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    seq = seq.translate(tab)
    return seq

def find_guides(seq):
    """
    function that finds all the NGG PAM sites and extends the guide sequence for that site
    """
    gg_pos = [i for i, s in enumerate(seq) if seq[i:i + 2] == "GG" and i >= 21]  # not 22, because 0 based
    guides = []
    for i in gg_pos:
        guides.append(seq[i - 20: i +2:])
    return guides

with open(args.i) as handle:
    content=handle.readlines()
    sequences= [content[l][:-1].upper() for l in range(1, len(content), 2)]
    #print(sequences)
    for s in sequences:
        all_guides.extend(find_guides(s))
        all_guides.extend(find_guides(rev_c(s)))
unique_guides= list(set(all_guides))
        

with open(args.o, "w") as outfile:
    for guideN in range(0, len(unique_guides)):
        outfile.writelines(f">guide{guideN}\n")
        outfile.writelines(unique_guides[guideN]+"\n")
