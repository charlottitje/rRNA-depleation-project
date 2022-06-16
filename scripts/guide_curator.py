import glob
import copy
import argparse


parser = argparse.ArgumentParser(description='combines barrnap rRNA fasta files from multiple species in a single file to streamline ribocutter input and define a name for the combination of genomes for futre refference')

parser.add_argument("--gff", type=str, help='path to dictionary of gff files')
parser.add_argument("--sam",type=str, help='path to dictionary of sam alingment files')
parser.add_argument("--set", type=str, help='desired name of the guide set')

args = parser.parse_args()



######
#function definitions

def short_species_name(organism:str):
    """
    get the shortend name of a species
    aspegillus_niger -> a. niger
    
    input
    -----
    organism:str
        name of species
        
    output
    ------
    str
        shortend species name
    """
    return organism.replace(str(organism[1:organism.index("_")+1]), ".")

def location_finder(chromosome, place, error_range=20):
    """
    sam files have a chomosome number and starting place of the alingment, 
    use these to find the gene region that this alingment is in and output this region
    """
    location=[]
    for gene in genes_list:
        #gene[0] = chromsosme
        #geme[1] = begin
        #gen[2] = eind
        if chromosome == gene[0] and int(gene[1])-error_range <= int(place) and int(gene[2])+error_range >= int(place):
                location.append(gene)
                
    return location

def sort_dict(dictionary):
    """
    sorts a dictionary on the lenght of the list in the value from shortest to longest
    """
    return (sorted(dictionary, key = lambda key: len(dictionary[key])))

def order_seq(sequence):
    return "TTCTAATACGACTCACTATAG"+sequence+"GTTTTAGAGCTAGA"

#####
#MAIN
    

    ## ditionary off gene locations and there corosponding name of the gene and the species it came from
gene_names={}

##list of all the gene locations
genes_list=[]

#all the locations that a guide is located in
guide_locations={}

#the sequences of the guides
guide_seq={}

#list of all species in the set
species_list=[]

### read out the gff files
for file in glob.glob(f"{args.gff}/*"):

    species= short_species_name(file[file.index("/")+1:-4])
    with open(file, "r") as gff:

        for l in gff.readlines()[1:]:
            line=l.split("\t")
            chromosome= line[0]
            positionA=line[3]
            positionB=line[4]
            name=line[8][5:line[8].index(";")]

            genes_list.append([chromosome, positionA, positionB]) #add the found gene region in the list of all the genes
            gene_names[str([chromosome, positionA, positionB])]=f"{species} {name}" #add the name in the gene_names dictionary under the region as key


### read out the alingment files
#BTfile="bowtie_out/penicillium_rubens.candida_auris.sam"
for BTfile in glob.glob(f"{args.sam}/*"):


    with open(BTfile, "r") as handle:
        long_guide_species= BTfile[BTfile.index(".")+1:-4]
        guide_species= short_species_name(BTfile[BTfile.index(".")+1:-4])
        species_list.append(BTfile[BTfile.index(".")+1:-4])
        for l in [i for i in handle.readlines() if not i.startswith("@")]:
            line=l.split("\t")
            #print(line)

            #get all guides that have no mismatches and fined the gene they are located on
            if line[11] == "XA:i:0":
                guide_name=f"{guide_species} {line[0]}"
                
                chromosome_code= line[2]
                place=line[3]
                digest_site=location_finder(chromosome_code, place)
                
                for site in digest_site:
                    
                    #exception case when a list has not yet been initiated
                    if guide_name in guide_locations.keys():
                        guide_locations[guide_name].append(site)
                    else:
                        guide_locations[guide_name]=[site]
                        with open(f"guide_finder/{long_guide_species}.fasta", "r")as fastafile:
                            content=fastafile.readlines()
                            for rowN in range(0,len(content)):
                            
                                if content[rowN] == f">{line[0]}\n":
                                    guide_seq[guide_name]= content[rowN+1][:-1]
                                    break


#copy the information to be able to alter it

cut_sites= copy.deepcopy(guide_locations)
to_digest= copy.deepcopy(genes_list)
curated_guides= []

####main curating algorithm
### find the guide that has the most amount of cut sites

for g in range(0, len(cut_sites)):
    most=sort_dict(cut_sites)[-1]
    if len(cut_sites[most]) == 0:
        
        break
    curated_guides.append(most)

    #loop over all the cutsites of the guide to remove them from the to_digest list and the dictionary with all the cutsites of the guides
    for i in cut_sites[most]:
        if i in to_digest:
            to_digest.remove(i)

        for k in cut_sites: #loop over all guides in the dictionary and get the list of there cut sites
            if i in cut_sites[k]: #if the curent cut site is in this list remove from posible cut sites
                cut_sites[k].remove(i)

    del cut_sites[most]
    if len(to_digest) == 0: #stop finding new guides when al the genes in to digest are being targeted
        break

## build output tsv file
with open(f"{args.set}.tsv", "w") as tsvfile:
    tsvfile.writelines(["\t".join(["included species:"] + list(set(species_list))), "\n"])
    tsvfile.writelines(f"number of guides in final set: {len(curated_guides)}\n")
    for g in curated_guides:
        location_line=["; ".join([gene_names[str(i)], ": ".join(i)]) for i in guide_locations[g]]
        guideline=[g, guide_seq[g], str(len(guide_locations[g]))]+ location_line
        tsvfile.writelines(["\t".join(guideline), "\n"])
        
    tsvfile.writelines(f"could not digest: {len(to_digest)}/{len(genes_list)}\n")
    cndline=["; ".join([gene_names[str(i)], ": ".join(i)]) for i in to_digest]
    tsvfile.writelines("\n".join(cndline))
    
## build order.txt file
with open(f"{args.set}.order.txt", "w") as orderfile:
    for guide in curated_guides:
        orderfile.writelines([order_seq(guide_seq[guide][:-2]), "\n"])

