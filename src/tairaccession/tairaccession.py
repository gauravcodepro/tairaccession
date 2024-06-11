#!/usr/bin/python3

def uniprotAgi(agi_file, ids_file):
    """
    _summary_
    this functions takes as input a uniprot_agi
    file and a file of ids and returns the converted ids
    or the links between the agi and the uniprot
    Arguments:
        agi_file -- _description_
        uniprot_agi file you can download from: 
        https://www.arabidopsis.org/download_files/Proteins/Id_conversions/Uniprot2AGI-JAN2023.txt
        ids_file -- _description_
        a list of the ids which you want to search for the association.
    Returns:
        _description_
        returns a list with a nested tuple which contains
        the uniprot id at the index 0 and the AGI as a list
    """
    uniprot_ag = {}
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
      for line in ids.readlines():
        agi_ids.append(line.strip())
        final_ids = [i for i in agi_ids if i!= ""]
    with open(os.path.abspath(os.path.join(os.getcwd(), agi_file)), "r") as uni:
      for line in uni.readlines():
        uniprot_ag[line.strip().split("\t")[0]] = [j for i in (list(map(lambda n: n.split(";"), \
                                                                line.strip().split("\t")[1:]))) for j in i]
      return [(k,v) for k,v in uniprot_ag.items() for i in final_ids if i in v]


def agiUniprot(agi_uniprot_file, ids_file):
    """
    _summary_
    this function takes a AGI2Uniprot and provides the conversions from
    the AGI to Uniprot. It returns a tuple with all the information on 
    the ids. You can provide a id file with the ids to be searched on a single line. 
    Arguments:
        agi_uniprot_file -- _description_
        Arabidopsis2Uniprot conversion file you can download from: 
        https://www.arabidopsis.org/download_files/Proteins/Id_conversions/AGI2uniprot-Jul2023.txt
        ids_file -- _description_
        a list of the ids which you want to search for the association.
    Returns:
        _description_
        returns a list with a nested tuple which contains
        the agi id at the index 0 and the uniprot as a list
    """
    ag_uniprot = {}
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
      for line in ids.readlines():
        agi_ids.append(line.strip())
        final_ids = [i for i in agi_ids if i!= ""]
    with open(os.path.abspath(os.path.join(os.getcwd(), agi_uniprot_file)), "r") as uni:
      for line in uni.readlines():
        ag_uniprot[line.strip().split("\t")[0]] = ''.join(line.strip().split("\t")[1:])
      return [(k,v) for k,v in ag_uniprot.items() for i in final_ids if i == k]


def uniprotTair(tair_uniprot_file, ids_file):
    """
    _summary_
    this function takes a TAIR2Uniprot mapping and returns a 
    nested tuple with the IDs and the Uniprot and thier locus
    ids. 
    Arguments:
        tair_uniprot_file -- _description_
        You can download this file from the TAIR from: 
        "https://www.arabidopsis.org/download_files/Proteins/Id_conversions/TAIR2UniprotMapping-Jul2023.txt"
        ids_file -- _description_
         a list of the ids which you want to search for the association.
    Returns:
        _description_
        returns a list with a nested tuple which contains
        the uniprot id at the index 0 and the locus and the AGI as a list
    """
    tair_uniprot = {}
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
      for line in ids.readlines():
        agi_ids.append(line.strip())
        final_ids = [i.split(".")[0] for i in [i.upper() for i in agi_ids if i!= ""]]
    with open(os.path.abspath(os.path.join(os.getcwd(), tair_uniprot_file)), "r") as uni:
      for line in uni.readlines():
        tair_uniprot[line.strip().split("\t")[0]] = line.strip().split("\t")[1:]
      return [(k,v) for k,v in tair_uniprot.items() for i in final_ids if i == v[1]]


def agiSpliceCoordinates(ids_file, gff_file, gene_type = None):
  """
  _summary_
  this functions takes a gff file and a specific type and returns the 
  coordinates of those ids including the splice variants. This function
  Arguments:
      ids_file -- _this takes the AGI ids from a text file_
      gff_file -- _this takes a gff3 file TAIR10_GFF3_genes.gff_
  Keyword Arguments:
      gene_type -- _description_ (default: {None})
      this is a keyworded argument, if gene is given then
      it will return the gene based coordinates, if the exon
      is given then it will give the exon based coordinates,
      if the three_prime_UTR is given then it will give the 
      three_prime_UTR based coordinates and if the five_prime_UTR
      is given then it will give the five_prime_UTR based arguments.
      A separate function getSpliceagiCoordinates is there if you
      want to have the information on the parent genes as well. 
      if you are searching for the gene then you should not list the
      splice variants and if you are searching for the splice variants
      then you should not list the genes. Splice variants starts with the 
      [geneid].1 or [gene].{0-9} on a regular expression. 
  Returns:
      _returns a nested list with the AGI and the coordinates_
  """
  if ids_file and gff_file and gene_type == "gene":
        tair = pd.read_csv(gff_file, sep = "\t")
        renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
        gene_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "gene").dropna()
        gene_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "gene").dropna()["Start"].to_list()))
        gene_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "gene").dropna()["End"].to_list()))
        gene_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "gene").dropna()["Strand"].to_list()))
        gene_type_gene_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]].where(renaming_tair["gene_type"] == "gene").dropna() 
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))
        gene_type_gene_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "gene").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: \
                                                               n.replace("ID=", "")))["Gene_ID"].to_list()
        arabidopsis_gene = []
        for i in range(len(gene_type_start)):
            arabidopsis_gene.append([gene_type_gene_ID_AGI[i],gene_type_start[i], gene_type_end[i], gene_type_strand[i]])
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        selected_gene = [i for i in arabidopsis_gene for j in final_ids if j in i[0]]
        return selected_gene
  if ids_file and gff_file and gene_type == "exon":
        tair = pd.read_csv(gff_file, sep = "\t")
        renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
        exon_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "exon").dropna()
        exon_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "exon").dropna()["Start"].to_list()))
        exon_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "exon").dropna()["End"].to_list()))
        exon_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "exon").dropna()["Strand"].to_list()))
        exon_type_exon_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))
        exon_type_exon_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()
        arabidopsis_exon = []
        for i in range(len(exon_type_start)):
            arabidopsis_exon.append([exon_type_exon_ID_AGI[i], exon_type_start[i], exon_type_end[i], exon_type_strand[i]])
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        selected_exon = [i for i in arabidopsis_exon for j in final_ids if j in i[0]]
        return selected_exon
  if ids_file and gff_file and gene_type == "three_prime_UTR":
        tair = pd.read_csv(gff_file, sep = "\t")
        renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
        three_prime_UTR_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "three_prime_UTR").dropna()
        three_prime_UTR_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "three_prime_UTR").dropna()["Start"].to_list()))
        three_prime_UTR_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "three_prime_UTR").dropna()["End"].to_list()))
        three_prime_UTR_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "three_prime_UTR").dropna()["Strand"].to_list()))
        three_prime_UTR_type_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "three_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))
        three_prime_UTR_type_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "three_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()
        arabidopsis_three_prime_UTR = []
        for i in range(len(three_prime_UTR_type_start)):
            arabidopsis_three_prime_UTR.append([three_prime_UTR_type_ID_AGI[i],\
                                                   three_prime_UTR_type_start[i],\
                                                   three_prime_UTR_type_end[i],\
                                                   three_prime_UTR_type_strand[i]])
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
        final_ids =[i.upper() for i in agi_ids if i!= ""]
        selected_three_prime_UTR = [i for i in arabidopsis_three_prime_UTR for j in final_ids if j in i[0]]
        return selected_three_prime_UTR
  if ids_file and gff_file and gene_type == "five_prime_UTR":
        tair = pd.read_csv(gff_file, sep = "\t")
        renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
        five_prime_UTR_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "five_prime_UTR").dropna()
        five_prime_UTR_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "five_prime_UTR").dropna()["Start"].to_list()))
        five_prime_UTR_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "five_prime_UTR").dropna()["End"].to_list()))
        five_prime_UTR_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "five_prime_UTR").dropna()["Strand"].to_list()))
        five_prime_UTR_type_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "five_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))
        five_prime_UTR_type_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "five_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()
        arabidopsis_five_prime_UTR = []
        for i in range(len(five_prime_UTR_type_start)):
            arabidopsis_five_prime_UTR.append([five_prime_UTR_type_ID_AGI[i],\
                                                   five_prime_UTR_type_start[i],\
                                                   five_prime_UTR_type_end[i],\
                                                   five_prime_UTR_type_strand[i]])
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
        final_ids =[i.upper() for i in agi_ids if i!= ""]
        selected_five_prime_UTR = [i for i in arabidopsis_five_prime_UTR for j in final_ids if j in i[0]]
        return selected_five_prime_UTR
  if ids_file and gff_file and gene_type == "cds":
        tair = pd.read_csv(gff_file, sep = "\t")
        renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
        cds_type = renaming_tair[["gene_type", "Start", "End", "Strand"]].where(renaming_tair["gene_type"] == "CDS").dropna()
        cds_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "CDS").dropna()["Start"].to_list()))
        cds_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "CDS").dropna()["End"].to_list()))
        cds_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "CDS").dropna()["Strand"].to_list()))
        cds_type_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "CDS").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.split(",")[0]).apply(lambda n: n.replace("Parent=", "")))
        cds_type_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "CDS").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.split(",")[0]).apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()  
        arabidopsis_cds = []
        for i in range(len(cds_type_start)):
            arabidopsis_cds.append([cds_type_ID_AGI[i], cds_type_start[i],cds_type_end[i], cds_type_strand[i]])
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        selected_cds = [i for i in arabidopsis_cds for j in final_ids if j in i[0]]
        return selected_cds

def agiGO(ids_file, association_file):
    """
  _summary_
  this provides the Go information for the AGI ids presented in 
  a file and presents them as tuples with the AGI at index 0 and 
  the GO at index 1. This takes the gene_association.tair as the 
  association file.
  Arguments:
      ids_file -- _a AGI ids file for the search_
      association_file -- _a association file for the search_
  Returns:
      _description_
      this provides a nested tuples of the association fetched
      from the association file for those agis. 
    """
    with open(os.path.join(os.getcwd(),association_file), "r") as goslim:
        with open(os.path.join(os.getcwd(),association_file + "name"), "w") as gofinal:
            for line in goslim.readlines():
                if line.startswith("!"): 
                    continue     
                gofinal.write(line)
    agi_ids = []          
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
        final_ids = [i.split(".")[0] for i in [i.upper() for i in agi_ids if i!= ""]]            
        go_data = pd.read_csv(os.path.join(os.getcwd(),association_file + "name"), sep = "\t")
        AGI = go_data.iloc[::,1].to_list()
        GO = go_data.iloc[::,4].to_list()
        return [j for i in final_ids for j in ([(i,j) for i,j in zip(AGI,GO)]) if j[0] == i]

def agiTAIR(ids_file, association_file):
  """
  _summary_
  this provides the Go information for the AGI ids presented in 
  a file and presents them as tuples with the AGI at index 0 and 
  the TAIR communication at index 1. This takes the gene_association.tair 
  as the association file.
  Arguments:
      ids_file -- _a AGI ids file for the search_
      association_file -- _a association file for the search_§  
  Returns:
      _description_
      a nested list of the AGI and the TAIR.
  """
  with open(os.path.join(os.getcwd(),association_file), "r") as goslim:
      with open(os.path.join(os.getcwd(),association_file + "name"), "w") as gofinal:
          for line in goslim.readlines():
              if line.startswith("!"): 
                  continue     
              gofinal.write(line)
  agi_ids = []          
  with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
    for line in ids.readlines():
        agi_ids.append(line.strip())
    final_ids = [i.split(".")[0] for i in [i.upper() for i in agi_ids if i!= ""]]            
    go_data = pd.read_csv(os.path.join(os.getcwd(),association_file + "name"), sep = "\t")
    AGI = go_data.iloc[::,1].to_list()
    tair_communication = go_data.iloc[::,5].apply(lambda n: n.split(":")[2]).apply(lambda n: n.replace("|PMID", "")).to_list()
    return set([j for i in final_ids for j in ([(i,j) for i,j in zip(AGI,tair_communication)]) if j[0] == i])

def agiDescription(ids_file, association_file):
  """
  _summary_
  this provides the Go information for the AGI ids presented in 
  a file and presents them as tuples with the AGI at index 0 and 
  the Description associated with the gene at index 1. This takes 
  the gene_association.tair as the association file.
  Arguments:
      ids_file -- _a AGI ids file for the search_
      association_file -- _a association file for the search_
  Returns:
      _description_
      a nested list of the AGI and the Description.
  """
  with open(os.path.abspath(os.path.join(os.getcwd(),association_file), "r")) as goslim:
    with open(os.path.abspath(os.path.join(os.getcwd(),association_file + "name"), "w")) as gofinal:
        for line in goslim.readlines():
          if line.startswith("!"): 
            continue     
          gofinal.write(line)
  agi_ids = []        
  with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
    for line in ids.readlines():
        agi_ids.append(line.strip())
    final_ids = [i.split(".")[0] for i in [i.upper() for i in agi_ids if i!= ""]]            
    go_data = pd.read_csv(os.path.join(os.getcwd(),association_file + "name"), sep = "\t")
    AGI = go_data.iloc[::,1].to_list()
    description = go_data.iloc[::,10].to_list()
    return set([j for i in final_ids for j in ([(i,j) for i,j in zip(AGI,description)]) if j[0] == i])


def phytozomePacID(gff_file, ids_file):
    """
    _summary_
    this function analyzes the phytozome file for the 
    search of the corresponding pacid and provided a 
    agi_id file, phytozome file, it will search for the 
    corresponding pacid. 
    Arguments:
        gff_file -- _description_
        phytozome gff file and you can download from the phytozome 
        "Athaliana_447_Araport11.gene_exons.gff3"
        ids_file -- _description_
        takes the agiIDs for the search
    Returns:
        _description_
        returns a nested list with agi listed at 
        list[0] and the pacid listed at list[1].
    """
    with open(os.path.abspath(os.path.join(os.getcwd(),gff_file)), "r") as phytozome:
        with open(os.path.abspath(os.path.join(os.getcwd(),gff_file + "name")), "w") as phytozomer:
            for line in phytozome.readlines():
                if line.startswith("!"): 
                    continue     
                phytozomer.write(line)
        phytozomedataframe = pd.read_csv(os.path.abspath(os.path.join(os.getcwd(),gff_file + "name")), sep = "\t")
        mRNA = phytozomedataframe.iloc[::,[2,8]]. \
                  where(phytozomedataframe.iloc[::,[2,8]]["gene"] == "mRNA").dropna()
        name = [i.split("=")[1] for i in ([j for i in ([i.split(";") \
                        for i in (mRNA.iloc[::,1].to_list())]) \
                                for j in i if j.startswith("Name=")])]
        pacid = [j for i in ([i.split(";") \
                          for i in (mRNA.iloc[::,1].to_list())]) \
                                    for j in i if j.startswith("pacid=")]
        agiPacID = [(i,j) for i,j in zip(name,pacid)]
        agi_ids = []
        final_ids = []
        with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
            for line in ids.readlines():
                agi_ids.append(line.strip())
            final_ids = [i.upper() for i in agi_ids if i!= ""]
            return [i for i in agiPacID for j in final_ids if j==i[0]]
        
def prepareFunctionalNamePhytozome(phytozome_file, output_file):
    """
    _summary_
    this function prepares the files for the pacID analysis from the phytozome
    gene names files. Example file: Athaliana_447_Araport11.geneName.txt 
    Arguments:
        phytozome_file -- _description_
        phytozome file which contains the gene names and the description.
    Return: return_description    
        output_file -- _description_
        a nested list with the pacid and the gene names along with the splice variants. 
    """
    functionalName = {}
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file)), "r") as file:
        for line in file.readlines():
            functionalName[line.strip().split("\t")[0]] = ''.join(j for i in \
                                                        ([line.strip().split("\t")[2:]]) for j in i)
    with open(os.path.abspath(os.path.join(os.getcwd(),output_file)), "w") as processed:
        print(f"functionalName={functionalName}", file=processed)

def preparegeneNamePhytozome(phytozome_file, output_file):
    """
    _summary_
    this function prepares the files for the pacID analysis from the phytozome
    functional names files. Example file: Athaliana_447_Araport11.defline.txt 
    Arguments:
        phytozome_file -- _description_
        phytozome file which contains the gene names and the description.
    Return: return_description    
        output_file -- _description_
        a nested list with the pacid and the gene names along with the splice variants. 
    """   
    geneName = {}
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file)), "r") as file:
        for line in file.readlines():
            geneName[line.strip().split("\t")[0]] = ''.join([j for i in \
                                                            ([line.strip().split("\t")[1]]) for j in i])
    with open(os.path.abspath(os.path.join(os.getcwd(),output_file)), "w") as processed:
        print(f"geneName={geneName}", file=processed) 

def FunctionalNames(ids_file):
    import tairaccession
    from tairaccession.functionalNames import functionalNames
    """
    _summary_
    this function takes an id file and scan
    against the phytozome and presents the 
    gene name and the functional name. 
    Arguments:
        ids_files -- _description_
        this contains the ids you want to 
        search for the functional name.
    Returns:
        _description_
        provides a nested list where the 
        list[0] is the agi provided in the
        file and list[1] is the functional
        name.
    """
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        return [(k,v) for k,v in tairaccession.functionalNames.functionalNames.items() for i in final_ids if i==k]


def GeneNames(ids_file):
    import tairaccession
    from tairaccession.geneNames import geneNames
    """
    _summary_
    this function takes an id file and scan
    against the phytozome and presents the 
    gene name and the gene name. 
    Arguments:
        ids_files -- _description_
        this contains the ids you want to 
        search for the gene name
    Returns:
        _description_
        provides a nested list where the 
        list[0] is the agi provided in the
        file and list[1] is the gene
        name.
    """
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        return [(k,v) for k,v in tairaccession.geneNames.geneNames.items() for i in final_ids if i==k]


def visualizeAgiCDS(ids_file, gff_file, genome_name, size, save_path):
    tair = pd.read_csv(gff_file, sep = "\t")
    renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
    cds_type = renaming_tair[["gene_type", "Start", "End", "Strand"]].where(renaming_tair["gene_type"] == "CDS").dropna()
    cds_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "CDS").dropna()["Start"].to_list()))
    cds_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "CDS").dropna()["End"].to_list()))
    cds_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "CDS").dropna()["Strand"].to_list()))
    cds_type_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "CDS").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.split(",")[0]).apply(lambda n: n.replace("Parent=", "")))
    cds_type_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "CDS").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.split(",")[0]).apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()  
    arabidopsis_cds = []
    for i in range(len(cds_type_start)):
        arabidopsis_cds.append([cds_type_ID_AGI[i], cds_type_start[i],cds_type_end[i],cds_type_strand[i]])
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
    final_ids = [i.upper() for i in agi_ids if i!= ""]
    global selected_cds
    selected_cds = [tuple([i[1],i[2],i[3]]) for i in arabidopsis_cds for j in final_ids if j in i[0]]
    from pygenomeviz import GenomeViz
    name, genome_size = genome_name, int(size)
    gv = GenomeViz()
    selected_list = [tuple([i[1],i[2],i[3]]) for i in arabidopsis_cds for j in final_ids if j in i[0]]
    track = gv.add_feature_track(name, genome_size)
    for idx, cds in enumerate(selected_list, 1):
        start, end, strand = cds
        track.add_feature(start, end, strand, label=f"CDS{idx:02d}")
    gv.savefig(os.path.abspath(os.path.join(os.getcwd(),save_path,"plot_coding.png")))
    return selected_cds    


def visualizeExons(ids_file, gff_file, genome_name, size,save_path):
     tair = pd.read_csv(gff_file, sep = "\t")
     renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
     exon_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "exon").dropna()
     exon_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "exon").dropna()["Start"].to_list()))
     exon_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "exon").dropna()["End"].to_list()))
     exon_type_strand = list(map(lambda n: 1 if n == "+" else -1,renaming_tair[["gene_type", "Start", "End", "Strand"]].
                                           where(renaming_tair["gene_type"] == "exon").dropna()["Strand"].to_list()))
     exon_type_exon_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))
     exon_type_exon_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n: n.replace("Parent=", "")))["Gene_ID"].to_list()
     arabidopsis_exon = []
     for i in range(len(exon_type_start)):
            arabidopsis_exon.append([exon_type_exon_ID_AGI[i], exon_type_start[i], exon_type_end[i], exon_type_strand[i]])
     agi_ids = []
     final_ids = []
     with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
     final_ids = [i.upper() for i in agi_ids if i!= ""]
     global selected_exon
     selected_exon = [i for i in arabidopsis_exon for j in final_ids if j in i[0]]
     from pygenomeviz import GenomeViz
     ame, genome_size = genome_name, int(size)
     gv = GenomeViz()
     selected_list = [tuple([i[1],i[2],i[3]]) for i in arabidopsis_exon for j in final_ids if j in i[0]]
     track = gv.add_feature_track(name, genome_size)
     for idx, cds in enumerate(selected_list, 1):
        start, end, strand = cds
        track.add_feature(start, end, strand, label=f"exon{idx:02d}")
     gv.savefig(os.path.abspath(os.path.join(os.getcwd(),save_path,"plot_exon.png")))
     return selected_exon
 
def readTairNCBI(input, output, arg_id):
    """
    _summary_
    this function establishes the missing link between the 
    tair and the pubmed, given the file ATH_GO_GOSLIM.txt
    and an output file name and the specific locus or the
    gene id, it will first fetch the corresponding pubmed id
    and then will fetch the corresponding abstract for those 
    pubmed id. Establishes a connecting link between the 
    gene and locus to language model training. 

    Arguments:
        input -- _description_
        ATH_GO_GOSLIM.txt from the TAIR release
        output -- _description_
        the output file name
        arg_id -- _description_
        gene or locus id. 

    Returns:
        _description_
        A nested list with the pubmed id and the corresponding 
        publication details.
    """
    with open(input, "r") as file_read:
        with open(output, "w") as file_out:
            for line in file_read.readlines():
                if line.startswith("!"):
                    continue
                file_out.write(line)
    with open(output, "r") as file_re_read:
        read_file = pd.read_csv(output, sep = "\t")
        correspondence = read_file.iloc[::,[1,12]].iloc[::,1]. \
                           apply(lambda n: n.split("|")).to_list()
        gene_id = list(map(lambda n: n.replace("gene:", ""), \
                            map(lambda n: n.replace("locus:", ""), \
                                         read_file.iloc[::,1].to_list())))
    pubmed = []
    for i in range(len(correspondence)):
        pubmed.append([gene_id[i],correspondence[i]])
    format_id = set([j.replace("PMID:", "") for i in ([i.split() \
                    for i in ([j for i in [pubmed[i][1] for i \
                                in range(len(pubmed)) if pubmed[i][0] == arg_id] \
                                          for j in i]) if "PMID" in i]) for j in i])
    format_id_links = list(format_id)
    format_check = [] 
    for i in range(len(format_id_links)):
        format_check.append(f"https://pubmed.ncbi.nlm.nih.gov/{format_id_links[i]}/")
    ncbi_derive_information = {}
    for i in range(len(format_check)):
        ncbi_derive_information[format_check[i]] = ''.join([i.get_text().strip() \
            for i in BeautifulSoup(urlopen(format_check[i]), \
                "html.parser").find_all("div", class_ = "abstract-content selected")])
    return [(k,v) for k,v in ncbi_derive_information.items() if k or v != ""]
