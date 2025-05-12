import glob
import os
from Bio import Entrez

def get_taxid(accession_numbers):
    Entrez.email = "your@email.com"
    handle = Entrez.esummary(db="nucleotide", id=accession_numbers, retmode="xml")
    records = Entrez.read(handle)
    taxid_dict = {}
    for record in records:
        accession_number = record["Caption"]
        try:
            taxid = int(record["TaxId"])
            taxid_dict[accession_number] = taxid
        except KeyError:
            print(f"TaxID not found for {accession_number}")
    handle.close()
    return taxid_dict


def read_run_fastani(fastani):
    ID_list = []
    with open(fastani, "r") as f_in:
        for line in f_in:
            if "/" in line.split("\t")[1]:
                accession = line.split("\t")[1].split("/")[2].split(".")[0]
            else:
                accession = line.split("\t")[1].split(".")[0]
            ID_list.append(accession)
    return ID_list

def read_abundance(abundance):
    ID_list = []
    with open(abundance, "r") as f_in:
        for line in f_in:
            accession = line.split(",")[0].split(".")[0]
            ID_list.append(accession)
    return ID_list

def create_source_list(source_file):
    source_list = []
    with open(source_file, "r") as s_in:
        for line in s_in:
            taxID = line.split(",")[1].strip()
            source_list.append(taxID)
    return source_list

def check_abundance(source_list, abundance):
    overlap_list = []
    count = 0
    mismatch = 0
    with open(abundance, "r") as a_in:
        for line in a_in:
            OTU = line.split(",")[1].strip()
            if OTU in source_list:
                if not OTU in overlap_list:
                    number = source_list.count(OTU)
                    if number > 1:
                        for x in range(0, number, 1):
                            overlap_list.append(OTU)
                            count = count + 1
                    else:
                        overlap_list.append(OTU)
                        count = count +1
            else:
                mismatch = mismatch + 1

    return count, overlap_list, mismatch

    file_location = "/home/andrew/Downloads/fastANI_MAGs_examplesFS/fastANI_MAGs_examplesFS/"

for directory in glob.glob(file_location + "*/"):
    print(directory)
    for file in glob.glob(directory + "*fastani.txt"):
        with open(file, "r") as f_in:
            total = len(f_in.readlines())
        with open(f"{os.path.splitext(file)[0]}_taxID.txt", "w") as f_out:
            taxid_dict = get_taxid(read_run_fastani(file))
            for accession_number, taxid in taxid_dict.items():
                f_out.write(f"{accession_number},{taxid}\n")
    for file in glob.glob(directory + "*names.csv"):
        with open(f"{os.path.splitext(file)[0]}_taxID.txt", "w") as f_out:
            taxid_dict = get_taxid(read_abundance(file))
            for accession_number, taxid in taxid_dict.items():
                f_out.write(f"{accession_number},{taxid}\n")
        source = create_source_list(f"{os.path.splitext(file)[0]}_taxID.txt")
        #print(source)
    for ID_file in glob.glob(directory + "*fastani_taxID.txt"):
        if "deamSim" in ID_file:
            deam_check = check_abundance(source, ID_file)
            print(f"{os.path.basename(ID_file)}:{deam_check[0]}")
            print(f"False positives: {deam_check[2]}")
        else:
            sim_check = check_abundance(source, ID_file)
            print(f"{os.path.basename(ID_file)}:{sim_check[0]}")
            print(f"False positives: {sim_check[2]}")
