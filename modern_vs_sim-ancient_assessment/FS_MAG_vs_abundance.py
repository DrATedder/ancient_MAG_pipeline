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
            genome_length = int(record["Length"])
        except KeyError:
            print(f"TaxID not found for {accession_number}")
    handle.close()
    return taxid, genome_length


def collect_abundance_data(directory):
    abundance_dict = {}
    tmp_dict = {}
    for file in sorted(glob.glob(directory + "*_abundance.txt")):
        sample = os.path.basename(file.split("-")[0])
        replicate = int(file.split("-")[1].split("_")[0])
        with open(file, "r") as f_in:
            for line in f_in:
                OTU = line.split("\t")[0].strip()
                abundance = float(line.split("\t")[1].strip())
                if replicate == 1:
                        ### hopefully deal with duplicates
                    if OTU in tmp_dict:
                        tmp_dict[OTU].append(abundance)
                    else:
                        tmp_dict[OTU] = [abundance]
                else:
                    tmp_dict[OTU].append(abundance)
            if replicate == 3:
                with open(f"{file.split('-')[0]}_abundance_taxid.txt", "w") as f_out:
                    tmp_list = []
                    for k, v in tmp_dict.items():
                        info = get_taxid(k)
                        f_out.write(f"{k},{info[0]},{info[1]},{(sum(v))/3}\n")
                        tmp_list.append((sum(v))/3)
                tmp_dict = {}
    return "file(s) created."


def read_run_fastani(fastani):
    ID_list = []
    with open(fastani, "r") as f_in:
        for line in f_in:
            if "/" in line.split("\t")[1]:
                accession = line.split("\t")[1].split("/")[2].split(".")[0]
                taxid = get_taxid(accession)
                ID_list.append(int(taxid[0]))
            else:
                accession = line.split("\t")[1].split(".")[0]
                taxid = get_taxid(accession)
                ID_list.append(int(taxid[0]))
    return ID_list


abundance_directory = "~/simulated_abundance_files/"
fastani_directory = "~/fastANI_reports_allFS/deamSim/"


collect_abundance_data(abundance_directory)

with open("~/combined_MAG_data.csv", "w") as f_out:
    f_out.write("accession,taxID,genome_size,mean_abundance,presence\n")
    for file in sorted(glob.glob(fastani_directory + "*fastani.txt")):
        fastani = read_run_fastani(file)
        run = file.split("run")[1].split("_")[0]
        match_file = f'{os.path.basename(file).split("run")[0]}UoBsim{file.split("run")[1].split("_")[0]}_abundance_taxid.txt'
        print(match_file)
        with open(abundance_directory + match_file, "r") as f_in:
            for line in f_in:
                print(line)

                taxid = int(line.split(",")[1])
                if taxid in fastani:
                    f_out.write(f'{line.strip()},1\n')
                else:
                    f_out.write(f'{line.strip()},0\n')


