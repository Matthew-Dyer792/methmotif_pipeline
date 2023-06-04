# file stream script to integrate read count ID with full lineage
import csv
from pathlib import Path
import argparse
from re import L
import json


def read_metascope_out(args):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    output_name = args.outDir / args.file.name.replace('.csv', '.tsv')

    # open csv files
    with open(args.file, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')

        # get all taxonomic IDs to make one Entrez call
        id_list = []
        for line in reader:
            id_list.append(line.get('TaxonomyID'))

    i = 0
    exit = -(len(id_list) // -10000)
    records = []

    while i < exit:
        start = 0 + (i * 10000)
        if 10000 + (i * 10000) < len(id_list):
            end = 10000 + (i * 10000)
        else: 
            end = len(id_list)

        # print(start, end)

        search_list = id_list[start:end]

        # use Entrez through biopython to get lineage information from NCBI's taxonomy database
        Entrez.email = "mjd443@mun.ca"
        # handle = Entrez.efetch(db="Taxonomy", id=id_list, retmode="xml")
        handle = Entrez.efetch(db="Taxonomy", id=search_list, retmode="xml")
        temp_records = Entrez.read(handle)

        records = records + temp_records
        i += 1
        time.sleep(10)

    # print(len(records))

    # open csv files
    with open(args.file, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')

        i = 0
        for line in reader:
            # get the current lines record
            record = records[i]
            # call the 'write_file' function
            write_file(output_name, line, record)

            i += 1
            # print(i, '/', len(records))


def write_file(output_name, line, record):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the lowest taxonomic level (species, strain, or serotype)
        if record['Rank'] not in ['no rank', 'subspecies', 'strain']:
            line[record['Rank']] = record['ScientificName']

        # access the expanded lineage information
        lineage = record["LineageEx"]

        # set the names of the fields to be written to the new file
        for rank in lineage:
            if rank.get('Rank') == 'superkingdom' and rank.get('ScientificName') == 'Bacteria':
                field_names = ['read_count', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            elif rank.get('Rank') == 'superkingdom' and rank.get('ScientificName') == 'Viruses':
                # field_names = ['read_count', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'serotype']
                field_names = ['read_count', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'serotype']

        # search the Entrez output and set the full lineage name
        for name in lineage:
            rank = name.get('Rank')
            if rank in field_names:
                line[rank] = name.get('ScientificName')

        # remove unused columns from the line dictionary
        entries_to_remove = ['Proportion', 'TaxonomyID', 'EMreads', 'EMProportion', 'Genome']

        for key in line.keys():
            if key not in field_names and key not in entries_to_remove:
                # print("THIS IS A BAD", key)
                entries_to_remove.append(key)

        for entry in entries_to_remove:
            line.pop(entry)

        # entries_to_remove = line.keys()
        # if entry not in field_names:

        # print(line.keys())

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')
        writer.writerow(line)


def main():
    # setup command line arguments to be passed to the script
    parser = argparse.ArgumentParser(
        prog='Template to input json file',
        description='A program to convert the encode template.json into a valid input.json for the chip-seq pipeline'
    )
    parser.add_argument('-f', '--file', nargs='?', required=True, help='Template.json file')
    parser.add_argument('-o', '--outDir', nargs='?', required=True, help='directory for the output file(s)')
    args = parser.parse_args()

    # directory containing the template.json file
    args.file = Path(args.file)

    read_metascope_out(args)


if __name__ == "__main__":
    main()