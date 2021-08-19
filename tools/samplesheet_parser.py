#!/usr/bin/env python3

# Parse covid samplesheet run by nextseq
# outputs a json file

from sample_sheet import SampleSheet
import json
import argparse
import os
import csv
import itertools

def arg():
    parser = argparse.ArgumentParser(prog="samplesheet_parser.py")
    parser.add_argument("-f", "--filepath", help="path to samplesheet file to parse")
    parser.add_argument("-r", "--run", help="runid")
    args = parser.parse_args()
    return args


# Parse input samplesheet
def sample_sheet(path,run):
    Sheet = SampleSheet(path)
    data = {}
    for sample in Sheet.samples:
        sample_name = sample['Sample_ID']
        description = sample['Description']

        # Skip controls
        if sample_name.startswith(('NegCtrl', 'PosCtrl', 'PosKon', 'NegKon')):
            continue
        else:
            data[sample['Sample_ID']] = []
            data[sample['Sample_ID']].append({
                'referensnummer': description.split("_")[0],
                'date': description.split("_")[1],
                'runtype': description.split("_")[2],
                'age': description.split("_")[3],
                'gender': description.split("_")[4],
                'lab_reference': description.split("_")[5],
                'postalcode': description.split("_")[6],
                'ct_value': description.split("_")[7]
            })
    #print metadata into json file
    with open(f"/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/{run}/metadata/{run}_metadata.json", 'w') as outfile:
        json.dump(data, outfile,indent=4)
    #print metadata into csv file
    csv_columns=['Sample_ID','referensnummer','date','runtype', 'age','gender','lab_reference','postalcode','ct_value']
    with open(f"/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/{run}/metadata/{run}_metadata.csv", "w") as csv_file:
        writer = csv.DictWriter(csv_file,csv_columns)
        writer.writeheader()
        for sample_id,metadata in data.items():
            row = {'Sample_ID': sample_id}
            row.update(metadata[0]) 
            writer.writerow(row)

def main():
    args = arg()
    path = args.filepath 
    run = args.run
    sample_sheet(path,run)


if __name__ == "__main__":
    main()
