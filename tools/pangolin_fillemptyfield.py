#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import datetime as dt
import fnmatch
import glob

def arg():
    parser = argparse.ArgumentParser(prog="pangolin_fillemptyfield.py")
    parser.add_argument("-f", "--filepath", help="path to pangolin file to parse")
    parser.add_argument("-e", "--eurofins", action="store_true", help="check for eurofins files automatically")
    parser.add_argument("-n", "--nextseq", action="store_true", help="check for nextseq files automatically")
    parser.add_argument("-g", "--gensam", action="store_true", help="check for nextseq files automatically that will be uploaded to GENSAM")
    args = parser.parse_args()
    return args


# Automatic check for files proveded by eurofins (every 24 hrs)
def check_files_eurofins():
    now = dt.datetime.now()
    ago = now-dt.timedelta(minutes=1440)

    path_list = []
    for path in glob.glob('/medstore/results/clinical/SARS-CoV-2-typing/eurofins_data/goteborg/2021*/*', recursive=True):
        st = os.stat(path)
        mtime = dt.datetime.fromtimestamp(st.st_ctime) # ctime for time of change of file
        if mtime > ago:
            #print('%s modified %s'%(path, mtime))
            path_list.append(path)
    return path_list


# Automatic check for files proveded by nextseq (every 24 hrs)
def check_files_nextseq():
    now = dt.datetime.now()
    ago = now-dt.timedelta(minutes=1440)

    path_list = []
    for path in glob.glob('/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/21*/lineage/*', recursive=True):
        st = os.stat(path)
        mtime = dt.datetime.fromtimestamp(st.st_ctime) # ctime for time of change of file
        if mtime > ago:
            #print('%s modified %s'%(path, mtime))
            path_list.append(path)
    return path_list


def automatic(path, args):
    # Specifics for nextseq data uploaded to GENSAM
    if args.gensam:
        for f in path:
            if fnmatch.fnmatch(os.path.basename(f), "*_lineage_report.txt"):
                print("updating: " + f)
                df = pd.DataFrame(pd.read_csv(f, sep=",")) # csv file input
                for i,row in df.iterrows():
                    df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4])) # change taxon names

                df.to_csv(os.path.dirname(os.path.abspath(f))+"/"+os.path.basename(os.path.dirname(f))+"_"+os.path.basename(f).replace(".txt","_gensam.txt"), index=None, header=True, sep="\t") # tab sep output

    # Specifics for nextseq data uploaded to HCP
    if args.nextseq:
        for f in path:
            if fnmatch.fnmatch(os.path.basename(f), "*_lineage_report.txt"):
                print("updating: " + f)
                df = pd.DataFrame(pd.read_csv(f, sep=",")).fillna(value = "NULL") # csv file input, fill empty with NULL
                for i,row in df.iterrows():
                    df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4])) # change taxon names

                df.to_csv(os.path.dirname(os.path.abspath(f))+"/"+os.path.basename(os.path.dirname(f))+"_"+os.path.basename(f).replace(".txt","_fillempty.txt"), index=None, header=True, sep="\t") # tab sep output
    else:
        # Specifics for other files
        for f in path:
            if fnmatch.fnmatch(os.path.basename(f), "*_pangolin_lineage_classification.txt"):
                print("updating: " + f)
                df = pd.DataFrame(pd.read_csv(f, sep="\t")).fillna(value = "NULL") # tab sep input
                df.to_csv(os.path.dirname(os.path.abspath(f))+"/"+os.path.basename(os.path.dirname(f))+"_"+os.path.basename(f).replace(".txt","_fillempty.txt"), index=None, header=True, sep="\t") # tab sep output


# If only one file is selected, not automatic.
def fill_empty_cells(args):
    if args.nextseq:
        df = pd.DataFrame(pd.read_csv(args.filepath, sep=",")).fillna(value = "NULL")
        for i,row in df.iterrows():
            df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4]))

        df.to_csv(os.path.dirname(os.path.abspath(args.filepath))+"/"+os.path.basename(args.filepath).replace(".txt","_fillempty.txt"), index=None, header=True, sep="\t")

    if args.gensam:
        df = pd.DataFrame(pd.read_csv(args.filepath, sep=","))
        for i,row in df.iterrows():
            df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4]))

        df.to_csv(os.path.dirname(os.path.abspath(args.filepath))+"/"+os.path.basename(args.filepath).replace(".txt","_gensam.txt"), index=None, header=True, sep="\t")

    else:
        df = pd.DataFrame(pd.read_csv(args.filepath, sep=",")).fillna(value = "NULL")
        for i,row in df.iterrows():
            df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split\
("_")[1:4])) # change taxon names

        df.to_csv(os.path.dirname(os.path.abspath(args.filepath))+"/"+os.path.basename(args.filepath).replace("\
.txt","_fillempty.txt"), index=None, header=True, sep="\t") # tab sep output with NULL

def main():
    args = arg()
    if args.eurofins:
        path = check_files_eurofins()
        automatic(path,args)

    if args.nextseq or args.gensam and args.filepath:
        fill_empty_cells(args)
    elif args.nextseq or args.gensam:
        path = check_files_nextseq()
        automatic(path,args)
    
    else:
        fill_empty_cells(args)


if __name__ == "__main__":
    main()
