#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import datetime 
import logging
import fnmatch
import glob
from NGPinterface.hcp import HCPManager
from tools.samplesheet_parser import sample_sheet
from tools.check_files import check_files
from tools import log
from tools.microReport import nextseq as microreport
from tools.clc_sync import clc
from tools.emailer import email_micro

def arg():
    parser = argparse.ArgumentParser(prog="direkttest_cronscript.py")
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument("-ep", "--endpoint",
                            help="endpoint url")
    requiredNamed.add_argument("-aki", "--aws_access_key_id",
                            help="aws access key id")
    requiredNamed.add_argument("-sak", "--aws_secret_access_key",
                            help="aws secret access key")
    requiredNamed.add_argument("-b", "--bucket",
                            help="bucket name")
    requiredNamed.add_argument("-r", "--run",
                            help="runid for nextseq")
    parser.add_argument("-p", "--password", 
                            help="CLC password")
    parser.add_argument("--sshkey", 
                            help="GENSAM upload sshkey- password")
    args = parser.parse_args()

    return args


def setup_logger(name, log_path=None):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    stream_handle = logging.StreamHandler()
    stream_handle.setLevel(logging.DEBUG)
    stream_handle.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(stream_handle)

    if log_path:
        file_handle = logging.FileHandler(log_path, 'a')
        file_handle.setLevel(logging.DEBUG)
        file_handle.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handle)

    return logger


@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
def pangolin(path_list):
    # Specifics for nextseq data uploaded to GENSAM
    for f in path_list:
        if fnmatch.fnmatch(os.path.basename(f), "*_lineage_report.txt"):
            print("updating: " + f)
            df = pd.DataFrame(pd.read_csv(f, sep=",")) # csv file input
            for i,row in df.iterrows():
                df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4])) # change taxon names

            df.to_csv(os.path.dirname(os.path.abspath(f))+"/"+os.path.basename(f).replace(".txt","_gensam.txt"), index=None, header=True, sep="\t") # tab sep output without NULL

    # Specifics for nextseq data uploaded to HCP
    for f in path_list:
        if fnmatch.fnmatch(os.path.basename(f), "*_lineage_report.txt"):
            print("updating: " + f)
            df = pd.DataFrame(pd.read_csv(f, sep=",")).fillna(value = "NULL") # csv file input, fill empty with NULL
            for i,row in df.iterrows():
                df["taxon"] = df["taxon"].replace(row["taxon"], "_".join(row["taxon"].split("_")[1:4])) # change taxon names

            df.to_csv(os.path.dirname(os.path.abspath(f))+"/"+os.path.basename(f).replace(".txt","_fillempty.txt"), index=None, header=True, sep="\t") # tab sep output


@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
# Send artic csv file and pangolin result file to micro sftp
def micro_report(run):
    nextseqdir = "/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/" 
    articdir = "/medstore/results/clinical/SARS-CoV-2-typing/artic_results/"
    syncdir = "/seqstore/remote/outbox/sarscov2-micro/shared/nextseq"
    syncedfiles = "/medstore/results/clinical/SARS-CoV-2-typing/microbiologySync/syncedFiles_nextseq.txt"
    logfile = "/medstore/logs/pipeline_logfiles/sars-cov-2-typing/microReport_nextseq.log"
    run_name = run

    synclist = microreport(nextseqdir, articdir, syncdir, syncedfiles, logfile) 

    # Notify Microbiology about new data, only if actually synced
    if len(synclist) > 0:
        email_subject = 'Results from Artic pipeline now on sFTP and CLC'
        email_body = f'Artic/pangolin results and virus fasta from the run {run_name} is now available on the sFTP and CLC, respectively.'
        email_micro(email_subject, email_body)

@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
# Parse samplesheet and put in metadata json, for HCP upload
def samplesheet_parser(samplesheet_path,run):
    sample_sheet(samplesheet_path,run)


@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
# Import consensus fasta files to CLC
def clc_sync(password, run):
    # Variables for CLC upload
    server = "medair.sahlgrenska.gu.se"
    port = 7777
    user = "cmduser"
    clc(password,run,server,port,user) 


@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
# Upload files and json to selected bucket on HCP.
def upload_fastq(hcp_paths,hcpm,logger):
    for file_pg in hcp_paths:
        if not hcpm.search_objects(file_pg):
            try:
                hcpm.upload_file(file_pg, "covid-wgs/"+os.path.basename(file_pg))
                logger.info(f"uploading: {file_pg}")
            except Exception as e:
                logger.error(e)
                continue
        else:
            logger.error(f"Object already exists {file_pg}")
            continue


@log.log_error("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/nextseqwrapper_cronjob.log")
# Upload fasta, fastq and pangolin files to GENSAM
def gensam_upload(sshkey,run):
    #logsub = open("~/inbox/gensam.log", "w") 
    #cmd = ["gensamupload/gensamupload.py", "-r", run, "--sshkey-password", args.sshkey]
    cmd = ["python", "/apps/bio/repos/sars-cov-2-typing/gensamupload/gensamupload.py", "-r", run, "--sshkey-password", sshkey, "--no-mail", "--no-upload"]
    print(cmd)
    process1 = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(process1.stdout)
    print("ERR")
    #print(process1.stderr)

    # while process.wait() is None:
    #     pass
    # process.stdout.close()
    logsub.close()

def main():
    args = arg()
    sshkey=args.sshkey
    # Get runID
    if args.run:
        run = args.run
    else:
        for d in check_files("/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/2*"):
            run = os.path.split(d)[-1]
            
    print(run)
    # Set up the logfile
    now = datetime.datetime.now()
    logfile = os.path.join("/medstore/logs/pipeline_logfiles/sars-cov-2-typing/HCP_upload/", "HCP_upload_nextseq" + now.strftime("%y%m%d_%H%M%S") + ".log")
    logger = setup_logger('hcp_log', logfile)

    # Fix pangolin files for HCP and GENSAM
    pangolin_path = check_files("/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/2*/lineage/*")
    if len(pangolin_path) < 1:
        sys.exit()
    else:
        pangolin(pangolin_path)

    # Sync pangolin and artic files to micro sftp
    micro_report(run)

    # Parse nextseq samplesheet for metadata
    samplesheet_path = f'/seqstore/instruments/nextseq_500175_gc/Demultiplexdir/{run}/SampleSheet.csv'
    if not os.path.exists(f"/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/{run}/metadata"):
        os.makedirs(f"/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/{run}/metadata")

    samplesheet_parser(samplesheet_path,run)

    # Import consensus fasta files to CLC
    clc_sync(args.password,run)

    # Upload files to HCP
    hcp_paths = check_files("/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/2*/*/*")
    upload_fastq(hcp_paths, hcpm, logger)
    
    # Upload files to GENSAM
    #gensam_upload(sshkey,run)



if __name__ == "__main__":
    main()
