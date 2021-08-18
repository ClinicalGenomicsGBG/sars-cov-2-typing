#!/usr/bin/env python

import click
import datetime
import os
import sys
import glob
from sample_sheet import SampleSheet
from collections import defaultdict
import pysftp
import smtplib
from email.message import EmailMessage
from shutil import copyfile
import csv
from config import *

@click.command()
@click.option('-r', '--runid', required=True,
              help='RUNID of the seqrun to upload')
@click.option('-d', '--demultiplexdir', required=True,
              default='/seqstore/instruments/nextseq_500175_gc/Demultiplexdir',
              help='Path to demultiplexdir')
@click.option('-i', '--inputdir', required=True,
              default='/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/files_to_gensam',
              help='Path to location of files to upload')
@click.option('-s', '--samplesheetname', required=True,
              default='SampleSheet.csv',
              help='Name of SampleSheet file')
@click.option('--regioncode', required=True,
              default='14',
              help='FOHM region code')
@click.option('--labcode', required=True,
              default='SE300',
              help='FOHM lab code')
@click.option('-l', '--logdir', required=True,
              default='/medstore/logs/pipeline_logfiles/sars-cov-2-typing/GENSAM-upload',
              help='Path to directory where logs should be created')
@click.option('--gensamhost', required=True,
              default='gensam-sftp.folkhalsomyndigheten.se',
              help='FOHM GENSAM hostname')
@click.option('--sftpusername', required=True,
              default='se300',
              help='Username to the GENSAM sFTP')
@click.option('--sshkey', required=True,
              default='~/.ssh/id_rsa',
              help='Path/to/private/sshkey')
@click.option('--sshkey-password', required=True,
              help='SSH key password')
@click.option('-g', '--gensamcsvdir', required=True,
              default='/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/gensam_upload',
              help='Path to dir where GENSAM upload csv file should be saved')
@click.option('--manualcsv',
              help='Manually specify a CSV file to upload to GENSAM with samples and info. Also used to specify which samples to upload')
@click.option('--uploadedsamples',
              default='/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data/gensam_upload/uploaded2gensam.txt',
              help='File containing already uploaded samples.')
@click.option('--no-mail', is_flag=True,
              help="Set if you do NOT want e-mails to be sent")
@click.option('--no-upload', is_flag=True,
              help="Set if you do NOT want to upload files to FOHM. Will still try to connect to the sFTP.")
def main(runid, demultiplexdir, logdir, inputdir, samplesheetname, regioncode, labcode, sshkey, 
         sshkey_password, gensamhost, sftpusername, gensamcsvdir, manualcsv, uploadedsamples, no_mail, no_upload):

    #Run checks on all given inputs
    checkinput(runid, demultiplexdir, inputdir, regioncode, labcode, logdir, 
               samplesheetname, gensamcsvdir, manualcsv, uploadedsamples)
    
    # Start the logging
    now = datetime.datetime.now()
    logfile = os.path.join(logdir, "GENSAM-upload_" + now.strftime("%y%m%d_%H%M%S") + ".log")
    logfile_sftp = os.path.join(logdir, "GENSAM-upload_" + now.strftime("%y%m%d_%H%M%S") + "_sFTP.log")
    #logfile = os.path.join(logdir, "GENSAM-upload_debug.log")
    #logfile_sftp = os.path.join(logdir, "GENSAM-upload_debug_sftp.log")

    log = open(logfile, "a")
    log.write("----\n")
    log.write(writelog("LOG", "Starting GENSAM upload workflow"))
 
    #Get the path to samplesheet
    sspath = os.path.join(demultiplexdir, runid, samplesheetname)
    #Read in all sampleIDs
    samples = sample_sheet(sspath)

    #Get a list of all fastq and fasta files to upload
    syncdict = defaultdict(lambda: defaultdict(dict))
    log.write(writelog("LOG", "Finding all files to upload."))
    for sample in samples:
        #Find all fastq files to upload.
        fastqpath = os.path.join(inputdir, runid, 'fastq')
        for fastqfile in glob.glob(fastqpath + "/" + sample + "*fastq.gz"):
            targetlink = os.readlink(fastqfile) #readlink
            if targetlink.endswith("R1_001.fastq.gz"): #Fastq files need to have this extension right now. It's ugly
                syncdict[sample]['fastq']['R1'] = targetlink
            elif targetlink.endswith("R2_001.fastq.gz"):
                syncdict[sample]['fastq']['R2'] = targetlink
            else:
                log.write(writelog("ERROR", "Found fastq file with ending other than R1(R2)_001.fastq.gz"))
                if not no_mail:
                    email_error(logfile, "FASTQ UPLOAD")
                sys.exit("ERROR: Found fastq file with ending other than R1(R2)_001.fastq.gz")
        #Find all fasta files to upload based on existing links.
        fastapath = os.path.join(inputdir, runid, 'fasta')
        for fastafile in glob.glob(fastapath + "/" + sample + "*consensus.fa"):
            targetlink = os.readlink(fastafile) #readlink
            #Store info in dict
            syncdict[sample]['fasta'] = targetlink

    #Check that all fastq files are paired
    #This code feels a bit clunky
    
    for sample in syncdict:
        if syncdict[sample]['fastq']['R1']:
            if not syncdict[sample]['fastq']['R2']:
                log.write(writelog("ERROR", "No R2 file found for " + syncdict[sample]['fastq']['R1'] + "."))
                if not no_mail:
                    email_error(logfile, "FASTQ PAIRING")
                sys.exit("ERROR: No R2 file found for " + syncdict[sample]['fastq']['R1'] + ".")
        if syncdict[sample]['fastq']['R2']:
            if not syncdict[sample]['fastq']['R1']:
                log.write(writelog("ERROR", "No R1 file found for " + syncdict[sample]['fastq']['R2'] + "."))
                if not no_mail:
                    email_error(logfile, "FASTQ PAIRING")
                sys.exit("ERROR: No R1 file found for " + syncdict[sample]['fastq']['R2'] + ".")

    #Check how manny files there is to upload
    for filetype in ['fastq', 'fasta']:
        numfiles = countkeys(syncdict, filetype)
        if filetype == 'fastq':
            log.write(writelog("LOG", "Found " + str(numfiles) + " " + filetype  + " pairs to upload."))
        else:
            log.write(writelog("LOG", "Found " + str(numfiles) + " " + filetype  + " files to upload."))

    #Make an csv file with FOHM info
    #Also store all samples from this file in a list
    gensamcsv_samples = []
    date_simple = now.strftime("%Y-%m-%d")
    if manualcsv:
        #check if the specified file is in correct place
        gensam_csv = os.path.join(gensamcsvdir, "_".join((regioncode,labcode,date_simple, "komplettering.csv")))
        if os.path.abspath(manualcsv) == gensam_csv:
            log.write(writelog("LOG", "The provided GENSAM csv file seem to be in correct place. Will use it as is."))
        else:
            log.write(writelog("LOG", "Copying the provided GENSAM csv file to " + gensamcsvdir))
            try:
                copyfile(manualcsv, gensamdir)
            except:
                log.write(writelog("ERROR", "Could not copy the provided GENSAM csv file to " + gensamcsvdir))
        
        gensam_csv = manualcsv

        #Load in the samples from the file
        with open(gensam_csv) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                #skip header
                if row[0].startswith('provnummer'):
                    continue

                #add sample names to list
                gensamcsv_samples.append(row[0])
                

    else:
        gensam_csv = os.path.join(gensamcsvdir, "_".join((regioncode,labcode,date_simple, "komplettering.csv")))
        csvout = open(gensam_csv, "w")
        csvout.write("provnummer,urvalskriterium,GISAID_accession\n")
        for sample, criteria in samples.items():
            if syncdict[sample]['fasta']:
                if criteria not in selection_criteria:
                    sys.exit("ERROR: selected criteria for sample" + sample + " is not valid.")
                else:
                    csvout.write(','.join((sample.split('_')[0], selection_criteria[criteria] , " \n")))
                    gensamcsv_samples.append(sample)
                    log.write(writelog("LOG", "Wrote GENSAM .csv file to : " + gensam_csv))
    
    #Open the connection to the sFTP
    if no_upload:
        log.write(writelog("LOG", "No-upload flag set. Will just try the sFTP connection."))
    else:
        log.write(writelog("LOG", "Starting sFTP upload."))

    try:
        sftp = pysftp.Connection(gensamhost, username=sftpusername, private_key=sshkey, private_key_pass=sshkey_password, log=logfile_sftp)
        sftp.chdir("till-fohm")
    except:
        log.write(writelog("ERROR", "Establishing sFTP connection failed. Check the sFTP log @ " + logfile_sftp))
        if not no_mail:
            email_error(logfile, "sFTP CONNECTION")
        sys.exit("ERROR: Establishing sFTP connection failed. Check the sFTP log @ " + logfile_sftp)

    #Upload all files to the FOHM FTP
    #Read in all files previously uploaded. Perhaps there is a faster, indexwed way?
    uploads_dict = {}
    with open(uploadedsamples) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
            #add sample names to dict
                uploads_dict[row[0]] = ''
                

    #Open file with previously uploaded files for writing
    uploads_samples = open(uploadedsamples, "a")
    #fastq and fasta file
    for sample in syncdict:
        #Skip sample if not in the gensamcsv file. Only really matters if there is a manually supplied csv file
        if not sample in gensamcsv_samples:
            log.write(writelog("WARNING", f'Sample {sample} not found in {os.path.basename(gensam_csv)}, skipping.'))
            continue

        #Is sample already uploaded (and no-upload flag not set)? If so, skip
        if sample in uploads_dict:
            log.write(writelog("WARNING", f'Sample {sample} found in file of previously uploaded samples @ {uploadedsamples}. Skipping it.'))
            continue

        #Get all fastq files and construct correct names
        if syncdict[sample]['fastq']['R1']: #Check if sample has fastq files to upload
            fastqR1_src = syncdict[sample]['fastq']['R1']
            fastqR2_src = syncdict[sample]['fastq']['R2']

            samplename_R1 = sample.replace("_", "-") + '_1.fastq.gz'
            fastqR1_trgt = '_'.join((regioncode, labcode, samplename_R1))
            samplename_R2 = sample.replace("_", "-") + '_2.fastq.gz'
            fastqR2_trgt = '_'.join((regioncode, labcode, samplename_R2))
            print(sample, fastqR1_trgt, sample.split('_')[0])
    
            #Upload to sFTP
            if not no_upload:
                sftp.put(fastqR1_src, fastqR1_trgt)
                sftp.put(fastqR2_src, fastqR2_trgt)
                                
        # Get all fasta files and construct correct names
        if syncdict[sample]['fasta']: #Check if sample has fastq files to upload
            fasta_src = syncdict[sample]['fasta']
            samplename_fasta = sample.replace("_", "-") + '.consensus.fasta'
            fasta_trgt = '_'.join((regioncode, labcode, samplename_fasta))
            #Upload to sFTP
            if not no_upload:
                sftp.put(fasta_src, fasta_trgt)

        #Save sample name in uploaded samples file unless
        if not no_upload:
            uploads_samples.write(f'{sample}\n')


    #Pangolin classification lineage file    
    lineagepath = os.path.join(inputdir, runid, 'lineage', runid + "_lineage_report_gensam.txt")
    lineage_trgt = '_'.join((regioncode, labcode, date_simple, "pangolin_classification.txt"))
    if not no_upload:
        sftp.put(lineagepath, lineage_trgt)

    #Close the uploaded samples file
    uploads_samples.close()

    #Close the sFTP connection
    sftp.close()

    if no_upload:
        log.write(writelog("LOG", "Completed test of sFTP connection."))
    else:
        log.write(writelog("LOG", "Finished the sFTP upload."))

    #Send an e-mail to FOHM (and clinicalgenomics) that upload has happened
    #csv file should be attached
    if no_upload:
        #email_fohm(gensam_csv)
        log.write(writelog("LOG", "No-upload flag set. Skipping mail to FOHM"))
    elif no_mail:
        log.write(writelog("LOG", "No-mail flag set. Skipping mail to FOHM"))
    else:
        email_fohm(gensam_csv)
        
    #Finished the workflow
    log.write(writelog("LOG", "Finished GENSAM upload workflow."))
    log.close()

def checkinput(runid, demultiplexdir, inputdir, regioncode, labcode, logdir, 
               samplesheetname, gensamcsvdir, manualcsv, uploadedsamples):
    #Make sure the samplesheet exists
    sspath = os.path.join(demultiplexdir, runid, samplesheetname)
    if not os.path.isfile(sspath):
        sys.exit("ERROR: Could not find SamleSheet @ " +  sspath)

    #If there is a manual GENSAM csv file. Check that it eists
    if manualcsv and not os.path.isfile(manualcsv):
        sys.exit("ERROR: The specified manual GENSAM CSV file does not seem to exist @ " +  os.path.abspath(manualcsv))

    #Check file containing previously uploaded samples.
    if not os.path.isfile(uploadedsamples):
        sys.exit("ERROR: Can't find file containing previously not uploaded sample @ " +  uploadedsamples)

    #Check for all fasta, fastq and lineage directories are in place
    for dirname in ("fastq", "fasta", "lineage"):
        dirpath = os.path.join(inputdir, runid, dirname)
        if not os.path.exists(dirpath):
            sys.exit("ERROR: No " + dirname + " directory found.")
        #Specifically check if the lineage file is in place
        if not os.path.exists(os.path.join(inputdir, runid, 'lineage', runid + "_lineage_report.txt")):
            sys.exit("ERROR: No lineage file for run " + runid + " found.")
    #Make sure region code is an accepted one

    if regioncode not in acceptedregions:
        sys.exit("ERROR: Region code " + regioncode + " is not an accepted one.")

    #Make sure labcode is an accepted one
    if labcode not in acceptedlabs:
        sys.exit("ERROR: Lab-code " + labcode + " is not an accepted one.")

    #Check that the logdir is there and accesible
    if not os.path.exists(logdir):
        sys.exit("ERROR: Can not find " + logdir + ". Perhaps you need to create it?")
    if not os.access(logdir, os.W_OK):
        sys.exit("ERROR: No write permissions in " + logdir + ".") 

    #Check that the gensam dir for sorting csv files is there and accesible
    if not os.path.exists(gensamcsvdir):
        sys.exit("ERROR: Can not find " + gensamcsvdir + ". Perhaps you need to create it?")
    if not os.access(gensamcsvdir, os.W_OK):
        sys.exit("ERROR: No write permissions in " + gensamcsvdir + ".") 
 

def sample_sheet(sspath):
    Sheet = SampleSheet(sspath)
    data = {}
    for sample in Sheet.samples:
        sample_name = sample['Sample_Name']
        sample_description = sample['Description']
        
        if len(sample_description.split("_")) != 15:
            sample_criteria= 7 #'Information saknas'
            print(sample, '----> No selection criteria given')
        else:
            sample_criteria = int(sample_description.split("_")[-1]) #eller [15]
            
        #Skip samples set to runType = 01 (desc. field 3)
        #These are samples which have been re-sequenced, so they have already been uploaded to GENSAM.
        runtype = sample_description.split("_")[2]
        if runtype == '01':
            continue

        #Skip controls
        if sample_name.startswith(('NegCtrl', 'PosCtrl', 'PosKon', 'NegKon')):
            continue
        else:
            #populate a dictionary with samples and their selection criteria value
            sample_id = sample_name.split('_')[0]
            data[sample_id] = sample_criteria
    return data
    
def writelog(logtype, message):
    now = datetime.datetime.now()
    logstring = "[" + now.strftime("%Y-%m-%d %H:%M:%S") + "] - " + logtype + " - " + message + "\n"
    return logstring

def countkeys(dictname, group):
    count = 0
    for keys in dictname:
        if dictname[keys][group]:
            count += 1
    return count

def email_error(logloc, errorstep):
    msg = EmailMessage()
    msg.set_content("Errors were encountered during the automatic upload of samples to FOHM GENSAM.\n\n" +
    "The error occured during: " + errorstep + "\n"
    "Please check the log file @  " + logloc)

    msg['Subject'] = "ERROR: GENSAM upload"
    msg['From'] = "clinicalgenomics@gu.se"
    msg['To'] = "clinicalgenomics@gu.se"

    #Send the messege
    s = smtplib.SMTP('smtp.gu.se')
    s.send_message(msg)
    s.quit()

def email_fohm(csvfile):
    print(csvfile)
    csv_filename = os.path.basename(csvfile)
    
    msg = EmailMessage()
    msg.set_content("Bifogat är en lista över uppladdade prover.\n\nMed vänliga hälsningar,\n Clinical Genomics Göteborg")

    msg['Subject'] = csv_filename
    msg['From'] = "clinicalgenomics@gu.se"
    #msg['To'] = "gensam@folkhalsomyndigheten.se"
    #msg['Cc'] = "clinicalgenomics@gu.se"
    msg['To'] = "sima.rahimi@gu.se"
    # Add the attachment
    with open(csvfile, 'r') as f:
        data = f.read()
        print(data)
        #msg.add_attachment(data, maintype='text', subtype='txt', filename=csv_filename)
        
    #Send the messege
    #s = smtplib.SMTP('smtp.gu.se')
    #s.send_message(msg)
    #s.quit()
    #with smtplib.SMTP('smtp.gu.se') as s:
    #    s.send_message(msg)

if __name__ == '__main__':
    main()
