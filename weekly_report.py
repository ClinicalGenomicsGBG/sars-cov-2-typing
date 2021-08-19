#!/user/bin/env python

import click
import datetime
import logging
import os
import re
from collections import defaultdict
import csv
import glob
from tools.emailer import email_general

@click.command()
@click.option('--logdir', required=True,
              default='/medstore/logs/pipeline_logfiles/sars-cov-2-typing/weekly_report',
              help='Path/to/logfile/directory, uses default if path not specified')
@click.option('--nextseqdir', required=True,
              default='/medstore/results/clinical/SARS-CoV-2-typing/nextseq_data',
              help='Path/to/nextseq_data/directory, uses default if path not specified')
@click.option('--eurofinsdir', required=True,
              default='/medstore/results/clinical/SARS-CoV-2-typing/eurofins_data/goteborg',
              help='Path/to/eurofind_data/directory, uses default if path not specified')
@click.option('-o', '--outfile', required=True,
              help='Path/to/output/file')
@click.option('-e', '--send-email', is_flag=True,
              help='Send results as e-mail')
def main (logdir, nextseqdir, eurofinsdir, outfile, send_email):
    #Set up the logfile and start loggin
    now = datetime.datetime.now()
    logfile = os.path.join(logdir, "weekly_" + now.strftime("%y%m%d_%H%M%S") + ".log")
    logger = setup_logger('weekly', logfile)
    logger.info('Starting the weekly report workflow.')

    #Find all nextseq runs done
    nextseq_runs = find_nextseqruns(nextseqdir)
    
    logger.info(f'The following NextSeq runs were found: {" ".join(nextseq_runs)}')

    #Find number of samples for nextseq run
    logger.info('Finding all samples from in-house sequencing runs')
    nextseq_dict = defaultdict(lambda: defaultdict(dict))
    for run in nextseq_runs:
        #Find which week the run belongs to
        rundate = run.split("_")[0]
        runweek = datetime.datetime.strptime(rundate, '%y%m%d').isocalendar()[1]
        
        #Get all the fasta samples
        fastadir = os.path.join(nextseqdir, run, 'fasta')
        #Check if there is a fastadir, if not ignore.
        #This could be due to the pipeline having started 
        #for this sample, but haven't finished
        if os.path.exists(fastadir):
            run_num_samples = num_fasta(fastadir)
            nextseq_dict[runweek][run]['fastas'] = run_num_samples
            if run_num_samples == 0:
                logger.warning(f'Could not find any fasta files for {run}')
        else:
            logger.warning(f'No fasta directory found for run: {run}. Pipeline still not finished?')
            continue

        #Find all pangolin types for nextseq runs
        #First check that the pangolin dir exists
        lineagepath = os.path.join(nextseqdir, run, 'lineage',run + "_lineage_report.txt")
        if os.path.exists(lineagepath):
            nextseq_dict[runweek][run]['lineages'] = pangolin_types(lineagepath, ',')
        else:
            logger.warning(f'No lineage dir found for run: {run}. Pipeline still not finished?')
        
    #Find all eurofins samples which have been sequenced
    eurofins_batches = find_eurofinsruns(eurofinsdir)
    logger.info(f'Found {len(eurofins_batches)} Eurofins batches')
    
    #Number of samples
    eurofins_dict = defaultdict(lambda: defaultdict(dict))
    for batch in eurofins_batches:
        batchdate = "-".join(batch.split("-")[0:3])
        batchweek = datetime.datetime.strptime(batchdate, '%Y-%m-%d').isocalendar()[1]
        batchpath = os.path.join(eurofinsdir, batch)
     
        #Count all the fasta samples
        run_num_samples = len(glob.glob1(batchpath, "*consensus.fasta"))
        if run_num_samples == 0:
            logger.warning(f'Could not find any fasta files for {batch}')
            continue
        else:
            eurofins_dict[batchweek][batch]['fastas'] = run_num_samples
        
        #Number of all different pangolin types
        lineagefiles = glob.glob1(batchpath, "*lineage_classification.txt")

        #Trim the list of lineages files if there are multiples 
        if len(lineagefiles) == 0:
            logger.warning(f'Could not find any pangolin types for {batch}')
            continue
        elif len(lineagefiles) > 1:
            lineagefiles = trimlineages(lineagefiles)

        #Get all strains
        lineagepath = os.path.join(eurofinsdir, batch, lineagefiles[0])
        eurofins_dict[batchweek][batch]['lineages'] = pangolin_types(lineagepath, '\t')

        #Check if same number of fastas as strains
        sum_fastas = eurofins_dict[batchweek][batch]['fastas']
        sum_lineages = sum(eurofins_dict[batchweek][batch]['lineages'].values()) 
        if sum_fastas != sum_lineages:
            logger.warning(f'Sum of fasta files and strains do not match in {batch}')

    #Make an output file (csv + excel)
    #Probably keep appending to the same old file
    logger.info(f'Writing data to {outfile}')

    #Open the output file
    try:
        outf = open(outfile, "w")
    except:
        logger.error(f'Could not open {outfile} for writing output.')

    #Write the data from nextseq
    write_nextseq(nextseq_dict, outf, logger)
    # Write the data from eurofins
    write_eurofins(eurofins_dict, outf)

    #Close the output file in case it should be sent as e-mail
    outf.close()

    #Send an e-mail if flagged to do
    if send_email:
        recipients = ["anders.lind.cgg@gu.se"]

        logger.info(f'Sending e-mail to {",".join(recipients)}')
        #Set up an e-mail
        msg_to = recipients
        msg_from = "clinicalgenomics@gu.se"
        message_subject = 'COVIDSeq in-house sequencing results'
        message_body = 'The latest pangolin types from the COVIDseq can be found in the attached file.'
        attachment = outfile
        email_general(msg_to, msg_from, message_subject, message_body, attachment)

def write_eurofins(eurofins_dict, outf):
    #Print header
    outf.write('Eurofins sequencing\n')
    outf.write('Week\tBatches\tSequenced Genomes\t')
    #Find all strains sequenced so far
    all_strains = sorted(liststrains(eurofins_dict))

    outf.write("\t".join(all_strains) + "\n")

    #Print eurofins for all weeks
    for week in eurofins_dict:
        num_runs = len(eurofins_dict[week])
        num_fastas = 0
        for batch in eurofins_dict[week]:
            num_fastas += eurofins_dict[week][batch]['fastas']

        num_strains = strain_nums(all_strains, eurofins_dict, week)

        outf.write(f'{week}\t{num_runs}\t{num_fastas}')
        for strain in sorted(num_strains):
            outf.write(f'\t{num_strains[strain]}')
        outf.write("\n")

def write_nextseq(nextseq_dict, outf, logger):
    #Print header
    outf.write('In-House sequencing\n')
    outf.write('Week\tRuns\tSequenced Genomes\t')
    #Find all strains sequenced so far
    all_strains = sorted(liststrains(nextseq_dict))
    outf.write("\t".join(all_strains) + "\n")

    #Print nextdata for all weeks
    for week in nextseq_dict:
        num_runs = len(nextseq_dict[week])
        num_fastas = 0
        for run in nextseq_dict[week]:
            num_fastas += nextseq_dict[week][run]['fastas']

        num_strains = strain_nums(all_strains, nextseq_dict, week)

        #Check if number of fastas and lineages match
        #Don't print weeks where it does not.
        #print(f'{week}: {num_fastas} {sum(num_strains.values())}')
        if num_fastas != sum(num_strains.values()):
            logger.info(f'Num fastas and num strains not matching for week {week}, skipping writing output.')
            continue
        #     logger.warning(f'Number of fastas in {week} is not same as num of strains. '
        #                    f'Fastas: {num_fastas}, strains: {sum(num_strains.values())}')

        outf.write(f'{week}\t{num_runs}\t{num_fastas}')
        for strain in sorted(num_strains):
            outf.write(f'\t{num_strains[strain]}')
        outf.write("\n")

def trimlineages(lineage_list):
    for element in lineage_list:
        if element.startswith('20'):
            lineage_list.remove(element)

    return lineage_list

def strain_nums(strainlist, nextseq_dict, week):
    #Build a dict with all seen strains
    straindict = {}
    for strain in strainlist:
        straindict[strain] = 0

    #for week in nextseq_dict:
    for run in nextseq_dict[week]:
        if 'lineages' in nextseq_dict[week][run].keys():
            for strain in nextseq_dict[week][run]['lineages']:
                straindict[strain] += nextseq_dict[week][run]['lineages'][strain]


    return straindict

def liststrains(nextseq_dict):
    strains = []
    for week in nextseq_dict:
        for run in nextseq_dict[week]:
            if 'lineages' in nextseq_dict[week][run].keys():
                for strain in nextseq_dict[week][run]['lineages']:
                    if not strain in strains:
                        strains.append(strain)
    return strains

def pangolin_types (lineagepath, delim):
    pango_dict = {}
    #Collect number of each found strain
    with open(lineagepath) as csv_file:
        next(csv_file) #Skip header
        csv_reader = csv.reader(csv_file, delimiter=delim)
        for row in csv_reader:
            taxon = row[0]
            strain = row[1]
            #Skip negative controls
            if taxon.lower().startswith('consensus_neg'):
                continue
            #Skip positive controls
            elif taxon.lower().startswith('consensus_pos'):
                continue
            elif strain.lower().startswith('none'):
                continue
            else:
                if strain in pango_dict:
                    pango_dict[strain] += 1
                else:
                    pango_dict[strain] = 1

    return pango_dict
    
def num_fasta (fastadir):
    num_fasta = 0
    for fasta in os.listdir(fastadir):
        if not fasta.lower().startswith('neg') \
                and not fasta.lower().startswith('pos')\
                and fasta.lower().endswith('.fa'):
            num_fasta += 1
    return num_fasta

def find_nextseqruns (nextseqdir):
    dir_list = []
    #Find all runs in folder
    for dirname in next(os.walk(nextseqdir))[1]:
        if re.match('^2[1-4]', dirname) and len(dirname.split("_")) == 4:
            dirdate = dirname.split("_")[0]
            dir_list.append(dirname)

            #This is a remnant from when this function only returned runs from last week
            #d = datetime.datetime.strptime(dirdate, "%y%m%d")
            #if ((d - now).days) >= -7:
            #    dir_list.append(dirname)

    return dir_list

def find_eurofinsruns (nextseqdir):
    dir_list = []
    #Find all runs in folder
    for dirname in next(os.walk(nextseqdir))[1]:
        dir_list.append(dirname)

    return dir_list

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

if __name__ == '__main__':
    main()    
