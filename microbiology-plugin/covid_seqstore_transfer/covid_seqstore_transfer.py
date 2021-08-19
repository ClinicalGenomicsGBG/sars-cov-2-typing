import re
import os
import csv
import json
import glob
import shutil
import subprocess

import config

from ion.plugin import IonPlugin, PluginCLI, RunType, RunLevel
from ion.utils import blockprocessing


class ProgressRender:
    def __init__(self):
        self.lines = []
        self.output = 'progress_block.html'  # Relative path to plugin working directory

    def add_subheader(self, line):
        subheader = '<h3>{}</h3>'.format(line)
        self.lines.append(subheader)

    def add_line(self, line):
        """Add line to be rendered in .html."""
        self.lines.append(line)

    def clear(self):
        """Remove previous lines."""
        self.lines = []

    def render(self):
        """Render progress by writing to .html."""
        html_pre = '<html><body>'
        html_post = '</body></html>'

        joined_lines = '<br>'.join(self.lines)
        html = html_pre + joined_lines + html_post

        with open(self.output, 'w') as out:
            out.write(html)

        return html


class SampleCollection:
    """
    Class for handling all samples part of a run.

    Barcode-samplename mapping is done with barcode as the main access point.
    In theory, a sample may be spread out over several barcodes.
    Plugin results and output names are generated on a barcode basis.
    """
    def __init__(self):
        self.sample_info = {}

    def add_sample(self, barcode, sample_name):
        self.sample_info[barcode] = sample_name

    def remove_sample(self, barcode):
        return self.sample_info.pop(barcode)

    @property
    def barcodes(self):
        return self.sample_info.keys()

    @property
    def sample_names(self):
        return self.sample_info.values()

    def __iter__(self):
        return self.sample_info


def clean_samplename(sample_name):
    """Return clean and formatted samplename in the event of typos."""
    for character in [' ']:
        sample_name.replace(character, '')

    for character in ['_']:
        sample_name.replace(character, '-')

    return sample_name


def is_covid_sample(sample_name):
    return re.search('D[A-Z]2[0-9]-[0-9]+', sample_name, re.IGNORECASE)


def find_plugin_outputs(plugin_name, root_path):
    plugin_output_paths = glob.glob(os.path.join(root_path, '{}_*'.format(plugin_name)))

    plugin_outputs = {}
    for plugin_output_path in plugin_output_paths:
        index = plugin_output_path.split('.')[-1]
        plugin_outputs[index] = plugin_output_path

    return plugin_outputs


def parse_pangolin_csv(pangolin_csv_path):
    """Parse csv into more fun to read format."""
    collected_info = {}
    with open(pangolin_csv_path, 'r') as inp:
        csv_handle = csv.DictReader(inp)

        for row in csv_handle:
            barcode = row.pop('Barcode')
            collected_info[barcode] = row

    return collected_info


def read_fastas(fasta_path):
    fasta_collection = {}

    with open(fasta_path, 'r') as fa:
        fasta_name = ''
        sequence = ''

        for line in fa.readlines():
            if line.startswith('>'):
                if fasta_name:
                    fasta_collection[fasta_name] = sequence
                    sequence = ''

                # >IonCode_0101_DC21-45858
                fasta_name = line.strip().replace('>', '').rsplit('_', 1)[0]
            else:
                sequence += line.strip()

    return fasta_collection


def read_vcf_paths(plugin_path):
    vcf_collection = {}

    all_paths = glob.glob(os.path.join(plugin_path, '*'))

    regex = 'IonCode_0[0-9]{3}$'
    sample_outputs = [path for path in all_paths if re.search(regex, path)]

    for sample_output in sample_outputs:
        barcode = os.path.basename(sample_output)
        vcf_path = os.path.join(sample_output, 'TSVC_variants_{}.vcf.gz'.format(barcode))
        vcf_collection[barcode] = vcf_path

    return vcf_collection


class covid_seqstore_transfer(IonPlugin):
    version = "0.1.2.1"
    runtypes = [RunType.COMPOSITE]
    runlevel = [RunLevel.LAST]
    depends = [config.pangolin_plugin_name, config.variant_caller_name]

    def launch(self):
        with open('barcodes.json', 'r') as fh:
            self.barcodes = json.load(fh)  # NOTE: Not available in memory like startplugin

        root_report_path = self.startplugin['runinfo']['report_root_dir']

        run_name = self.startplugin['expmeta']['run_name']
        result_name = self.startplugin['expmeta']['results_name']
        # instrument_name = self.startplugin['expmeta']['instrument']
        # experiment_id = self.startplugin['expmeta']['results_name'].split('_')[-1]

        # Divide found samples into covid & non-covid, add barcode information
        covid_samples = SampleCollection()
        non_covid_samples = SampleCollection()
        for sample_name, barcode_info in self.startplugin['plan']['barcodedSamples'].items():
            sample_name = clean_samplename(sample_name)
            sample_barcodes = barcode_info['barcodes']

            if not is_covid_sample(sample_name):
                for barcode in sample_barcodes:
                    non_covid_samples.add_sample(barcode, sample_name)

            for barcode in sample_barcodes:
                covid_samples.add_sample(barcode, sample_name)

        # Filter away samples without matching metadata
        no_metadata_samples = SampleCollection()
        for barcode, sample_name in covid_samples.sample_info.copy().items():  # NOTE: Copy as we remove what we iterate through otherwise
            try:
                self.startplugin['pluginconfig']['input_metadata'][sample_name]  # NOTE: This might cause an issue with misspellings since we clean names previously
            except KeyError:
                no_metadata_samples.add_sample(barcode, sample_name)
                covid_samples.remove_sample(barcode)

        # Parse for plugin outputs
        root_plugin_output_path = os.path.join(root_report_path, 'plugin_out')

        # Pangolin
        plugin_name = config.pangolin_plugin_name
        plugin_outputs = find_plugin_outputs(plugin_name, root_plugin_output_path)
        latest_plugin_output_id = max(plugin_outputs, key=plugin_outputs.get)
        latest_plugin_output_path = plugin_outputs[latest_plugin_output_id]

        pangolin_csv_path = os.path.join(latest_plugin_output_path, '{}.xls'.format(run_name))  #NOTE: It's actually a csv
        pangolin_csv_info = parse_pangolin_csv(pangolin_csv_path)

        # Filter away samples that failed pangolin QC
        qc_failed_samples = SampleCollection()
        for barcode, pangolin_result in pangolin_csv_info.items():
            if barcode not in covid_samples.barcodes:  # Already filtered away
                continue

            if pangolin_result['status'] != 'passed_qc' or pangolin_result['passes'] != 'Passed':
                sample_id = covid_samples.remove_sample(barcode)
                qc_failed_samples.add_sample(barcode, sample_id)
                continue

        pangolin_fasta_path = os.path.join(latest_plugin_output_path, '{}.fasta'.format(run_name))
        # Read fastas into memory on barcode keys
        sample_fastas = read_fastas(pangolin_fasta_path)

        # Variant caller
        plugin_name = config.variant_caller_name
        plugin_outputs = find_plugin_outputs(plugin_name, root_plugin_output_path)
        latest_plugin_output_id = max(plugin_outputs, key=plugin_outputs.get)
        latest_plugin_output_path = plugin_outputs[latest_plugin_output_id]

        # Read vcf paths into memory on barcode keys
        sample_vcfs = read_vcf_paths(latest_plugin_output_path)

        # Dump data at designated location, under result name
        output_path = os.path.join(config.root_dump_path, result_name)
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        transfer_failed_samples = SampleCollection()
        for sample_barcode, sample_name in covid_samples.sample_info.items():
            try:
                sample_fasta_sequence = sample_fastas[sample_barcode]
                sample_fasta_output_path = os.path.join(output_path, '{}.fa'.format(sample_name))

                # FASTA
                with open(sample_fasta_output_path, 'w') as out:
                    out.write('>{}\n'.format(sample_barcode))
                    out.write(sample_fasta_sequence)

                # VCF
                sample_vcf_path = sample_vcfs[sample_barcode]
                sample_vcf_name = '{}.vcf.gz'.format(sample_name)
                sample_vcf_output_path = os.path.join(output_path, sample_vcf_name)
                shutil.copyfile(sample_vcf_path, sample_vcf_output_path)

                # FASTQ
                sample_bam_path = self.barcodes[sample_barcode]['bam_filepath']
                sample_fastq_name = '{}.fastq'.format(sample_name)
                sample_fastq_output_path = os.path.join(output_path, sample_fastq_name)
                if not sample_fastq_output_path:
                    command = blockprocessing.bam2fastq_command(sample_bam_path, sample_fastq_output_path)
                    subprocess.check_call(command, shell=True)  # NOTE Security issue. Enables shell commands when user defines sample names
                    subprocess.check_call(['gzip', sample_fastq_output_path])  # Edits in place

            except Exception as e:
                self.log.error(e)
                transfer_failed_samples.add_sample(sample_barcode, sample_name)
                covid_samples.remove_sample(sample_barcode)

        # Dump pangolin result
        pangolin_csv_output_path = os.path.join(output_path, 'pangolin.csv')
        shutil.copyfile(pangolin_csv_path, pangolin_csv_output_path)

        # Dump metadata
        metadata = self.startplugin['pluginconfig']['input_metadata']
        # Remove failed prior to dump
        for _, sample_name in qc_failed_samples.sample_info.items():
            try:
                metadata.pop(sample_name)
            except KeyError:  # Not all failed may have metadata
                continue
        metadata_output_path = os.path.join(output_path, 'metadata.json')
        with open(metadata_output_path, 'w') as out:
            json.dump(metadata, out, indent=4)

        # Setup class for rendering html to user on plugin page
        progress_renderer = ProgressRender()

        progress_renderer.add_subheader('Transferred Samples:')
        for sample_barcode, sample_name in covid_samples.sample_info.items():
            metadata = self.startplugin['pluginconfig']['input_metadata'][sample_name]
            progress_renderer.add_line('\t'.join([sample_barcode, sample_name, metadata]))

        progress_renderer.add_subheader('Transfer Failed Samples:')
        for sample_barcode, sample_name in transfer_failed_samples.sample_info.items():
            progress_renderer.add_line('\t'.join([sample_barcode, sample_name]))

        progress_renderer.add_subheader('QC Failed Samples:')
        for sample_barcode, sample_name in qc_failed_samples.sample_info.items():
            progress_renderer.add_line('\t'.join([sample_barcode, sample_name]))

        progress_renderer.add_subheader('No Metadata Samples:')
        for sample_barcode, sample_name in no_metadata_samples.sample_info.items():
            progress_renderer.add_line('\t'.join([sample_barcode, sample_name]))

        progress_renderer.render()

        return True


if __name__ == "__main__":
    PluginCLI()
