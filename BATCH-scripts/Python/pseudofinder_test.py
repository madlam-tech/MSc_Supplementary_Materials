#!/usr/bin/env python3

from .common import get_args, bold, print_with_time
import os
import sys
from os.path import dirname
import re
import subprocess
from collections import OrderedDict
from time import strftime


def manage_folders(path: str):
    """Creates a folder to store test results and returns the name."""

    time = strftime("%Y%m%d")
    folder_str = time + "_results_"
    folder_num = 1

    try:  # simplest solution
        os.makedirs(path + folder_str + str(folder_num))

    except FileExistsError:  # unless the folder already exists.
        current_folders = os.listdir(path)
        results_folders = [folder for folder in current_folders if re.match(time, folder)]
        folder_numbers = [int(x.replace(folder_str, '')) for x in results_folders]
        biggest_folder_number = sorted(folder_numbers, key=int, reverse=True)[0]
        folder_num = biggest_folder_number + 1
        os.makedirs(path + folder_str + str(folder_num))

    return path + folder_str + str(folder_num)


def test_command(command_name: str, full_command: str):
    """
    Tests the given pseudofinder command to make sure that the command runs without an error.
    """
    print_with_time(bold('Testing the %s command.' % command_name))
    print_with_time('Full shell command: %s' % full_command)

    try:
        subprocess.run(full_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print_with_time('Command failure: ' + command_name)
        exit()


def main():
    # Inputs for all the test commands
    args = get_args('test')

    python = sys.executable
    path_to_master = dirname(dirname((__file__)))
    path_to_test_data = path_to_master + "/test/"
    path_to_pseudofinder = path_to_master + "/pseudofinder.py"

    if args.outfolder:
        if not args.outfolder.endswith("/"):
            args.outfolder = args.outfolder + "/"
        folder_name = manage_folders(args.outfolder)
    else:
        folder_name = manage_folders(path_to_test_data)

    output_prefix = folder_name + "/test"
    reannotate_output_prefix = output_prefix + "_reannotate"
    blastp_file = output_prefix + "_proteome.faa.blastP_output.tsv"
    blastx_file = output_prefix + "_intergenic.fasta.blastX_output.tsv"
    log_file = output_prefix + "_log.txt"

    if args.genome is None:
        args.genome = path_to_test_data + "candidatus_tremblaya_princeps_PCIT.gbf"
        print_with_time('No genome specified, will use default test genome: %s' % args.genome)
    if args.diamond:
        diamond_param = ' --diamond'
    else:
        diamond_param = ''

    ctl_full_path = path_to_master + "/codeml-2.ctl"
    ref_pep = path_to_test_data + "Mycobacterium_tuberculosis_H37Rv-subset.faa"
    ref_nuc = path_to_test_data + "Mycobacterium_tuberculosis_H37Rv-subset.ffn"
    pep = path_to_test_data + "Mycobacterium_leprae-subset.faa"
    nuc = path_to_test_data + "Mycobacterium_leprae-subset.ffn"
    dndsOutput = "dnds_output"

    # A dictionary to store the names and shell commands for each section of pseudofinder
    command_dict = OrderedDict()
    command_dict['Annotate'] = f"{python} {path_to_pseudofinder} annotate -g {args.genome} -db {args.database} " \
                               f"-op {output_prefix} -t {args.threads} {diamond_param}"

    command_dict['Reannotate'] = f"{python} {path_to_pseudofinder} reannotate -g {args.genome} -log {log_file} " \
                                 f"-op {reannotate_output_prefix}"

    command_dict['Visualize'] = f"{python} {path_to_pseudofinder} visualize -g {args.genome} -op {output_prefix} " \
                                f"-p {blastp_file} -x {blastx_file} -log {log_file}"

    for command in command_dict:
        test_command(command_name=command, full_command=command_dict[command])


if __name__ == '__main__':
    main()

