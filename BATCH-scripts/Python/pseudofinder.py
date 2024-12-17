#!/usr/bin/env python3

"""
pseudofinder.py: A script to find pseudogene candidates in annotated genbank files.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.

Installation requirements: Please see github repo for detailed explanation of requirements.
 """
__author__ = "Mitch Syberg-Olsen, Arkadiy Garber & Filip Husnik"
__version__ = "1.1.0"
__maintainer__ = "Filip Husnik"
__email__ = "filip.husnik@gmail.com"

try:
    from sys import argv, stderr
    from modules import *


except ModuleNotFoundError as err:
    stderr.write(f"Pseudofinder encountered an error: {str(err)}.\n"
                 f"If you are a conda user, please enable the Pseudofinder environment with the following command:\n"
                 f"conda activate pseudofinder\n"
                 f"If you are not a conda user, please check that you have installed the required dependencies.\n")
    exit()

except ImportError as err:
    stderr.write(f"Pseudofinder encountered an error: {str(err)}.\n"
                 f"Please check that no files are missing from the pseudofinder/modules folder.\n")
    exit()

menu_string = "pseudofinder.py [ annotate | reannotate | visualize | sleuth | breaker | test | help | version | citation ]"
errorMessage = f"Options: {common.bold(menu_string)}\n"

try:
    module = argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

try:
    if module == "annotate":
        annotate.main()

    elif module == "reannotate":
        reannotate.main()

    elif module == "visualize":
        visualize.main()

    elif module == "sleuth":
        sleuth.main()

    elif module == "breaker":
        breaker.main()

    elif module == "test":
        pseudofinder_test.main()

    elif module == "help":
        stderr.write("\tpseudofinder.py annotate: Flags candidate pseudogenes.\n"
                     "\tpseudofinder.py reannotate: Begins the annotate pipeline post-BLAST.\n"
                     "\tpseudofinder.py visualize: Generates a 3D plot to visualize different combinations of "
                     "settings.\n"
                     "\tpseudofinder.py sleuth: pairwise comparison against a reference genome."
                     "Pseudogenes inferred from relaxed selection.\n"
                     "\tpseudofinder.py breaker: generates random potentially pseudogene-inducing mutations in a"
                     "subset of genes, and generates a new set of contigs with the mutated genes\n"
                     "\tpseudofinder.py test: Runs all commands on a test dataset and checks that the outputs"
                     "are as expected.\n")

    elif module == "version":
        print(f"Pseudofinder version {__version__}")

    elif module == "citation":
        print(f"Syberg-Olsen MJ*, Graber AI*, Keeling PJ, McCutcheon JP, Husnik F. Pseudofinder: detection of "
              f"pseudogenes in prokaryotic genomes, bioRxiv 2021, doi: https://doi.org/10.1101/2021.10.07.463580. "
              f"GitHub repository: https://github.com/filip-husnik/pseudofinder/.")

    # calling directly is only for dev purposes. For users, it is always called within annotate module.
    elif module == "interactive":
        interactive.main()

    else:
        stderr.write(errorMessage)
        exit()

except KeyboardInterrupt:
    print("\n")
    common.print_with_time("Pseudofinder process cancelled by user. Exiting now.")
    exit()
