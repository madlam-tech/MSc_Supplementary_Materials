#!/usr/bin/env python3
__package__ = 'modules'
import sys
from . import common
import os
import statistics
import sys
from copy import deepcopy
from typing import NamedTuple, List, OrderedDict
from collections import defaultdict
from Bio import SeqIO
from enum import Enum

# All possible types of regions
RegionType = Enum('RegionType', ['ORF',
                                 'intergenic',
                                 'shortpseudo',
                                 'fragmentedpseudo',
                                 'intergenicpseudo'])

RegionInfo = NamedTuple('RegionInfo', [('contig', str),
                                       ('query', str),
                                       ('genbank_locus_tags', list),
                                       ('pseudo_locus_tag', str),
                                       ('start', int),
                                       ('end', int),
                                       ('strand', str),
                                       ('hits', List[BlastHit]),
                                       ('note', str),
                                       ('region_type', RegionType)])

# A collection of regions (ORFs and intergenic regions) on the same contig.
Contig = NamedTuple('Contig', [('regions', List[RegionInfo]),
                               ('name', str),
                               ('number', int)])




class RegionInfo_WIP:  # TODO: This is the start of changing data types into classes. Perhaps long term project, not urgent
    def __init__(self):
        self.contig = None
        self.genbank_locus_tag = None
        self.genbank_locus_tag_list = []
        self.pseudo_locus_tag = None
        self.start = None
        self.end = None
        self.strand = None
        self.hits = []
        self.source = None
        self.note = None
        self.region_type = None

    def __str__(self):
        return False #stub

    def nucleotide_length(self):
        #TODO: Make sure this is uniform across blastp/blastx and all that
        return self.end - self.start

    def ratio_gene_length_to_avg_hit_length(self):
        hit_lengths = [hit.nucleotide_length for hit in self.hits]
        avg_hit_length = sum(hit_lengths) / len(hit_lengths)
        return self.nucleotide_length() / avg_hit_length

    def lastItem(ls):
        x = ''
        for i in ls:
            if i != "":
                x = i
        return x

    def allButTheLast(iterable, delim):
        x = ''
        length = len(iterable.split(delim))
        for i in range(0, length - 1):
            x += iterable.split(delim)[i]
            x += delim
        return x[0:len(x) - 1]

    def pseudogene_reasoning(self):
        """
        Cases:
            1. Reason: Predicted fragmentation of a single gene.
            2. Reason: ORF is %s%% of the average length of hits to this gene.
            3. Reason: Intergenic region with %s blast hits.
        """
        if self.region_type == RegionType.fragmentedpseudo:
            return "Reason: Predicted fragmentation of a single gene."
        elif self.region_type == RegionType.shortpseudo:
            return "Reason: ORF is %s%% of the average length of hits to this gene." % self.ratio_gene_length_to_avg_hit_length()
        elif self.region_type == RegionType.intergenicpseudo:
            return "Reason: Intergenic region with %s blast hits." % len(self.hits)

    def write_gff_note(self):
        note = "Note=pseudogene candidate. %s" % self.pseudogene_reasoning()
        colour = "colour=229 204 255"
        locus_tag = "locus_tag=%s" % self.pseudo_locus_tag
        genbank_locus_tags = "gbk_locus_tags=%s" % ",".join(self.genbank_locus_tag_list)
        return ";".join([note, colour, locus_tag, genbank_locus_tags])

    def gff_entry(self):  # gff-version 3 compliant entry
        seqid = "gnl|Prokka|%s" % self.contig
        source = "pseudofinder"
        type = "gene"
        start = self.start
        end = self.end
        score = "."
        strand = self.strand
        phase = "."
        attributes = self.write_gff_note()

        return "\t".join([seqid, source, type, start, end, score, strand, phase, attributes])


with common.suppress_output_to_console():
    print('hello 1')
    os.system("echo hello 2")


def write_test_genome_output(file_dict, genome):
    with open(file_dict['base_filename']+"test_genome.gbk", "w") as output_handle:
        SeqIO.write(genome, output_handle, "genbank")



def get_CDSs(gbk: str, out_fasta: str) -> None:
    """Parse genbank input file for coding sequences (CDSs) and write the
    nucleotide sequences to the output file with coordinates."""

    with open(gbk, "r") as input_handle:
        with open(out_fasta, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.extract(seq_record.seq)))

    print('%s\tCDS extracted from:\t\t\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), gbk, out_fasta,)),
    sys.stdout.flush()


def get_proteome(gbk: str, out_faa: str) -> None:
    """Parse genbank input file for coding sequences (CDSs) and write them to the output file with coordinates."""

    with open(gbk, "r") as input_handle:
        with open(out_faa, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        assert len(seq_feature.qualifiers['translation']) == 1
                        output_handle.write(">%s %s %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],
                                                                 seq_record.name,
                                                                 seq_feature.location,
                                                                 seq_feature.qualifiers['translation'][0]))

    print('%s\tProteome extracted from:\t\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), gbk, out_faa,)),
    sys.stdout.flush()


def get_intergenic_regions(args, out_fasta: str) -> None:
    """Parse genbank input file for intergenic regions and write them to the output file with coordinates.

    Copied/modified from "get_interregions" by Iddo Friedberg & Ian MC Fleming
    Released under Biopython license. http://www.biopython.org/DIST/LICENSE
    The original code extracts all regions strand-dependently, even if there is a gene on the other strand
    Such strand information is not needed here, so I arbitrarily select plus strand sequence."""

    # Resets 'fasta' if it contains content already
    open(out_fasta, 'w').close()

    # Parse all contigs in the multicontig genbank
    for contig in SeqIO.parse(args.genome, "genbank"):  # contig = all information for an entire contig
        gene_list = []  # List of coding regions extracted from genbank file.
        intergenic_records = []  # List of intergenic regions that has been extracted from in between coding regions.

        for feature in contig.features:  # Loop over the contig, get the gene features on each of the strands
            if feature.type == "gene":  # Only present if prokka was run with --compliant flag
                start_position = feature.location._start.position
                end_position = feature.location._end.position
                gene_list.append((start_position, end_position))

        if args.contig_ends is True:
            # Put 'gene' at the start of the contig (position 0). This will force the next 'for loop'
            # to consider intergenic space between position '0' and the beginning of the first gene.
            gene_list.insert(0, (0, 0))

            contig_end = len(contig.seq)  # Apparently this is the fastest way to retrieve the end of a contig
            # Append a 'gene' the end of the contig. This will force the next 'for loop' to consider
            # intergenic space between the last gene and the end of the contig.
            gene_list.append((contig_end, contig_end))

        for i, gene in enumerate(gene_list):  # Compare current start position to previous end position
            last_end = gene_list[i - 1][1]
            this_start = gene_list[i][0]

            if this_start - last_end >= args.intergenic_length:  # Default 30bp.

                intergenic_region = SeqRecord(seq=contig.seq[last_end:this_start],        # Nucleotide sequence in range
                                              id="%s_ign_%d" % (contig.name, i),          # Individual ID
                                              description="%s %d-%d %s" % (contig.name,   # Description including name,
                                                                           last_end + 1,  # start position
                                                                           this_start,    # end position
                                                                           "+"))          # strand (default +)

                intergenic_records.append(intergenic_region)

        # Write to the intergenic records file
        SeqIO.write(intergenic_records, open(out_fasta, "a"), "fasta")

    print('%s\tIntergenic regions extracted from:\t%s\n'
          '\t\t\tWritten to file:\t\t\t%s.' % (current_time(), args.genome, out_fasta,)),
    sys.stdout.flush()


def parse_blast(fasta_file: str, blast_file: str, blast_format: str) -> List[RegionInfo]:
    """This function needs to take a blast query and extract the relevant information (RegionInfo)."""
    print('%s\tExtracting information from %s file.' % (current_time(), blast_format)),
    sys.stdout.flush()

    query_dict = {}  # Dictionary of information relating to each query

    # Read fasta file, build query dict without knowledge of blasthits
    with open(fasta_file, 'r') as fasta:
        lines = [line for line in fasta.readlines() if re.match("^>", line)]
        if blast_format == 'blastp':
            StatisticsDict['ProteomeOrfs'] = len(lines)
        for line in lines:
            fields_in_line = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))
            query = fields_in_line[0][1:]
            query_dict[query] = {'contig': fields_in_line[1],
                                 'query': query,
                                 'start': int(fields_in_line[2]),
                                 'end': int(fields_in_line[3]),
                                 'strand': fields_in_line[4],
                                 'hits': []}

    # Iterate through the tsv file, adding information from each line to the appropriate query
    with open(blast_file, 'r') as tsv:
        lines = tsv.readlines()
        for line in lines:
            fields_in_line = list(filter(None, re.split("\s|(?<=[0-9])-|\[|\]|:|\(|\)", line)))
            line_query = fields_in_line[0]
            blast_hit = BlastHit(query=line_query,
                                 subject_accession=fields_in_line[1],
                                 percent_ident=float(fields_in_line[2]),
                                 aa_length=int(fields_in_line[3])*3, # convert aa to nucleotide
                                 mismatch=int(fields_in_line[4]),
                                 gapopen=int(fields_in_line[5]),
                                 q_start=int(fields_in_line[6]),
                                 q_end=int(fields_in_line[7]),
                                 s_start=int(fields_in_line[8]),
                                 s_end=int(fields_in_line[9]),
                                 evalue=float(fields_in_line[10]),
                                 bitscore=float(fields_in_line[11]),
                                 stitle=" ".join(fields_in_line[12:]))

            query_dict[line_query]['hits'].append(blast_hit)

    # Write all of the information to a RegionInfo, which is then added to a list of RegionInfo
    if blast_format == 'blastp':
        note = "From BlastP;colour=51 153 102"
        region_type = RegionType.ORF
    elif blast_format == 'blastx':
        note = "From BlastX"
        region_type = RegionType.intergenic

    region_list = []
    for query in query_dict:
        region_list.append(RegionInfo(contig=query_dict[query]['contig'],
                                      query=query,
                                      genbank_locus_tags=[query],
                                      pseudo_locus_tag="",
                                      start=query_dict[query]['start'],
                                      end=query_dict[query]['end'],
                                      strand=query_dict[query]['strand'],
                                      hits=query_dict[query]['hits'],
                                      note=note,
                                      region_type=region_type))
    return region_list


def split_regions_into_contigs(lori: List[RegionInfo]) -> List[Contig]:
    """Takes a list of regions and splits them based on which contig it belongs to.
    Contig is defined above as 'List[RegionInfo]', so 'List[Contig]' is a list of lists."""

    # collects all contig names. Doesn't store duplicates, so a contig name will not be stored more than once.
    contig_names = list(set([ri.contig for ri in lori]))
    StatisticsDict['NumberOfContigs'] = len(contig_names)
    contig_list = []  # this will store the output

    for contig_name in contig_names:
        # Finds all numbers in the contig name (ie. '15' in EOKKIDHA_15) and returns them as a single integer
        contig_number = int("".join(re.findall('\d', str(contig_name))))
        regions_on_contig = []  # stores the List[RegionInfo] to be contained on a contig

        for ri in lori:
            if ri.contig == contig_name:  # if the region's contig name matches, it is added to that contig
                regions_on_contig.append(ri)

        # once all regions have been added, that list of regions is appended as a 'Contig' to the list of contigs.
        contig_list.append(Contig(regions=regions_on_contig, name=contig_name, number=contig_number))

    return contig_list


def annotate_pseudos(args, contig: Contig) -> Contig:
    """
    This function will take input blast files and return a list of all pseudogene candidates.
    """

    # 1: Look through list of regions and find individual ORFs that could be pseudogenes.
    individual_pseudos, intergenic_pseudos = check_individual_ORFs(args=args, lori=contig.regions)

    # 2: Update list of regions with any pseudogenes that were flagged from step #1.
    updated_list = replace_pseudos_in_list(pseudos=individual_pseudos+intergenic_pseudos, regions=contig.regions)

    # 3: Check adjacent regions to see if they could be pseudogene fragments.
    #   This function returns two lists: [0] = Individual pseudogenes
    #                                    [1] = Merged pseudogenes
    all_pseudos = check_adjacent_regions(args=args, lori=updated_list)

    final_regions = add_locus_tags(lori=(all_pseudos[0] + all_pseudos[1]), contig=contig.name)

    # returns both individual and merged pseudogenes as a single list, with locus tags added.
    return Contig(regions=final_regions, name=contig.name, number=contig.number)


def check_individual_ORFs(args, lori: List[RegionInfo]) -> tuple:
    """This function will take an input of regions and return two lists:
    [0]: a list of individual ORFs that could be pseudogenes.
    [1]: a list of intergenic regions that could be pseudogenes."""

    initial_blastp_list = []  # This list will contain all ORFs that have enough blast hits to be considered.
    blastp_pseudos = []      # This list will contain the resulting pseudogenes
    blastx_pseudos = []


    for region in lori:
        # Only include regions that were already as genes from whichever
        # annotation software, and that have at least 2 blast hits.
        if region.region_type == RegionType.ORF and len(region.hits) > 2:
            initial_blastp_list.append(region)

        # Include a blastx hit if it meets the minimum criteria defined by args.intergenic_threshold
        # For example, if a blastx region has 5 blast hits, the blast hitcap is 15 hits, and the threshold is 0.20,
        # the region will pass. ( 5/15 >= 0.2 ) is true.
        elif region.region_type == RegionType.intergenic and len(region.hits)/int(args.hitcap) >= args.intergenic_threshold:
            pseudo = convert_region_to_pseudo(region=region,
                                              ratio=None,  # this value is only used for BlastP-derived pseudos (below)
                                              number_of_hits=len(region.hits))
            blastx_pseudos.append(pseudo)

    for region in initial_blastp_list:

        # Retrieves lengths of genes that this region has blasted against
        list_of_database_lengths = [hit.nucleotide_length for hit in region.hits]

        # Calculates the average length of genes that this region has blasted against
        average_database_length = sum(list_of_database_lengths) / len(list_of_database_lengths)

        # Calculates the length of this region
        region_length = (region.end - region.start)

        # ratio of the region's length to the average length of hits.
        ratio = (region_length/average_database_length)

        if ratio < args.length_pseudo:
            pseudo = convert_region_to_pseudo(region=region,
                                              ratio=ratio*100,      # Multiplied by 100 to convert to percentage
                                              number_of_hits=None)  # Not important for BlastP hits
            blastp_pseudos.append(pseudo)

    return blastp_pseudos, blastx_pseudos


def convert_region_to_pseudo(region: RegionInfo, ratio: float, number_of_hits: int) -> RegionInfo:
    """Flags a region as a pseudogene by adding a note, that will appear in the GFF file.
    Regions must be explicitly rewritten because NamedTuples are immutable."""

    if region.region_type == RegionType.ORF:
        message = 'Note=pseudogene candidate. ' \
                  'Reason: ORF is %s%% of the average length of hits to this gene.;' \
                  'colour=229 204 255' % (round(ratio, 1))  # 'colour=' makes this region appear coloured in Artemis.
        pseudo_type = RegionType.shortpseudo

    elif region.region_type == RegionType.intergenic:
        message = 'Note=pseudogene candidate. ' \
                  'Reason: Intergenic region with %s blast hits.;' \
                  'colour=229 204 255' % number_of_hits  # 'colour=' makes this region appear coloured in Artemis.
        pseudo_type = RegionType.intergenicpseudo

    pseudogene = RegionInfo(contig=region.contig,
                            query=region.query,
                            genbank_locus_tags=region.genbank_locus_tags,
                            pseudo_locus_tag="",
                            start=region.start,
                            end=region.end,
                            strand=region.strand,
                            hits=region.hits,
                            note=message,
                            region_type=pseudo_type)
    return pseudogene


def replace_pseudos_in_list(pseudos: List[RegionInfo], regions: List[RegionInfo]) -> List[RegionInfo]:
    """This function prevents duplicates of regions that would occur if a gene was
    labelled and pseudogene and the original gene was not removed from the list."""

    final_list = []

    for region in regions:
        if pseudo_present(region, pseudos)[0]:  # if a pseudogene is present at the same position, write the pseudo
            final_list.append(pseudo_present(region, pseudos)[1])
        else:
            final_list.append(region)  # if it is not, write the gene
    return final_list


def pseudo_present(region: RegionInfo, pseudos: List[RegionInfo]) -> tuple:
    """Takes a particular gene and checks if that gene has been flagged as a pseudogene.
    Returns two pieces of information.
    0. If a pseudogene has been annotated at this location (bool)
    1. The identity of the (pseudo)gene at this location (RegionInfo)"""

    for pseudo in pseudos:
        if pseudo.start == region.start:
            return True, pseudo
        else:
            pass
    return False, region


def check_adjacent_regions(args, lori: List[RegionInfo]) -> tuple:
    """This function will take input blast files and return a list of all pseudogene candidates.

    lori: List of regions you want to run through.
    contig_number: the position of the contig in a list of contigs. Used for printing information.
    cutoff: refer to arg.shared_hits. Percentage of hits shared between two regions to consider joining them."""

    sorted_lori = sorted(lori, key=lambda r: r.start)
    merged_list = []  # List of merged pseudogenes stored as RegionInfo
    individual_list = []  # List of individual pseudogenes stored as RegionInfo
    i = 0   # Iterator

    while i < len(sorted_lori)-1 and len(sorted_lori) > 1:
        new_pseudo_made = False
        try:
            # compare_regions() checks that the two regions pass certain criteria
            if compare_regions(args, r1=sorted_lori[i], r2=sorted_lori[i + 1]) is True:
                new_pseudo_made = True    # this bool will be important later on in this function
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 1])   # if they pass, create a pseudogene

                # this is to keep track of overall statistics. If the regions are plain ORFs or ORFs annotated
                # as short pseudos, the counter will increase by 1 for each of them.
                for region in [sorted_lori[i], sorted_lori[i + 1]]:
                    if region.region_type == RegionType.ORF or region.region_type == RegionType.shortpseudo:
                        StatisticsDict['FragmentedOrfs'] += 1

                del sorted_lori[i + 1]  # remove items that were joined together
                del sorted_lori[i]

            # If regions [i] and [i+1] fail to join (above), look at regions [i] and [i+2].
            elif compare_regions(args=args, r1=sorted_lori[i], r2=sorted_lori[i + 2]) is True:

                new_pseudo_made = True  # this boolean will be important later on in this function
                pseudo = join_regions(sorted_lori[i], sorted_lori[i + 2])  # if they pass, create a pseudogene

                # same as above ^
                for region in [sorted_lori[i], sorted_lori[i + 1], sorted_lori[i + 2]]:
                    if region.region_type == RegionType.ORF or region.region_type == RegionType.shortpseudo:
                        StatisticsDict['FragmentedOrfs'] += 1

                del sorted_lori[i + 2]  # remove items that were joined together, and [i+1] because it's in between them
                del sorted_lori[i + 1]
                del sorted_lori[i]

            # If the pieces were not assembled but one of them is an 'individual pseudogene',
            # it is added to the individual_list
            elif sorted_lori[i].region_type == RegionType.shortpseudo or sorted_lori[i].region_type == RegionType.intergenicpseudo:
                pseudo = sorted_lori[i]
                # Deletes an item in individual_list if it has the same start position as an individual pseudo.
                individual_list[:] = [item for item in individual_list if item.start is not pseudo.start]
                individual_list.append(pseudo)

            # If the region in question fits none of the critera, move on.
            else:
                pass

        except IndexError:  # This will be triggered when 'i' equals the length of the list of pseudos
            pass

        # 'new_pseudo_made' resets to false every loop
        # so it will only be 'True' if two regions have just been merged together
        if new_pseudo_made is True:
            # Deletes an item in merged_list if that item has the same start position as the pseudogene.
            # It works like:
            #   merged_list[:] = a new version of merged_list, that contains items from merged_list,
            #   unless that item is nested within the new pseudogene.
            merged_list[:] = [item for item in merged_list if (item.start is not pseudo.start) and (item.end is not pseudo.end)]

            # Adds the merged region to a list to keep track of all merged regions
            merged_list.append(pseudo)
            # Adds the merged region to the original list so that it will continue to be considered
            sorted_lori.append(pseudo)

            # Re-sorts the list, because two regions will have been removed and one new one added (see just above).
            sorted_lori = sorted(sorted_lori, key=lambda r: r.start)

            i = i - 1  # Resets the iterator so that new region can be tested by join_regions()

        # If new_pseudo_made is False, then the iterator moves forward in the list to keep checking new regions.
        else:
            i = i + 1

        # This will remove rare cases where a pseudogene isnt handled correctly and remains in the individual_list
        # despite being a part of a merged pseudogene in merged_list. No touchy.
        individual_list[:] = [item for item in individual_list if item.start not in [pseudo.start for pseudo in merged_list]]

    # Once the loop finishes, add all statistics to StatisticsDict for reporting in the log file.
    StatisticsDict['PseudogenesTotal'] += len(individual_list) + len(merged_list)
    StatisticsDict['PseudogenesShort'] += len([item for item in individual_list if item.region_type == RegionType.shortpseudo])
    StatisticsDict['PseudogenesIntergenic'] += len([item for item in individual_list if item.region_type == RegionType.intergenicpseudo])
    StatisticsDict['PseudogenesFragmented'] += len(merged_list)

    return individual_list, merged_list


def compare_regions(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    """Takes two regions and decides if they are similar enough to join together."""

    # A list of conditions that must be met in order for two regions to be joined
    if (
        region_proximity(r1, r2) < args.distance and      # Closer than cutoff default (1000bp)
        matching_hit_critera(args, r1, r2) is True and    # Have enough matching blast hits
        r1.strand == r2.strand and                        # Same strand
        not (r1.region_type == RegionType.intergenic and r2.region_type == RegionType.intergenic)  # They are not both intergenic regions
    ):
        return True

    else:
        return False


def region_proximity(r1: RegionInfo, r2: RegionInfo) -> int:
    """Takes two regions and returns their distance from each other in # of nucleotides."""
    # sorts the two regions by starting point, so the math will always be consistent.
    sorted_by_start = sorted([r1, r2], key=lambda r: r.start)

    # substracts the end position of the first from the start position of the second
    # this value can actually be negative if a gene starts before the previous one finishes
    return sorted_by_start[1].start - sorted_by_start[0].end


def matching_hit_critera(args, r1: RegionInfo, r2: RegionInfo) -> bool:
    """This function determines if two regions meet the minimum blast hit criteria to be joined together."""

    if len(r1.hits) != 0 and len(r2.hits) != 0:
        # sorts the two regions based on number of blast hits.
        s = sorted([r1, r2], key=lambda r: len(r.hits))

        # math: (Number of shared hits) / (Total number of hits from the region with the least hits) >= cutoff value.
        if number_of_matching_hits(r1, r2)/len(s[0].hits) >= args.shared_hits:
            return True
        else:
            return False
    else:
        return False


def number_of_matching_hits(r1: RegionInfo, r2: RegionInfo) -> int:
    """This function returns the number of blast hits that two regions have in common."""

    r1_accessions = set([blasthit.subject_accession for blasthit in r1.hits])
    r2_accessions = set([blasthit.subject_accession for blasthit in r2.hits])

    return len(set(r1_accessions) & set(r2_accessions))


def join_regions(r1: RegionInfo, r2: RegionInfo) -> RegionInfo:
    """This function needs to take two regions and merge their locations."""

    # concatenates hits from both regions, discards any duplicates, and sorts them by e-value.
    merged_hits = sort_hits_by_eval(list(set(r1.hits + r2.hits)))

    merged_region = RegionInfo(contig=r1.contig,
                               query=r1.query+","+r2.query+",",
                               genbank_locus_tags=r1.genbank_locus_tags + r2.genbank_locus_tags,
                               pseudo_locus_tag="",
                               start=min([r1.start, r2.start]),
                               end=max([r1.end, r2.end]),
                               strand=r1.strand,
                               hits=merged_hits,
                               note='Note=pseudogene candidate. Reason: Predicted fragmentation of a single gene.;'
                                    'colour=229 204 255',  # 'colour=' makes this region appear coloured in Artemis.
                               region_type=RegionType.fragmentedpseudo)
    return merged_region


def sort_hits_by_eval(lobh: List[BlastHit]) -> List[BlastHit]:
    """Sorts a list of blasthits by e-value from low to high (returning the hit with the lowest evalue first)."""

    sorted_list = sorted(lobh, key=lambda r: r.evalue)
    return sorted_list


def sort_contigs(loc: List[Contig]) -> List[Contig]:
    """Takes a list of contigs and sorts it numerically."""

    sortedlist = sorted(loc, key=lambda c: c.number)

    return sortedlist


def add_locus_tags(lori: List[RegionInfo], contig: str) -> List[RegionInfo]:
    """Adds numerically increasing locus tags to a list of regions."""

    sorted_by_start = sorted(lori, key=lambda r: r.start)

    final_list = []

    for counter, region in enumerate(sorted_by_start):
        tagged_region = RegionInfo(contig=region.contig,
                                   query=region.query,
                                   genbank_locus_tags=region.genbank_locus_tags,
                                   pseudo_locus_tag=str("%s_%04d" % (contig, counter+1)),
                                   start=region.start,
                                   end=region.end,
                                   strand=region.strand,
                                   hits=region.hits,
                                   note=region.note,
                                   region_type=region.region_type)

        final_list.append(tagged_region)

    return final_list


def write_genes_to_gff(args, lopg: List[RegionInfo], gff: str) -> None:
    """Takes an input list of genes and writes them to a GFF file in proper format."""

    with open(gff, 'w') as gff_output_handle:
        # write header
        gff_output_handle.write("##gff-version 3\n#!annotation-date\t%s\n" % (current_time()))  # first line
        for i, seq_record in enumerate(SeqIO.parse(args.genome, "genbank")):  # writes one line for each contig
            entry_elements = ["##sequence-region",                # Necessary to comply with GFF3 formatting
                              "gnl|Prokka|%s" % seq_record.name,    # contig seqid
                              1,                                  # contig start
                              len(seq_record)]

            gff_output_handle.write(' '.join(map(str, entry_elements))+'\n')

        # write genes
        for region in lopg:
            locus_tag = "locus_tag=%s" % region.pseudo_locus_tag
            genbank_locus_tags = "gbk_locus_tags=%s" % ",".join(region.genbank_locus_tags)
            attributes = ";".join([region.note, locus_tag, genbank_locus_tags])
            entry_elements = ["gnl|Prokka|%s" % region.contig,
                              "pseudofinder",
                              "gene",
                              region.start,
                              region.end,
                              '.',
                              region.strand,
                              '.',
                              attributes]

            gff_output_handle.write('\t'.join(map(str, entry_elements))+'\n')


def get_intact_genes(contig: Contig, pseudos: List[RegionInfo]) -> Contig:
    """"Inspects a contig for genes that have not been annotated as pseudogenes, and returns them."""

    # All regions on a contig, sorted by start position
    region_list = sorted(contig.regions, key=lambda r: r.start)
    # Begin with all regions on the contig
    intact_list = sorted(contig.regions, key=lambda r: r.start)

    # Iterate through regions on the contig
    for region in region_list:
        for pseudo in pseudos:
            # This will be true if the region is nested within a pseudogene
            if region.start >= pseudo.start and region.end <= pseudo.end:
                intact_list.remove(region)  # Remove that region from the final list
                break  # this will speed up the function by ending the 'for pseudo' loop if a match is successful
            else:
                pass

    intact_genes = Contig(regions=intact_list, name=contig.name, number=contig.number)

    return intact_genes



#TODO: FINISH THIS BOY
def write_intact_to_fasta(infile: str, outfile: str, contigs: List[Contig]) -> None:
    """Parses a multifasta file for regions and returns them in a list."""
    #
    # parsed = SeqIO.parse(infile, 'fasta')
    # print([fasta for fasta in parsed][0].description)
    # print([fasta for fasta in parsed][0].description)
    # exit()
    with open(infile, 'r') as infile, open(outfile, 'a') as output:
        lines = infile.readlines()
        for contig in contigs:
            region_index = 0
            for line_number, line in enumerate(lines):
                try:
                    if re.match("^>%s %s" % (regions[region_index].query, contig), line):
                        output.write("%s\n%s" % (line, lines[line_number+1]))
                        region_index += 1
                except IndexError:
                    pass


def write_pseudos_to_fasta(args, pseudofinder_regions: List[RegionInfo], outfile: str) -> None:
    """Parse genbank input file for the regions provided and write them to the output file in fasta format."""

    fasta_list = []

    for contig in SeqIO.parse(args.genome, "genbank"):
        try:
            pseudofinder_regions_on_contig = [region for region in pseudofinder_regions if region.contig == contig.name]
        except IndexError:
            continue

        coord_list = [(region.start, region.end) for region in pseudofinder_regions_on_contig]
        for counter, coordinate in enumerate(coord_list):
            fasta_list.append(SeqRecord(seq=contig.seq[coordinate[0]:coordinate[1]], id="%s_%04d" % (contig.name, counter + 1),
                                        description="%s-%s +" % (coordinate[0], coordinate[1])))

    SeqIO.write(fasta_list, open(outfile, "w"), "fasta")


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def fastaReader(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def fastaReader2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def integrate_dnds(func_gff: str, pseudo_gff: str, dnds_out: str, func_faa: str, func_ffn: str, proteome: str,
                   cds: str, pseudo_fasta: str, max_ds: float, min_ds: float, max_dnds: float):
    # args = common.get_args('dnds')
    funcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    func = open(func_gff)  # reading original Intact gff outfile to dictionary
    for i in func:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            funcDict[locus]["list"] = ls
            funcDict[locus]["locus"] = locus
        else:
            pass

    pseudoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    pseudo = open(pseudo_gff)  # reading original pseudos gff outfile to dictionary
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            for j in locus.split(","):
                pseudoDict[j]["list"] = ls
                pseudoDict[j]["locus"] = locus
        else:
            pass

    # processing dnds output
    newPseudosDict = defaultdict(list)
    dndsfile = open(dnds_out)
    for i in dndsfile:
        ls = i.rstrip().split(",")
        locus = ls[0]
        if ls[2] != "dN":
            dn = float(ls[2])
            ds = float(ls[3])
            DNDS = float(ls[4])
            ref = ls[1]
            if ds < float(max_ds) and ds > float(min_ds):
                if DNDS > float(max_dnds):
                    if locus in funcDict.keys():
                        gffList = (funcDict[locus]["list"])
                        originalLocus = (funcDict[locus]["locus"])
                        commentsList = (gffList[8].split(";"))
                        comment = "Note=pseudogene candidate. Reason: dN/dS value between this gene and reference (" + ref + ") is " + str(
                            round(DNDS, 3)) + "."
                        newCommentsList = ([comment, commentsList[1], commentsList[2], commentsList[3]])
                        newCommentsLine = ";".join(newCommentsList)
                        gffLine = (("\t".join(gffList)))
                        newgffLine = (allButTheLast(gffLine, "\t") + "\t" + newCommentsLine + "\t" + str(DNDS))

                        newPseudosDict[originalLocus].append(newgffLine)
                    elif locus in pseudoDict.keys():

                        gffList = (pseudoDict[locus]["list"])
                        originalLocus = (pseudoDict[locus]["locus"])
                        commentsList = (gffList[8].split(";"))
                        newComment = (commentsList[
                                          0] + " Another Reason: dN/dS value between this gene and reference (" + ref + ") is " + str(
                            round(DNDS, 3)) + ".")
                        newCommentsList = [newComment, commentsList[1], commentsList[2], commentsList[3]]
                        newCommentsLine = ";".join(newCommentsList)
                        gffLine = (("\t".join(gffList)))
                        newgffLine = (allButTheLast(gffLine, "\t") + "\t" + newCommentsLine + "\t" + str(DNDS))

                        newPseudosDict[originalLocus].append(newgffLine)

    newPseudosDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in newPseudosDict.keys():
        dndsList = []
        for j in newPseudosDict[i]:
            dndsList.append(float(j.split("\t")[9]))
        ls = j.split("\t")
        DNDS = (statistics.mean(dndsList))
        commentsList = (ls[8].split(";"))
        comment = (commentsList[0])
        newComment = (comment.split(") is ")[0] + ") is " + str(round(DNDS, 3)))
        newCommentsList = [newComment, commentsList[1], commentsList[2], commentsList[3]]
        newCommentsLine = ";".join(newCommentsList)
        newgffLine = ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[5] + "\t" + ls[
            6] + "\t" + ls[7] + "\t" + newCommentsLine
        newPseudosDict2[i] = newgffLine

    # re-writing Intact GFF output file
    func = open(func_gff)
    newFunctionalOut = func_gff.split(".")[0] + "-new.gff"
    out = open(newFunctionalOut, "w")
    for i in func:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            if locus not in newPseudosDict2.keys():
                out.write(i.rstrip() + "\n")
        else:
            out.write(i.rstrip() + "\n")
    out.close()

    # re-writing pseudos GFF output file
    pseudo = open(pseudo_gff)
    newPseudoOut = pseudo_gff.split(".")[0] + "-new.gff"
    out = open(newPseudoOut, "w")
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = (lastItem(ls[8].split(";")))
            locus = locus.split("=")[1]
            if locus in newPseudosDict2.keys():
                out.write(newPseudosDict2[locus] + "\n")
                newPseudosDict2.pop(locus, None)
            else:
                out.write(i.rstrip() + "\n")
        else:
            out.write(i.rstrip() + "\n")

    StatisticsDict["dnds"] = len(newPseudosDict2.keys())
    for i in newPseudosDict2.keys():
        out.write(newPseudosDict2[i] + "\n")

    out.close()

    os.system("mv %s %s" % (newPseudoOut, pseudo_gff))
    os.system("mv %s %s" % (newFunctionalOut, func_gff))

    # resorting the new GFF output file
    sortDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    pseudo = open(pseudo_gff)
    out = open(newPseudoOut, "w")
    for i in pseudo:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            start = int(ls[3])
            sortDict[start] = i.rstrip()
        else:
            out.write(i.rstrip() + "\n")

    for i in sorted(sortDict.keys()):
        out.write(sortDict[i] + "\n")

    out.close()

    os.system("mv %s %s" % (newPseudoOut, pseudo_gff))

    # write intact genes to FASTA files
    faa = open(proteome)
    faa = fastaReader(faa)
    ffn = open(cds)
    ffn = fastaReader(ffn)
    ffn2 = open(cds)
    ffn2 = fastaReader2(ffn2)

    out = open(func_faa, "w")
    for i in faa.keys():
        locus = i.split(" ")[0]
        if locus in funcDict.keys():
            out.write(">" + i + "\n")
            out.write(faa[i] + "\n")
    out.close()

    out = open(func_ffn, "w")
    for i in ffn.keys():
        locus = i.split(" ")[0]
        if locus in funcDict.keys() and locus not in newPseudosDict2.keys():
            out.write(">" + i + "\n")
            out.write(ffn[i] + "\n")
    out.close()

    pseudos = open(pseudo_fasta)
    newPseudoSeqs = pseudo_fasta + "-new.fasta"
    out = open(newPseudoSeqs, "w")
    for i in pseudos:
        out.write(i.rstrip() + "\n")

    for i in newPseudosDict2.keys():
        ls = (newPseudosDict2[i].split("\t"))
        contig = ls[0].split("|")[2]
        header = contig + " " + ls[3] + "-" + ls[4] + " " + ls[6]
        seq = (ffn2[i])
        out.write(">" + header + "\n")
        out.write(seq + "\n")

    os.system("mv %s %s" % (newPseudoSeqs, pseudo_fasta))
    return len(newPseudosDict2.keys())


def add_dnds_info_to_genome(args, genome, dnds_folder):
    """Adds dN/dS info from analysis to genome"""

    dnds_list = convert_csv_to_dnds(dnds_folder+"/dnds-summary.csv")
    relevant_features = extract_features_from_genome(args, genome, 'CDS')
    feature_dict = {feature.qualifiers['locus_tag'][0]: feature for feature in relevant_features} # list comprehension
    for item in dnds_list:
        feature_dict[item.locus_tag].qualifiers['dnds'].append(item)



def convert_csv_to_dnds(dnds_file):
    """
    Reads values from the dNdS csv file and converts each row into a dnds_data entry.
    Returns a list of dnds_data
    """
    dnds_list = []
    with open(dnds_file, 'r') as csv:
        next(csv)
        for line in csv.readlines():
            fields = [common.literal_eval(x) for x in re.split(",", line.rstrip("\n"))]
            dnds_list.append(dnds_data(*fields[:-1])) # unpacking
    return dnds_list



def assign_pseudotype(feature, pseudotype):
    """
    Fully replaceable pseudotypes:
        - Input.general
        - NotPseudo.intact

    Fully irreplaceable pseudotypes:
        - NotPseudo.consumed
        - MultiIssue.sleuth

    Pseudotypes that can be combined to create a MultiIssue.general pseudotype:
        - Input.indel
        - Input.internalstop
        - any single Blast pseudotype

    Blast PseudoType logic:
        - Blast.fragment will replace any of the following:
            - truncated
            - truncated + short alignment
            - truncated + intergenic
            - intergenic
            -


    Pseudotypes that will create MultiIssue.sleuth pseudotype if they are combined with anything else:
        - Any sleuth pseudotype
    """
    current_pseudotype = feature.qualifiers['pseudo_type']
    assigned_pseudotype = pseudotype

    # If the feature is assigned to no longer be a pseudogene, that takes absolute precendence
    if type(assigned_pseudotype) is PseudoType.NotPseudo:
        feature.qualifiers['pseudo_type'] = assigned_pseudotype

    # If the feature has been consumed, do not act on it
    elif current_pseudotype is PseudoType.NotPseudo.consumed:
        pass

    # Input pseudo refinement
    elif type(current_pseudotype) is PseudoType.Input:
        if current_pseudotype is PseudoType.Input.general:
            feature.qualifiers['pseudo_type'] = assigned_pseudotype
        else:
            feature.qualifiers['pseudo_type'] = PseudoType.MultiIssue.general

    return False

def adjacent_fragments_proximity(f1, f2):
    f1_hits = set(f1.qualifiers['hits'])
    f2_hits = set(f2.qualifiers['hits'])
    common_hits = list(f1_hits & f2_hits)
    mean_common_length = sum([blasthit_length(hit, 'nt') for hit in common_hits]) / len(common_hits)

    if len(common_hits) == 0:
        return False

    else:
        f1_set = set([hit.subject_accession for hit in f1_hits])
        f2_set = set([hit.subject_accession for hit in f2_hits])

    return False # stub