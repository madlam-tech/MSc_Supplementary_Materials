# Data to be written to the TSV file
data = """
GCA_000013425.1_ASM1342v1\tStaphylococcus aureus subsp. aureus NCTC 8325 reference genome ASM1342v1
GCA_000709415.1_ASM70941v1\tStaphylococcus xylosus reference genome ASM70941v1
GCA_001611955.1_ASM161195v1\tStaphylococcus haemolyticus reference genome ASM161195v1
GCA_002208805.2_ASM220880v2\tStaphylococcus pettenkoferi reference genome ASM220880v2
GCA_003571725.1_ASM357172v1\tStaphylococcus warneri reference genome ASM357172v1
GCA_003812505.1_ASM381250v1\tStaphylococcus hominis reference genome ASM381250v1
GCA_006094375.1_ASM609437v1\tStaphylococcus epidermidis reference genome ASM609437v1
GCA_016599795.1_ASM1659979v1\tStaphylococcus pasteuri reference genome ASM1659979v1
GCA_025272815.1_ASM2527281v1\tStaphylococcus capitis subsp. capitis reference genome ASM2527281v1
GCF_000006765.1_ASM676v1\tPseudomonas aeruginosa PAO1 reference genome ASM676v1
GCF_000010125.1_ASM1012v1\tStaphylococcus saprophyticus subsp. saprophyticus ATCC 15305 = NCTC 7292 genome assembly ASM1012v1
GCF_000063585.1_ASM6358v1\tClostridium botulinum A str. ATCC 3502 reference genome ASM6358v1
GCF_000195955.2_ASM19595v2\tMycobacterium tuberculosis H37Rv reference genome ASM19595v2
GCF_000246755.1_ASM24675v1\tTreponema pallidum subsp. pertenue str. SamoaD reference genome
GCF_000709415.1_ASM70941v1\tStaphylococcus xylosus reference genome ASM70941v1
GCF_003387165.1_ASM338716v1\tJeotgalicoccus halotolerans reference genome ASM338716v1
GCF_004359515.1_ASM435951v1\tMacrococcus carouselicus reference genome ASM435951v1
GCF_007814115.1_ASM781411v1\tStaphylococcus saprophyticus reference genome ASM781411v1
GCF_025272815.1_ASM2527281v1\tStaphylococcus capitis subsp. capitis reference genome ASM2527281v1
GCF_900458255.1_44343_D01\tStaphylococcus cohnii reference genome 44343_D01
"""

# Write the data to species_names.tsv
with open('species_names.tsv', 'w') as file:
    file.write(data.strip())

print("species_names.tsv file created successfully.")

