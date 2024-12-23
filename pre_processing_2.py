import pandas as pd

refrance_genome_file_path = 'EPI_ISL_402124-ORF1ab.fasta'
processed_bam_file_path = 'processed_bam_data.csv'
df = pd.read_csv(processed_bam_file_path)

with open(refrance_genome_file_path, 'r') as file:
    # Read each line in the file
    lines = [line for line in file]
    refrance_gnome = ''
    for i in range(1,len(lines)):
        refrance_gnome += lines[i]
        # print(line)

print(refrance_gnome)