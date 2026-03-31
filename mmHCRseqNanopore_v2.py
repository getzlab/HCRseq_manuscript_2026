import pandas as pd
import subprocess
import os

def list_directory(path = "/content/", detailed = False):
    """
    List directory contents similar to 'ls -la'.

    """
    if detailed:
        cmd = ["ls", "-la", path]
    else:
        cmd = ["ls", path]

    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    print(result.stdout)

def extract_amplicon(input_fastq, output_fastq = "matched.fastq"):
    """
    Run cutadapt to extract amplicon sequences.

    """
    # Fixed parameters
    fq_out = output_fastq
    FWDPRIMER = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTGACAACCACTACCTGAG"
    RCREVPRIMER = "GCATGGACGAGCTGTACAAGTGA"
    umi_len = 15
    bc_len = 10

    cmd = [
        "cutadapt",
        "-a", f"N{{{bc_len}}}{FWDPRIMER}...{RCREVPRIMER}N{{{umi_len}}}",
        "--discard-untrimmed",
        "--action", "retain",
        "--rc",
        "--json", "cutadapt_metrics.json",
        "-o", fq_out,
        input_fastq
    ]

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def align_fastq(reference_fasta: str, input_fastq: str, output_bam: str):
  """
    Aligned the sequence processed by cutadapt to reference sequence using minimap2.

    Returns:
        The name of the .bam file

  """
  # Create intermediate file names
  mmi_index = f"{reference_fasta}.mmi"
  sam_file = f"{output_bam.replace('.bam', '')}.sam"

  # 1. Create minimap2 index
  print(f"Creating index: {mmi_index}")
  subprocess.run(["minimap2", "-d", mmi_index, reference_fasta], check=True)

  # 2. Align with minimap2
  print(f"Aligning {input_fastq} to {reference_fasta}")
  with open(sam_file, 'w') as sam_out:
      subprocess.run(
          ["minimap2", "-ax", "map-ont", reference_fasta, input_fastq],
          stdout=sam_out,
          check=True
      )

  # 3. Convert SAM to sorted BAM
  print(f"Converting SAM to sorted BAM: {output_bam}")
  view_process = subprocess.Popen(
      ["samtools", "view", "-bS", sam_file],
      stdout=subprocess.PIPE
  )
  sort_process = subprocess.Popen(
      ["samtools", "sort", "-o", output_bam],
      stdin=view_process.stdout
  )
  view_process.stdout.close()
  sort_process.communicate()

  print(f"Alignment complete!")
  return output_bam

def index_bam(bam_file: str, bai_file: str):
  """
    Create a BAM index (.bai) file using samtools.

    Returns:
        The created .bai file
    """
  print(f"Indexing {bam_file} -> {bai_file}")
  subprocess.run(["samtools", "index", bam_file, bai_file], check=True)

  print(f"Index created: {bai_file}")
  return bai_file

def align_and_index(reference_fasta: str, input_fastq: str, output_bam: str) -> tuple:
    """
    Align FASTQ to reference and create BAM + BAI files.

    Returns:
        Tuple of (bam_file, bai_file)
    """
    # Align and get BAM (using previous function)
    bam_file = align_fastq(reference_fasta, input_fastq, output_bam)

    # Create index
    bai_file = index_bam(bam_file, bam_file.replace(".bam", ".bai"))

def bam_to_fasta(
    bam_file: str,
    fasta_file: str
):
    """
    Convert BAM file to FASTA using samtools.

    Args:
        bam_file: Input BAM file
        fasta_file: Output FASTA file

    Returns:
        Path to the created FASTA file
    """
    print(f"Converting {bam_file} to FASTA: {fasta_file}")

    cmd = f"samtools fasta {bam_file} > {fasta_file}"

    # Run the command directly using shell
    result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    return fasta_file

def calculate_marker_proportion(
    fasta_file: str,
    barcode: str,
    marker: str
) -> dict:
    """
    Calculate proportion of reads containing a specific marker.

    Args:
        fasta_file: Input FASTA file
        barcode: Barcode sequence to filter reads
        marker: Marker sequence to count

    Returns:
        Dictionary with results
    """
    # Get total reads with barcode
    total_cmd = f"seqkit grep -s -p '{barcode}' {fasta_file} | seqkit stats -T | cut -f 4 | tail -n 1"
    total_reads = int(subprocess.check_output(total_cmd, shell=True, text=True).strip())

    # Get reads with both barcode and marker
    matching_cmd = f"seqkit grep -s -p '{barcode}' {fasta_file} | seqkit grep -s -p '{marker}' | seqkit stats -T | cut -f 4 | tail -n 1"
    matching_reads = int(subprocess.check_output(matching_cmd, shell=True, text=True).strip())

    # Calculate proportion
    proportion = matching_reads / total_reads if total_reads > 0 else 0

    # Print report
    print("\n=== Sequence Match Report ===")
    print(f"Total reads with barcode '{barcode}':    {total_reads}")
    print(f"Reads with marker '{marker}':           {matching_reads}")
    print(f"Proportion:                           {proportion:.4f}")

def merge_bam_files(output_bam: str, bam1: str, bam2: str):
    """
    Merge two BAM files using samtools merge.
    """
    # Build the samtools merge command
    cmd = ["samtools", "merge", output_bam, bam1, bam2]

    # Execute the command
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    print(f"Successfully merged BAM files to {output_bam}")

def process_bam_to_crosstab(bam_file, output_file="extracted_sample_info.txt"):
    """
    Process BAM file and return crosstab of first 10 bases vs contigs.
    """
    # Run samtools and awk command
    cmd = f"samtools view {bam_file} | awk '{{print substr($10, 1, 10) \"\\t\" $3}}' > {output_file}" #extract first 10 nucleotide of the sequence, and get the reporter name
    subprocess.run(cmd, shell=True, check=True)

    # Read and parse the output file
    sequences, reporters = [], [] # sample is the
    with open(output_file) as f:
        for line in f:
            sequence_10nt, reporter = line.split()
            sequences.append(sequence_10nt)
            reporters.append(reporter)

    # Create crosstab
    return pd.crosstab(sequences, reporters)

def bam_to_summary(bam_file: str, sample_info_file: str):
    """
    Process BAM file and return a data frame summarizing the number of reads in each sample (rows) and reporter (column).
    sample info should contain at least two columns with these exact names:
        SampleIfo (the names of the sample),
        IndexBarcode2 (the sample barcodes),
    """
    base_name = os.path.splitext(bam_file)[0]
    extracted_sample_barcode = f"{base_name}_extracted_sample_info.txt" # this variable will contain the extracted sample barcode and the reporter name, separated by tab. e.g. "ACAGTTCCAG	GFP_BAR_003"
    # Run samtools and awk command
    cmd = f"samtools view {bam_file} | awk '{{print substr($10, 1, 10) \"\\t\" $3}}' > {extracted_sample_barcode}" #extract first 10 nucleotide of the sequence, and get the reporter name
    subprocess.run(cmd, shell=True, check=True)

    # Read and parse the output file
    sequences, reporters = [], [] # sample is the
    with open(extracted_sample_barcode) as f:
        for line in f:
            sequence_10nt, reporter = line.split()
            sequences.append(sequence_10nt)
            reporters.append(reporter)

    # Create crosstab
    read_summary_df = pd.crosstab(sequences, reporters)

    # Get sample info and map to sample barcodes.
    sample_info_df = pd.read_csv(sample_info_file)
    read_summary_df = read_summary_df.loc[sample_info_df['IndexBarcode2']]
    read_summary_df.index = sample_info_df.set_index('IndexBarcode2').loc[read_summary_df.index]['SampleInfo'].values
    return read_summary_df

def preprocess_nanopore_data(raw_fastq: str, reference_sequence: str):
    """
    Full preprocessing pipeline for nanopore data.
    Align reads to reference and create the .bam and .bai file
    """
    base_name = os.path.splitext(raw_fastq)[0]
    extracted_fastq = f"{base_name}_matched.fastq"
    output_bam_name = f"{base_name}.bam" 
    # Step 1: Extract amplicons
    extract_amplicon(input_fastq = raw_fastq, 
                     output_fastq = extracted_fastq)

    # Step 2: Align and index
    align_and_index(reference_fasta = reference_sequence, 
                    input_fastq = extracted_fastq,
                    output_bam = output_bam_name)

### Example usage:
# extract_amplicon(input_fastq = "PJ85B5_1_rerun_0050_0061_1.fastq",
#                  output_fastq = "matched_1.fastq")
# align_and_index(reference_fasta = "experiment7.fa",
#                 input_fastq = "matched_1.fastq",
#                 output_bam = "matched_1.bam")

# extract_amplicon(input_fastq = "PJ85B5_2_rerun_0062_0072_2.fastq",
#                  output_fastq = "matched_2.fastq")
# align_and_index(reference_fasta = "experiment7.fa",
#                 input_fastq = "matched_2.fastq",
#                 output_bam = "matched_2.bam")

# merge_bam_files(output_bam = "merged.bam", bam1 = "matched_1.bam", bam2 = "matched_2.bam")

# index_bam(bam_file = "merged.bam", bai_file = "merged.bai")

# list_directory(detailed=True)

# bam_to_summary(bam_file = "merged.bam", sample_info_file = "WUS_SampleInfo_KOCellLines_11042025.csv")

# calculate_marker_proportion("matched_test_RT.fasta", "CTTATGGAAT", "AGTCTCGAGACT")

# calculate_marker_proportion("matched_test_RT.fasta", "C", "AGTCTCGAGACT")

# calculate_marker_proportion("matched_test_RT.fasta", "CTTATGGAAT", "AAGTGC")

