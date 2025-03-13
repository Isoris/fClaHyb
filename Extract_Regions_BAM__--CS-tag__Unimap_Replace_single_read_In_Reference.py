import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Replace regions in a reference genome using alignments in a BAM file.")
parser.add_argument("--fasta", required=True, help="Input hybrid genome FASTA file.")
parser.add_argument("--bam", required=True, help="Input BAM file with alignments.")
parser.add_argument("--contigs", required=True, help="Input contig FASTA file.")
parser.add_argument("--bed", required=True, help="BED file with regions to replace.")
parser.add_argument("--output", required=True, help="Output modified FASTA file.")
args = parser.parse_args()

# Open BAM file and load hybrid genome
bamfile = pysam.AlignmentFile(args.bam, "rb")
hybrid_genome = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# Convert sequences to MutableSeq for updates
for chrom in hybrid_genome:
    hybrid_genome[chrom].seq = MutableSeq(str(hybrid_genome[chrom].seq))

# Load contig sequences from FASTA into a dictionary
contig_sequences = SeqIO.to_dict(SeqIO.parse(args.contigs, "fasta"))

# Initialize shift tracker for each chromosome
coordinate_shifts = {chrom: 0 for chrom in hybrid_genome}

# Read BED file
with open(args.bed) as bed:
    for line in bed:
        chrom, start, end = line.strip().split()[:3]
        start, end = int(start), int(end)

        # Apply cumulative shift to coordinates
        shift = coordinate_shifts.get(chrom, 0)
        adj_start = start + shift
        adj_end = end + shift

        print(f"Processing region: {chrom}:{start}-{end} (adjusted to {adj_start}-{adj_end})")

        # Check if chromosome exists in hybrid genome
        if chrom not in hybrid_genome:
            print(f"Chromosome {chrom} not found in hybrid genome.")
            continue

        # Fetch alignments overlapping the original region
        overlapping_reads = bamfile.fetch(chrom, start, end)
        replacement_seq = None

        # Find a contig that spans the entire region
        for read in overlapping_reads:
            if read.reference_start <= start and read.reference_end >= end:
                contig_name = read.query_name

                # Get the contig sequence
                if contig_name not in contig_sequences:
                    print(f"Contig {contig_name} not found in contigs FASTA.")
                    continue

                contig_seq_full = str(contig_sequences[contig_name].seq)

                # Handle reverse-complemented alignments
                if read.is_reverse:
                    contig_seq_full = str(Seq(contig_seq_full).reverse_complement())

                # Map the BED region to contig positions using aligned pairs
                aligned_pairs = read.get_aligned_pairs(matches_only=True)
                query_positions = []

                for query_pos, ref_pos in aligned_pairs:
                    if ref_pos is not None and start <= ref_pos < end:
                        query_positions.append(query_pos)

                if not query_positions:
                    print(f"No matching positions found in contig {contig_name} for region {chrom}:{start}-{end}")
                    continue

                # Extract the corresponding sequence from the contig
                contig_bed_start = min(query_positions)
                contig_bed_end = max(query_positions) + 1  # +1 because slicing is end-exclusive
                replacement_seq = contig_seq_full[contig_bed_start:contig_bed_end]

                print(f"Contig {contig_name} provides replacement sequence for {chrom}:{start}-{end}")
                print(f"Replacement sequence length: {len(replacement_seq)}")
                break  # Use the first suitable contig

        if replacement_seq is None:
            print(f"No suitable replacement sequence found for region {chrom}:{start}-{end}.")
            continue

        # Calculate length difference
        original_length = end - start
        new_length = len(replacement_seq)
        length_difference = new_length - original_length

        # Replace the region in the genome using adjusted coordinates
        seq = hybrid_genome[chrom].seq
        seq[adj_start:adj_end] = replacement_seq
        hybrid_genome[chrom].seq = seq

        print(f"Region {chrom}:{start}-{end} replaced with contig {contig_name}")
        print(f"Length difference: {length_difference}")

        # Update cumulative shift
        coordinate_shifts[chrom] = shift + length_difference
        print(f"Cumulative shift for {chrom}: {coordinate_shifts[chrom]}")

# Write the modified genome to a new FASTA file
with open(args.output, "w") as out_fasta:
    SeqIO.write(hybrid_genome.values(), out_fasta, "fasta")

print(f"Modified hybrid genome saved to {args.output}")