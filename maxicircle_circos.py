# ─────────────────────────────────────────────
#  MAXICIRCLE CIRCOS PLOT
#  Requires: biopython, pycirclize
# ─────────────────────────────────────────────

import subprocess
import sys

def install(pkg):
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])

install("biopython")
install("pycirclize")

import os
import numpy as np
from Bio import SeqIO
from pycirclize import Circos
from pycirclize.parser import Gff
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
#  DETECT ENVIRONMENT (Colab or local)
# ─────────────────────────────────────────────
def is_colab():
    try:
        import google.colab
        return True
    except ImportError:
        return False

# ─────────────────────────────────────────────
#  VISUALIZATION PARAMETERS
#  Adjust these values to customize the plot
# ─────────────────────────────────────────────
window_size          = 50       # Window size for GC content and GC skew calculation (bp)
xtick_interval       = 5000     # Interval between axis tick marks (bp)
track1               = 40       # Inner radius of GC content track
track2               = 50       # Inner radius of GC skew track
track3               = 85       # Inner radius of CDS track
track4               = 90       # Boundary between forward and reverse gene tracks
track5               = 95       # Outer radius of CDS track
track_color_forward  = '#FFE082'  # Color for forward strand genes
track_color_reverse  = '#FFE082'  # Color for reverse strand genes
track_color_divergent = '#FFAB40' # Color for divergent region
plotstyle_cds        = "box"    # Gene plot style: "box" or "arrow"

# ─────────────────────────────────────────────
#  REQUEST INPUT FILES FROM USER
# ─────────────────────────────────────────────
def request_file(description, extensions):
    """Ask the user for a file path and verify it exists."""
    while True:
        filename = input(f"\nEnter the {description} file (e.g. file{extensions[0]}): ").strip()
        if not filename:
            print("[WARN] No filename entered, please try again.")
            continue
        if not os.path.isfile(filename):
            print(f"[ERROR] File '{filename}' not found in current directory.")
            print(f"        Current directory: {os.getcwd()}")
            retry = input("        Try another filename? (y/n): ").strip().lower()
            if retry != "y":
                sys.exit("[INFO] Exiting.")
        else:
            print(f"[OK] File found: {filename}")
            return filename


def upload_files_colab():
    """In Colab, allow uploading files from the local machine."""
    from google.colab import files

    print("\n[INFO] Upload your FASTA file:")
    uploaded = files.upload()
    fasta_file = list(uploaded.keys())[0]
    print(f"[OK] FASTA file: {fasta_file}")

    print("\n[INFO] Upload your GFF annotation file:")
    uploaded = files.upload()
    gff_file = list(uploaded.keys())[0]
    print(f"[OK] GFF file: {gff_file}")

    return fasta_file, gff_file


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────
def main():

    # Get input files depending on environment
    if is_colab():
        fasta_file, gff_file = upload_files_colab()
    else:
        fasta_file = request_file("FASTA", [".fasta", ".fa"])
        gff_file   = request_file("GFF annotation", [".gff", ".gff3"])

    # Read FASTA
    print("\n[INFO] Reading FASTA file...")
    genome_record = SeqIO.read(fasta_file, "fasta")
    genome_seq    = str(genome_record.seq)
    genome_length = len(genome_seq)
    print(f"[INFO] Genome length: {genome_length} bp")

    # Extract species name from FASTA header for the plot title
    # Expected format: >Genus_species ... (first two words used as species name)
    # If your header has a different format, edit the title_text variable below
    fasta_header = genome_record.description
    species_name = " ".join(fasta_header.split()[:2]).replace("_", " ")
    title_text   = r"$\it{" + species_name.replace(" ", r"\ ") + r"}$" + "\nmaxicircle"

    # Calculate GC content and GC skew
    print("[INFO] Calculating GC content and GC skew...")
    gc_contents = []
    gc_skews    = []
    pos_list_gc = np.arange(0, genome_length - window_size + 1, window_size) + window_size / 2

    for i in pos_list_gc:
        start  = int(i - window_size / 2)
        end    = int(i + window_size / 2)
        window = genome_seq[start:end]
        g = window.count('G')
        c = window.count('C')
        gc_contents.append((g + c) / window_size)
        gc_skews.append((g - c) / window_size)

    # Load GFF annotation
    print("[INFO] Loading GFF annotation...")
    gff = Gff(gff_file)

    # Create Circos object with genome length from FASTA
    circos = Circos(sectors={gff.name: genome_length})
    circos.text(title_text, size=12)

    # Create sector and main track
    sector    = circos.get_sector(gff.name)
    cds_track = sector.add_track((track3, track5))
    cds_track.axis(fc="#FFFF", ec="none")

    # Extract forward and reverse strand genes
    f_cds_feats = gff.extract_features("gene", target_strand=+1)
    r_cds_feats = gff.extract_features("gene", target_strand=-1)

    # Plot genes on forward and reverse strands
    if track4 != track5:
        cds_track.genomic_features(f_cds_feats, plotstyle=plotstyle_cds,
                                   r_lim=(track4, track5), fc=track_color_forward)
        cds_track.genomic_features(r_cds_feats, plotstyle=plotstyle_cds,
                                   r_lim=(track3, track4), fc=track_color_reverse)
    else:
        cds_track.genomic_features(f_cds_feats, plotstyle=plotstyle_cds,
                                   r_lim=(track3, track5), fc=track_color_forward)
        cds_track.genomic_features(r_cds_feats, plotstyle=plotstyle_cds,
                                   r_lim=(track3, track4), fc=track_color_reverse)

    # Plot divergent region
    d_cds_feats = gff.extract_features("divergent_region", target_strand=-1)
    cds_track.genomic_features(d_cds_feats, plotstyle="box",
                               r_lim=(87, 93), fc=track_color_divergent)

    # Add divider line between forward and reverse tracks
    divider_track = sector.add_track((track4 + 0.1, track4 + 0.2))
    divider_track.axis(fc=track_color_divergent, ec="grey", lw=2)

    # Add gene labels
    labels, label_pos_list = [], []
    for feat in gff.extract_features("gene"):
        start      = int(feat.location.start)
        end        = int(feat.location.end)
        label_pos  = (start + end) / 2
        gene_name  = feat.qualifiers.get("Name", [None])[0]
        if gene_name:
            labels.append(gene_name)
            label_pos_list.append(label_pos)

    cds_track.xticks(label_pos_list, labels, label_size=10, label_orientation="vertical")

    # GC content track (inner)
    gc_content_track = sector.add_track((track1, track2 - 5))
    gc_content_track.fill_between(pos_list_gc, gc_contents, [0] * len(gc_contents),
                                  vmin=0, vmax=max(gc_contents), color="lightblue")

    # GC skew track
    gc_skew_track = sector.add_track((track2, track3 - 5))
    gc_skew_track.fill_between(pos_list_gc, gc_skews, [0] * len(gc_skews),
                               vmin=min(gc_skews), vmax=max(gc_skews), color="violet")

    # Genomic position axis
    cds_track.xticks_by_interval(xtick_interval, outer=False,
                                 label_formatter=lambda v: f"{v/1000:.1f} Kb")

    # Save and show figure
    output_file = "maxicircle_circos.png"
    circos.savefig(output_file, dpi=1000)
    fig = circos.plotfig()
    fig.show()
    print(f"\n[OK] Figure saved as: {output_file}")

    # Download in Colab
    if is_colab():
        from google.colab import files
        files.download(output_file)
        print("[OK] Download started.")


if __name__ == "__main__":
    main()
