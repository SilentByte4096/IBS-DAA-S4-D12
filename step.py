#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
De Bruijn graph assembly with optional long-read scaffolding, FASTQC analysis from raw data,
optional BUSCO alignment, and quality metrics.
"""

import os
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
import networkx as nx
import numpy as np
import subprocess
import shutil
import zipfile

# =============================================================================
# Constants (configurable)
# =============================================================================

# 1. Paths to input/output files
INPUT_FILE = r"C:\Users\nehar\OneDrive\Desktop\Done\Hope\code\input.fasta"  
REF_FILE = r"C:\Users\nehar\OneDrive\Desktop\Done\Hope\code\final.fasta"  
OUTPUT_FILE = "assembled_genome.fasta"           # Output FASTA of assembled scaffold(s)
REPORT_FILE = "assembly_report.pdf"              # Optional: if generating a report
PLOTS_DIR = "plots"                              # Directory to save coverage or k-mer plots

# 2. Short-read assembly config (ignored if long-read only)
K_SIZE_SHORT = 3                                 # ✅ 3 is fine for testing, use 21-55 for real data
MIN_COVERAGE_SHORT = 1                           # ✅ 1 = include all k-mers (no filtering)
READ_THRESHOLD = 5                               # ✅ Ignore reads shorter than this (good default)

# 3. Optional tools/config
FASTQC_DIR = "fastqc_output"                     # Dir to store FastQC results (if used)
RUN_BUSCO = "no"                                 # "yes" or "no" to run BUSCO analysis
ASSEMBLY_TYPE = "short"                           # "short", "long", or "hybrid"


# Ensure plots directory exists
if not os.path.exists(PLOTS_DIR):
    os.makedirs(PLOTS_DIR)

# =============================================================================
# 1. Read and Classify FASTA
# =============================================================================

def read_and_classify_fasta(fasta_file):
    """Read FASTA file and classify reads into short and long based on length."""
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Input file {fasta_file} not found.")
    
    short_reads = []
    long_reads = []
    current_seq = ""
    current_header = ""

    try:
        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        seq_length = len(current_seq)
                        if seq_length < READ_THRESHOLD:
                            short_reads.append(current_seq)
                        else:
                            long_reads.append(current_seq)
                    current_header = line[1:]
                    current_seq = ""
                else:
                    current_seq += line.upper()
            if current_seq:
                seq_length = len(current_seq)
                if seq_length < READ_THRESHOLD:
                    short_reads.append(current_seq)
                else:
                    long_reads.append(current_seq)

        print(f"Step 1: Read {len(short_reads) + len(long_reads)} sequences from {fasta_file}. Classified {len(short_reads)} short reads and {len(long_reads)} long reads completed.")
        return short_reads, long_reads

    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return [], []

# =============================================================================
# 2. De Bruijn Graph for Short Reads
# =============================================================================

def cut_kmer(sequence, kmer_size):
    """Cut a sequence into a k-mer iterator."""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def build_kmer_dict(reads, kmer_size):
    """Return a dict of kmer counts from reads."""
    if not reads:
        return {}
    kmer_count = {}
    for sequence in reads:
        kmers = cut_kmer(sequence, kmer_size)
        for kmer in kmers:
            kmer_count[kmer] = kmer_count.get(kmer, 0) + 1
    print(f"Step 2: Built k-mer dictionary with {len(kmer_count)} unique k-mers for short reads completed.")
    return kmer_count

def build_graph(kmer_count):
    """Build a networkx.DiGraph from a dict of kmer counts."""
    if not kmer_count:
        return nx.DiGraph()
    graph = nx.DiGraph()
    for kmer, count in kmer_count.items():
        if count >= MIN_COVERAGE_SHORT:
            node_1 = kmer[:-1]
            node_2 = kmer[1:]
            graph.add_edge(node_1, node_2, weight=count)
    print(f"Step 3: De Bruijn graph built with {len(graph.nodes)} nodes and {len(graph.edges)} edges completed.")
    return graph

def get_starting_nodes(graph):
    """Get the list of starting nodes."""
    return [node for node in graph.nodes() if not list(graph.predecessors(node))]

def get_sink_nodes(graph):
    """Get the list of sink nodes."""
    return [node for node in graph.nodes() if not list(graph.successors(node))]

def get_contigs(graph, starting_nodes, sink_nodes):
    """Return a list of (contigs, len(contigs)) tuples."""
    contigs = []
    for start_node in starting_nodes:
        for sink_node in sink_nodes:
            try:
                path = nx.shortest_path(graph, start_node, sink_node)
                contig = "".join([node[0] for node in path[:-1]] + [path[-1]])
                contigs.append((contig, len(contig)))
            except (nx.NodeNotFound, nx.NetworkXNoPath):
                continue
    print(f"Step 4: Generated {len(contigs)} contigs from short reads completed.")
    return contigs

# =============================================================================
# 3. Long Read Scaffolding 
# =============================================================================

def find_best_overlap(sequence1, sequence2, min_overlap=2):
    """Find the best overlap between two sequences."""
    best_overlap = 0
    for overlap_len in range(min(len(sequence1), len(sequence2)), min_overlap - 1, -1):
        if sequence1[-overlap_len:] == sequence2[:overlap_len]:
            best_overlap = overlap_len
            break
    return best_overlap

def scaffold_contigs(contigs, long_reads, min_overlap=2):
    """Scaffold short-read contigs or long reads directly using overlaps."""
    
    if not long_reads:
        return contigs if contigs else []

    if not contigs:
        if not long_reads:
            return []

        overlap_graph = nx.DiGraph()
        for i, read in enumerate(long_reads):
            overlap_graph.add_node(i, sequence=read)

        for i, read1 in enumerate(long_reads):
            for j, read2 in enumerate(long_reads):
                if i != j:
                    overlap_len = find_best_overlap(read1, read2, min_overlap)
                    if overlap_len >= min_overlap:
                        overlap_graph.add_edge(i, j, overlap=overlap_len)


        starting_nodes = [n for n in overlap_graph.nodes() if overlap_graph.in_degree(n) == 0]
        scaffolds = []

        for start in starting_nodes:
            path = [start]
            current = start
            while True:
                next_nodes = [n for n in overlap_graph.successors(current) if n not in path]
                if not next_nodes:
                    break
                next_node = max(next_nodes, key=lambda x: overlap_graph[current][x]['overlap'], default=None)
                if next_node is None:
                    break
                path.append(next_node)
                current = next_node

            scaffold_seq = long_reads[path[0]]
            for i in range(1, len(path)):
                prev_read = long_reads[path[i-1]]
                curr_read = long_reads[path[i]]
                overlap_len = overlap_graph[path[i-1]][path[i]]['overlap']
                scaffold_seq += curr_read[overlap_len:]

            scaffolds.append((scaffold_seq, len(scaffold_seq)))

        unique_scaffolds = []
        seen = set()
        for scaffold, length in scaffolds:
            if scaffold not in seen:
                seen.add(scaffold)
                unique_scaffolds.append((scaffold, length))

        unique_scaffolds.sort(key=lambda x: x[1], reverse=True)

        print(f"Step 5: Scaffolding completed with {len(unique_scaffolds)} unique scaffolds using graph-based long read assembly.")
        return unique_scaffolds

    else:
        return contigs

# =============================================================================
# 4. Graph Simplification
# =============================================================================

def path_average_weight(graph, path):
    """Compute the average weight on a path."""
    total = 0
    for node_1, node_2 in zip(path[:-1], path[1:]):
        total += graph[node_1][node_2].get("weight", 0)
    return total / (len(path) - 1) if total else 0

def remove_paths(graph, paths, delete_entry_node, delete_sink_node):
    """Remove specified paths from the graph."""
    for path in paths:
        graph.remove_nodes_from(path[(not delete_entry_node):(None if delete_sink_node else -1)])
    return graph

def select_best_path(graph, paths, path_lengths, avg_path_weights, delete_entry_node=False, delete_sink_node=False):
    """Select the best path based on weight and length."""
    random.seed(9001)
    if not avg_path_weights:
        return graph
    max_weight = max(avg_path_weights)
    best_weight_indexes = [i for i, w in enumerate(avg_path_weights) if w == max_weight]
    best_lengths = [path_lengths[i] for i in best_weight_indexes]
    best_path_index = best_weight_indexes[best_lengths.index(max(best_lengths))]
    return remove_paths(graph, paths[:best_path_index] + paths[best_path_index+1:], delete_entry_node, delete_sink_node)

def simplify_bubbles(graph):
    """Remove bubbles from the graph."""
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    for start in starting_nodes:
        for sink in sink_nodes:
            paths = list(nx.all_simple_paths(graph, start, sink))
            if len(paths) > 1:
                path_lengths = [len(p) for p in paths]
                avg_path_weights = [path_average_weight(graph, p) for p in paths]
                graph = select_best_path(graph, paths, path_lengths, avg_path_weights)
    print(f"Step 6a: Bubble simplification completed. Graph now has {len(graph.edges)} edges.")
    return graph

def solve_entry_tips(graph):
    """Remove entry tips."""
    starting_nodes = get_starting_nodes(graph)
    tips = []
    for start in starting_nodes:
        path = [start]
        successors = list(graph.successors(start))
        while len(successors) == 1 and not list(graph.predecessors(successors[0])):
            path.append(successors[0])
            successors = list(graph.successors(successors[0]))
        if len(path) > 1:
            tips.append(path)
    if tips:
        path_lengths = [len(p) for p in tips]
        avg_path_weights = [path_average_weight(graph, p) for p in tips]
        graph = select_best_path(graph, tips, path_lengths, avg_path_weights, delete_entry_node=True)
    print(f"Step 6b: Entry tip removal completed. Graph now has {len(graph.edges)} edges.")
    return graph

def solve_out_tips(graph):
    """Remove out tips."""
    sink_nodes = get_sink_nodes(graph)
    tips = []
    for sink in sink_nodes:
        path = [sink]
        predecessors = list(graph.predecessors(sink))
        while len(predecessors) == 1 and not list(graph.successors(predecessors[0])):
            path.append(predecessors[0])
            predecessors = list(graph.predecessors(predecessors[0]))
        if len(path) > 1:
            tips.append(path[::-1])
    if tips:
        path_lengths = [len(p) for p in tips]
        avg_path_weights = [path_average_weight(graph, p) for p in tips]
        graph = select_best_path(graph, tips, path_lengths, avg_path_weights, delete_sink_node=True)
    print(f"Step 6c: Out tip removal completed. Graph now has {len(graph.edges)} edges.")
    return graph

# =============================================================================
# 5. Quality Metrics and Analysis
# =============================================================================

def GC(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def run_fastqc(fasta_file):
    """Run FASTQC on the input file and parse the report (skipped for simplicity)."""
    print("Step 7: FASTQC analysis skipped for small dataset.")
    return "FASTQC skipped: No analysis performed."

def calculate_quality_metrics(scaffolds, reference=None):
    """Calculate N50, G20, alignment score, and accuracy."""
    if not scaffolds:
        return 0, 0, 0, 0, 0

    lengths = sorted([length for _, length in scaffolds], reverse=True)
    total_length = sum(lengths)
    # Correct N50 calculation: length where 50% of total length is covered
    cumulative_length = 0
    n50 = 0
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            n50 = length
            break

    g20 = 100.0 if any(l >= 20 for l in lengths) else 0.00  # 100% if any contig >= 20 bp
    gc_content = [GC(seq) for seq, _ in scaffolds]
    avg_gc = np.mean(gc_content) if gc_content else 0

    alignment_score = 0
    accuracy = 0
    if reference and os.path.exists(REF_FILE):
        try:
            with open("temp_scaffolds.fasta", "w") as f:
                for i, (seq, _) in enumerate(scaffolds):
                    f.write(f">scaffold_{i}\n{seq}\n")
            cmd = ["minimap2", "-x", "map-ont", REF_FILE, "temp_scaffolds.fasta", "-o", "alignments.paf"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                aligned_bases = 0
                correct_bases = 0
                with open("alignments.paf", "r") as paf:
                    for line in paf:
                        fields = line.split()
                        aligned_bases += int(fields[10])  
                        correct_bases += int(fields[9])   
                alignment_score = (aligned_bases / total_length) * 100 if total_length else 0
                accuracy = (correct_bases / aligned_bases) * 100 if aligned_bases else 0
            else:
                print(f"minimap2 failed: {result.stderr}")
                alignment_score = 0
                accuracy = 0
        except Exception as e:
            print(f"Error in alignment calculation: {e}")
            alignment_score = 0
            accuracy = 0
        finally:
            for temp_file in ["temp_scaffolds.fasta", "alignments.paf"]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
    else:
        alignment_score = 0
        accuracy = 0

    print(f"Step 8: Quality metrics calculation completed (N50: {n50}, G20: {g20:.2f}%).")
    return n50, g20, avg_gc, alignment_score, accuracy

def run_busco(scaffolds):
    """Run BUSCO analysis if chosen and installed."""
    if RUN_BUSCO.lower() != "yes" or not scaffolds:
        print("Step 9: BUSCO analysis skipped (not requested or no scaffolds).")
        return "BUSCO skipped: Not performed."

    try:
        with open("temp_scaffolds.fasta", "w") as f:
            for i, (seq, _) in enumerate(scaffolds):
                f.write(f">scaffold_{i}\n{seq}\n")
        
        cmd = ["busco", "-i", "temp_scaffolds.fasta", "-o", "busco_output", "-l", "bacteria_odb10", "-m", "genome", "--offline"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        os.remove("temp_scaffolds.fasta")
        print("Step 9: BUSCO analysis completed.")
        return result.stdout if result.returncode == 0 else f"BUSCO failed: Not performed."
    
    except Exception as e:
        print(f"Error running BUSCO: {e}")
        return f"BUSCO failed: Not performed."

# =============================================================================
# 6. Plots
# =============================================================================

def plot_debruijn_graph(graph):
    """Plot the De Bruijn graph (simplified)."""
    if not graph or not graph.nodes():
        print("Step 10a: De Bruijn graph plot not performed.")
        return False
    pos = nx.spring_layout(graph)
    plt.figure(figsize=(12, 8))
    nx.draw(graph, pos, node_size=50, node_color='lightblue', with_labels=False)
    plt.title("De Bruijn Graph")
    plt.savefig(os.path.join(PLOTS_DIR, "debruijn_graph.png"))
    plt.close()
    print(f"Step 10a: De Bruijn graph plot saved to {PLOTS_DIR}/debruijn_graph.png.")
    return True

def plot_circos(scaffolds):
    """Create and save a Circos plot."""
    if not scaffolds:
        print("Step 10b: Circos plot not performed.")
        return False
    print(f"Debug: Scaffolds content = {scaffolds}")  # Debug output
    plt.figure(figsize=(10, 10))
    for i, (_, length) in enumerate(scaffolds[:min(10, len(scaffolds))]):
        circle = Circle((0, 0), i+1, fill=False)
        plt.gca().add_artist(circle)
        plt.text(0, i+1, f'Scaffold {i+1} (len={length})', va='center')
    plt.axis('equal')
    plt.title("Circos Plot of Scaffolds")
    plt.savefig(os.path.join(PLOTS_DIR, "circos_plot.png"))
    plt.close()
    print(f"Step 10b: Circos plot saved to {PLOTS_DIR}/circos_plot.png.")
    return True

def plot_gc_content(gc_contents):
    """Create and save a GC content histogram."""
    if not gc_contents:
        print("Step 10c: GC content plot not performed.")
        return False
    plt.figure(figsize=(10, 6))
    plt.hist(gc_contents, bins=20, color='blue', alpha=0.7)
    plt.title("GC Content Distribution")
    plt.xlabel("GC Percentage")
    plt.ylabel("Frequency")
    plt.savefig(os.path.join(PLOTS_DIR, "gc_content_plot.png"))
    plt.close()
    print(f"Step 10c: GC content plot saved to {PLOTS_DIR}/gc_content_plot.png.")
    return True

def plot_g20(g20):
    """Plot G20 score."""
    if g20 == 0:
        print("Step 10d: G20 plot not performed.")
        return False
    plt.figure(figsize=(6, 6))
    plt.bar(["G20"], [g20], color='green')
    plt.title("G20 Score")
    plt.ylabel("Percentage")
    plt.savefig(os.path.join(PLOTS_DIR, "g20_plot.png"))
    plt.close()
    print(f"Step 10d: G20 plot saved to {PLOTS_DIR}/g20_plot.png.")
    return True

# =============================================================================
# 7. Output Functions
# =============================================================================

def save_contigs(contig_tuples, output_filename):
    """Save contigs/scaffolds to a FASTA file."""
    if not contig_tuples:
        print("Step 11: Saving contigs/scaffolds not performed.")
        return
    with open(output_filename, "w") as f:
        for i, (contig, length) in enumerate(contig_tuples):
            f.write(f">contig_{i}_len={length}\n")
            f.write("\n".join(contig[i:i+60] for i in range(0, len(contig), 60)) + "\n")
    print(f"Step 11: Saved {len(contig_tuples)} contigs/scaffolds to {output_filename} completed.")

def generate_pdf_report(n50, g20, avg_gc, alignment_score, accuracy, contig_tuples, busco_result, fastqc_result):
    """Generate a PDF report with all metrics."""
    if not contig_tuples:
        print("Step 12: Report generation not performed due to no data.")
        return
    doc = SimpleDocTemplate(REPORT_FILE, pagesize=letter)
    styles = getSampleStyleSheet()
    elements = []

    elements.append(Paragraph("Assembly Quality Report", styles['Heading1']))
    elements.append(Spacer(1, 12))

    # Mark as "Not Performed" if values are zero or default
    alignment_score_str = f"{alignment_score:.2f}%" if alignment_score > 0 else "Not Performed"
    accuracy_str = f"{accuracy:.2f}%" if accuracy > 0 else "Not Performed"
    busco_result_str = busco_result if "completed" in busco_result.lower() else "Not Performed"
    fastqc_result_str = fastqc_result if "performed" not in fastqc_result.lower() else "Not Performed"
    g20_str = f"{g20:.2f}%" if g20 > 0 else "Not Performed"
    gc_content_str = f"{avg_gc:.2f}%" if avg_gc > 0 else "Not Performed"

    data = [
        ["Metric", "Value"],
        ["Number of Contigs/Scaffolds", str(len(contig_tuples)) if contig_tuples else "Not Performed"],
        ["Total Length", f"{sum(length for _, length in contig_tuples)} bp" if contig_tuples else "Not Performed"],
        ["Longest Contig", f"{max((length for _, length in contig_tuples), default=0)} bp" if contig_tuples else "Not Performed"],
        ["N50", f"{n50} bp" if n50 > 0 else "Not Performed"],
        ["G20", g20_str],
        ["GC Content", gc_content_str],
        ["Alignment Score", alignment_score_str],
        ["Accuracy", accuracy_str],
        ["BUSCO Result", busco_result_str],
        ["FASTQC Summary", fastqc_result_str]
    ]

    table = Table(data)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 14),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 12),
    ]))
    elements.append(table)

    doc.build(elements)
    print(f"Step 12: PDF report generation {'completed' if contig_tuples else 'not performed'} and saved to {REPORT_FILE}.")

# =============================================================================
# Main Function
# =============================================================================

def main():
    print(f"Starting assembly pipeline for {INPUT_FILE}...")
    print(f"Using reference file: {REF_FILE}")
    
    # Step 1: Read and classify reads
    short_reads, long_reads = read_and_classify_fasta(INPUT_FILE)

    if not short_reads and not long_reads:
        print("Step 1: Reading and classifying not performed due to no sequences found.")
        return

    # Initialize graph and contigs/scaffolds
    graph = None
    contigs = []

    # Step 2-4: De Bruijn assembly for short reads (if chosen)
    if ASSEMBLY_TYPE in ["short", "both"]:
        kmer_count = build_kmer_dict(short_reads, K_SIZE_SHORT)
        graph = build_graph(kmer_count)
        if graph and graph.nodes():
            graph = simplify_bubbles(graph)
            graph = solve_entry_tips(graph)
            graph = solve_out_tips(graph)
            starting_nodes = get_starting_nodes(graph)
            sink_nodes = get_sink_nodes(graph)
            contigs = get_contigs(graph, starting_nodes, sink_nodes)
            print(f"Step 4: Short read assembly completed with {len(contigs)} contigs.")
        else:
            print("Step 4: Short read assembly not performed due to no valid graph.")

    # Step 5: Optional long read scaffolding
    scaffolds = contigs  # Default to short read contigs
    if ASSEMBLY_TYPE in ["long", "both"] and long_reads:
        scaffolds = scaffold_contigs(contigs, long_reads)
        if not scaffolds:
            print("Step 5: Scaffolding not performed due to no valid scaffolds generated.")
    else:
        print("Step 5: Scaffolding not performed due to no long reads or incorrect assembly type.")

    # Step 6-9: Quality analysis
    reference = open(REF_FILE, "r").read().replace("\n", "").upper() if os.path.exists(REF_FILE) else None
    fastqc_result = run_fastqc(INPUT_FILE)
    n50, g20, avg_gc, alignment_score, accuracy = calculate_quality_metrics(scaffolds, reference)
    busco_result = run_busco(scaffolds) if RUN_BUSCO.lower() == "yes" else "BUSCO skipped: Not performed."

    # Step 10: Generate plots
    if scaffolds:
        debruijn_plotted = plot_debruijn_graph(graph)
        circos_plotted = plot_circos(scaffolds)
        gc_plotted = plot_gc_content([GC(seq) for seq, _ in scaffolds] if scaffolds else [])
        g20_plotted = plot_g20(g20)
    else:
        print("Step 10: Plot generation not performed due to no scaffolds.")

    # Step 11: Save contigs/scaffolds
    save_contigs(scaffolds, OUTPUT_FILE)

    # Step 12: Generate report
    generate_pdf_report(n50, g20, avg_gc, alignment_score, accuracy, scaffolds, busco_result, fastqc_result)

    # Clean up
    if os.path.exists(FASTQC_DIR):
        shutil.rmtree(FASTQC_DIR)
    for temp_file in ["temp_scaffolds.fasta", "alignments.paf"]:
        if os.path.exists(temp_file):
            os.remove(temp_file)

if __name__ == "__main__":
    main()