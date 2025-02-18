# IBS-DAA-S4-D12

**Title:**  
"Genome Sequencing Using De Bruijn Graph"

**Implementation Details:**  
This project focuses on the application of De Bruijn graphs in genome sequencing. De Bruijn graphs are widely used in bioinformatics to efficiently assemble short DNA sequences (reads) into long contiguous sequences (contigs). The project aims to implement a method for genome sequencing that uses a De Bruijn graph-based algorithm. The steps of the implementation include:

  1. **Input Data:** Raw sequencing reads from high-throughput sequencing technologies such as Illumina.
  2. **K-mer Construction:** First, the reads are broken down into k-mers (substrings of length k). These k-mers will serve as the vertices of the De Bruijn graph.
  3. **Graph Construction:** The De Bruijn graph is constructed by linking k-mers that overlap by (k-1) bases. Each edge of the graph represents a possible overlap between two k-mers.
  4. **Graph Traversal:** Algorithms like Eulerian path or Hamiltonian path can be used to find the optimal traversal of the graph and reconstruct the genome sequence.
  5. **Error Correction and Contig Assembly:** The graph will be analyzed to identify errors, and the genome will be assembled into longer contigs.
  6. **Post-Processing:** The contigs are further processed to remove redundancy and errors to produce a final genome sequence.

**Abstract:**  
Genome sequencing is the process of determining the complete DNA sequence of an organism's genome. In recent years, high-throughput sequencing technologies have produced vast amounts of short DNA fragments, requiring efficient computational techniques for sequence assembly. This project explores the use of De Bruijn graphs, a data structure that has proven to be highly effective in genome assembly. The project implements a genome sequencing pipeline that leverages De Bruijn graph construction and traversal algorithms. By breaking down DNA sequences into k-mers, a De Bruijn graph is built to represent the overlaps between the k-mers. The resulting graph is then analyzed to assemble the genome through a combination of graph traversal and error-correction techniques. The project aims to provide an efficient method for genome assembly, applicable to both small-scale and large-scale sequencing efforts in the fields of bioinformatics and computational biology. 

**Subject:**  
- Data Analysis and Algorithms (DAA)  
- Intelligence of Biological Systems (IBS-2)
