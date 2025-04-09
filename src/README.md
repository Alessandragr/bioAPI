#Overview
The bioAPI is a bioinformatics tool designed to process, analyze, and manipulate genomic data. It supports operations such as sequence extraction, GC content calculation, k-mer indexing, FASTA/FASTQ file handling, and more. The program is implemented in C++ and provides a command-line interface for interacting with genomic data.

#Features

##FASTA/FASTQ File Handling:

Load and process genomic sequences from FASTA and FASTQ files.
Count sequences in MULTIFASTA files.
Extract sequences between start and stop codons.
Filter FASTQ files based on quality.

##Sequence Analysis:

Calculate GC content.
Identify degenerate bases.
Calculate average sequence quality.

##Sequence Manipulation:

Generate complementary and reverse-complementary sequences.
Remove poly-A/poly-T tails.
Remove sequence prefixes.

##Indexing and Pattern Matching:

Build suffix and LCP (Longest Common Prefix) tables.
Search for patterns in sequences.
Find repeated patterns.
Map reads to the genome using k-mers.

##Command-Line Interface:

Supports multiple options for processing and analyzing genomic data.
Provides detailed usage instructions.

##Dependencies
The program uses the following libraries and packages:

##Standard Libraries
<iostream>                          For input/output operations.
<fstream>:                          For file handling.
<string>:                           For string manipulation.
<sstream>:                          For string stream operations.
<algorithm>:                        For algorithms like sorting and transformations.
<stdexcept>:                        For exception handling.
<optional>:                         For optional values.
<vector>:                           For dynamic arrays.
<tuple>:                            For tuple operations.
<utility>:                          For utility functions like std::pair.
<regex>:                            For regular expressions.
<cmath>:                            For mathematical operations.
<numeric>:                          For numeric operations like std::iota.

##Third-Party Libraries
<zlib.h>:                          For handling GZIP file decompression.
<archive.h> and <archive_entry.h>  For handling TAR file decompression.
<minizip/unzip.h>                  For handling ZIP file decompression.


##Platform-Specific Libraries
<windows.h>: For Windows-specific console encoding (UTF-8).

##Installation
1. Clone the Repository:

git clone https://github.com/your-repo/bioAPI.git
cd bioAPI/src

2. Install Dependencies: Ensure the required libraries (e.g., zlib, libarchive, minizip) are installed on your system.


3. Build the Program: Use a C++ compiler (e.g., g++) to compile the program:

g++ -o bioAPI MainView.cpp Controller/GeneticMaterialController.cpp Model/GeneticMaterial.cpp -lz -larchive -lminizip

4. Run the Program:
./bioAPI --file <file_path> [options]

##Usage
The program provides a variety of command-line options for genomic data processing. Below is a list of available options:

###Command-Line Options
Option	Description
--file <file_path>	             Specify the input file (FASTA/FASTQ).
--count-fasta	                   Count the number of sequences in a MULTIFASTA file.
--complement	                   Calculate the complementary sequence.
--reverse-complement	           Calculate the reverse-complementary sequence.
--start-stop	                   Extract sequences between start and stop codons.
--extract	                       Extract subsequences based on start, length, and sequence ID.
--description	                   Get descriptions of sequences (e.g., ID, species, length).
--gc-content	                   Calculate the GC content of sequences.
--average-quality	               Calculate the average quality of a sequence.
--remove-prefix	                 Remove a prefix from sequences.
--trim-polyAT	                   Remove poly-A/poly-T tails from sequences.
--degenerate-bases	             Check for degenerate bases in sequences.
--filter-fastq	                 Filter FASTQ files based on quality.
--process-clean	                 Process and clean sequences (e.g., remove poly-A/T, filter by quality).
--suffix-table	                 Build and display the suffix and LCP tables.
--search-pattern	               Search for a specific pattern in the sequence.
--find-repeated	                 Find repeated patterns in the sequence.
--map-reads	                     Map reads to the genome using k-mers.
-h, --help	                     Display the help message.

##Examples

1. Count FASTA Sequences
./bioAPI --file genome.fasta --count-fasta

2. Calculate GC Content
./bioAPI --file genome.fasta --gc-content

3. Extract Sequences Between Start and Stop Codons
./bioAPI --file genome.fasta --start-stop ATG,TAA,TAG,TGA

4. Filter FASTQ by Quality
./bioAPI --file reads.fastq --filter-fastq 30

5. Build and Display Suffix Table
./bioAPI --file genome.fasta --suffix-table

##Error Handling
The program provides detailed error messages for common issues, such as:

Missing or invalid input files.
Unsupported file formats.
Invalid command-line arguments.

Example:
Error: Could not load the file content: genome.fasta

##Contributing
Contributions are welcome! Please follow these steps:

##Fork the repository.
Create a new branch for your feature or bug fix.
Submit a pull request with a detailed description of your changes.
License
This project is licensed under the MIT License. See the LICENSE file for details.

