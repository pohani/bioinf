# HIRGC - Genome Compression Tool

A tool for compressing and decompressing genome sequences.

## Prerequisites

- C++17 compatible compiler (g++)
- For Linux: Basic build tools

## Installation

### Windows

1. Navigate to the project directory
2. Compile the source code:
   ```bash
   g++ -std=c++17 hirgc.cpp -O3 -o hirgc
   ```

### Linux

1. Copy the following files to your Linux system:
   - `hirgc.cpp` (source code)
   - All files from the `./data/` directory (`ref.fa`, `target.fa`)

2. Install required dependencies:
   ```bash
   sudo apt-get update
   sudo apt-get install g++
   ```

3. Compile the source code:
   ```bash
   g++ -std=c++17 hirgc.cpp -O3 -o hirgc -Wall -Wextra
   ```

## Usage

### Compression

To compress a genome sequence:
```bash
./hirgc compress <reference_file> <target_file> <output_file>
```

Example:
```bash
./hirgc compress ./data/ref.fa ./data/target.fa ./data/compressed.bin
```

### Decompression

To decompress a compressed genome:
```bash
./hirgc decompress <reference_file> <compressed_file> <output_file>
```

Example:
```bash
./hirgc decompress ./data/ref.fa ./data/compressed.bin ./data/decompressed.fa
```

## File Descriptions

- `ref.fa`: Reference genome file
- `target.fa`: Target genome file to be compressed
- `compressed.bin`: Output compressed file
- `decompressed.fa`: Output decompressed file

