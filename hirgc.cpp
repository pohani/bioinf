/**
 * High-Information Reference-based Genome Compressor (HIRGC)
 * A specialized compression algorithm for genomic sequences that leverages reference-based compression
 * to achieve high compression ratios while maintaining perfect reconstruction.
 *
 * Written by Karlo Pohajda
 */

#include <iostream>
#include <chrono>
#include <string>
#include <cctype>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>

// Global storage for reference genome sequence
std::vector<char> referenceSequence;
// Mapping of nucleotides to their binary codes (A=0, C=1, G=2, T=3)
char nucleotideCodeMap[] = {0, 1, 2, 3};

// Metadata and sequence storage for target genome
std::string sequenceMetadata;
std::vector<char> targetSequence;
std::vector<char> decompressedSequence;

// Configuration parameters
int basesPerLine = 80;

struct AmbiguousInterval
{
    int start;  // 0‐based position in the original FASTA (counting only valid bases)
    int length; // how many consecutive ambiguous bases
    char base;  // exactly which IUPAC letter (‘N’, ‘R’, ‘Y’, …) to restore
};
static std::vector<AmbiguousInterval> ambiguousBaseIntervals;

// Hash table parameters for k-mer indexing
static constexpr int MAX_HASH_SIZE = 1 << 20;
static std::vector<int> hashTable;
static std::vector<int> nextPosition;

// Length of k-mers used for matching
static constexpr int K_MER_LENGTH = 11;

// Mapping of binary codes back to nucleotide characters
static const char nucleotideMap[4] = {'A', 'C', 'G', 'T'};

/**
 * Converts a nucleotide character to its corresponding index in the nucleotide map
 * @param nucleotide The nucleotide character to convert
 * @return Index in the nucleotide map (0-3) or -1 if invalid
 *
 * Written by Karlo Pohajda
 */
int getNucleotideIndex(char nucleotide)
{
    switch (nucleotide)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return -1;
    }
}

/**
 * Loads and processes the reference genome from a FASTA file
 * Converts nucleotides to binary codes and stores them in referenceSequence
 * @param referenceFile Path to the reference genome FASTA file
 *
 * Written by Karlo Pohajda
 */
void loadReferenceGenome(const std::string &referenceFile)
{
    std::cout << "Reading reference file: " << referenceFile << std::endl;

    std::ifstream in(referenceFile);
    if (!in.is_open())
    {
        std::cerr << "Error: could not open file " << referenceFile << std::endl;
        return;
    }

    referenceSequence.clear();
    std::string line;
    if (!std::getline(in, line))
    {
        std::cerr << "Error: file " << referenceFile << " is empty." << std::endl;
        return;
    }

    while (std::getline(in, line))
    {
        for (char nucleotide : line)
        {
            char upperNucleotide = static_cast<char>(std::toupper(static_cast<unsigned char>(nucleotide)));
            int index = getNucleotideIndex(upperNucleotide);
            if (index >= 0)
            {
                referenceSequence.push_back(nucleotideCodeMap[index]);
            }
        }
    }

    int referenceLength = static_cast<int>(referenceSequence.size());
    std::cout << "referenceLength: " << referenceLength << "\n";

    int previewLength = std::min(referenceLength, 50);
    std::cout << "First " << previewLength << " values of referenceSequence: ";
    for (int i = 0; i < previewLength; ++i)
    {
        std::cout << +referenceSequence[i];
    }
    std::cout << std::endl;
}

/**
 * Loads and processes the target genome from a FASTA file
 * Stores metadata and converts nucleotides to binary codes
 * @param targetFile Path to the target genome FASTA file
 *
 * Written by Karlo Pohajda
 */
void loadTargetGenome(const std::string &targetFile)
{
    std::cout << "Reading target file: " << targetFile << std::endl;
    std::ifstream in(targetFile);
    if (!in.is_open())
    {
        std::cerr << "Error: could not open file " << targetFile << std::endl;
        return;
    }

    std::string line;
    std::getline(in, sequenceMetadata);
    targetSequence.clear();

    while (std::getline(in, line))
    {
        int lineLength = (int)line.size();
        for (int i = 0; i < lineLength; ++i)
        {
            char nucleotide = line[i];
            if (std::islower(nucleotide))
            {
                nucleotide = (char)std::toupper(nucleotide);
            }

            int index = getNucleotideIndex(nucleotide);
            if (index >= 0)
            {
                targetSequence.push_back(nucleotideCodeMap[index]);
            }
        }
    }

    int targetLength = static_cast<int>(targetSequence.size());
    std::cout << "targetLength: " << targetLength << "\n";

    int previewLength = std::min(targetLength, 50);
    std::cout << "First " << previewLength << " values of targetSequence: ";
    for (int i = 0; i < previewLength; ++i)
    {
        std::cout << +targetSequence[i];
    }
    std::cout << std::endl;
}

/**
 * Scans the target genome for ambiguous bases (N) and records their positions
 * Used to preserve ambiguous regions during compression/decompression
 * @param targetFile Path to the target genome FASTA file
 *
 * Written by Karlo Pohajda
 */
void scanAmbiguousBases(const std::string &targetFile)
{
    ambiguousBaseIntervals.clear();

    std::ifstream in(targetFile);
    if (!in.is_open())
    {
        std::cerr << "Error: could not open target file for N‐scanning: "
                  << targetFile << std::endl;
        return;
    }

    std::string line;
    if (!std::getline(in, line))
        return;

    int position = 0;
    bool inAmbiguousRun = false;
    int runStart = 0;
    int runLength = 0;
    char runChar = '\0'; // which IUPAC letter (‘N’, ‘R’, ‘Y’, …) is in this run

    while (std::getline(in, line))
    {
        for (char nucleotide : line)
        {
            char upperNucleotide = static_cast<char>(std::toupper(static_cast<unsigned char>(nucleotide)));

            // Only A, C, G, T, plus any of the IUPAC ambiguity codes:
            if (upperNucleotide == 'A' || upperNucleotide == 'C' ||
                upperNucleotide == 'G' || upperNucleotide == 'T' ||
                upperNucleotide == 'N' || upperNucleotide == 'R' ||
                upperNucleotide == 'Y' || upperNucleotide == 'S' ||
                upperNucleotide == 'W' || upperNucleotide == 'K' ||
                upperNucleotide == 'M' || upperNucleotide == 'B' ||
                upperNucleotide == 'D' || upperNucleotide == 'H' ||
                upperNucleotide == 'V')
            {
                // Is this letter one of the “ambiguous” codes?
                bool isAmbiguous = (upperNucleotide == 'N' ||
                                    upperNucleotide == 'R' ||
                                    upperNucleotide == 'Y' ||
                                    upperNucleotide == 'S' ||
                                    upperNucleotide == 'W' ||
                                    upperNucleotide == 'K' ||
                                    upperNucleotide == 'M' ||
                                    upperNucleotide == 'B' ||
                                    upperNucleotide == 'D' ||
                                    upperNucleotide == 'H' ||
                                    upperNucleotide == 'V');

                if (isAmbiguous)
                {
                    if (!inAmbiguousRun)
                    {
                        // start a new run of ambiguous bases
                        inAmbiguousRun = true;
                        runStart = position;
                        runLength = 1;
                        runChar = upperNucleotide;
                    }
                    else
                    {
                        // we are already in a run—check if the letter is the same
                        if (upperNucleotide == runChar)
                        {
                            // continuing the same ambiguous‐letter run
                            ++runLength;
                        }
                        else
                        {
                            // the ambiguous letter changed, so close the old run…
                            ambiguousBaseIntervals.push_back(
                                AmbiguousInterval{runStart, runLength, runChar});
                            // …and start a new run with this new letter
                            runStart = position;
                            runLength = 1;
                            runChar = upperNucleotide;
                        }
                    }
                }
                else
                {
                    // it’s A/C/G/T. If we were in an ambiguous run, close it now.
                    if (inAmbiguousRun)
                    {
                        ambiguousBaseIntervals.push_back(
                            AmbiguousInterval{runStart, runLength, runChar});
                        inAmbiguousRun = false;
                    }
                    // (nothing else to do for A/C/G/T except advance position)
                }

                ++position; // count this base toward the “genome coordinate”
            }
            // else—some other character (e.g. FASTA header, whitespace)—skip entirely
        }
    }

    // If file ended while still in a run of ambiguous bases, close it out:
    if (inAmbiguousRun)
    {
        ambiguousBaseIntervals.push_back(
            AmbiguousInterval{runStart, runLength, runChar});
    }

    std::cout << "Found " << ambiguousBaseIntervals.size()
              << " ambiguous‐base intervals" << std::endl;
    if (!ambiguousBaseIntervals.empty())
    {
        const auto &first = ambiguousBaseIntervals[0];
        std::cout << "First interval: start=" << first.start
                  << ", length=" << first.length
                  << ", base=" << first.base << std::endl;
    }
}

/**
 * Builds a hash table index of k-mers from the reference sequence
 * Enables efficient matching of k-mers between target and reference
 *
 * Written by Karlo Pohajda
 */
void buildReferenceIndex()
{
    int referenceLength = static_cast<int>(referenceSequence.size());
    int windowCount = referenceLength - K_MER_LENGTH + 1;

    if (windowCount <= 0)
    {
        std::cerr << "Error: Reference sequence too short for K_MER_LENGTH="
                  << K_MER_LENGTH << std::endl;
        return;
    }

    int shiftBits = K_MER_LENGTH * 2 - 2;

    hashTable.assign(MAX_HASH_SIZE, -1);
    nextPosition.assign(windowCount, -1);

    uint64_t hashValue = 0;
    for (int k = K_MER_LENGTH - 1; k >= 0; --k)
    {
        hashValue = (hashValue << 2) | static_cast<uint8_t>(referenceSequence[k]);
    }
    int hashIndex = static_cast<int>(hashValue & (MAX_HASH_SIZE - 1));
    nextPosition[0] = hashTable[hashIndex];
    hashTable[hashIndex] = 0;

    for (int i = 1; i < windowCount; ++i)
    {
        hashValue >>= 2;
        hashValue |= (uint64_t(referenceSequence[i + K_MER_LENGTH - 1]) << shiftBits);

        hashIndex = static_cast<int>(hashValue & (MAX_HASH_SIZE - 1));
        nextPosition[i] = hashTable[hashIndex];
        hashTable[hashIndex] = i;
    }
}

std::ofstream debugLog;
bool enableDebugLogging = true;

/**
 * Initializes debug logging functionality
 * @param logFile Path to the debug log file
 *
 * Written by Karlo Pohajda
 */
void initializeDebugLogging(const std::string &logFile = "compression_debug.log")
{
    if (enableDebugLogging)
    {
        debugLog.open(logFile);
        if (!debugLog)
        {
            std::cerr << "Warning: Could not open debug log file. Continuing without logging." << std::endl;
            enableDebugLogging = false;
        }
        else
        {
            debugLog << "=== Compression Debug Log ===\n";
        }
    }
}

/**
 * Core compression function that finds matches between target and reference sequences
 * Uses k-mer matching and extension to identify long matches
 * Writes compressed data in a binary format
 * @param outputFile Path to write the compressed output
 *
 * Written by Karlo Pohajda
 */
void findMatches(const std::string &outputFile)
{
    initializeDebugLogging();

    if (enableDebugLogging)
    {
        debugLog << "Starting match search process with K_MER_LENGTH=" << K_MER_LENGTH << "\n";
        debugLog << "Target sequence length: " << targetSequence.size() << " bases\n";
        debugLog << "Reference sequence length: " << referenceSequence.size() << " bases\n\n";
    }

    buildReferenceIndex();

    std::ofstream out(outputFile,
                      std::ios::out | std::ios::trunc | std::ios::binary);
    if (!out)
    {
        std::cerr << "Error: cannot open " << outputFile << " for writing\n";
        return;
    }

    // Store sequence metadata
    uint32_t metadataLength = static_cast<uint32_t>(sequenceMetadata.length());
    out.write(reinterpret_cast<char *>(&metadataLength), 4);
    out.write(sequenceMetadata.c_str(), metadataLength);

    uint32_t ambiguousCount = static_cast<uint32_t>(ambiguousBaseIntervals.size());
    out.write(reinterpret_cast<const char *>(&ambiguousCount), 4);

    for (const auto &interval : ambiguousBaseIntervals)
    {
        uint32_t start = static_cast<uint32_t>(interval.start);
        uint32_t length = static_cast<uint32_t>(interval.length);
        char base = interval.base; // ASCII letter (‘N’, ‘R’, etc.)

        out.write(reinterpret_cast<const char *>(&start), 4);
        out.write(reinterpret_cast<const char *>(&length), 4);
        out.write(&base, 1);
    }

    int targetLength = static_cast<int>(targetSequence.size());
    int lastLiteralPos = 0;
    int lastRefPos = 0;
    int windowCount = targetLength - K_MER_LENGTH + 1;

    if (windowCount <= 0)
    {
        std::cerr << "Warning: Target sequence too short for K_MER_LENGTH="
                  << K_MER_LENGTH << ", emitting as literal" << std::endl;

        if (enableDebugLogging)
        {
            debugLog << "Target sequence too short. Emitting entire sequence as literal.\n";
            debugLog.close();
        }

        uint8_t tag = 0x00;
        uint32_t length = static_cast<uint32_t>(targetLength);
        out.write(reinterpret_cast<char *>(&tag), 1);
        out.write(reinterpret_cast<char *>(&length), 4);
        out.write(reinterpret_cast<char *>(&targetSequence[0]), targetLength);
        return;
    }

    int matchCount = 0;
    int literalCount = 0;
    int totalMatchLength = 0;
    int totalLiteralLength = 0;

    for (int i = 0; i < windowCount; ++i)
    {
        uint64_t hashValue = 0;
        for (int k = K_MER_LENGTH - 1; k >= 0; --k)
        {
            hashValue = (hashValue << 2) | static_cast<uint8_t>(targetSequence[i + k]);
        }
        int hashIndex = hashTable[hashValue & (MAX_HASH_SIZE - 1)];
        if (hashIndex < 0)
        {
            if (enableDebugLogging && i % 10000 == 0)
            {
                debugLog << "Position " << i << ": No hash match found\n";
            }
            continue;
        }

        int bestRefPos = -1;
        size_t bestMatchLength = 0;
        int chainPositionsChecked = 0;

        for (int k = hashIndex; k != -1; k = nextPosition[k])
        {
            chainPositionsChecked++;

            bool initialMatchVerified = true;
            for (int m = 0; m < K_MER_LENGTH; m++)
            {
                if (referenceSequence[k + m] != targetSequence[i + m])
                {
                    initialMatchVerified = false;
                    if (enableDebugLogging)
                    {
                        debugLog << "WARNING: Hash collision detected at position " << i
                                 << " (reference pos " << k << "). Initial k-mer doesn't actually match!\n";
                    }
                    break;
                }
            }

            if (!initialMatchVerified)
            {
                continue;
            }

            size_t refPos = k + K_MER_LENGTH, targetPos = i + K_MER_LENGTH, matchLength = K_MER_LENGTH;

            while (refPos < referenceSequence.size() && targetPos < static_cast<size_t>(targetLength) && referenceSequence[refPos] == targetSequence[targetPos])
            {
                if (enableDebugLogging && matchLength - K_MER_LENGTH < 10)
                {
                    debugLog << "  Extending match: pos=" << (matchLength - K_MER_LENGTH)
                             << " tar[" << targetPos << "]=" << +targetSequence[targetPos] << "(" << nucleotideMap[static_cast<unsigned char>(targetSequence[targetPos])] << ")"
                             << " ref[" << refPos << "]=" << +referenceSequence[refPos] << "(" << nucleotideMap[static_cast<unsigned char>(referenceSequence[refPos])] << ")\n";
                }
                ++refPos;
                ++targetPos;
                ++matchLength;
            }

            if (matchLength > bestMatchLength)
            {
                bestMatchLength = matchLength;
                bestRefPos = k;
            }
        }

        if (bestRefPos < 0)
        {
            if (enableDebugLogging && i % 10000 == 0)
            {
                debugLog << "Position " << i << ": Hash match but no valid extension\n";
            }
            continue;
        }

        if (enableDebugLogging)
        {
            std::string matchContext = "";
            int contextSize = 5;
            size_t matchStart = static_cast<size_t>(i);
            size_t refStart = static_cast<size_t>(bestRefPos);

            for (size_t j = std::max(size_t(0), matchStart - contextSize);
                 j < std::min(static_cast<size_t>(targetLength), matchStart + bestMatchLength + contextSize); j++)
            {
                char base = nucleotideMap[static_cast<unsigned char>(targetSequence[j])];
                if (j >= matchStart && j < matchStart + bestMatchLength)
                {
                    matchContext += base;
                }
                else
                {
                    matchContext += std::tolower(base);
                }
            }

            std::string refContext = "";
            for (size_t j = std::max(size_t(0), refStart - contextSize);
                 j < std::min(referenceSequence.size(), static_cast<size_t>(refStart) + bestMatchLength + contextSize); j++)
            {
                char base = nucleotideMap[static_cast<unsigned char>(referenceSequence[j])];
                if (j >= refStart && j < refStart + bestMatchLength)
                {
                    refContext += base;
                }
                else
                {
                    refContext += std::tolower(base);
                }
            }

            debugLog << "Found match at target pos " << i << ", ref pos " << bestRefPos << ":\n";
            debugLog << "  • Match length: " << bestMatchLength << " bases\n";
            debugLog << "  • Offset from previous ref position: " << bestRefPos - lastRefPos << "\n";
            debugLog << "  • Chain positions checked: " << chainPositionsChecked << "\n";
            debugLog << "  • Target context: " << matchContext << "\n";
            debugLog << "  • Reference context: " << refContext << "\n";

            bool matchVerified = true;
            std::string mismatchDetails = "";
            for (size_t v = 0; v < bestMatchLength; v++)
            {
                if (targetSequence[i + v] != referenceSequence[bestRefPos + v])
                {
                    matchVerified = false;
                    mismatchDetails += "    Mismatch at position " + std::to_string(v) + ": " +
                                       "target=" + nucleotideMap[static_cast<unsigned char>(targetSequence[i + v])] + ", " +
                                       "ref=" + nucleotideMap[static_cast<unsigned char>(referenceSequence[bestRefPos + v])] + "\n";
                }
            }

            if (!matchVerified)
            {
                debugLog << "  • WARNING! Match verification failed. Mismatches found:\n"
                         << mismatchDetails;
            }
            else
            {
                debugLog << "  • Match verified: All " << bestMatchLength << " characters match exactly\n";
            }

            if (i - lastLiteralPos > 0)
            {
                debugLog << "  • Emitting literal run of " << (i - lastLiteralPos) << " bases before this match\n";
            }
            debugLog << "\n";
        }

        int literalLength = i - lastLiteralPos;
        if (literalLength > 0)
        {
            uint8_t tag = 0x00;
            uint32_t length = static_cast<uint32_t>(literalLength);
            out.write(reinterpret_cast<char *>(&tag), 1);
            out.write(reinterpret_cast<char *>(&length), 4);
            out.write(reinterpret_cast<char *>(&targetSequence[lastLiteralPos]),
                      literalLength);

            literalCount++;
            totalLiteralLength += literalLength;
        }

        {
            uint8_t tag = 0x01;
            int32_t offset = static_cast<int32_t>(bestRefPos - lastRefPos);
            uint32_t matchLength = static_cast<uint32_t>(bestMatchLength - K_MER_LENGTH);
            out.write(reinterpret_cast<char *>(&tag), 1);
            out.write(reinterpret_cast<char *>(&offset), 4);
            out.write(reinterpret_cast<char *>(&matchLength), 4);

            matchCount++;
            totalMatchLength += bestMatchLength;
        }

        lastRefPos = bestRefPos + bestMatchLength;
        lastLiteralPos = i + bestMatchLength;
        i = lastLiteralPos - 1;
    }

    if (lastLiteralPos < targetLength)
    {
        uint8_t tag = 0x00;
        uint32_t length = static_cast<uint32_t>(targetLength - lastLiteralPos);
        out.write(reinterpret_cast<char *>(&tag), 1);
        out.write(reinterpret_cast<char *>(&length), 4);
        out.write(reinterpret_cast<char *>(&targetSequence[lastLiteralPos]),
                  length);

        literalCount++;
        totalLiteralLength += (targetLength - lastLiteralPos);

        if (enableDebugLogging)
        {
            debugLog << "Emitting final literal run of " << (targetLength - lastLiteralPos) << " bases\n";
        }
    }

    std::cout << "Compression stats: " << matchCount << " matches, "
              << literalCount << " literals\n";
    std::cout << "Match bytes: " << totalMatchLength
              << ", Literal bytes: " << totalLiteralLength << "\n";

    if (enableDebugLogging)
    {
        debugLog << "\n=== Compression Summary ===\n";
        debugLog << "Total matches: " << matchCount << "\n";
        debugLog << "Total literals: " << literalCount << "\n";
        debugLog << "Total match bytes: " << totalMatchLength << "\n";
        debugLog << "Total literal bytes: " << totalLiteralLength << "\n";
        double matchPercent = 100.0 * totalMatchLength / (totalMatchLength + totalLiteralLength);
        debugLog << "Match percentage: " << matchPercent << "%\n";
        debugLog << "Average match length: " << (matchCount > 0 ? totalMatchLength / matchCount : 0) << " bases\n";
        debugLog << "Average literal length: " << (literalCount > 0 ? totalLiteralLength / literalCount : 0) << " bases\n";
        debugLog.close();
        std::cout << "Debug log written to compression_debug.log" << std::endl;
    }
}

/**
 * Main compression function that orchestrates the entire compression process
 * @param referenceFile Path to the reference genome
 * @param targetFile Path to the target genome
 * @param outputFile Path to write the compressed output
 * @param debugMode Whether to enable detailed debug logging
 *
 * Written by Karlo Pohajda
 */
void compressGenome(const std::string &referenceFile, const std::string &targetFile, const std::string &outputFile, bool debugMode = true)
{
    std::cout << "Starting compression process..." << std::endl;

    enableDebugLogging = debugMode;

    loadReferenceGenome(referenceFile);
    loadTargetGenome(targetFile);
    scanAmbiguousBases(targetFile);
    findMatches(outputFile);

    std::ifstream origFile(targetFile, std::ios::binary | std::ios::ate);
    std::ifstream compFile(outputFile, std::ios::binary | std::ios::ate);

    if (origFile && compFile)
    {
        std::streamsize origSize = origFile.tellg();
        std::streamsize compSize = compFile.tellg();

        std::cout << "Original size: " << origSize << " bytes\n";
        std::cout << "Compressed size: " << compSize << " bytes\n";

        if (origSize > 0)
        {
            double ratio = static_cast<double>(compSize) / static_cast<double>(origSize);
            std::cout << "Compression ratio: " << ratio << " ("
                      << (ratio * 100) << "%)\n";
        }
    }

    std::cout << "Compression process finished." << std::endl;
}

/**
 * Validates the decompressed genome against the original target
 * Checks for perfect reconstruction and reports any mismatches
 * @param targetFile Path to the original target genome
 * @param decompressedFile Path to the decompressed genome
 *
 * Written by Karlo Pohajda
 */
void validateDecompressedGenome(const std::string &targetFile, const std::string &decompressedFile)
{
    std::string targetWithAmbiguous;
    std::ifstream targetFileStream(targetFile);

    if (!targetFileStream)
    {
        std::cerr << "Cannot open original target file for validation: " << targetFile << "\n";
        return;
    }

    std::string line;
    std::getline(targetFileStream, line);

    while (std::getline(targetFileStream, line))
    {
        for (char nucleotide : line)
        {
            if (std::isalpha(nucleotide))
            {
                targetWithAmbiguous += std::toupper(nucleotide);
            }
        }
    }
    targetFileStream.close();

    std::string decompressedWithAmbiguous;
    std::ifstream decompressedFileStream(decompressedFile);

    if (!decompressedFileStream)
    {
        std::cerr << "Cannot open decompressed file for validation: " << decompressedFile << "\n";
        return;
    }

    std::getline(decompressedFileStream, line);

    while (std::getline(decompressedFileStream, line))
    {
        for (char nucleotide : line)
        {
            if (std::isalpha(nucleotide))
            {
                decompressedWithAmbiguous += std::toupper(nucleotide);
            }
        }
    }
    decompressedFileStream.close();

    if (targetWithAmbiguous.size() != decompressedWithAmbiguous.size())
    {
        std::cerr << "Validation FAILED: Size mismatch!\n";
        std::cerr << "Target file size: " << targetWithAmbiguous.size() << " characters\n";
        std::cerr << "Decompressed file size: " << decompressedWithAmbiguous.size() << " characters\n";
        return;
    }

    int mismatchCount = 0;
    int ambiguousCount = 0;
    int firstMismatchPos = -1;

    for (size_t i = 0; i < targetWithAmbiguous.size(); ++i)
    {
        if (targetWithAmbiguous[i] == 'N')
        {
            ambiguousCount++;
        }

        if (targetWithAmbiguous[i] != decompressedWithAmbiguous[i])
        {
            mismatchCount++;
            if (firstMismatchPos == -1)
            {
                firstMismatchPos = i;
            }

            if (mismatchCount <= 10)
            {
                std::cerr << "Mismatch at position " << i << ": "
                          << "target=" << targetWithAmbiguous[i] << ", "
                          << "decompressed=" << decompressedWithAmbiguous[i] << "\n";
            }
        }
    }

    if (mismatchCount == 0)
    {
        std::cout << "Validation PASSED: Decompressed sequence matches target perfectly!\n";
        std::cout << "The output contains " << ambiguousCount << " 'N' characters.\n";
    }
    else
    {
        std::cerr << "Validation FAILED: " << mismatchCount << " mismatches found ("
                  << (static_cast<double>(mismatchCount) * 100.0 / targetWithAmbiguous.size())
                  << "% error rate)\n";

        if (firstMismatchPos >= 0)
        {
            int start = std::max(0, firstMismatchPos - 10);
            int end = std::min(static_cast<int>(targetWithAmbiguous.size()), firstMismatchPos + 11);

            std::cerr << "Context around first mismatch (position " << firstMismatchPos << "):\n";
            std::cerr << "Target:      ";
            for (int i = start; i < end; i++)
            {
                if (i == firstMismatchPos)
                    std::cerr << "[";
                std::cerr << targetWithAmbiguous[i];
                if (i == firstMismatchPos)
                    std::cerr << "]";
            }
            std::cerr << "\nDecompressed: ";
            for (int i = start; i < end; i++)
            {
                if (i == firstMismatchPos)
                    std::cerr << "[";
                std::cerr << decompressedWithAmbiguous[i];
                if (i == firstMismatchPos)
                    std::cerr << "]";
            }
            std::cerr << "\n";
        }
    }
}

/**
 * Decompresses a compressed genome using the reference sequence
 * Reconstructs the original sequence and restores ambiguous bases
 * @param referenceFile Path to the reference genome
 * @param inputFile Path to the compressed input file
 * @param outputFile Path to write the decompressed output
 *
 * Written by Karlo Pohajda
 */
void decompressGenome(const std::string &referenceFile,
                      const std::string &inputFile,
                      const std::string &outputFile)
{
    loadReferenceGenome(referenceFile);

    std::ifstream in(inputFile, std::ios::binary);
    if (!in)
    {
        std::cerr << "Error: cannot open " << inputFile << "\n";
        return;
    }

    uint32_t metadataLength;
    in.read(reinterpret_cast<char *>(&metadataLength), 4);
    std::vector<char> metadataBuffer(metadataLength);
    in.read(metadataBuffer.data(), metadataLength);
    std::string decompressedMetadata(metadataBuffer.begin(), metadataBuffer.end());

    ambiguousBaseIntervals.clear();
    uint32_t ambiguousCount;
    in.read(reinterpret_cast<char *>(&ambiguousCount), 4);

    for (uint32_t i = 0; i < ambiguousCount; ++i)
    {
        uint32_t start, length;
        char baseChar;

        in.read(reinterpret_cast<char *>(&start), 4);
        in.read(reinterpret_cast<char *>(&length), 4);
        in.read(&baseChar, 1); // single ASCII letter: ’N’, ’R’, etc.

        ambiguousBaseIntervals.push_back(
            AmbiguousInterval{static_cast<int>(start),
                              static_cast<int>(length),
                              baseChar});
    }

    decompressedSequence.clear();
    int32_t currentPos = 0;
    int32_t lastRefEnd = 0;
    int matchCount = 0;
    int literalCount = 0;
    int totalMatchLen = 0;
    int totalLitLen = 0;

    while (true)
    {
        uint8_t tag;
        if (!in.read(reinterpret_cast<char *>(&tag), 1))
            break;

        if (tag == 0x00)
        {
            uint32_t length;
            in.read(reinterpret_cast<char *>(&length), 4);
            std::vector<char> buffer(length);
            in.read(buffer.data(), length);
            decompressedSequence.insert(
                decompressedSequence.end(),
                buffer.begin(), buffer.end());
            currentPos += static_cast<int32_t>(length);
            literalCount += 1;
            totalLitLen += length;
        }
        else if (tag == 0x01)
        {
            int32_t offset;
            uint32_t mLen;
            in.read(reinterpret_cast<char *>(&offset), 4);
            in.read(reinterpret_cast<char *>(&mLen), 4);

            int32_t fullMatchLen = static_cast<int32_t>(mLen) + K_MER_LENGTH;
            int32_t refPos = lastRefEnd + offset;

            if (refPos < 0 ||
                refPos + fullMatchLen > static_cast<int32_t>(referenceSequence.size()))
            {
                std::cerr << "Decompress error: invalid ref_pos="
                          << refPos << " match_len=" << fullMatchLen
                          << " (ref size=" << referenceSequence.size() << ")\n";
                std::exit(1);
            }

            for (int32_t k = 0; k < fullMatchLen; ++k)
            {
                decompressedSequence.push_back(
                    referenceSequence[refPos + k]);
            }
            currentPos += fullMatchLen;
            lastRefEnd = refPos + fullMatchLen;
            matchCount += 1;
            totalMatchLen += fullMatchLen;
        }
        else
        {
            std::cerr << "Error: unknown tag byte " << int(tag) << "\n";
            return;
        }
    }

    std::cout << "Decompression stats: " << matchCount
              << " matches, " << literalCount << " literals\n";
    std::cout << "Match bytes: " << totalMatchLen
              << ", Literal bytes: " << totalLitLen << "\n";

    int totalAmbiguousBases = 0;
    for (const auto &interval : ambiguousBaseIntervals)
        totalAmbiguousBases += interval.length;

    std::string textSequence;
    textSequence.reserve(decompressedSequence.size() + totalAmbiguousBases);

    size_t currentIndex = 0;
    size_t insertedAmbiguous = 0;

    for (const auto &interval : ambiguousBaseIntervals)
    {
        int originalStart = interval.start;
        int length = interval.length;
        char baseChar = interval.base; // e.g. 'R', 'Y', …

        int adjustedStart = originalStart - static_cast<int>(insertedAmbiguous);

        while (currentIndex < decompressedSequence.size() &&
               static_cast<int>(currentIndex) < adjustedStart)
        {
            char letter = nucleotideMap[static_cast<unsigned char>(decompressedSequence[currentIndex])];
            textSequence.push_back(letter);
            ++currentIndex;
        }

        for (int i = 0; i < length; ++i)
        {
            textSequence.push_back(baseChar);
        }
        insertedAmbiguous += length;
    }

    while (currentIndex < decompressedSequence.size())
    {
        char letter = nucleotideMap[static_cast<unsigned char>(decompressedSequence[currentIndex])];
        textSequence.push_back(letter);
        ++currentIndex;
    }

    std::ofstream out(outputFile);
    if (!out)
    {
        std::cerr << "Error: cannot write " << outputFile << "\n";
        return;
    }
    out << decompressedMetadata << "\n";

    int column = 0;
    for (size_t i = 0; i < textSequence.size(); ++i)
    {
        out << textSequence[i];
        if (++column == basesPerLine && i < textSequence.size() - 1)
        {
            out << '\n';
            column = 0;
        }
    }

    std::cout << "[decompress] wrote "
              << textSequence.size() << " bases (including ambiguous) to "
              << outputFile << "\n";
}

/**
 * Main entry point for the HIRGC compressor
 * Supports three modes: compress, decompress, or both
 * Usage:
 *   compress: <program> compress <reference_file> <target_file> <output_file>
 *   decompress: <program> decompress <reference_file> <compressed_file> <output_decompressed_file>
 *   both: <program> both <reference_file> <target_file> <output_compressed_file>
 *
 * Written by Karlo Pohajda
 */
int main(int argc, char *argv[])
{
    std::string mode = "both";
    std::string referenceFile = "./data/ref.fa";
    std::string targetFile = "./data/target.fa";
    std::string compressedFile = "./data/out.hrgc";
    std::string decompressedFile = "./data/decompressed.hrgc";

    if (argc > 1)
    {
        mode = argv[1];
        if (mode != "compress" && mode != "decompress" && mode != "both")
        {
            std::cerr << "Error: First argument must be 'compress', 'decompress', or 'both'\n";
            return 1;
        }

        if (mode == "compress" && argc != 5)
        {
            std::cerr << "Error: For compress mode, all arguments are required:\n"
                      << "Usage: " << argv[0] << " compress <reference_file> <target_file> <output_compressed_file>\n";
            return 1;
        }
        else if (mode == "decompress" && argc != 5)
        {
            std::cerr << "Error: For decompress mode, all arguments are required:\n"
                      << "Usage: " << argv[0] << " decompress <reference_file> <compressed_file> <output_decompressed_file>\n";
            return 1;
        }
        else if (mode == "both" && argc != 5)
        {
            std::cerr << "Error: For both mode, all arguments are required:\n"
                      << "Usage: " << argv[0] << " both <reference_file> <target_file> <output_compressed_file>\n";
            return 1;
        }

        referenceFile = argv[2];
        if (mode == "compress" || mode == "both")
        {
            targetFile = argv[3];
            compressedFile = argv[4];
        }
        else if (mode == "decompress")
        {
            compressedFile = argv[3];
            decompressedFile = argv[4];
        }
    }

    auto startTime = std::chrono::steady_clock::now();

    if (mode == "compress")
    {
        if (std::remove(compressedFile.c_str()) != 0)
            perror(("Error deleting " + compressedFile).c_str());
        if (std::remove(decompressedFile.c_str()) != 0)
            perror(("Error deleting " + decompressedFile).c_str());

        compressGenome(referenceFile, targetFile, compressedFile);
    }
    else if (mode == "decompress")
    {
        if (std::remove(decompressedFile.c_str()) != 0)
            perror(("Error deleting " + decompressedFile).c_str());

        decompressGenome(referenceFile, compressedFile, decompressedFile);
    }
    else if (mode == "both")
    {
        if (std::remove(compressedFile.c_str()) != 0)
            perror(("Error deleting " + compressedFile).c_str());
        if (std::remove(decompressedFile.c_str()) != 0)
            perror(("Error deleting " + decompressedFile).c_str());

        compressGenome(referenceFile, targetFile, compressedFile);

        // loadTargetGenome(targetFile);

        decompressGenome(referenceFile, compressedFile, decompressedFile);

        validateDecompressedGenome(targetFile, decompressedFile);
    }

    auto endTime = std::chrono::steady_clock::now();
    double elapsedSeconds = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "Operation took " << elapsedSeconds << " seconds." << std::endl;

    return 0;
}
