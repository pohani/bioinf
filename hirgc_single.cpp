// hirgc_single_modified.cpp
// Single-file implementation of a referential genome compressor/decompressor.
// Written by Karlo Pohajda.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// -------------------------------------------------------------------------------------------------
// GenomeCommand
// Represents a single instruction: either copy from reference or emit a literal.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
struct GenomeCommand
{
    enum class Kind
    {
        Copy = 0,
        Literal = 1
    } kind;                 ///< Type of this command.
    size_t refPosition;     ///< Position in reference to copy from.
    size_t length;          ///< Number of bases to copy or length of literals.
    std::string literalSeq; ///< Sequence to output if Kind::Literal.
};

// -------------------------------------------------------------------------------------------------
// EntropyProcessor
// Handles binary serialization and deserialization of GenomeCommand streams.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
class EntropyProcessor
{
public:
    // -------------------------------------------------------------------------------------------------
    // encode
    // Serialize a sequence of commands into a binary stream.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    void encode(const std::vector<GenomeCommand> &cmds, std::ostream &os)
    {
        uint64_t total = cmds.size();
        os.write(reinterpret_cast<const char *>(&total), sizeof(total));

        for (const auto &cmd : cmds)
        {
            uint8_t t = static_cast<uint8_t>(cmd.kind);
            os.write(reinterpret_cast<const char *>(&t), sizeof(t));
            os.write(reinterpret_cast<const char *>(&cmd.refPosition), sizeof(cmd.refPosition));
            os.write(reinterpret_cast<const char *>(&cmd.length), sizeof(cmd.length));
            if (cmd.kind == GenomeCommand::Kind::Literal)
            {
                os.write(cmd.literalSeq.data(), cmd.literalSeq.size());
            }
        }
    }

    // -------------------------------------------------------------------------------------------------
    // decode
    // Read back a sequence of commands from a binary stream.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    void decode(std::istream &is, std::vector<GenomeCommand> &cmds)
    {
        uint64_t total = 0;
        is.read(reinterpret_cast<char *>(&total), sizeof(total));
        cmds.clear();
        cmds.reserve(total);

        for (uint64_t i = 0; i < total; ++i)
        {
            uint8_t t = 0;
            is.read(reinterpret_cast<char *>(&t), sizeof(t));

            GenomeCommand cmd;
            cmd.kind = static_cast<GenomeCommand::Kind>(t);
            is.read(reinterpret_cast<char *>(&cmd.refPosition), sizeof(cmd.refPosition));
            is.read(reinterpret_cast<char *>(&cmd.length), sizeof(cmd.length));

            if (cmd.kind == GenomeCommand::Kind::Literal)
            {
                cmd.literalSeq.resize(cmd.length);
                is.read(&cmd.literalSeq[0], cmd.length);
            }
            cmds.push_back(std::move(cmd));
        }
    }
};

// -------------------------------------------------------------------------------------------------
// ReferenceIndexer
// Loads a reference FASTA and builds a k-mer to positions index.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
class ReferenceIndexer
{
public:
    // -------------------------------------------------------------------------------------------------
    // Constructor
    // pathToFasta: file path of the reference FASTA
    // kmerSize: length of k-mers for indexing
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    ReferenceIndexer(std::string pathToFasta, size_t kmerSize = 28)
        : m_fastaPath(std::move(pathToFasta)), m_kmerSize(kmerSize)
    {
    }

    // -------------------------------------------------------------------------------------------------
    // buildIndex
    // Loads the FASTA, concatenates sequences, and builds the k-mer index.
    // Returns true on success.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    bool buildIndex()
    {
        std::ifstream in(m_fastaPath);
        if (!in)
        {
            std::cerr << "Error: cannot open reference file '" << m_fastaPath << "'\n";
            return false;
        }

        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '>')
                continue;
            m_referenceSeq += line;
        }
        if (m_referenceSeq.empty())
        {
            std::cerr << "Error: no sequence data found in '" << m_fastaPath << "'\n";
            return false;
        }

        // Build k-mer index
        size_t N = m_referenceSeq.size();
        for (size_t i = 0; i + m_kmerSize <= N; ++i)
        {
            // Replace string_view with substr for compatibility
            std::string kmer = m_referenceSeq.substr(i, m_kmerSize);
            m_index[kmer].push_back(i);
        }
        return true;
    }

    // -------------------------------------------------------------------------------------------------
    // queryPositions
    // Returns all start positions in the reference matching the given k-mer.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    const std::vector<size_t> &queryPositions(const std::string &kmer) const
    {
        static const std::vector<size_t> empty;
        auto it = m_index.find(kmer);
        return (it == m_index.end() ? empty : it->second);
    }

    // -------------------------------------------------------------------------------------------------
    // getKmerSize
    // Returns the k-mer length used for indexing.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    size_t getKmerSize() const { return m_kmerSize; }

    // -------------------------------------------------------------------------------------------------
    // getReference
    // Returns the full concatenated reference sequence.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    const std::string &getReference() const { return m_referenceSeq; }

private:
    std::string m_fastaPath;
    size_t m_kmerSize;
    std::string m_referenceSeq;
    std::unordered_map<std::string, std::vector<size_t>> m_index;
};

// -------------------------------------------------------------------------------------------------
// GenomeCompressor
// Produces a compressed binary representation of a target genome using a reference.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
class GenomeCompressor
{
public:
    // -------------------------------------------------------------------------------------------------
    // Constructor
    // refIndexer: loaded reference index
    // minMatch: minimal match length to emit as COPY
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    GenomeCompressor(const ReferenceIndexer &refIndexer, size_t minMatch = 28)
        : m_refIndexer(refIndexer), m_minMatch(minMatch)
    {
    }

    // -------------------------------------------------------------------------------------------------
    // compress
    // Reads target FASTA, generates copy/literal commands, and writes binary to outFile.
    // Returns true on success.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    bool compress(const std::string &targetFasta, const std::string &outFile)
    {
        std::ifstream in(targetFasta);
        if (!in)
        {
            std::cerr << "Error: cannot open target file '" << targetFasta << "'\n";
            return false;
        }

        std::string targetSeq, line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '>')
                continue;
            targetSeq += line;
        }
        if (targetSeq.empty())
        {
            std::cerr << "Error: no sequence data in target '" << targetFasta << "'\n";
            return false;
        }

        std::vector<GenomeCommand> cmds;
        size_t N = targetSeq.size();
        size_t i = 0;
        const auto &refSeq = m_refIndexer.getReference();
        size_t k = m_refIndexer.getKmerSize();

        while (i < N)
        {
            size_t bestLen = 0;
            size_t bestPos = 0;

            if (i + k <= N)
            {
                // Replace string_view with substr for compatibility
                std::string kmer = targetSeq.substr(i, k);
                for (size_t pos : m_refIndexer.queryPositions(kmer))
                {
                    size_t len = 0;
                    while (i + len < N && pos + len < refSeq.size() && targetSeq[i + len] == refSeq[pos + len])
                    {
                        ++len;
                    }
                    if (len > bestLen)
                    {
                        bestLen = len;
                        bestPos = pos;
                    }
                }
            }

            if (bestLen >= m_minMatch)
            {
                cmds.push_back({GenomeCommand::Kind::Copy, bestPos, bestLen, ""});
                i += bestLen;
            }
            else
            {
                cmds.push_back({GenomeCommand::Kind::Literal, 0, 1, std::string(1, targetSeq[i])});
                ++i;
            }
        }

        std::ofstream out(outFile, std::ios::binary);
        if (!out)
        {
            std::cerr << "Error: cannot write to '" << outFile << "'\n";
            return false;
        }

        EntropyProcessor coder;
        coder.encode(cmds, out);
        return true;
    }

private:
    const ReferenceIndexer &m_refIndexer;
    size_t m_minMatch;
};

// -------------------------------------------------------------------------------------------------
// GenomeDecompressor
// Reconstructs the target genome fasta from compressed binary and reference.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
class GenomeDecompressor
{
public:
    // -------------------------------------------------------------------------------------------------
    // Constructor
    // refIndexer: loaded reference index
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    GenomeDecompressor(const ReferenceIndexer &refIndexer)
        : m_refIndexer(refIndexer)
    {
    }

    // -------------------------------------------------------------------------------------------------
    // decompress
    // Reads binary commands, reconstructs sequence, and writes FASTA to outFasta.
    // Returns true on success.
    // Written by Karlo Pohajda.
    // -------------------------------------------------------------------------------------------------
    bool decompress(const std::string &inFile, const std::string &outFasta)
    {
        std::ifstream in(inFile, std::ios::binary);
        if (!in)
        {
            std::cerr << "Error: cannot open compressed file '" << inFile << "'\n";
            return false;
        }

        std::vector<GenomeCommand> cmds;
        EntropyProcessor coder;
        coder.decode(in, cmds);

        std::string rebuilt;
        const auto &refSeq = m_refIndexer.getReference();
        for (const auto &cmd : cmds)
        {
            if (cmd.kind == GenomeCommand::Kind::Copy)
            {
                rebuilt.append(refSeq, cmd.refPosition, cmd.length);
            }
            else
            {
                rebuilt += cmd.literalSeq;
            }
        }

        std::ofstream out(outFasta);
        if (!out)
        {
            std::cerr << "Error: cannot write FASTA to '" << outFasta << "'\n";
            return false;
        }

        out << ">decompressed_sequence\n";
        for (size_t pos = 0; pos < rebuilt.size(); pos += 80)
        {
            out << rebuilt.substr(pos, 80) << "\n";
        }
        return true;
    }

private:
    const ReferenceIndexer &m_refIndexer;
};

// -------------------------------------------------------------------------------------------------
// main
// Parses command-line arguments and dispatches compression or decompression.
// Written by Karlo Pohajda.
// -------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cerr << "Usage:\n  " << argv[0]
                  << " compress <reference.fa> <target.fa> <output.bin>\n"
                     "  "
                  << argv[0]
                  << " decompress <reference.fa> <input.bin> <output.fa>\n";
        return 1;
    }

    std::string mode = argv[1];
    std::string refPath = argv[2];

    ReferenceIndexer indexer(refPath);
    if (!indexer.buildIndex())
    {
        return 1;
    }

    if (mode == "compress")
    {
        GenomeCompressor compressor(indexer);
        return compressor.compress(argv[3], argv[4]) ? 0 : 1;
    }
    else if (mode == "decompress")
    {
        GenomeDecompressor decompressor(indexer);
        return decompressor.decompress(argv[3], argv[4]) ? 0 : 1;
    }
    else
    {
        std::cerr << "Unknown mode: '" << mode << "'. Use 'compress' or 'decompress'.\n";
        return 1;
    }
}
