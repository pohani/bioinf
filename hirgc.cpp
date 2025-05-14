// hirgc_single.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

struct Command
{
    enum Type
    {
        COPY = 0,
        LITERAL = 1
    } type;
    size_t pos;
    size_t length;
    std::string literals;
};

class EntropyCoder
{
public:
    void encode(const std::vector<Command> &commands, std::ostream &os)
    {
        uint64_t count = commands.size();
        os.write(reinterpret_cast<const char *>(&count), sizeof(count));
        for (const auto &cmd : commands)
        {
            uint8_t t = static_cast<uint8_t>(cmd.type);
            os.write(reinterpret_cast<const char *>(&t), 1);
            os.write(reinterpret_cast<const char *>(&cmd.pos), sizeof(cmd.pos));
            os.write(reinterpret_cast<const char *>(&cmd.length), sizeof(cmd.length));
            if (cmd.type == Command::LITERAL)
            {
                os.write(cmd.literals.data(), cmd.literals.size());
            }
        }
    }

    void decode(std::istream &is, std::vector<Command> &commands)
    {
        // implemented earlier
    }
};

class ReferenceIndex
{
    // implemented earlier
};

class Compressor
{
    const ReferenceIndex &refIndex_;
    size_t min_match_;

public:
    Compressor(const ReferenceIndex &refIndex, size_t min_match = 28)
        : refIndex_(refIndex), min_match_(min_match) {}

    bool compress(const std::string &target_fasta, const std::string &output_file)
    {
        std::ifstream in(target_fasta);
        std::string line, target;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '>')
                continue;
            target += line;
        }
        if (target.empty())
            return false;

        std::vector<Command> commands;
        size_t n = target.size();
        size_t i = 0;
        while (i < n)
        {
            size_t best_len = 0, best_pos = 0;
            if (i + refIndex_.kmer_size() <= n)
            {
                std::string kmer = target.substr(i, refIndex_.kmer_size());
                for (size_t pos : refIndex_.query(kmer))
                {
                    size_t len = 0;
                    while (i + len < n &&
                           pos + len < refIndex_.sequence().size() &&
                           target[i + len] == refIndex_.sequence()[pos + len])
                    {
                        ++len;
                    }
                    if (len > best_len)
                    {
                        best_len = len;
                        best_pos = pos;
                    }
                }
            }

            if (best_len >= min_match_)
            {
                commands.push_back({Command::COPY, best_pos, best_len, ""});
                i += best_len;
            }
            else
            {
                commands.push_back({Command::LITERAL, 0, 1, std::string(1, target[i])});
                ++i;
            }
        }

        std::ofstream out(output_file, std::ios::binary);
        if (!out)
            return false;
        EntropyCoder coder;
        coder.encode(commands, out);
        return true;
    }
};

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " compress <ref.fa> <target.fa> <out.bin>\n"
                  << "  " << argv[0] << " decompress <ref.fa> <in.bin> <out.fa>\n";
        return 1;
    }

    std::string mode = argv[1];
    std::string ref_path = argv[2];
    ReferenceIndex idx(ref_path);
    if (!idx.load())
    {
        std::cerr << "Failed to load reference: " << ref_path << std::endl;
        return 1;
    }

    if (mode == "compress")
    {
        Compressor comp(idx);
        if (!comp.compress(argv[3], argv[4]))
        {
            std::cerr << "Compression failed" << std::endl;
            return 1;
        }
    }
    else if (mode == "decompress")
    {
        // decompress logic to be implemented
    }
    else
    {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }

    return 0;
}
