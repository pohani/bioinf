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
    void encode(const std::vector<Command> &commands, std::ostream &os) { /* … */ }
    void decode(std::istream &is, std::vector<Command> &commands) { /* … */ }
};

class ReferenceIndex
{
    // as before
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
        // to be implemented
        return false;
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
