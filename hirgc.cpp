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
    std::string fasta_path_;
    size_t k_;
    std::string reference_sequence_;
    std::unordered_map<std::string, std::vector<size_t>> index_;

public:
    ReferenceIndex(const std::string &fasta_path, size_t k = 28)
        : fasta_path_(fasta_path), k_(k) {}

    bool load()
    {
        std::ifstream in(fasta_path_);
        if (!in)
            return false;
        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '>')
                continue;
            std::istringstream iss(line);
            std::string chunk;
            while (iss >> chunk)
                reference_sequence_ += chunk;
        }
        if (reference_sequence_.empty())
            return false;

        size_t n = reference_sequence_.size();
        for (size_t i = 0; i + k_ <= n; ++i)
        {
            std::string kmer = reference_sequence_.substr(i, k_);
            index_[kmer].push_back(i);
        }
        return true;
    }

    const std::vector<size_t> &query(const std::string &kmer) const
    {
        static const std::vector<size_t> empty;
        auto it = index_.find(kmer);
        return it == index_.end() ? empty : it->second;
    }

    size_t kmer_size() const { return k_; }
    const std::string &sequence() const { return reference_sequence_; }
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
    // remaining logic to be implemented

    return 0;
}
