// hirgc_single.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

struct Command {
    enum Type { COPY = 0, LITERAL = 1 } type;
    size_t pos;
    size_t length;
    std::string literals;
};

class EntropyCoder {
public:
    void encode(const std::vector<Command>& commands, std::ostream& os) { /* … */ }
    void decode(std::istream& is, std::vector<Command>& commands) { /* … */ }
};

class ReferenceIndex {
    std::string fasta_path_;
    size_t k_;
    std::string reference_sequence_;
    std::unordered_map<std::string, std::vector<size_t>> index_;

public:
    ReferenceIndex(const std::string& fasta_path, size_t k = 28)
        : fasta_path_(fasta_path), k_(k) {}

    bool load() {
        // to be implemented
        return false;
    }

    const std::vector<size_t>& query(const std::string& kmer) const {
        static const std::vector<size_t> empty;
        return empty;
    }

    size_t kmer_size() const { return k_; }
    const std::string& sequence() const { return reference_sequence_; }
};

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " compress <ref.fa> <target.fa> <out.bin>\n"
                  << "  " << argv[0] << " decompress <ref.fa> <in.bin> <out.fa>\n";
        return 1;
    }

    std::string mode = argv[1];
    // reference loading and compress/decompress logic to be implemented

    return 0;
}
