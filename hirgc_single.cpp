// hirgc_single.cpp
// Single-file implementation of HiRGC referential genome compressor/decompressor

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// --- Command & EntropyCoder ---
struct Command {
    enum Type { COPY = 0, LITERAL = 1 } type;
    size_t pos;
    size_t length;
    std::string literals;
};

class EntropyCoder {
public:
    void encode(const std::vector<Command>& commands, std::ostream& os) {
        uint64_t count = commands.size();
        os.write(reinterpret_cast<const char*>(&count), sizeof(count));
        for (const auto& cmd : commands) {
            uint8_t t = static_cast<uint8_t>(cmd.type);
            os.write(reinterpret_cast<const char*>(&t), 1);
            os.write(reinterpret_cast<const char*>(&cmd.pos), sizeof(cmd.pos));
            os.write(reinterpret_cast<const char*>(&cmd.length), sizeof(cmd.length));
            if (cmd.type == Command::LITERAL) {
                os.write(cmd.literals.data(), cmd.literals.size());
            }
        }
    }

    void decode(std::istream& is, std::vector<Command>& commands) {
        uint64_t count;
        is.read(reinterpret_cast<char*>(&count), sizeof(count));
        commands.clear();
        commands.reserve(count);
        for (uint64_t i = 0; i < count; ++i) {
            uint8_t t;
            is.read(reinterpret_cast<char*>(&t), 1);
            Command cmd;
            cmd.type = static_cast<Command::Type>(t);
            is.read(reinterpret_cast<char*>(&cmd.pos), sizeof(cmd.pos));
            is.read(reinterpret_cast<char*>(&cmd.length), sizeof(cmd.length));
            if (cmd.type == Command::LITERAL) {
                cmd.literals.resize(cmd.length);
                is.read(&cmd.literals[0], cmd.length);
            }
            commands.push_back(cmd);
        }
    }
};

// --- ReferenceIndex ---
class ReferenceIndex {
    std::string fasta_path_;
    size_t k_;
    std::string reference_sequence_;
    std::unordered_map<std::string, std::vector<size_t>> index_;

public:
    ReferenceIndex(const std::string& fasta_path, size_t k = 28)
        : fasta_path_(fasta_path), k_(k) {}

    bool load() {
        std::ifstream in(fasta_path_);
        if (!in) return false;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '>') continue;
            std::istringstream iss(line);
            std::string chunk;
            while (iss >> chunk) reference_sequence_ += chunk;
        }
        if (reference_sequence_.empty()) return false;

        // build k-mer index
        size_t n = reference_sequence_.size();
        for (size_t i = 0; i + k_ <= n; ++i) {
            std::string kmer = reference_sequence_.substr(i, k_);
            index_[kmer].push_back(i);
        }
        return true;
    }

    const std::vector<size_t>& query(const std::string& kmer) const {
        static const std::vector<size_t> empty;
        auto it = index_.find(kmer);
        return it == index_.end() ? empty : it->second;
    }

    size_t kmer_size() const { return k_; }
    const std::string& sequence() const { return reference_sequence_; }
};

// --- Compressor ---
class Compressor {
    const ReferenceIndex& refIndex_;
    size_t min_match_;

public:
    Compressor(const ReferenceIndex& refIndex, size_t min_match = 28)
        : refIndex_(refIndex), min_match_(min_match) {}

    bool compress(const std::string& target_fasta, const std::string& output_file) {
        std::ifstream in(target_fasta);
        std::string line, target;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '>') continue;
            target += line;
        }
        if (target.empty()) return false;

        std::vector<Command> commands;
        size_t n = target.size();
        size_t i = 0;
        while (i < n) {
            size_t best_len = 0;
            size_t best_pos = 0;

            if (i + refIndex_.kmer_size() <= n) {
                std::string kmer = target.substr(i, refIndex_.kmer_size());
                for (size_t pos : refIndex_.query(kmer)) {
                    size_t len = 0;
                    while (i + len < n && pos + len < refIndex_.sequence().size() &&
                           target[i + len] == refIndex_.sequence()[pos + len]) {
                        ++len;
                    }
                    if (len > best_len) {
                        best_len = len;
                        best_pos = pos;
                    }
                }
            }

            if (best_len >= min_match_) {
                commands.push_back({Command::COPY, best_pos, best_len, ""});
                i += best_len;
            } else {
                commands.push_back({Command::LITERAL, 0, 1, std::string(1, target[i])});
                ++i;
            }
        }

        std::ofstream out(output_file, std::ios::binary);
        if (!out) return false;
        EntropyCoder coder;
        coder.encode(commands, out);
        return true;
    }
};

// --- Decompressor ---
class Decompressor {
    const ReferenceIndex& refIndex_;

public:
    Decompressor(const ReferenceIndex& refIndex) : refIndex_(refIndex) {}

    bool decompress(const std::string& compressed_file, const std::string& output_fasta) {
        std::ifstream in(compressed_file, std::ios::binary);
        if (!in) return false;

        std::vector<Command> commands;
        EntropyCoder coder;
        coder.decode(in, commands);

        std::string out_seq;
        for (const auto& cmd : commands) {
            if (cmd.type == Command::COPY) {
                out_seq += refIndex_.sequence().substr(cmd.pos, cmd.length);
            } else {
                out_seq += cmd.literals;
            }
        }

        std::ofstream out(output_fasta);
        if (!out) return false;
        out << ">decompressed\n";
        for (size_t i = 0; i < out_seq.size(); i += 80)
            out << out_seq.substr(i, 80) << '\n';
        return true;
    }
};

// --- main ---
int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage:\n  " << argv[0] << " compress <ref.fa> <target.fa> <out.bin>\n";
        std::cerr << "  " << argv[0] << " decompress <ref.fa> <in.bin> <out.fa>\n";
        return 1;
    }

    std::string mode = argv[1];
    std::string ref_path = argv[2];
    ReferenceIndex idx(ref_path);
    if (!idx.load()) {
        std::cerr << "Failed to load reference: " << ref_path << std::endl;
        return 1;
    }

    if (mode == "compress") {
        Compressor comp(idx);
        if (!comp.compress(argv[3], argv[4])) {
            std::cerr << "Compression failed" << std::endl;
            return 1;
        }
    } else if (mode == "decompress") {
        Decompressor decomp(idx);
        if (!decomp.decompress(argv[3], argv[4])) {
            std::cerr << "Decompression failed" << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }

    return 0;
}
