// hirgc_single.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct Command {
    enum Type { COPY = 0, LITERAL = 1 } type;
    size_t pos;
    size_t length;
    std::string literals;
};

class EntropyCoder {
public:
    void encode(const std::vector<Command>& commands, std::ostream& os) {
        // to be implemented
    }

    void decode(std::istream& is, std::vector<Command>& commands) {
        // to be implemented
    }
};

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " compress <ref.fa> <target.fa> <out.bin>\n"
                  << "  " << argv[0] << " decompress <ref.fa> <in.bin> <out.fa>\n";
        return 1;
    }

    std::string mode = argv[1];
    if (mode == "compress") {
        // compress logic to be implemented
    } else if (mode == "decompress") {
        // decompress logic to be implemented
    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }

    return 0;
}
