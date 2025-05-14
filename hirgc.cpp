// hirgc_single.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

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
