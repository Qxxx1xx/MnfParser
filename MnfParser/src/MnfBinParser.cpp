#include <MnfParser/MnfBinParser.h>

namespace MnfParser {
MnfBinParser::MnfBinParser(const std::string& file_path)
    : mnf_file{ std::ifstream(file_path, std::ios::binary) }
{
    if (!mnf_file) {
        throw std::runtime_error("Could not open file " + file_path);
    }
    checkMnfVersion();
}
} // namespace MnfParser