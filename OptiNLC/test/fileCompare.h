#include <fstream>
#include <string>

inline bool filesAreEqual(const std::string& path1, const std::string& path2) {
    std::ifstream file1(path1, std::ios::binary);
    std::ifstream file2(path2, std::ios::binary);

    if (!file1.is_open() || !file2.is_open()) {
        return false;
    }

    file1.seekg(0, std::ios::end);
    file2.seekg(0, std::ios::end);
    if (file1.tellg() != file2.tellg()) {
        return false;
    }

    file1.seekg(0);
    file2.seekg(0);

    constexpr std::size_t bufferSize = 4096;
    char buffer1[bufferSize];
    char buffer2[bufferSize];

    while (file1 && file2) {
        file1.read(buffer1, bufferSize);
        file2.read(buffer2, bufferSize);

        if (file1.gcount() != file2.gcount()) return false;
        if (std::memcmp(buffer1, buffer2, file1.gcount()) != 0) return false;
    }

    return true;
}

