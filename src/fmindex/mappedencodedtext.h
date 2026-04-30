/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Lore Depuydt <lore.depuydt@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef MAPPEDENCODEDTEXT_H
#define MAPPEDENCODEDTEXT_H

#include "../alphabet.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

template <size_t S> class MappedEncodedText {
  private:
#if !defined(__clang__)
    static
#endif
        const uint64_t B = ceil(log2(S));
#if !defined(__clang__)
    static
#endif
        const uint64_t bitmask = (-1ull) ^ (-1ull >> B);
#if !defined(__clang__)
    static
#endif
        const uint64_t hasOverflowBits =
            ~(-1ull << (B - 1));
    const static std::array<uint64_t, 2> overflowMasks;

    size_t tSize = 0;
    size_t vectorSize = 0;
    const uint64_t* encodedText = nullptr;
    Alphabet<S> sigma;

#ifdef _WIN32
    HANDLE hFile = INVALID_HANDLE_VALUE;
    HANDLE hMapFile = NULL;
#else
    int fd = -1;
#endif
    void* mapped = nullptr;
    size_t mappedSize = 0;

#if !defined(__clang__)
    static
#endif
        uint64_t
        hasOverflow(uint64_t index)
#if defined(__clang__)
            const
#endif
    {
        assert(index < 64);
        uint64_t maskIndex =
            ((1ull << (63 - index)) & hasOverflowBits) >> (63 - index);
        return overflowMasks[maskIndex];
    }

    void closeHandles() {
#ifdef _WIN32
        if (mapped != nullptr) {
            UnmapViewOfFile(mapped);
            mapped = nullptr;
        }
        if (hMapFile != NULL) {
            CloseHandle(hMapFile);
            hMapFile = NULL;
        }
        if (hFile != INVALID_HANDLE_VALUE) {
            CloseHandle(hFile);
            hFile = INVALID_HANDLE_VALUE;
        }
#else
        if (mapped != nullptr) {
            munmap(mapped, mappedSize);
            mapped = nullptr;
        }
        if (fd != -1) {
            close(fd);
            fd = -1;
        }
#endif
        encodedText = nullptr;
        tSize = 0;
        vectorSize = 0;
        mappedSize = 0;
    }

  public:
    MappedEncodedText() {
        std::vector<length_t> charCounts(NUM_CHAR, 0);
        charCounts[static_cast<unsigned char>('$')] = 1;
        charCounts[static_cast<unsigned char>('A')] = 1;
        charCounts[static_cast<unsigned char>('C')] = 1;
        charCounts[static_cast<unsigned char>('G')] = 1;
        charCounts[static_cast<unsigned char>('T')] = 1;
        sigma = Alphabet<S>(charCounts);
    }

    ~MappedEncodedText() {
        closeHandles();
    }

    MappedEncodedText(const MappedEncodedText&) = delete;
    MappedEncodedText& operator=(const MappedEncodedText&) = delete;

    MappedEncodedText(MappedEncodedText&& other) noexcept
        : tSize(other.tSize), vectorSize(other.vectorSize),
          encodedText(other.encodedText), sigma(other.sigma),
#ifdef _WIN32
          hFile(other.hFile), hMapFile(other.hMapFile),
#else
          fd(other.fd),
#endif
          mapped(other.mapped), mappedSize(other.mappedSize) {
#ifdef _WIN32
        other.hFile = INVALID_HANDLE_VALUE;
        other.hMapFile = NULL;
#else
        other.fd = -1;
#endif
        other.mapped = nullptr;
        other.encodedText = nullptr;
        other.tSize = 0;
        other.vectorSize = 0;
        other.mappedSize = 0;
    }

    MappedEncodedText& operator=(MappedEncodedText&& other) noexcept {
        if (this != &other) {
            closeHandles();
            tSize = other.tSize;
            vectorSize = other.vectorSize;
            encodedText = other.encodedText;
            sigma = other.sigma;
#ifdef _WIN32
            hFile = other.hFile;
            hMapFile = other.hMapFile;
            other.hFile = INVALID_HANDLE_VALUE;
            other.hMapFile = NULL;
#else
            fd = other.fd;
            other.fd = -1;
#endif
            mapped = other.mapped;
            mappedSize = other.mappedSize;
            other.mapped = nullptr;
            other.encodedText = nullptr;
            other.tSize = 0;
            other.vectorSize = 0;
            other.mappedSize = 0;
        }
        return *this;
    }

    bool load(const std::string& filename) {
        closeHandles();

#ifdef _WIN32
        hFile = CreateFile(filename.c_str(), GENERIC_READ, FILE_SHARE_READ,
                           NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (hFile == INVALID_HANDLE_VALUE) {
            return false;
        }

        LARGE_INTEGER fileSize;
        if (!GetFileSizeEx(hFile, &fileSize) || fileSize.QuadPart <= 0) {
            closeHandles();
            return false;
        }
        mappedSize = static_cast<size_t>(fileSize.QuadPart);

        hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
        if (hMapFile == NULL) {
            closeHandles();
            return false;
        }

        mapped = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, mappedSize);
        if (mapped == NULL) {
            closeHandles();
            return false;
        }
#else
        fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            return false;
        }

        struct stat sb;
        if (fstat(fd, &sb) == -1 || sb.st_size <= 0) {
            closeHandles();
            return false;
        }
        mappedSize = static_cast<size_t>(sb.st_size);

        mapped = mmap(NULL, mappedSize, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped == MAP_FAILED) {
            mapped = nullptr;
            closeHandles();
            return false;
        }
#endif

        if (mappedSize < (2 * sizeof(size_t))) {
            closeHandles();
            return false;
        }

        const char* bytes = static_cast<const char*>(mapped);
        std::memcpy(&tSize, bytes, sizeof(size_t));
        std::memcpy(&vectorSize, bytes + sizeof(size_t), sizeof(size_t));

        const size_t expectedSize =
            (2 * sizeof(size_t)) + (vectorSize * sizeof(uint64_t));
        if (mappedSize != expectedSize) {
            closeHandles();
            return false;
        }

        encodedText = reinterpret_cast<const uint64_t*>(bytes + (2 * sizeof(size_t)));
        return true;
    }

    size_t size() const {
        return tSize;
    }

    uint64_t getEncodedLetter(const uint64_t index) const {
        if (encodedText == nullptr || index >= tSize) {
            throw std::out_of_range("MappedEncodedText index out of bounds");
        }

        uint64_t w = (index * B) / 64;
        uint64_t b = (index * B) % 64;
        uint64_t bits = (encodedText[w] & (bitmask >> b)) << b;
        uint64_t mask = hasOverflow(b) & ((-1ull) ^ (-1ull >> (B - (64 - b))));
        uint64_t bitsNext = encodedText[w + 1] & mask;

        return (bits >> (64 - B)) +
               (bitsNext >> (64 - __builtin_popcountll(mask)));
    }

    char decodeLetter(const uint64_t index) const {
        return sigma.i2c(getEncodedLetter(index));
    }

    std::string decodeSubstring(size_t start, size_t length) const {
        if (start > tSize || start + length > tSize) {
            throw std::out_of_range("MappedEncodedText substring out of bounds");
        }

        std::string substring(length, '\0');
        for (size_t i = 0; i < length; ++i) {
            substring[i] = decodeLetter(start + i);
        }
        return substring;
    }
};

template <size_t S>
const std::array<uint64_t, 2> MappedEncodedText<S>::overflowMasks = {0ull, -1ull};

#endif
