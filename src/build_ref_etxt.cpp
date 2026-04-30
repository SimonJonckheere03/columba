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

#include "alphabet.h"
#include "fmindex/encodedtext.h"
#include "seqfile.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace std;

namespace {

char replaceNonACGT(char original, std::minstd_rand& gen,
                    const std::string& seed, size_t& seedIndex) {
    (void)seed;
    (void)seedIndex;
    static const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        std::uniform_int_distribution<size_t> distribution(
            0, validChars.length() - 1);
        return validChars[distribution(gen)];
    }
    return original;
}

char replaceNonACGTWithSeed(char original, std::minstd_rand& gen,
                            const std::string& seed, size_t& seedIndex) {
    (void)gen;
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        char replacement = seed[seedIndex];
        seedIndex = (seedIndex + 1) % seed.size();
        return replacement;
    }
    seedIndex = 0;
    return original;
}

bool isGzipped(FileType fileType, const string& filename) {
    if (fileType == FileType::FASTA_GZ) {
#ifndef HAVE_ZLIB
        throw runtime_error("Error: " + filename +
                            " is gzipped, but Columba was built without zlib.");
#else
        return true;
#endif
    }
    return false;
}

SeqFile openFastaFile(const string& fastaFile) {
    FileType fileType;
    tie(fileType, ignore) = getFileType(fastaFile);

    if (fileType != FileType::FASTA && fileType != FileType::FASTA_GZ) {
        throw runtime_error("Error: " + fastaFile + " is not a FASTA file.");
    }

    SeqFile file(isGzipped(fileType, fastaFile));
    file.open(fastaFile);
    return file;
}

void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::function<char(char, size_t&)> replaceFunc) {
    SeqFile file = openFastaFile(fastaFile);

    std::string sequence;
    std::string line;
    size_t seedIndex = 0;

    while (file.good()) {
        file.getLine(line);
        if (line.empty() || (line.size() == 1 && line[0] == '\n')) {
            continue;
        }
        if (line.back() == '\n') {
            line.pop_back();
        }

        if (line[0] == '>') {
            if (!sequence.empty()) {
                for (char& c : sequence) {
                    c = replaceFunc(static_cast<char>(std::toupper(
                                        static_cast<unsigned char>(c))),
                                    seedIndex);
                    concatenation += c;
                }
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }

    for (char& c : sequence) {
        c = replaceFunc(static_cast<char>(std::toupper(static_cast<unsigned char>(c))),
                        seedIndex);
        concatenation += c;
    }

    file.close();
}

void writeReferencePrefixEncodedText(const string& fastaFile,
                                     const string& outputFile,
                                     size_t seedLength) {
    std::minstd_rand gen(42);
    std::string seed;
    size_t seedIndex = 0;
    for (size_t i = 0; i < seedLength; ++i) {
        seed += replaceNonACGT('N', gen, seed, seedIndex);
    }

    std::function<char(char, size_t&)> replaceFunc;
    if (seedLength == 0) {
        replaceFunc = [&gen, &seed](char c, size_t& index) -> char {
            return replaceNonACGT(c, gen, seed, index);
        };
    } else {
        replaceFunc = [&gen, &seed](char c, size_t& index) -> char {
            return replaceNonACGTWithSeed(c, gen, seed, index);
        };
    }

    std::string text;
    concatenateAndTransform(fastaFile, text, replaceFunc);

    std::vector<length_t> charCounts(NUM_CHAR, 0);
    charCounts[static_cast<unsigned char>('$')] = 1;
    charCounts[static_cast<unsigned char>('A')] = 1;
    charCounts[static_cast<unsigned char>('C')] = 1;
    charCounts[static_cast<unsigned char>('G')] = 1;
    charCounts[static_cast<unsigned char>('T')] = 1;
    Alphabet<ALPHABET> sigma(charCounts);
    EncodedText<ALPHABET> encodedText(sigma, text);
    encodedText.write(outputFile);
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <reference.fa[.gz]> <output.ref.etxt> <seedLength>\n";
        return 1;
    }

    try {
        const string fastaFile = argv[1];
        const string outputFile = argv[2];
        const size_t seedLength = static_cast<size_t>(stoull(argv[3]));
        writeReferencePrefixEncodedText(fastaFile, outputFile, seedLength);
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << '\n';
        return 1;
    }

    return 0;
}
