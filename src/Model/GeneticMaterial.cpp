#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <optional>
#include <zlib.h>
#include <archive.h>
#include <archive_entry.h>
#include <vector>
#include <minizip/unzip.h>
#include <stdio.h>
#include <cctype>
#include <variant>
#include <numeric> 


#include <optional>
#include <utility>
#include <regex>
#include <cmath>
#include <tuple>


#include <vector>
#include <string>
#include <optional>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>


#include <optional>
#include <vector>
#include <string>


#include <vector>
#include <string>
#include <algorithm>
#include <optional>




#include <utility> // Para std::pair

#include "GeneticMaterial.h"



GeneticMaterial::GeneticMaterial(const std::string& filePath) : filePath(filePath) {
    loadGenomeSequence();
}



std::string GeneticMaterial::getFilePath() const {
    return filePath;
}


bool GeneticMaterial::startsWith(const std::string& str, const std::string& prefix) {
    if (prefix.size() > str.size()) return false;
    return str.compare(0, prefix.size(), prefix) == 0;
}


bool GeneticMaterial::endsWith(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}
   
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Décompression de fichiers ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Função que permite fazer a descompressão do arquivo de extensão .gz
// std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>>
// GeneticMaterial::decompressGzipFile(const std::string& filePath) {
//     gzFile gzFile = gzopen(filePath.c_str(), "rb");
//     if (!gzFile) {
//         throw std::runtime_error("Não foi possível abrir o arquivo GZIP: " + filePath);
//     }

//     constexpr size_t bufferSize = 8192; // 8 KB
//     char buffer[bufferSize];
//     std::ostringstream decompressedContent;

//     int bytesRead;
//     while ((bytesRead = gzread(gzFile, buffer, sizeof(buffer))) > 0) {
//         decompressedContent.write(buffer, bytesRead);
//     }

//     if (bytesRead < 0) { // Verifica erro do gzread
//         int errNum;
//         const char* errorMsg = gzerror(gzFile, &errNum);
//         gzclose(gzFile);
//         throw std::runtime_error("Erro ao ler o arquivo GZIP: " + std::string(errorMsg));
//     }

//     gzclose(gzFile);

//     if (decompressedContent.str().empty()) {
//         throw std::runtime_error("Falha ao descompactar o arquivo GZIP: " + filePath);
//     }

//     // Processa o conteúdo para dividir em FASTA/FASTQ
//     std::vector<GeneticMaterial::FileContent> fastaFiles;
//     std::vector<GeneticMaterial::FileContent> fastqFiles;

//     std::istringstream contentStream(decompressedContent.str());
//     std::string line;
//     GeneticMaterial::FileContent currentFile;

//     while (std::getline(contentStream, line)) {
//         if (line.empty()) continue;
//         if (startsWith(line, ">")) { // Início de um arquivo FASTA
//             if (currentFile.content.has_value()) {
//                 if (currentFile.fileType == "FASTA") {
//                     fastaFiles.push_back(currentFile);
//                 } else if (currentFile.fileType == "FASTQ") {
//                     fastqFiles.push_back(currentFile);
//                 }
//             }
//             currentFile = GeneticMaterial::FileContent();
//             currentFile.fileType = "FASTA";
//             currentFile.content = line + "\n";
//         } else if (startsWith(line, "@")) { // Início de um arquivo FASTQ
//             if (currentFile.content.has_value()) {
//                 if (currentFile.fileType == "FASTA") {
//                     fastaFiles.push_back(currentFile);
//                 } else if (currentFile.fileType == "FASTQ") {
//                     fastqFiles.push_back(currentFile);
//                 }
//             }
//             currentFile = GeneticMaterial::FileContent();
//             currentFile.fileType = "FASTQ";
//             currentFile.content = line + "\n";
//         } else {
//             if (!currentFile.content.has_value()) {
//                 currentFile.content = "";
//             }
//             currentFile.content.value() += line + "\n";
//         }
//     }

//     // Adiciona o último arquivo processado
//     if (currentFile.content.has_value()) {
//         if (currentFile.fileType == "FASTA") {
//             fastaFiles.push_back(currentFile);
//         } else if (currentFile.fileType == "FASTQ") {
//             fastqFiles.push_back(currentFile);
//         }
//     }

//     // Normaliza e valida os arquivos após o processamento
//     for (auto& file : fastaFiles) {
//         if (file.content) {
//             file.content = normalizeNewlines(*file.content);
//             if (!verifySequence(*file.content, "FASTA")) {
//                 std::cerr << "Erro: Arquivo FASTA inválido após descompactação: " << file.fileName << std::endl;
//             }
//         }
//     }

//     for (auto& file : fastqFiles) {
//         if (file.content) {
//             file.content = normalizeNewlines(*file.content);
//             if (!verifySequence(*file.content, "FASTQ")) {
//                 std::cerr << "Erro: Arquivo FASTQ inválido após descompactação: " << file.fileName << std::endl;
//             }
//         }
//     }

//     return {fastaFiles, fastqFiles};
// }

std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>>
GeneticMaterial::decompressGzipFile(const std::string& filePath) {
    std::cout << "Starting decompression of GZIP file: " << filePath << std::endl;

    gzFile gzFile = gzopen(filePath.c_str(), "rb");
    if (!gzFile) {
        throw std::runtime_error("Could not open the GZIP file: " + filePath);
    }

    constexpr size_t bufferSize = 8192; // 8 KB
    char buffer[bufferSize];
    std::ostringstream decompressedContent;

    int bytesRead;
    while ((bytesRead = gzread(gzFile, buffer, sizeof(buffer))) > 0) {
        decompressedContent.write(buffer, bytesRead);
    }

    if (bytesRead < 0) {
        int errNum;
        const char* errorMsg = gzerror(gzFile, &errNum);
        gzclose(gzFile);
        throw std::runtime_error("Error reading the GZIP file: " + std::string(errorMsg));
    }

    gzclose(gzFile);

    if (decompressedContent.str().empty()) {
        throw std::runtime_error("Failed to decompress the GZIP file: " + filePath);
    }

    std::cout << "Decompression completed. Processing content..." << std::endl;

    std::vector<GeneticMaterial::FileContent> fastaFiles;
    std::vector<GeneticMaterial::FileContent> fastqFiles;

    std::istringstream contentStream(decompressedContent.str());
    std::string line;
    GeneticMaterial::FileContent currentFile;

    while (std::getline(contentStream, line)) {
        if (line.empty()) continue;

        if (startsWith(line, ">")) {
            if (currentFile.content.has_value()) {
                fastaFiles.push_back(currentFile);
            }
            currentFile = GeneticMaterial::FileContent();
            currentFile.fileType = "FASTA";
            currentFile.fileName = filePath + "_decompressed.fasta";
            currentFile.content = line + "\n";
        } else if (startsWith(line, "@")) {
            if (currentFile.content.has_value()) {
                fastqFiles.push_back(currentFile);
            }
            currentFile = GeneticMaterial::FileContent();
            currentFile.fileType = "FASTQ";
            currentFile.fileName = filePath + "_decompressed.fastq";
            currentFile.content = line + "\n";
        } else {
            if (!currentFile.content.has_value()) {
                currentFile.content = "";
            }
            currentFile.content.value() += line + "\n";
        }
    }

    if (currentFile.content.has_value()) {
        if (currentFile.fileType == "FASTA") {
            fastaFiles.push_back(currentFile);
        } else if (currentFile.fileType == "FASTQ") {
            fastqFiles.push_back(currentFile);
        }
    }

    std::cout << "Processing completed. Returning decompressed files." << std::endl;
    return {fastaFiles, fastqFiles};
}








std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> 
GeneticMaterial::descompressTarFile(const std::string& filePath) {
    struct archive* tarFile = archive_read_new();
    if (!tarFile) {
        throw std::runtime_error("Failed to allocate memory for the TAR handler.");
    }

    struct archive_entry* entry;
    std::vector<GeneticMaterial::FileContent> fastaFiles;
    std::vector<GeneticMaterial::FileContent> fastqFiles;

    // Configuring libarchive to handle TAR files
    archive_read_support_filter_all(tarFile);
    archive_read_support_format_all(tarFile);

    if (archive_read_open_filename(tarFile, filePath.c_str(), 10240) != ARCHIVE_OK) {
        archive_read_free(tarFile);
        throw std::runtime_error("Could not open the TAR file: " + filePath);
    }

    while (archive_read_next_header(tarFile, &entry) == ARCHIVE_OK) {
        std::string currentFileName = archive_entry_pathname(entry);
        std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);

        if (endsWith(currentFileName, ".fasta") || endsWith(currentFileName, ".fastq")) {
            std::ostringstream fileContent;
            const void* buffer;
            size_t size;
            int64_t offset;

            while (true) {
                int readBytes = archive_read_data_block(tarFile, &buffer, &size, &offset);
                if (readBytes == ARCHIVE_EOF) break;
                if (readBytes < ARCHIVE_OK) {
                    archive_read_free(tarFile);
                    throw std::runtime_error("Error reading the file: " + currentFileName);
                }
                fileContent.write(static_cast<const char*>(buffer), size);
            }

            if (!fileContent.str().empty()) {
                GeneticMaterial::FileContent fileData;
                fileData.fileName = currentFileName;
                fileData.fileType = endsWith(currentFileName, ".fasta") ? "FASTA" : "FASTQ";
                fileData.content = fileContent.str();

                if (fileData.fileType == "FASTA") {
                    fastaFiles.push_back(std::move(fileData));
                } else {
                    fastqFiles.push_back(std::move(fileData));
                }
            }
        }
        archive_read_data_skip(tarFile);
    }

    archive_read_free(tarFile);
    return {std::move(fastaFiles), std::move(fastqFiles)};
}








std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> 
GeneticMaterial::decompressZipFile(const std::string& filePath) {
    unzFile zipFile = unzOpen(filePath.c_str());
    if (!zipFile) {
        throw std::runtime_error("Could not open the ZIP file: " + filePath);
    }

    if (unzGoToFirstFile(zipFile) != UNZ_OK) {
        unzClose(zipFile);
        throw std::runtime_error("Error accessing the first file in the ZIP: " + filePath);
    }

    std::vector<GeneticMaterial::FileContent> fastaFiles;
    std::vector<GeneticMaterial::FileContent> fastqFiles;

    do {
        char fileName[256];
        unz_file_info fileInfo;

        if (unzGetCurrentFileInfo(zipFile, &fileInfo, fileName, sizeof(fileName), nullptr, 0, nullptr, 0) != UNZ_OK) {
            unzClose(zipFile);
            throw std::runtime_error("Error retrieving file information in the ZIP.");
        }

        std::string currentFileName(fileName);
        std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);

        if (endsWith(currentFileName, ".fasta") || endsWith(currentFileName, ".fastq")) {
            if (unzOpenCurrentFile(zipFile) != UNZ_OK) {
                unzClose(zipFile);
                throw std::runtime_error("Error opening file inside the ZIP: " + currentFileName);
            }

            std::ostringstream fileContent;
            char buffer[4096];
            int bytesRead;

            while ((bytesRead = unzReadCurrentFile(zipFile, buffer, sizeof(buffer))) > 0) {
                fileContent.write(buffer, bytesRead);
            }

            if (bytesRead < 0) {
                unzCloseCurrentFile(zipFile);
                unzClose(zipFile);
                throw std::runtime_error("Error reading the file inside the ZIP: " + currentFileName);
            }

            unzCloseCurrentFile(zipFile);

            if (!fileContent.str().empty()) {
                GeneticMaterial::FileContent fileData;
                fileData.fileName = currentFileName;
                fileData.fileType = endsWith(currentFileName, ".fasta") ? "FASTA" : "FASTQ";
                fileData.content = fileContent.str();

                if (fileData.fileType == "FASTA") {
                    fastaFiles.push_back(std::move(fileData));
                } else {
                    fastqFiles.push_back(std::move(fileData));
                }
            }
        }
    } while (unzGoToNextFile(zipFile) == UNZ_OK);

    unzClose(zipFile);
    return {std::move(fastaFiles), std::move(fastqFiles)};
}






///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Lecture des données ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> GeneticMaterial::openFile(const std::string& filePath) {
//     std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> fileData;

//     size_t position = filePath.find_last_of('.');
//     if (position == std::string::npos) {
//         throw std::runtime_error("Arquivo sem extensão, insira um arquivo com uma extensão válida: FASTA, FASTQ, ZIP ou TAR.");
//     }

//     std::string fileExtension = filePath.substr(position + 1);
//     std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

//     try {
//         if (fileExtension == "gz") {
//             if (filePath.find(".fasta.gz") != std::string::npos || filePath.find(".fastq.gz") != std::string::npos) {
//                 fileData = decompressGzipFile(filePath);
//             } else {
//                 throw std::runtime_error("Extensão não suportada para arquivos GZIP: " + filePath);
//             }
//         } else if (fileExtension == "tar") {
//             fileData = descompressTarFile(filePath);
//         } else if (fileExtension == "zip") {
//             fileData = decompressZipFile(filePath);
//         } else if (fileExtension == "fasta" || fileExtension == "fastq") {
//             // Chama a função loadFileContent e verifica se o retorno é válido
//             auto result = loadFileContent(filePath);
//             if (result) {
//                 const auto& [fastaFiles, fastqFiles] = *result;
//                 // Verifica se as listas não são nullopt e atribui às listas de fileData
//                 if (fastaFiles) {
//                     fileData.first = *fastaFiles;
//                 }
//                 if (fastqFiles) {
//                     fileData.second = *fastqFiles;
//                 }
//         } else {
//             throw std::runtime_error("Falha ao carregar conteúdo do arquivo: " + filePath);
//         }
//     } else {
//         throw std::runtime_error("Extensão não suportada: " + fileExtension +
//             "\nPor favor, insira arquivos do tipo: FASTA, FASTQ, ZIP ou TAR.");
//         }
//     } catch (const std::exception& e) {
//         throw std::runtime_error("Erro ao processar o arquivo: " + std::string(e.what()));
//     }

//     return fileData;
// }
std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> GeneticMaterial::openFile(const std::string& filePath) {
    std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> fileData;

    size_t position = filePath.find_last_of('.');
    if (position == std::string::npos) {
        throw std::runtime_error("File without extension, please provide a file with a valid extension: FASTA, FASTQ, ZIP, or TAR.");
    }

    std::string fileExtension = filePath.substr(position + 1);
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

    try {
        if (fileExtension == "gz") {
            if (filePath.find(".fasta.gz") != std::string::npos || filePath.find(".fastq.gz") != std::string::npos) {
                fileData = decompressGzipFile(filePath);
            } else {
                throw std::runtime_error("Unsupported extension for GZIP files: " + filePath);
            }
        } else if (fileExtension == "tar") {
            fileData = descompressTarFile(filePath);
        } else if (fileExtension == "zip") {
            fileData = decompressZipFile(filePath);
        } else if (fileExtension == "fasta" || fileExtension == "fastq") {
            // Calls the loadFileContent function and checks if the return is valid
            auto result = loadFileContent(filePath);
            if (result) {
                const auto& [fastaFiles, fastqFiles] = *result;
                // Checks if the lists are not nullopt and assigns them to fileData lists
                if (fastaFiles) {
                    fileData.first = *fastaFiles;
                }
                if (fastqFiles) {
                    fileData.second = *fastqFiles;
                }
            } else {
                throw std::runtime_error("Failed to load file content: " + filePath);
            }
        } else {
            throw std::runtime_error("Unsupported extension: " + fileExtension +
                "\nPlease provide files of type: FASTA, FASTQ, ZIP, or TAR.");
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Error processing the file: " + std::string(e.what()));
    }

    return fileData;
}




std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> 
GeneticMaterial::loadFileContent(const std::variant<std::string, std::vector<std::string>>& input) {
    std::vector<std::string> filePaths;

    // Checks if the input is a string or a vector of strings
    if (std::holds_alternative<std::string>(input)) {
        filePaths.push_back(std::get<std::string>(input)); // Adds the single path to the list
    } else {
        filePaths = std::get<std::vector<std::string>>(input); // Directly converts to the list
    }

    // Vectors for FASTA and FASTQ files
    std::optional<std::vector<GeneticMaterial::FileContent>> fastaFiles;
    std::optional<std::vector<GeneticMaterial::FileContent>> fastqFiles;

    // Processes each file in the list
    for (const auto& filePath : filePaths) {

        std::ifstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open the file: " << filePath << std::endl;
            throw std::runtime_error("Could not open the file: " + filePath);
        }

        // Process the file in chunks
        constexpr size_t bufferSize = 4096; // 4 KB
        char buffer[bufferSize];
        std::stringstream fileContent;

        while (file.read(buffer, bufferSize)) {
            fileContent.write(buffer, file.gcount());
        }

        // Checks if there are remaining bytes at the end of the file
        if (file.gcount() > 0) {
            fileContent.write(buffer, file.gcount());
        }

        file.close();

        std::cout << "File successfully read: " << filePath << std::endl;

        // Normalizes the file content
        std::string rawContent = fileContent.str();

        rawContent.erase(
            std::remove_if(rawContent.begin(), rawContent.end(),
                [](char c) { 
                    return !std::isprint(static_cast<unsigned char>(c)) && 
                        !std::isspace(static_cast<unsigned char>(c)); 
                }),
            rawContent.end()
        );

        std::string normalizedContent;
        normalizedContent.reserve(rawContent.size());
        for (size_t i = 0; i < rawContent.size(); ++i) {
            if (rawContent[i] == '\r') {
                if (i + 1 < rawContent.size() && rawContent[i + 1] == '\n') {
                    normalizedContent += '\n';
                    ++i;
                } else {
                    normalizedContent += '\n';
                }
            } else {
                normalizedContent += rawContent[i];
            }
        }

        // Creates the FileContent object with the normalized content
        GeneticMaterial::FileContent fileData;
        fileData.fileName = filePath;
        fileData.content = normalizedContent;

        // Checks the file type and adds it to the appropriate vector
        if (endsWith(filePath, ".fasta")) {
            std::cout << "File identified as FASTA: " << filePath << std::endl;
            if (!fastaFiles) {
                fastaFiles = std::vector<GeneticMaterial::FileContent>();
            }
            fileData.fileType = "FASTA";

            if (!verifySequence(fileData.content.value(), "FASTA")) {
                std::cerr << "Error: Invalid FASTA file: " << filePath << std::endl;
                continue;
            }

            fastaFiles->push_back(fileData);
        } else if (endsWith(filePath, ".fastq")) {
            std::cout << "File identified as FASTQ: " << filePath << std::endl;
            if (!fastqFiles) {
                fastqFiles = std::vector<GeneticMaterial::FileContent>();
            }
            fileData.fileType = "FASTQ";

            if (!verifySequence(fileData.content.value(), "FASTQ")) {
                std::cerr << "Error: Invalid FASTQ file: " << filePath << std::endl;
                continue;
            }

            fastqFiles->push_back(fileData);
        } else {
            std::cerr << "Error: Unsupported file type: " << filePath << std::endl;
            throw std::runtime_error("Unsupported file type: " + filePath);
        }
    }

    return std::make_optional(std::make_pair(std::move(fastaFiles), std::move(fastqFiles)));
}






// std::string normalizeNewlines(const std::string& input) {
//     std::string output;
//     output.reserve(input.size());

//     for (size_t i = 0; i < input.size(); ++i) {
//         if (input[i] == '\r') {
//             // Ignora \r seguido de \n (Windows) ou mantém \r sozinho (Mac antigo)
//             if (i + 1 < input.size() && input[i+1] == '\n') {
//                 output += '\n';
//                 ++i; // Pula o próximo caractere
//             } else {
//                 output += '\n'; // Converte \r solitário para \n
//             }
//         } else {
//             output += input[i];
//         }
//     }
//     return output;
// }

std::string GeneticMaterial::normalizeNewlines(const std::string& input) {
    std::string output;
    output.reserve(input.size());

    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] == '\r') {
            // Ignores \r followed by \n (Windows) or keeps \r alone (old Mac)
            if (i + 1 < input.size() && input[i + 1] == '\n') {
                output += '\n';
                ++i; // Skips the next character
            } else {
                output += '\n'; // Converts lone \r to \n
            }
        } else {
            output += input[i];
        }
    }
    return output;
}










bool GeneticMaterial::verifySequence(const std::string& content, const std::string& fileType) {
    std::cout << "Starting validation of file type: " << fileType << std::endl;

    bool errorFound = false;
    std::string line;
    bool isHeader = true;
    std::string sequence;
    std::string quality;
    int lineNumber = 0;

    std::istringstream stream(content);

    while (std::getline(stream, line)) {
        if (line.empty()) continue;

        if (fileType == "FASTA") {
            if (line[0] == '>') {
                isHeader = true;
            } else {
                // Check if the line contains only valid nucleotides
                for (char c : line) {
                    if (!isValidNucleotide(c)) {
                        std::cerr << "Error: Invalid character in FASTA sequence: " << c << std::endl;
                        errorFound = true;
                    }
                }
                isHeader = false;
            }
        } else if (fileType == "FASTQ") {
            switch (lineNumber % 4) {
                case 0:
                    if (line[0] != '@') {
                        std::cerr << "Error: Invalid FASTQ header. Must start with '@'." << std::endl;
                        return false;
                    }
                    break;
                case 1:
                    sequence = line;
                    for (char c : sequence) {
                        if (!isValidNucleotide(c)) {
                            std::cerr << "Error: Invalid character in FASTQ sequence: " << c << std::endl;
                            errorFound = true;
                        }
                    }
                    break;
                case 3:
                    quality = line;
                    if (quality.length() != sequence.length()) {
                        std::cerr << "Error: FASTQ quality line must have the same length as the sequence." << std::endl;
                        return false;
                    }
                    for (char c : quality) {
                        if (!isValidQualityChar(c)) {
                            std::cerr << "Error: Invalid character in FASTQ quality line: " << c << std::endl;
                            errorFound = true;
                        }
                    }
                    break;
            }
            lineNumber++;
        }
    }

    std::cout << "Validation completed. Errors found: " << (errorFound ? "Yes" : "No") << std::endl;
    return !errorFound;
}




// Auxiliary method to check if a character is a valid nucleotide
bool GeneticMaterial::isValidNucleotide(char c) {
    c = toupper(c);
    bool isValid = (c == 'A' || c == 'T' || c == 'C' || c == 'G' || c == 'U' || 
                    c == 'N' || c == 'R' || c == 'Y' || c == 'W' || c == 'S' || 
                    c == 'K' || c == 'M' || c == 'B' || c == 'D' || c == 'H' || c == 'V');
    if (!isValid) {
        std::cerr << "Invalid character found: " << c << std::endl;
    }
    return isValid;
}

// Auxiliary method to check if a character is valid in a FASTQ quality line
bool GeneticMaterial::isValidQualityChar(char c) {
    // Quality characters are usually in the ASCII range from '!' (33) to '~' (126)
    return (c >= '!' && c <= '~');
}




std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> 
GeneticMaterial::getLoadedFileData() const {
    std::vector<FileContent> fastaFiles;
    std::vector<FileContent> fastqFiles;

    for (const auto& file : m_loadedFiles) {
        if (file.fileType == "FASTA") {
            fastaFiles.push_back(file);
        } else if (file.fileType == "FASTQ") {
            fastqFiles.push_back(file);
        }
    }

    return {fastaFiles, fastqFiles};
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Traitement des données ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


// Function implementation
std::vector<GeneticMaterial::SequenceDescription> GeneticMaterial::getSequenceDescriptions(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
                                  std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles) {
    std::vector<GeneticMaterial::SequenceDescription> descriptions;

    if (!loadedFiles) {
        return descriptions;
    }

    auto extractSpeciesName = [](const std::string& header) -> std::optional<std::string> {
        // Pattern 1: >gi|number|ref|ID|Species_Name|...
        std::regex pattern1(R"((?:.*\|){3}([^|]+)\|)");
        // Pattern 2: >ID Species_Name additional_information
        std::regex pattern2(R"(^\S+\s+([A-Za-z]+_[a-z]+)\b)");
        
        std::smatch matches;
        
        if (std::regex_search(header, matches, pattern1) && matches.size() > 1) {
            return matches[1].str();
        } else if (std::regex_search(header, matches, pattern2) && matches.size() > 1) {
            return matches[1].str();
        }
        
        return std::nullopt;
    };

    auto [fastaFiles, fastqFiles] = *loadedFiles;

    // Open the file to save descriptions
    std::ofstream outFile("descriptions.txt");
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not create the file descriptions.txt");
    }

    // Process FASTA files
    if (fastaFiles) {
        for (const auto& file : *fastaFiles) {
            if (!file.content) continue;

            std::stringstream ss(*file.content);
            std::string line;
            std::string currentHeader;
            std::string currentSequence;

            while (std::getline(ss, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!currentHeader.empty()) {
                        GeneticMaterial::SequenceDescription desc;
                        desc.fileType = "FASTA";
                        desc.sequenceId = currentHeader.substr(0, currentHeader.find(' '));
                        desc.speciesName = extractSpeciesName(currentHeader);
                        desc.sequenceLength = currentSequence.length();
                        
                        if (currentHeader.find(' ') != std::string::npos) {
                            desc.additionalInfo = currentHeader.substr(currentHeader.find(' ') + 1);
                        }
                        
                        descriptions.push_back(desc);

                        // Save the description to the file
                        outFile << "FileType: " << desc.fileType << "\n"
                            << "SequenceID: " << desc.sequenceId << "\n"
                            << "SpeciesName: " << (desc.speciesName ? *desc.speciesName : "N/A") << "\n"
                            << "SequenceLength: " 
                            << (desc.sequenceLength.has_value() ? std::to_string(desc.sequenceLength.value()) : "N/A") << "\n"
                            << "AdditionalInfo: " << (desc.additionalInfo ? *desc.additionalInfo : "N/A") << "\n\n";

                        currentSequence.clear();
                    }
                    currentHeader = line.substr(1);
                } else {
                    currentSequence += line;
                }
            }

            if (!currentHeader.empty()) {
                GeneticMaterial::SequenceDescription desc;
                desc.fileType = "FASTA";
                desc.sequenceId = currentHeader.substr(0, currentHeader.find(' '));
                desc.speciesName = extractSpeciesName(currentHeader);
                desc.sequenceLength = currentSequence.length();
                
                if (currentHeader.find(' ') != std::string::npos) {
                    desc.additionalInfo = currentHeader.substr(currentHeader.find(' ') + 1);
                }
                
                descriptions.push_back(desc);

                // Save the description to the file
                outFile << "FileType: " << desc.fileType << "\n"
                    << "SequenceID: " << desc.sequenceId << "\n"
                    << "SpeciesName: " << (desc.speciesName ? *desc.speciesName : "N/A") << "\n"
                    << "SequenceLength: " 
                    << (desc.sequenceLength.has_value() ? std::to_string(desc.sequenceLength.value()) : "N/A") << "\n"
                    << "AdditionalInfo: " << (desc.additionalInfo ? *desc.additionalInfo : "N/A") << "\n\n";
            }
        }
    }

    if (fastqFiles) {
        for (const auto& file : *fastqFiles) {
            if (!file.content) continue;

            std::stringstream ss(*file.content);
            std::string line;
            int lineCount = 0;
            std::string currentHeader;
            std::string currentSequence;

            while (std::getline(ss, line)) {
                if (line.empty()) continue;

                lineCount++;
                
                if (lineCount == 1) {
                    if (line[0] == '@') {
                        currentHeader = line.substr(1);
                    }
                } else if (lineCount == 2) {
                    currentSequence = line;
                } else if (lineCount == 4) {
                    GeneticMaterial::SequenceDescription desc;
                    desc.fileType = "FASTQ";
                    desc.sequenceId = currentHeader.substr(0, currentHeader.find(' '));
                    desc.speciesName = extractSpeciesName(currentHeader);
                    desc.sequenceLength = currentSequence.length();
                    
                    if (currentHeader.find(' ') != std::string::npos) {
                        desc.additionalInfo = currentHeader.substr(currentHeader.find(' ') + 1);
                    }
                    
                    descriptions.push_back(desc);

                    // Save the description to the file
                    outFile << "FileType: " << desc.fileType << "\n"
                        << "SequenceID: " << desc.sequenceId << "\n"
                        << "SpeciesName: " << (desc.speciesName ? *desc.speciesName : "N/A") << "\n"
                        << "SequenceLength: " 
                        << (desc.sequenceLength.has_value() ? std::to_string(desc.sequenceLength.value()) : "N/A") << "\n"
                        << "AdditionalInfo: " << (desc.additionalInfo ? *desc.additionalInfo : "N/A") << "\n\n";
                    lineCount = 0;
                    currentHeader.clear();
                    currentSequence.clear();
                }
            }
        }
    }

    outFile.close();
    std::cout << "Descriptions saved in descriptions.txt" << std::endl;

    return descriptions;
}





// Função auxiliar para contar sequências em um conteúdo MULTIFASTA
size_t GeneticMaterial::countSequencesInMultifasta(const std::vector<FileContent>& fastaFiles) const {
    size_t totalSequences = 0;

    for (const auto& file : fastaFiles) {
        if (file.content) {

            std::istringstream stream(*file.content);
            std::string line;

            while (std::getline(stream, line)) {
                if (!line.empty() && line[0] == '>') {
                    std::cout << "Cabeçalho encontrado: " << line << std::endl; // Depuração: Exibe cabeçalhos encontrados
                    totalSequences++;
                }
            }
        } else {
            std::cerr << "Erro: O conteúdo do arquivo " << file.fileName << " está vazio ou não foi carregado." << std::endl;
        }
    }

    std::cout << "Total de sequências encontradas: " << totalSequences << std::endl;
    return totalSequences;
}









// Auxiliary function to get the complement of a nucleotide
char GeneticMaterial::getComplement(char nucleotide) {
    switch (toupper(nucleotide)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A'; // For RNA
        case 'R': return 'Y'; // R = A/G, Y = C/T
        case 'Y': return 'R';
        case 'S': return 'S'; // S = C/G
        case 'W': return 'W'; // W = A/T
        case 'K': return 'M'; // K = G/T, M = A/C
        case 'M': return 'K';
        case 'B': return 'V'; // B = C/G/T, V = A/C/G
        case 'V': return 'B';
        case 'D': return 'H'; // D = A/G/T, H = A/C/T
        case 'H': return 'D';
        case 'N': return 'N'; // N = any base
        default: return nucleotide; // Keeps unrecognized characters
    }
}

// Function to calculate the complementary sequence and save it to a file
std::string GeneticMaterial::getComplementSequenceAndSave(const std::string& sequence, const std::string& outputFilePath) {
    // Calculate the complementary sequence
    std::string complement;
    complement.reserve(sequence.size());
    
    for (char c : sequence) {
        complement.push_back(GeneticMaterial::getComplement(c));
    }

    // Open the file for writing
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not create the file: " + outputFilePath);
    }

    // Write the complementary sequence to the file
    outFile << complement;

    // Close the file
    outFile.close();

    std::cout << "Complementary sequence saved in: " << outputFilePath << std::endl;

    // Return the complementary sequence
    return complement;
}




std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> 
GeneticMaterial::calculateComplementSequences(
    const std::optional<std::pair<std::optional<std::vector<FileContent>>, std::optional<std::vector<FileContent>>>>& loadedFiles) {
    if (!loadedFiles.has_value()) {
        throw std::runtime_error("GM: No data loaded to calculate complementary sequences.");
    }

    const auto& [fastaFiles, fastqFiles] = *loadedFiles;

    std::optional<std::vector<FileContent>> complementFastaFiles;
    std::optional<std::vector<FileContent>> complementFastqFiles;

    if (fastaFiles) {
        complementFastaFiles = std::vector<FileContent>();
        for (const auto& file : *fastaFiles) {
            // Save the complement to a file
            std::string outputFilePath = file.fileName + "_complement.fasta";

            std::string complement = getComplementSequenceAndSave(file.content.value(), outputFilePath);

            std::ofstream outFile(outputFilePath);
            if (!outFile.is_open()) {
                throw std::runtime_error("Error opening the file to save: " + outputFilePath);
            }
            outFile << ">Complement of " << file.fileName << "\n";
            outFile << complement << "\n";
            outFile.close();

            std::cout << "Complementary sequence saved in: " << outputFilePath << std::endl;

            // Add the complement to the vector
            complementFastaFiles->push_back({outputFilePath, "FASTA", complement});
        }
    }

    if (fastqFiles) {
        complementFastqFiles = std::vector<FileContent>();
        for (const auto& file : *fastqFiles) {
            // Save the complement to a file
            std::string outputFilePath = file.fileName + "_complement.fastq";
            std::string complement = getComplementSequenceAndSave(file.content.value(), outputFilePath);

            std::ofstream outFile(outputFilePath);
            if (!outFile.is_open()) {
                throw std::runtime_error("Error opening the file to save: " + outputFilePath);
            }
            outFile << "@Complement of " << file.fileName << "\n";
            outFile << complement << "\n";
            outFile.close();

            std::cout << "Complementary sequence saved in: " << outputFilePath << std::endl;

            // Add the complement to the vector
            complementFastqFiles->push_back({outputFilePath, "FASTQ", complement});
        }
    }

    return std::make_optional(std::make_pair(std::move(complementFastaFiles), std::move(complementFastqFiles)));
}


// Function that calculates the reverse complement
std::string GeneticMaterial::getReverseComplement(const std::string& sequence) {
    std::string complement;
    complement.reserve(sequence.size());
    
    // First calculate the complement
    for (char c : sequence) {
        complement.push_back(GeneticMaterial::getComplement(c));
    }
    
    // Then reverse the sequence
    std::reverse(complement.begin(), complement.end());
    
    return complement;
}



// Updated main function to calculate the reverse complement
std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> 
GeneticMaterial::calculateReverseComplementSequences(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles) {
    if (!loadedFiles) {
        return std::nullopt;
    }

    auto [fastaFiles, fastqFiles] = *loadedFiles;

    // Process FASTA files
    if (fastaFiles) {
        for (auto& file : *fastaFiles) {
            if (file.content) {
                std::string& content = *file.content;
                std::stringstream ss(content);
                std::string line;
                std::stringstream newContent;

                // Output file path
                std::string outputFilePath = file.fileName + "_reverse_complement.fasta";
                std::ofstream outFile(outputFilePath);
                if (!outFile.is_open()) {
                    throw std::runtime_error("Could not create the file: " + outputFilePath);
                }

                while (std::getline(ss, line)) {
                    if (line.empty()) continue;
                    
                    if (line[0] == '>') {
                        // Keep the header but add a note
                        newContent << line << " [reverse complement]" << '\n';
                        outFile << line << " [reverse complement]" << '\n';
                    } else {
                        // Process the sequence to get the reverse complement
                        std::string revComp = GeneticMaterial::getReverseComplement(line);
                        newContent << revComp << '\n';
                        outFile << revComp << '\n';
                    }
                }

                file.content = newContent.str();
                outFile.close();
                std::cout << "Reverse complement sequence saved in: " << outputFilePath << std::endl;
            }
        }
    }

    // Process FASTQ files
    if (fastqFiles) {
        for (auto& file : *fastqFiles) {
            if (file.content) {
                std::string& content = *file.content;
                std::stringstream ss(content);
                std::string line;
                std::stringstream newContent;
                int lineCount = 0;
                std::string qualityLine;

                // Output file path
                std::string outputFilePath = file.fileName + "_reverse_complement.fastq";
                std::ofstream outFile(outputFilePath);
                if (!outFile.is_open()) {
                    throw std::runtime_error("Could not create the file: " + outputFilePath);
                }

                while (std::getline(ss, line)) {
                    if (line.empty()) continue;
                    
                    lineCount++;
                    
                    if (lineCount == 2) {
                        // Sequence line - calculate reverse complement
                        std::string revComp = GeneticMaterial::getReverseComplement(line);
                        newContent << revComp << '\n';
                        outFile << revComp << '\n';
                    } else if (lineCount == 4) {
                        // Quality line - needs to be reversed as well
                        std::reverse(qualityLine.begin(), qualityLine.end());
                        newContent << qualityLine << '\n';
                        outFile << qualityLine << '\n';
                        lineCount = 0; // Reset for the next record
                    } else {
                        // Keep other lines (header, + line)
                        if (lineCount == 1) {
                            // Add note to the header
                            newContent << line << " [reverse complement]" << '\n';
                            outFile << line << " [reverse complement]" << '\n';
                        } else {
                            newContent << line << '\n';
                            outFile << line << '\n';
                        }
                    }
                }

                file.content = newContent.str();
                outFile.close();
                std::cout << "Reverse complement sequence saved in: " << outputFilePath << std::endl;
            }
        }
    }

    return std::make_optional(std::make_pair(std::move(fastaFiles), std::move(fastqFiles)));
}




std::vector<std::string> GeneticMaterial::getSequenceBetweenStartAndStop(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
                                  std::optional<std::vector<GeneticMaterial::FileContent>>>>& fileData,
    const std::string& startCodon,
    const std::vector<std::string>& stopCodons) 
{
    std::vector<std::string> extractedSequences;

    if (!fileData || !fileData->first) {
        throw std::runtime_error("No valid FASTA file loaded.");
    }

    const auto& fastaFiles = fileData->first.value();

    for (const auto& file : fastaFiles) {
        if (!file.content.has_value()) continue;

        std::string raw = file.content.value();
        std::stringstream ss(raw);
        std::string line, sequence, header;

        // Remove headers (lines starting with '>') and store the header
        while (std::getline(ss, line)) {
            if (!line.empty() && line[0] == '>') {
                header = line; // Save the original header
            } else if (!line.empty()) {
                sequence += line;
            }
        }

        // Search for start and stop codons
        size_t startPos = sequence.find(startCodon);
        if (startPos == std::string::npos) {
            continue;
        }

        startPos += startCodon.length();
        size_t stopPos = std::string::npos;
        for (const auto& stopCodon : stopCodons) {
            size_t pos = sequence.find(stopCodon, startPos);
            if (pos != std::string::npos && (stopPos == std::string::npos || pos < stopPos)) {
                stopPos = pos;
            }
        }

        if (stopPos == std::string::npos) {
            continue;
        }

        std::string codingRegion = sequence.substr(startPos, stopPos - startPos);
        extractedSequences.push_back(codingRegion);

        // Save the extracted sequence to a FASTA file
        std::string outputFileName = file.fileName + "_extracted.fasta";
        std::ofstream outFile(outputFileName, std::ios::app);
        if (!outFile.is_open()) {
            throw std::runtime_error("Could not create the file: " + outputFileName);
        }
        outFile << header << " [extracted between " << startCodon << " and stop codons]\n";
        outFile << codingRegion << "\n";
        outFile.close();
    }

    return extractedSequences;
}










// Extract a sous-sequence
std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>
GeneticMaterial::extractSubsequences(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
    size_t start, size_t length,
    const std::optional<std::string>& sequenceId) 
{
    if (!loadedFiles) {
        return std::nullopt;
    }

    auto [fastaFiles, fastqFiles] = *loadedFiles;

    auto processSequence = [&](const std::string& content, const std::string& fileType, const std::string& fileName) -> std::string {
        std::stringstream ss(content);
        std::string line;
        std::stringstream newContent;
        bool inTargetSequence = false;
        bool sequenceFound = false;
        std::string currentSequence;
        std::string header;

        // Caminho do arquivo de saída
        std::string outputFileName = fileName + "_subsequences." + (fileType == "FASTA" ? "fasta" : "fastq");
        std::ofstream outFile(outputFileName, std::ios::app);
        if (!outFile.is_open()) {
            throw std::runtime_error("Could not create the file " + outputFileName);
        }

        while (std::getline(ss, line)) {
            if (line.empty()) continue;

            if (fileType == "FASTA" && line[0] == '>') {
                if (!currentSequence.empty()) {
                    std::string subseq = GeneticMaterial::extractFromSequence(currentSequence, start, length);
                    newContent << subseq << "\n";
                    outFile << header << " [subsequence " << start << "-" << (start + length) << "]\n";
                    outFile << subseq << "\n";
                    currentSequence.clear();
                }

                inTargetSequence = (!sequenceId || line.find(*sequenceId) != std::string::npos);
                header = line; // Salva o cabeçalho original
                newContent << line << "\n";
                sequenceFound = sequenceFound || inTargetSequence;
            } else if (fileType == "FASTQ" && line[0] == '@') {
                // Lógica similar para FASTQ
            } else if (inTargetSequence) {
                currentSequence += line;
            }
        }

        if (!currentSequence.empty()) {
            std::string subseq = GeneticMaterial::extractFromSequence(currentSequence, start, length);
            newContent << subseq << "\n";
            outFile << header << " [subsequence " << start << "-" << (start + length) << "]\n";
            outFile << subseq << "\n";
        }

        if (sequenceId && !sequenceFound) {
            throw std::runtime_error("Sequence ID not found: " + *sequenceId);
        }

        outFile.close();
        return newContent.str();
    };

    // Processa arquivos FASTA
    if (fastaFiles) {
        for (auto& file : *fastaFiles) {
            if (file.content) {
                file.content = processSequence(*file.content, "FASTA", file.fileName);
            }
        }
    }

    // Processa arquivos FASTQ
    if (fastqFiles) {
        for (auto& file : *fastqFiles) {
            if (file.content) {
                file.content = processSequence(*file.content, "FASTQ", file.fileName);
            }
        }
    }

    return std::make_optional(std::make_pair(std::move(fastaFiles), std::move(fastqFiles)));
}






std::string GeneticMaterial::extractFromSequence(const std::string& sequence, size_t start, size_t length) {
    if (start >= sequence.size()) {
        throw std::out_of_range("Start position is out of sequence bounds");
    }

    size_t end = std::min(start + length, sequence.size());
    return sequence.substr(start, end - start);
}










double GeneticMaterial::calculateSequenceGC(const std::string& sequence, 
    const std::optional<std::pair<size_t, size_t>>& range) 
{
    size_t start = 0;
    size_t end = sequence.length();
    
    if (range) {
        start = std::min(range->first, end);
        end = std::min(range->second, end);
        if (start >= end) return 0.0;
    }

    size_t gcCount = 0;
    size_t validBases = 0;
    
    for (size_t i = start; i < end; ++i) {
        char c = toupper(sequence[i]);
        if (c == 'G' || c == 'C') {
            gcCount++;
            validBases++;
        } else if (c == 'A' || c == 'T' || c == 'U' || c == 'N') {
            validBases++;
        }
        // Ignores others characteres and spaces
    }

    return (validBases > 0) ? (static_cast<double>(gcCount) / validBases * 100.0) : 0.0;
}










std::optional<std::vector<std::pair<std::string, double>>> 
GeneticMaterial::calculateGCContent(
    const std::optional<std::pair<std::optional<std::vector<FileContent>>,
    std::optional<std::vector<FileContent>>>>& loadedFiles,
    std::optional<std::pair<size_t, size_t>> range)
{
    std::vector<std::pair<std::string, double>> gcContents;

    if (!loadedFiles) {
        return std::nullopt;
    }

    auto [fastaFiles, fastqFiles] = *loadedFiles;

    auto processSequence = [&](const std::string& content, const std::string& fileType) {
        std::stringstream ss(content);
        std::string line;
        std::string currentHeader;
        std::string currentSequence;
        bool inSequence = false;

        while (std::getline(ss, line)) {
            if (line.empty()) continue;

            if ((fileType == "FASTA" && line[0] == '>') || 
                (fileType == "FASTQ" && line[0] == '@')) {
                
                // Process the previous sequence
                if (!currentSequence.empty()) {
                    double gc = GeneticMaterial::calculateSequenceGC(currentSequence, range);

                    // Extract the ID up to the comma
                    std::string trimmedHeader = currentHeader;
                    size_t commaPos = currentHeader.find(',');
                    if (commaPos != std::string::npos) {
                        trimmedHeader = currentHeader.substr(0, commaPos);
                    }

                    gcContents.emplace_back(trimmedHeader, gc);
                    currentSequence.clear();
                }
                
                currentHeader = line.substr(1); // Remove '>' or '@'
                inSequence = (fileType == "FASTA");
                
            } else if (fileType == "FASTQ") {
                // Specific logic for FASTQ (lines 2 and 4)
                static int lineCount = 0;
                lineCount++;
                
                if (lineCount == 2) { // Sequence line
                    currentSequence = line;
                } else if (lineCount == 4) { // End of record
                    double gc = GeneticMaterial::calculateSequenceGC(currentSequence, range);

                    // Extract the ID up to the comma
                    std::string trimmedHeader = currentHeader;
                    size_t commaPos = currentHeader.find(',');
                    if (commaPos != std::string::npos) {
                        trimmedHeader = currentHeader.substr(0, commaPos);
                    }

                    gcContents.emplace_back(trimmedHeader, gc);
                    currentSequence.clear();
                    lineCount = 0;
                }
            } else if (inSequence) {
                currentSequence += line;
            }
        }

        // Process the last sequence in the file
        if (!currentSequence.empty()) {
            double gc = GeneticMaterial::calculateSequenceGC(currentSequence, range);

            // Extract the ID up to the comma
            std::string trimmedHeader = currentHeader;
            size_t commaPos = currentHeader.find(',');
            if (commaPos != std::string::npos) {
                trimmedHeader = currentHeader.substr(0, commaPos);
            }

            gcContents.emplace_back(trimmedHeader, gc);
        }
    };

    // Process FASTA files
    if (fastaFiles) {
        for (const auto& file : *fastaFiles) {
            if (file.content) {
                processSequence(*file.content, "FASTA");
            }
        }
    }

    // Process FASTQ files
    if (fastqFiles) {
        for (const auto& file : *fastqFiles) {
            if (file.content) {
                processSequence(*file.content, "FASTQ");
            }
        }
    }

    return gcContents;
}











std::optional<std::vector<GeneticMaterial::FileContent>> GeneticMaterial::filterFastqByQuality(
    const std::optional<std::vector<GeneticMaterial::FileContent>>& fastqFiles,
    const std::string& qualityStr,
    const std::string& outputDir) {

    if (!fastqFiles || fastqFiles->empty()) {
        throw std::runtime_error("No FASTQ files loaded or invalid data.");
    }

    int minQuality;
    try {
        minQuality = std::stoi(qualityStr);
    } catch (...) {
        throw std::invalid_argument("Invalid minimum quality: " + qualityStr);
    }

    const auto& files = *fastqFiles;
    std::vector<GeneticMaterial::FileContent> filteredFiles;

    for (const auto& fileData : files) {
        if (!fileData.content) continue;

        std::istringstream iss(*fileData.content);
        std::ostringstream oss;
        std::string line1, line2, line3, line4;
        int totalReads = 0, passedReads = 0;

        while (std::getline(iss, line1) && 
               std::getline(iss, line2) && 
               std::getline(iss, line3) && 
               std::getline(iss, line4)) {
            
            totalReads++;
            if (calculateAverageQuality(line4) >= minQuality) {
                passedReads++;
                oss << line1 << '\n' << line2 << '\n' << line3 << '\n' << line4 << '\n';
            }
        }

        if (passedReads > 0) {
            // Extract only the file name (without the path)
            std::string fileName = fileData.fileName;
            size_t lastSlash = fileName.find_last_of("/\\");
            if (lastSlash != std::string::npos) {
                fileName = fileName.substr(lastSlash + 1);
            }

            // Save the filtered file in the current directory
            std::string outputFileName = "filtered_" + fileName;
            std::ofstream outFile(outputFileName);
            if (!outFile.is_open()) {
                throw std::runtime_error("Could not create the file: " + outputFileName);
            }
            outFile << oss.str();
            outFile.close();

            filteredFiles.push_back({
                outputFileName,
                "FASTQ",
                oss.str()
            });
        }

        std::cout << "File: " << fileData.fileName 
                  << " - Total reads: " << totalReads 
                  << ", Passed reads: " << passedReads 
                  << " (" << (totalReads > 0 ? (100 * passedReads / totalReads) : 0) 
                  << "%)\n";
    }

    return filteredFiles.empty() ? std::nullopt : std::optional<std::vector<GeneticMaterial::FileContent>>(filteredFiles);
}







// Calculates the average quality of a FASTQ sequence
double GeneticMaterial::calculateAverageQuality(const std::string& qualityStr) {
    if (qualityStr.empty()) return 0.0;
    double sum = 0.0;
    for (char c : qualityStr) {
        sum += static_cast<int>(c) - 33; // Phred+33
    }
    return sum / qualityStr.length();
}



// Removes a specific prefix from a sequence
void GeneticMaterial::removePrefix(std::string& sequence, std::string& quality, const std::string& prefix) {
    if (sequence.find(prefix) == 0) {
        sequence.erase(0, prefix.length());
        quality.erase(0, prefix.length());
    }
}



// Removes poly-A/poly-T tails from a sequence
void GeneticMaterial::trimPolyAT(std::string& sequence, std::string& quality) {
    while (!sequence.empty() && (sequence.back() == 'A' || sequence.back() == 'T' ||
                                sequence.back() == 'a' || sequence.back() == 't')) {
        sequence.pop_back();
        if (!quality.empty()) quality.pop_back();
    }
}



// Checks if a sequence contains degenerate characters
bool GeneticMaterial::containsDegenerateBases(const std::string& sequence) {
    const std::string validBases = "ACGTacgt";
    return std::any_of(sequence.begin(), sequence.end(),
        [&validBases](char c) { return validBases.find(c) == std::string::npos; });
}



std::string GeneticMaterial::toString(const ProcessedSequence& sequence) {
    std::ostringstream oss;
    oss << "Header: " << sequence.header << "\n"
        << "Sequence: " << sequence.sequence << "\n"
        << "Quality: " << sequence.quality << "\n"
        << "Average Quality: " << sequence.avgQuality;
    return oss.str();
}





std::optional<std::vector<GeneticMaterial::ProcessedSequence>> GeneticMaterial::processAndCleanSequences(
    const std::optional<std::vector<FileContent>>& loadedFiles,
    const std::optional<std::string>& prefixToRemove,
    bool removePolyAT,
    double minQuality,
    bool removeDegenerate) 
{
    std::vector<ProcessedSequence> cleanedSequences;

    if (!loadedFiles || loadedFiles->empty()) {
        std::cout << "No files loaded for processing.\n";
        return std::nullopt;
    }

    for (const auto& fileData : *loadedFiles) {
        if (!fileData.content) {
            std::cout << "File without content: " << fileData.fileName << "\n";
            continue;
        }

        // Generate the output file name
        std::string outputFileName = fileData.fileName + "_processed_cleaned.fastq";

        // Open the output file
        std::ofstream outFile(outputFileName);
        if (!outFile.is_open()) {
            throw std::runtime_error("Could not create the output file: " + outputFileName);
        }

        std::istringstream fileStream(*fileData.content);
        std::string line1, line2, line3, line4;

        while (std::getline(fileStream, line1) && 
               std::getline(fileStream, line2) && 
               std::getline(fileStream, line3) && 
               std::getline(fileStream, line4)) {

            // Create the processed sequence structure
            GeneticMaterial::ProcessedSequence seq{line1, line2, line4, GeneticMaterial::calculateAverageQuality(line4)};

            // Remove the prefix, if specified
            if (prefixToRemove && !prefixToRemove->empty()) {
                std::cout << "Removing prefix: " << *prefixToRemove << "\n";
                GeneticMaterial::removePrefix(seq.sequence, seq.quality, *prefixToRemove);
            }

            // Remove poly-A/poly-T tails, if enabled
            if (removePolyAT) {
                std::cout << "Removing poly-A/poly-T tails\n";
                GeneticMaterial::trimPolyAT(seq.sequence, seq.quality);
            }

            // Skip sequences with average quality below the threshold
            if (seq.avgQuality < minQuality) {
                std::cout << "Sequence skipped due to low quality: " << seq.avgQuality << "\n";
                continue;
            }

            // Skip sequences with degenerate characters, if enabled
            if (removeDegenerate && GeneticMaterial::containsDegenerateBases(seq.sequence)) {
                std::cout << "Sequence skipped due to degenerate bases.\n";
                continue;
            }

            // Add the cleaned sequence to the list, if not empty
            if (!seq.sequence.empty()) {
                cleanedSequences.push_back(seq);

                // Write the cleaned sequence to the output file
                outFile << seq.header << "\n";
                outFile << seq.sequence << "\n";
                outFile << "+\n";
                outFile << seq.quality << "\n";
            }
        }

        outFile.close();

        // Add the generated file name to the output structure
        if (!cleanedSequences.empty()) {
            std::cout << "File generated: " << outputFileName << "\n";
        }
    }

    if (!cleanedSequences.empty()) {
        std::cout << "Processing completed. " << cleanedSequences.size() << " cleaned sequences processed.\n";
        return cleanedSequences;
    } else {
        std::cout << "No sequences met the cleaning criteria.\n";
        return std::nullopt;
    }
}








///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Suffix Table /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<int> GeneticMaterial::buildSuffixArray(const std::string& sequence) {
    if (sequence.empty()) {
        throw std::runtime_error("The sequence is empty. It is not possible to build the suffix array.");
    }

    int n = sequence.size();
    std::vector<int> suffixArray(n);
    std::iota(suffixArray.begin(), suffixArray.end(), 0); // Initializes with 0, 1, ..., n-1

    std::sort(suffixArray.begin(), suffixArray.end(), [&sequence](int i, int j) {
        return sequence.substr(i) < sequence.substr(j);
    });

    return suffixArray;
}

std::vector<int> GeneticMaterial::buildLCPArray(const std::string& sequence, const std::vector<int>& suffixArray) {
    int n = suffixArray.size();
    std::vector<int> lcpArray(n, 0);

    for (int i = 1; i < n; ++i) {
        int lcp = 0;
        int x = suffixArray[i - 1];
        int y = suffixArray[i];

        while (x + lcp < sequence.size() && y + lcp < sequence.size() && sequence[x + lcp] == sequence[y + lcp]) {
            lcp++;
        }

        lcpArray[i] = lcp;
    }

    return lcpArray;
}



int GeneticMaterial::getSequenceLength() const {
    return sequence.size();
}

void GeneticMaterial::setSequence(const std::string& seq) {
    if (seq.empty()) {
        throw std::runtime_error("The provided sequence is empty.");
    }
    sequence = seq;
}

const std::vector<int>& GeneticMaterial::getSuffixArray() const {
    return suffixArray;
}

std::string GeneticMaterial::getFactor(int i, int k) const {
    if (suffixArray.empty()) {
        throw std::runtime_error("Suffix array has not been initialized.");
    }
    if (sequence.empty()) {
        throw std::runtime_error("Sequence has not been initialized.");
    }
    if (i < 0 || i >= suffixArray.size() || suffixArray[i] + k > sequence.size()) {
        throw std::out_of_range("Invalid index or length.");
    }
    return sequence.substr(suffixArray[i], k);
}

std::vector<int> GeneticMaterial::findRepeatedFactors(int k) const {
    if (lcpArray.empty()) {
        throw std::runtime_error("LCP array has not been initialized.");
    }

    std::vector<int> repeatedFactors;
    for (int i = 1; i < lcpArray.size(); ++i) {
        if (lcpArray[i] >= k) {
            repeatedFactors.push_back(i);
        }
    }
    return repeatedFactors;
}

const std::vector<int>& GeneticMaterial::getLCPArray() const {
    return lcpArray;
}



void GeneticMaterial::buildIndex() {
    if (sequence.empty()) {
        throw std::runtime_error("The sequence is empty. It is not possible to build the suffix table.");
    }

    // Build the suffix table
    suffixArray.resize(sequence.size());
    for (size_t i = 0; i < sequence.size(); ++i) {
        suffixArray[i] = i;
    }
    std::sort(suffixArray.begin(), suffixArray.end(), [this](int a, int b) {
        return sequence.substr(a) < sequence.substr(b);
    });

    // Build the LCP (Longest Common Prefix) table
    lcpArray.resize(sequence.size());
    for (size_t i = 1; i < suffixArray.size(); ++i) {
        size_t len = 0;
        while (suffixArray[i] + len < sequence.size() &&
               suffixArray[i - 1] + len < sequence.size() &&
               sequence[suffixArray[i] + len] == sequence[suffixArray[i - 1] + len]) {
            ++len;
        }
        lcpArray[i] = len;
    }
}



std::string GeneticMaterial::extractSequenceFromLoadedData(
    const std::optional<std::pair<
        std::vector<FileContent>,
        std::vector<FileContent>>>& loadedData) const {
    if (!loadedData.has_value()) {
        throw std::runtime_error("No data loaded to extract the sequence.");
    }

    const auto& [fastaFiles, fastqFiles] = *loadedData;

    // Prioritize FASTA files
    if (!fastaFiles.empty()) {
        std::string sequence;
        for (const auto& file : fastaFiles) {
            if (file.content.has_value()) {
                sequence += file.content.value(); // Concatenate all FASTA sequences
            }
        }
        return sequence;
    }

    // If no FASTA files, try FASTQ files
    if (!fastqFiles.empty()) {
        std::string sequence;
        for (const auto& file : fastqFiles) {
            if (file.content.has_value()) {
                sequence += file.content.value(); // Concatenate all FASTQ sequences
            }
        }
        return sequence;
    }

    throw std::runtime_error("No valid sequence found in the loaded data.");
}




void GeneticMaterial::buildAndSaveSuffixTable() {
    // Open the FASTA file
    std::ifstream fastaFile(filePath);
    if (!fastaFile.is_open()) {
        throw std::runtime_error("Error opening the FASTA file: " + filePath);
    }

    // Extract the sequence name (first line of the file)
    std::string sequenceName;
    std::getline(fastaFile, sequenceName);
    if (sequenceName.empty() || sequenceName[0] != '>') {
        throw std::runtime_error("Invalid format: the file does not contain a valid FASTA header.");
    }

    // Remove the '>' character from the beginning of the sequence name
    sequenceName = sequenceName.substr(1);

    // Read the complete sequence (concatenate all subsequent lines)
    std::string sequence;
    std::string line;
    while (std::getline(fastaFile, line)) {
        if (!line.empty() && line[0] != '>') {
            sequence += line;
        }
    }
    fastaFile.close();

    // Build the suffix table
    suffixArray.resize(sequence.size());
    for (size_t i = 0; i < sequence.size(); ++i) {
        suffixArray[i] = i;
    }
    std::sort(suffixArray.begin(), suffixArray.end(), [&sequence](int a, int b) {
        return sequence.substr(a) < sequence.substr(b);
    });

    // Save the suffix table to the file "suffix-table.txt"
    std::ofstream suffixTableFile("suffix-table.txt");
    if (!suffixTableFile.is_open()) {
        throw std::runtime_error("Error opening the file suffix-table.txt for writing.");
    }

    // Write the sequence name at the beginning of the file
    suffixTableFile << ">Sequence: " << sequenceName << "\n";
    suffixTableFile << "Suffix Table:\n";

    // Write the suffix table
    for (size_t i = 0; i < suffixArray.size(); ++i) {
        suffixTableFile << i << ": " << sequence.substr(suffixArray[i]) << "\n";
    }

    suffixTableFile.close();
    std::cout << "Suffix table saved in 'suffix-table.txt'.\n";
}





int GeneticMaterial::searchPattern(const std::string& pattern) {
    if (sequence.empty()) {
        throw std::runtime_error("The sequence is empty. It is not possible to perform the search.");
    }

    if (suffixArray.empty()) {
        throw std::runtime_error("The suffix array has not been built.");
    }

    // Performs binary search using the suffix array
    auto it = std::lower_bound(suffixArray.begin(), suffixArray.end(), pattern,
        [this](int index, const std::string& pat) {
            return sequence.substr(index) < pat;
        });

    if (it != suffixArray.end() && sequence.substr(*it, pattern.size()) == pattern) {
        return *it; // Returns the position of the found pattern
    }

    return -1; // Returns -1 if the pattern is not found
}




std::vector<int> GeneticMaterial::findRepeatedFactors(int minLength) {
    if (sequence.empty()) {
        throw std::runtime_error("The sequence is empty. It is not possible to find repeated patterns.");
    }

    if (suffixArray.empty() || lcpArray.empty()) {
        throw std::runtime_error("The suffix array or the LCP array has not been built.");
    }

    std::vector<int> repeatedFactors;

    // Iterates through the LCP array to find repeated patterns with length >= minLength
    for (size_t i = 1; i < lcpArray.size(); ++i) {
        if (lcpArray[i] >= minLength) {
            repeatedFactors.push_back(suffixArray[i]);
        }
    }

    return repeatedFactors;
}








///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Mapping //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void GeneticMaterial::loadGenomeSequence() {
    // Uses the openFile function to load data from the file
    auto loadedData = openFile(filePath);

    // Extracts the loaded FASTA and FASTQ sequences
    auto [fastaFiles, fastqFiles] = loadedData;

    // Concatenates all FASTA and FASTQ sequences into a single genomic sequence
    genomeSequence.clear();

    for (const auto& file : fastaFiles) {
        if (file.content.has_value()) {
            genomeSequence += file.content.value();
        }
    }

    for (const auto& file : fastqFiles) {
        if (file.content.has_value()) {
            genomeSequence += file.content.value();
        }
    }

    if (genomeSequence.empty()) {
        throw std::runtime_error("No valid sequence found in the loaded files.");
    }

    std::cout << "Genomic sequence successfully loaded. Size: " << genomeSequence.size() << " bases.\n";
}

std::string GeneticMaterial::getReverseComplementKmer(const std::string& kmer) const {
    std::string complement;
    for (char c : kmer) {
        switch (c) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            default: complement += 'N'; break; // N represents unknown nucleotides
        }
    }
    std::reverse(complement.begin(), complement.end());
    return complement;
}

void GeneticMaterial::buildKmerIndex(int k) {
    if (genomeSequence.empty()) {
        throw std::runtime_error("The genomic sequence is empty. Ensure the data has been loaded correctly.");
    }

    kmerIndex.clear();
    for (size_t i = 0; i <= genomeSequence.size() - k; ++i) {
        std::string kmer = genomeSequence.substr(i, k);
        kmerIndex[kmer].push_back(i);
    }

    std::cout << "K-mer index successfully built. Total k-mers: " << kmerIndex.size() << "\n";
}



std::vector<GeneticMaterial::KmerMatch> GeneticMaterial::mapReadToGenome(const std::string& read, int k) const {
    std::vector<KmerMatch> matches;

    for (size_t i = 0; i <= read.size() - k; ++i) {
        std::string kmer = read.substr(i, k);

        // Checks the k-mer in the index
        auto it = kmerIndex.find(kmer);
        if (it != kmerIndex.end()) {
            for (int pos : it->second) {
                matches.push_back({pos, false});
            }
        }

        // Checks the reverse complement
        std::string reverseComplement = getReverseComplementKmer(kmer);
        it = kmerIndex.find(reverseComplement);
        if (it != kmerIndex.end()) {
            for (int pos : it->second) {
                matches.push_back({pos, true});
            }
        }
    }

    return matches;
}



// Function to get the support of a k-mer
int GeneticMaterial::getKmerSupport(const std::string& kmer) const {
    auto it = kmerIndex.find(kmer);
    if (it != kmerIndex.end()) {
        return it->second.size();
    }
    return 0;
}



// Function to diagnose variations in a read
std::optional<std::string> GeneticMaterial::diagnoseVariation(const std::string& read, const std::vector<KmerMatch>& matches, int k) const {
    if (matches.empty()) {
        return std::nullopt;
    }

    // Checks if the k-mers are correctly aligned
    int firstPosition = matches.front().position;
    for (size_t i = 1; i < matches.size(); ++i) {
        int expectedPosition = firstPosition + static_cast<int>(i * k);
        if (matches[i].position != expectedPosition) {
            return "Possible variation detected (insertion, deletion, or substitution).";
        }
    }

    return "Read successfully mapped without variations.";
}



void GeneticMaterial::mapReadsAndReport(const std::vector<std::string>& reads, int k) {
    // Builds the k-mer index
    buildKmerIndex(k);

    for (const auto& read : reads) {
        std::cout << "Mapping read: " << read << "\n";

        // Maps the read to the genome
        auto matches = mapReadToGenome(read, k);

        if (matches.empty()) {
            std::cout << "No mapping found for the read.\n";
            continue;
        }

        // Displays the total number of occurrences of the read
        std::cout << "Total occurrences of the read: " << matches.size() << "\n";

        // Variation diagnosis
        auto diagnosis = diagnoseVariation(read, matches, k);
        if (diagnosis) {
            std::cout << "Diagnosis: " << *diagnosis << "\n";
        }

        // Displays the positions and orientation (direct or reverse complement)
        for (const auto& match : matches) {
            std::cout << "Position: " << match.position
                      << ", Reverse complement: " << (match.isReverseComplement ? "Yes" : "No") << "\n";
        }

        // Displays the support of each k-mer in the read
        std::cout << "Support of k-mers:\n";
        for (size_t i = 0; i <= read.size() - k; ++i) {
            std::string kmer = read.substr(i, k);
            int support = getKmerSupport(kmer);
            std::cout << "  k-mer: " << kmer << ", Support: " << support << "\n";
        }

        std::cout << "----------------------------------------\n";
    }
}



/////////////////////////////////////////////////




// Reads the sequence and returns how many alignments can be performed with it
// Function to calculate the binomial coefficient: C(n, k) = n! / (k! * (n-k)!)
uint64_t binomialCoefficient(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;

    uint64_t result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }
    return result;
}


// Function to calculate the number of possible global alignments
uint64_t calculateGlobalAlignments(int length1, int length2) {
    int totalLength = length1 + length2;
    return binomialCoefficient(totalLength, length1);
}

// Function to calculate the number of possible local alignments
uint64_t calculateLocalAlignments(int length1, int length2) {
    return static_cast<uint64_t>((length1 + 1) * (length2 + 1) * (length1 + length2) / 2);
}

















