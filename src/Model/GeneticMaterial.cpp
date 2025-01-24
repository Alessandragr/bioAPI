#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <zlib.h>
#include <tar.h>

class GeneticMaterial {

    protected:
        std::string filePath;
        std::string fileContent;
        std::string fileExtension;

    public:
        GeneticMaterial(const std::string& path) : filePath(path) {}

        std::string verifyAndLoadFile() {
            size_t position = filePath.find_last_of('.');
                if(position==std::string::npos) {
                    throw std::runtime_error("Arquivo sem extensão, insira um arquivo com uma extensão válida: FASTA, FASTQ, ZIP or TAR.");
                }

            fileExtension = filePath.substr(position + 1);
            std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

            if(fileExtension=="gz") {
                descompressGzipFile();
            } else if(fileExtension=="tar") {
                descompressTarFile();
            } else if(fileExtension=="zip") {
                decompressZipFile();
            } else if (fileExtension=="fasta" || fileExtension=="fastq") {
                loadFileContent();
            } else {
                throw std::runtime_error("Extensão não suportada: " + fileExtension + 
                "\n" + "Please insert files of type: FASTA, FASTQ, ZIP or TAR");
            }
                
        }

        void loadFileContent() {
            std::ifstream file(filePath);
            if(!file.is_open()) {
                throw std::runtime_error("Não foi possível abrir o arquivo: " + filePath);
            }

            std::ostringstream buffer;
            buffer << file.rdbuf();
            fileContent = buffer.str();
            file.close();
        }

        void descompressGzipFile() {
            gzFile gzFile = gzopen(filePath.c_str, "rb");
            if(!gzFile) {
                throw std::runtime_error("Não foi possível abrir o arquivo GZIP: " + filePath);
            }

            char buffer[1024];
            std::ostringstream decompressedContent;

            int bytesRead;
            while((bytesRead = gzread(gzFile, buffer, sizeof(buffer))) > 0) {
                decompressedContent.write(buffer, bytesRead);
            }

            gzclose(gzFile);

            if(decompressedContent.str().empty()) {
                throw std::runtime_error("Falha ao descompactar o arquivo GZIP: " + filePath);
            }

            fileContent = decompressedContent.str();
        }

        void descompressTarFile() {}

        void decompressZipFile() {}




};