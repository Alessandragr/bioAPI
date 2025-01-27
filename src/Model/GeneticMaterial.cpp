#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <optional>
#include <zlib.h>
#include <tar.h>
#include <archive.h>
#include <archive_entry.h>
#include <vector>
#include <unzip.h>

class GeneticMaterial {

    struct FileContent {
    std::string fileName;
    std::optional<std::string> fileType; // Pode ser "FASTA", "FASTQ" ou std::nullopt
    std::optional<std::string> content; // Conteúdo do arquivo, ou std::nullopt se não carregado

    fileData.content = (decompressedContent.str().empty()) ? std::nullopt : std::make_optional(decompressedContent.str());

    }; 

    protected:
        std::string filePath;
        std::string fileContent;
        std::string fileExtension;

    public:
        GeneticMaterial(const std::string& path) : filePath(path) {}



        std::pair<std::vector<FileContent>, std::vector<FileContent>> verifyAndLoadFile(const std::string& filePath) {
            std::pair<std::vector<FileContent>, std::vector<FileContent>> fileData;
            
            size_t position = filePath.find_last_of('.');
                if(position==std::string::npos) {
                    throw std::runtime_error("Arquivo sem extensão, insira um arquivo com uma extensão válida: FASTA, FASTQ, ZIP or TAR.");
                }

            std::string fileExtension = filePath.substr(position + 1);
            std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);
            try {
                if(fileExtension=="gz") {
                    fileData = descompressGzipFile(filePath);
                } else if(fileExtension=="tar") {
                    fileData = descompressTarFile(filePath);
                } else if(fileExtension=="zip") {
                    fileData = decompressZipFile(filePath);
                } else if (fileExtension=="fasta" || fileExtension=="fastq") {
                fileData = loadFileContent(filePath);
                } else {
                    throw std::runtime_error("Extensão não suportada: " + fileExtension + 
                    "\n" + "Please insert files of type: FASTA, FASTQ, ZIP or TAR");
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Erro ao processar o arquivo: " + std::string(e.what()));
            }

            return fileData;
        }




        void loadFileContent(const std::variant<std::string, std::vector<std::string>>& input, 
            std::vector<FileContent>& fastaFiles, 
            std::vector<FileContent>& fastqFiles) {

            std::vector<std::string> filePaths;

            if (std::holds_alternative<std::string>(input)) {
                filePaths.push_back(std::get<std::string>(input)); // Adiciona o único caminho à lista
            } else {
                filePaths = std::get<std::vector<std::string>>(input); // Converte diretamente para a lista
            }

            // Processa cada arquivo na lista
            for (const auto& filePath : filePaths) {
                std::ifstream file(filePath, std::ios::binary);

                if(!file.is_open()) {
                    throw std::runtime_error("Não foi possível abrir o arquivo: " + filePath);
                }

                // Processar o arquivo em blocos
                constexpr size_t bufferSize = 4096; // 4 KB
                char buffer[bufferSize];
                std::stringstream fileContent;

                while (file.read(buffer, bufferSize)) {
                    fileContent.append(buffer, file.gcount());
                }

                // Verifica se há bytes restantes no final do arquivo
                if (file.gcount() > 0) {
                    fileContent.append(buffer, file.gcount());
                }

                file.close();


                FileContent fileData;
                fileData.fileName = filePath;
                fileData.content = fileContent.str();

                if (filePath.ends_with(".fasta")) {
                    fileData.fileType = "FASTA";
                    fastaFiles.push_back(fileData);
                } else if (filePath.ends_with(".fastq")) {
                    fileData.fileType = "FASTQ";
                    fastqFiles.push_back(fileData);
                } else {
                    throw std::runtime_error("Tipo de arquivo não suportado: " + filePath);
                }
            }
        }



        std::vector<FileContent> decompressGzipFile(const std::string& filePath) {
            gzFile gzFile = gzopen(filePath.c_str(), "rb");
            if (!gzFile) {
                throw std::runtime_error("Não foi possível abrir o arquivo GZIP: " + filePath);
            }

            constexpr size_t bufferSize = 8192; // 8 KB
            char buffer[bufferSize];
            std::ostringstream decompressedContent;

            int bytesRead;
            while ((bytesRead = gzread(gzFile, buffer, sizeof(buffer))) > 0) {
                decompressedContent.write(buffer, bytesRead);
            }

            if (bytesRead < 0) { // Verifica erro do gzread
                int errNum;
                const char* errorMsg = gzerror(gzFile, &errNum);
                gzclose(gzFile);
                throw std::runtime_error("Erro ao ler o arquivo GZIP: " + std::string(errorMsg));
            }

            gzclose(gzFile);

            if (decompressedContent.str().empty()) {
                throw std::runtime_error("Falha ao descompactar o arquivo GZIP: " + filePath);
            }

            // Processa o conteúdo para dividir em FASTA/FASTQ
            std::vector<FileContent> fileContents;
            std::string content = decompressedContent.str();

            std::istringstream contentStream(content);
            std::string line;
            FileContent currentFile;

            while (std::getline(contentStream, line)) {
                if (line.starts_with(">")) { // Início de um arquivo FASTA
                    if (!currentFile.fileName.empty()) {
                        fileContents.push_back(currentFile);
                        currentFile = FileContent(); // Reinicia o arquivo atual
                    }
                    currentFile.fileType = "FASTA";
                    currentFile.content += line + "\n";
                } else if (line.starts_with("@")) { // Início de um arquivo FASTQ
                    if (!currentFile.fileName.empty()) {
                        fileContents.push_back(currentFile);
                        currentFile = FileContent(); // Reinicia o arquivo atual
                    }
                    currentFile.fileType = "FASTQ";
                    currentFile.content += line + "\n";
                } else {
                    currentFile.content += line + "\n";
                }
            }

            // Adiciona o último arquivo processado
            if (!currentFile.content.empty()) {
                fileContents.push_back(currentFile);
            }

            return fileContents;
        }

        





        std::pair<std::vector<FileContent>, std::vector<FileContent>> descompressTarFile(const std::string& filePath) {        
        struct archive* tarFile = archive_read_new();
        struct archive_entry* entry;
        std::vector<FileContent> fastaFiles;
        std::vector<FileContent> fastqFiles; 

        // Configurando o libarchive para lidar com arquivos TAR
        archive_read_support_filter_all(tarFile);
        archive_read_support_format_all(tarFile);

        if (archive_read_open_filename(tarFile, filePath.c_str(), 10240) != ARCHIVE_OK) {
            throw std::runtime_error("Não foi possível abrir o arquivo TAR: " + filePath);
        }

        while (archive_read_next_header(tarFile, &entry) == ARCHIVE_OK) {
            std::string currentFileName = archive_entry_pathname(entry);

            // Converte o nome do arquivo para lowercase para evitar problemas de case
            std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);

            // Filtra apenas arquivos .fasta ou .fastq
            if (currentFileName.ends_with(".fasta") || currentFileName.ends_with(".fastq")) {
                std::cout << "Arquivo encontrado: " << currentFileName << std::endl;


                std::string fileType = (currentFileName.ends_with(".fasta")) ? "FASTA" : "FASTQ";

                // Lê o conteúdo do arquivo
                std::stringstream fileContent;
                const void* buffer;
                size_t size;
                int64_t offset;

                while (true) {
                    int readBytes = archive_read_data_block(tarFile, &buffer, &size, &offset);
                    if (readBytes == ARCHIVE_EOF) break;
                    if (readBytes < ARCHIVE_OK) {
                        throw std::runtime_error("Erro ao ler o arquivo: " + currentFileName);
                    }

                    // Adiciona o conteúdo ao stream
                    fileContent.write(static_cast<const char*>(buffer), size);
                }

                 // Cria o objeto FileContent
                FileContent fileData;
                fileData.fileName = currentFileName;
                ffileData.fileType = std::make_optional(fileType);
                fileData.content = fileContent.str();


                if (fileType == "FASTA") {
                    fastaFiles.push_back(fileData);  // Adiciona à lista de FASTA
                } else {
                    fastqFiles.push_back(fileData);  // Adiciona à lista de FASTQ
                }

            } else {
                std::cout << "Arquivo ignorado: " << currentFileName << std::endl;
            }
            // Avança para o próximo arquivo no TAR
            archive_read_data_skip(tarFile);
        }

        archive_read_free(tarFile);  // Libera o recurso do arquivo TAR
        return {fastaFiles, fastqFiles};  // Retorna as listas separadas de FASTA e FASTQ
    }





        std::pair<std::vector<FileContent>, std::vector<FileContent>> decompressZipFile(const std::string& filePath) {      
            unzFile zipFile = unzOpen(filePath.c_str());

            if (!zipFile) {
                throw std::runtime_error("Não foi possível abrir o arquivo ZIP: " + filePath);
            }


            if (unzGoToFirstFile(zipFile) != UNZ_OK) {
                unzClose(zipFile);
                throw std::runtime_error("Erro ao acessar o primeiro arquivo no ZIP: " + filePath);
            }

            std::vector<std::string> fastaFiles;
            std::vector<std::string> fastqFiles;

            do {
                // Obtém o nome do arquivo atual no ZIP
                char fileName[256];
                unz_file_info fileInfo;
                
                if (unzGetCurrentFileInfo(zipFile, &fileInfo, fileName, sizeof(fileName), nullptr, 0, nullptr, 0) != UNZ_OK) {
                    throw std::runtime_error("Erro ao obter informações do arquivo no ZIP.");
                }

                std::string currentFileName(fileName);
                std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);


                if (currentFileName.ends_with(".fasta") || currentFileName.ends_with(".fastq")) {
                    std::cout << "Arquivo encontrado: " << currentFileName << std::endl;

                    // Abre o arquivo atual dentro do ZIP
                    if (unzOpenCurrentFile(zipFile) != UNZ_OK) {
                        throw std::runtime_error("Erro ao abrir arquivo dentro do ZIP: " + currentFileName);
                    }

                    // Lê o conteúdo do arquivo
                    std::string fileContent;
                    char buffer[4096];
                    int bytesRead;
                    while ((bytesRead = unzReadCurrentFile(zipFile, buffer, sizeof(buffer))) > 0) {
                        fileContent.append(buffer, bytesRead);
                    }

                    if (bytesRead < 0) {
                        unzCloseCurrentFile(zipFile);
                        throw std::runtime_error("Erro ao ler o arquivo dentro do ZIP: " + currentFileName);
                    }

                    unzCloseCurrentFile(zipFile);

                    // Armazena o conteúdo no vetor apropriado
                    FileContent fileData = {currentFileName, fileContent};
                    if (currentFileName.ends_with(".fasta")) {
                        fastaFiles.push_back(fileData);
                    } else if (currentFileName.ends_with(".fastq")) {
                        fastqFiles.push_back(fileData);
                    }
                } else {
                    std::cout << "Arquivo ignorado: " << currentFileName << std::endl;
                }
            } while (unzGoToNextFile(zipFile) != UNZ_END_OF_LIST_OF_FILE);

            unzClose(zipFile);
            return {fastaFiles, fastqFiles};
            }

        virtual ~GeneticMaterial() {}
};
        // Verificações para saber se os arquivos fasta e fastq são válidos
        // Falta acerta essa questão dos arquivos passados como fasta ou fastq. Quantos arquivos podem ser passado de uma só vez?
