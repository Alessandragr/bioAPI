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
#include <unzip.h>
#include <stdio.h>
#include <cctype>
#include <variant>


bool startsWith(const std::string& str, const std::string& prefix) {
    if (prefix.size() > str.size()) return false;
    return str.compare(0, prefix.size(), prefix) == 0;
}


bool endsWith(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}



class GeneticMaterial {

    struct FileContent {
    std::string fileName;
    std::optional<std::string> fileType; // Pode ser "FASTA", "FASTQ" ou std::nullopt
    std::optional<std::string> content; // Conteúdo do arquivo, ou std::nullopt se não carregado

    

    };

    protected:
        std::string filePath;
        std::string fileContent;
        std::string fileExtension;
        FileContent fileData;

    public:
        GeneticMaterial(const std::string& path) : filePath(path) {}
              


        std::pair<std::vector<FileContent>, std::vector<FileContent>> decompressGzipFile(const std::string& filePath) {
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
            std::vector<FileContent> fastaFiles;
            std::vector<FileContent> fastqFiles;

            std::istringstream contentStream(decompressedContent.str());
            std::string line;
            FileContent currentFile;

            while (std::getline(contentStream, line)) {
                if (line.empty()) continue;
                if (startsWith(line, ">")) { // Início de um arquivo FASTA
                    if (!currentFile.content.has_value()) {
                        if (currentFile.fileType == "FASTA") {
                            fastaFiles.push_back(currentFile);
                        } else if (currentFile.fileType == "FASTQ") {
                            fastqFiles.push_back(currentFile);
                        }
                    }
                    currentFile = FileContent();
                    currentFile.fileType = "FASTA";
                    *currentFile.content += line + "\n"; // Acessa o valor do optional
                } else if (startsWith(line, "@")) { // Início de um arquivo FASTQ
                     if (currentFile.fileType == "FASTA") {
                        fastaFiles.push_back(currentFile);
                    } else if (currentFile.fileType == "FASTQ") {
                        fastqFiles.push_back(currentFile);
                    }
                    currentFile = FileContent();
                    currentFile.fileType = "FASTQ";
                    *currentFile.content += line + "\n"; // Acessa o valor do optional
                } else {
                    *currentFile.content += line + "\n"; // Acessa o valor do optional
                }
            }

            // Adiciona o último arquivo processado
            if (!currentFile.content.has_value()) {
                if (currentFile.fileType == "FASTA") {
                    fastaFiles.push_back(currentFile);
                } else if (currentFile.fileType == "FASTQ") {
                    fastqFiles.push_back(currentFile);
                }
            }

            return {fastaFiles, fastqFiles};
        }

        



        std::pair<std::vector<FileContent>, std::vector<FileContent>> descompressTarFile(const std::string& filePath) {
            struct archive* tarFile = archive_read_new();

            if (!tarFile) {
                throw std::runtime_error("Falha ao alocar memória para o manipulador de TAR.");
            }

            struct archive_entry* entry;
            std::vector<FileContent> fastaFiles;
            std::vector<FileContent> fastqFiles;

            // Configurando o libarchive para lidar com arquivos TAR
            archive_read_support_filter_all(tarFile);
            archive_read_support_format_all(tarFile);

            if (archive_read_open_filename(tarFile, filePath.c_str(), 10240) != ARCHIVE_OK) {
                archive_read_free(tarFile);  // Libera a memória antes de lançar a exceção e evita o vazamento de dados
                throw std::runtime_error("Não foi possível abrir o arquivo TAR: " + filePath);
            }

            while (archive_read_next_header(tarFile, &entry) == ARCHIVE_OK) {
                std::string currentFileName = archive_entry_pathname(entry);

                // Converte o nome do arquivo para lowercase para evitar problemas de case
                std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);

                // Filtra apenas arquivos .fasta ou .fastq
                if (endsWith(currentFileName, ".fasta") || endsWith(currentFileName, ".fastq")) {
                    std::cout << "Arquivo encontrado: " << currentFileName << std::endl;

                    std::string fileType = (endsWith(currentFileName, ".fasta")) ? "FASTA" : "FASTQ";

                    // Lê o conteúdo do arquivo
                    std::ostringstream fileContent;
                    const void* buffer;
                    size_t size;
                    int64_t offset;

                    while (true) {
                        int readBytes = archive_read_data_block(tarFile, &buffer, &size, &offset);
                        if (readBytes == ARCHIVE_EOF) break;
                        if (readBytes < ARCHIVE_OK) {
                            archive_read_free(tarFile);
                            throw std::runtime_error("Erro ao ler o arquivo: " + currentFileName);
                        }

                        if (buffer) {  // Verifica se buffer não é nullptr (caso extremo)
                            fileContent.write(static_cast<const char*>(buffer), size);  // Adiciona o conteúdo ao stream
                        }
                    }

                    // Cria o objeto FileContent
                    FileContent fileData;
                    fileData.fileName = currentFileName;
                    fileData.fileType = std::make_optional(fileType);
                    fileData.content = fileContent.str();

                    if (fileType == "FASTA") {
                        fastaFiles.push_back(std::move(fileData)); // Adiciona à lista de FASTA
                    } else {
                        fastqFiles.push_back(std::move(fileData));  // Adiciona à lista de FASTQ
                    }

                } else {
                    std::cout << "Arquivo ignorado: " << currentFileName << std::endl;
                }
                // Avança para o próximo arquivo no TAR
                archive_read_data_skip(tarFile);
            }

            archive_read_free(tarFile);  // Libera o recurso do arquivo TAR
            return {std::move(fastaFiles), std::move(fastqFiles)};  // Retorna as listas separadas de FASTA e FASTQ
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

        // Vetores para armazenar os arquivos FASTA e FASTQ
        std::vector<FileContent> fastaFiles;
        std::vector<FileContent> fastqFiles;

        do {
            // Obtém o nome do arquivo atual no ZIP
            char fileName[256];
            unz_file_info fileInfo;

            if (unzGetCurrentFileInfo(zipFile, &fileInfo, fileName, sizeof(fileName), nullptr, 0, nullptr, 0) != UNZ_OK) {
                unzClose(zipFile);
                throw std::runtime_error("Erro ao obter informações do arquivo no ZIP.");
            }

            std::string currentFileName(fileName);
            std::transform(currentFileName.begin(), currentFileName.end(), currentFileName.begin(), ::tolower);

            // Filtra apenas arquivos .fasta ou .fastq
            if (endsWith(currentFileName, ".fasta") || endsWith(currentFileName, ".fastq")) {
                std::cout << "Arquivo encontrado: " << currentFileName << std::endl;

                // Abre o arquivo atual dentro do ZIP
                if (unzOpenCurrentFile(zipFile) != UNZ_OK) {
                    unzClose(zipFile);
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
                    unzClose(zipFile);
                    throw std::runtime_error("Erro ao ler o arquivo dentro do ZIP: " + currentFileName);
                }

                unzCloseCurrentFile(zipFile);

                // Cria o objeto FileContent
                FileContent fileData;
                fileData.fileName = currentFileName;
                fileData.fileType = std::make_optional(endsWith(currentFileName, ".fasta") ? "FASTA" : "FASTQ");
                fileData.content = std::move(fileContent);

                // Armazena o conteúdo no vetor apropriado
                if (endsWith(currentFileName, ".fasta")) {
                    fastaFiles.push_back(std::move(fileData));
                } else if (endsWith(currentFileName, ".fastq")) {
                    fastqFiles.push_back(std::move(fileData));
                }
            } else {
                std::cout << "Arquivo ignorado: " << currentFileName << std::endl;
            }
        } while (unzGoToNextFile(zipFile) == UNZ_OK);

        unzClose(zipFile);
        return {std::move(fastaFiles), std::move(fastqFiles)};
    }




    std::pair<std::vector<FileContent>, std::vector<FileContent>> verifyAndLoadFile(const std::string& filePath) {
        std::pair<std::vector<FileContent>, std::vector<FileContent>> fileData;

        size_t position = filePath.find_last_of('.');
        if (position == std::string::npos) {
            throw std::runtime_error("Arquivo sem extensão, insira um arquivo com uma extensão válida: FASTA, FASTQ, ZIP ou TAR.");
        }

        std::string fileExtension = filePath.substr(position + 1);
        std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

        try {
            if (fileExtension == "gz") {
                fileData = decompressGzipFile(filePath);
            } else if (fileExtension == "tar") {
                fileData = descompressTarFile(filePath);
            } else if (fileExtension == "zip") {
                fileData = decompressZipFile(filePath);
            } else if (fileExtension == "fasta" || fileExtension == "fastq") {
                // Chama a função loadFileContent e verifica se o retorno é válido
                auto result = loadFileContent(filePath);
                if (result) {
                    const auto& [fastaFiles, fastqFiles] = *result;
                    // Verifica se as listas não são nullopt e atribui às listas de fileData
                    if (fastaFiles) {
                        fileData.first = *fastaFiles;
                    }
                    if (fastqFiles) {
                        fileData.second = *fastqFiles;
                    }
            } else {
                throw std::runtime_error("Falha ao carregar conteúdo do arquivo: " + filePath);
            }
        } else {
            throw std::runtime_error("Extensão não suportada: " + fileExtension +
                "\nPor favor, insira arquivos do tipo: FASTA, FASTQ, ZIP ou TAR.");
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Erro ao processar o arquivo: " + std::string(e.what()));
        }

        return fileData;
    }






    std::optional<std::pair<std::optional<std::vector<FileContent>>, std::optional<std::vector<FileContent>>>> loadFileContent(const std::variant<std::string, 
        std::vector<std::string>>& input) {

        std::vector<std::string> filePaths;

        // Verifica se o input é uma string ou um vetor de strings
        if (std::holds_alternative<std::string>(input)) {
            filePaths.push_back(std::get<std::string>(input)); // Adiciona o único caminho à lista
        } else {
            filePaths = std::get<std::vector<std::string>>(input); // Converte diretamente para a lista
        }

        // Vetores de arquivos FASTA e FASTQ
        std::optional<std::vector<FileContent>> fastaFiles;
        std::optional<std::vector<FileContent>> fastqFiles;

        // Processa cada arquivo na lista
        for (const auto& filePath : filePaths) {
            std::ifstream file(filePath, std::ios::binary);

            if (!file.is_open()) {
                throw std::runtime_error("Não foi possível abrir o arquivo: " + filePath);
            }

            // Processar o arquivo em blocos
            constexpr size_t bufferSize = 4096; // 4 KB
            char buffer[bufferSize];
            std::stringstream fileContent;

            while (file.read(buffer, bufferSize)) {
                fileContent.write(buffer, file.gcount()); // Usa write em vez de append
            }

            // Verifica se há bytes restantes no final do arquivo
            if (file.gcount() > 0) {
                fileContent.write(buffer, file.gcount()); // Usa write em vez de append
            }

            file.close();

            // Cria o objeto FileContent
            FileContent fileData;
            fileData.fileName = filePath;
            fileData.content = fileContent.str();

            // Verifica o tipo de arquivo e adiciona ao vetor apropriado
            if (endsWith(filePath, ".fasta")) {
                if (!fastaFiles) {
                fastaFiles = std::vector<FileContent>();  // Inicializa o vetor se necessário
            }
            fileData.fileType = "FASTA";
            fastaFiles->push_back(fileData);
            } else if (endsWith(filePath, ".fastq")) {
                if (!fastqFiles) {
                    fastqFiles = std::vector<FileContent>();  // Inicializa o vetor se necessário
                }
                fileData.fileType = "FASTQ";
                fastqFiles->push_back(fileData);
            } else {
                throw std::runtime_error("Tipo de arquivo não suportado: " + filePath);
            }
        }
        // Retorna a lista de arquivos FASTA e FASTQ como opcional
        return std::make_optional(std::make_pair(std::move(fastaFiles), std::move(fastqFiles)));
    }
        





    bool verifySequence(const std::string& fastaContent) {
        bool erroEncontrado = false;
        std::string line;
        bool isHeader = true;

        // Percorre cada caractere da string
        for (size_t i = 0; i < fastaContent.size(); ++i) {
            char c = fastaContent[i];

            // Verifica se é uma nova linha (fim de sequência ou cabeçalho)
            if (c == '\n' || c == '\r') {
                // Se a linha estava com cabeçalho (iniciada por ">", ignorar)
                if (line.empty() || line[0] == '>') {
                    // Linha de cabeçalho encontrada, ignora e reinicia
                    line.clear();
                    isHeader = true;
                } else {
                    // Verifica os caracteres da sequência
                    for (char seqChar : line) {
                        if (!isalpha(seqChar) && seqChar != 'N' && seqChar != 'n') {
                            std::cerr << "Caractere inválido encontrado na sequência: " << seqChar << std::endl;
                            erroEncontrado = true;
                        }
                    }
                    line.clear();
                    isHeader = false;
                }
            } else {
                // Acumula os caracteres da sequência ou cabeçalho
                if (isHeader || c == '>') {
                    line.push_back(c);
                } else {
                    line.push_back(c); // Acumula caracteres da sequência
                }
            }
        }

        return !erroEncontrado;  // Retorna true se não houver caracteres inválidos
    }




    // Lê a sequência e me retorna quantos alinhamentos é possível fazer com ela
    // Função para calcular o coeficiente binomial: C(n, k) = n! / (k! * (n-k)!)
    uint64_t binomialCoefficient(int n, int k) {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;

        uint64_t result = 1;
        for (int i = 1; i <= k; ++i) {
            result = result * (n - i + 1) / i;
        }
        return result;
    }

    // Função para calcular o número de alinhamentos globais possíveis
    uint64_t calculateGlobalAlignments(int length1, int length2) {
        int totalLength = length1 + length2;
        return binomialCoefficient(totalLength, length1);
    }

    // Função para calcular o número de alinhamentos locais possíveis
    uint64_t calculateLocalAlignments(int length1, int length2) {
        return static_cast<uint64_t>((length1 + 1) * (length2 + 1) * (length1 + length2) / 2);
    }

    virtual ~GeneticMaterial() {}
};
        // Verificações para saber se os arquivos fasta e fastq são válidos.
