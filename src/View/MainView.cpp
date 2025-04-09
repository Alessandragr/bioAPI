#include "MainView.h"
#include "../Controller/GeneticMaterialController.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <locale>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>

#ifdef _WIN32
#include <windows.h>
#endif

// Implementação para Windows
#ifdef _WIN32
void MainView::initWindowsConsole() {
    SetConsoleOutputCP(65001); // UTF-8
}
#endif

void MainView::displayMessage(const std::string& message) {
    std::cout << message << std::endl;
}

void MainView::showUsage() {
    std::cout << "Usage: program --file <file_path> [option]\n"
              << "Available options:\n"
              << "  --count-fasta                 Count FASTA sequences\n"
              << "  --complement                  Calculate complementary sequence\n"
              << "  --reverse-complement          Calculate reverse-complement sequence\n"
              << "  --start-stop                  Extract sequence between codons\n"
              << "  --extract                     Extract subsequences\n"
              << "  --description                 Get sequence descriptions\n"
              << "  --gc-content                  Calculate GC content\n"
              << "  --average-quality             Calculate average sequence quality\n"
              << "  --remove-prefix               Remove sequence prefix\n"
              << "  --trim-polyAT                 Remove poly-A/poly-T tails\n"
              << "  --degenerate-bases            Check for degenerate bases\n"
              << "  --filter-fastq                Filter FASTQ by quality\n"
              << "  --process-clean               Process and clean sequences\n"
              << "  --suffix-table                Build and display suffix and LCP table\n"
              << "  --search-pattern              Search for a pattern in the sequence\n"
              << "  --find-repeated               Find repeated patterns\n"
              << "  --map-reads                   Map reads to the genome using k-mers\n"
              << "  -h, --help                    Display this help message\n";
}

void MainView::processCommandLine(int argc, char* argv[]) {
    // Configura encoding do console
    #ifdef _WIN32
    initWindowsConsole();
    #endif
    
    // Variáveis para os parâmetros
    struct ProcessConfig {
        std::optional<std::string> prefixToRemove = std::nullopt;
        bool removePolyAT = true;
        double minQuality = 20.0;
        bool removeDegenerate = true;
    };

    std::string filePath;
    bool hasFile = false;
    
    ProcessConfig config;
    bool processClean = false;
    
    int opt;
    int option_index = 0;
    GeneticMaterialController* controller = nullptr;
   
    static struct option long_options[] = {
        {"file", required_argument, nullptr, 'f'},               // Especificar o arquivo de entrada
        {"count-fasta", no_argument, nullptr, 'n'},              // Contar sequências FASTA
        {"complement", no_argument, nullptr, 'c'},               // Calcular sequência complementar
        {"reverse-complement", no_argument, nullptr, 'r'},       // Calcular sequência reverso-complementar
        {"start-stop", required_argument, nullptr, 's'},         // Extrair sequência entre códons
        {"extract", required_argument, nullptr, 'x'},            // Extrair subsequências
        {"description", no_argument, nullptr, 'd'},              // Obter descrição das sequências
        {"gc-content", no_argument, nullptr, 'g'},               // Calcular conteúdo GC
        {"average-quality", required_argument, nullptr, 'q'},    // Calcular qualidade média
        {"remove-prefix", required_argument, nullptr, 'p'},      // Remover prefixo da sequência
        {"trim-polyAT", required_argument, nullptr, 't'},        // Remover caudas poly-A/poly-T
        {"degenerate-bases", required_argument, nullptr, 'b'},   // Verificar bases degeneradas
        {"filter-fastq", required_argument, nullptr, 'F'},       // Filtrar FASTQ por qualidade
        {"process-clean", no_argument, nullptr, 'P'},            // Processar e limpar sequências
        {"help", no_argument, nullptr, 'h'},                     // Exibir ajuda
        {"prefix-to-remove", required_argument, nullptr, 'k'},   // Prefixo a ser removido
        {"remove-polyAT", required_argument, nullptr, 'a'},      // Remover caudas poly-A/poly-T
        {"min-quality", required_argument, nullptr, 'm'},        // Qualidade mínima
        {"suffix-table", no_argument, nullptr, 'S'},             // Construir e exibir a tabela de sufixos e LCP
        {"search-pattern", required_argument, nullptr, 'B'},     // Buscar padrão na sequência
        {"find-repeated", no_argument, nullptr, 'R'},            // Encontrar padrões repetidos
        {"map-reads", required_argument, nullptr, 'M'},          // Mapear reads usando k-mer
        {nullptr, 0, nullptr, 0}                                 // Finalizador
    };

    try {
        while ((opt = getopt_long(argc, argv, "f:ncrs:x:d:gqptbFPh", long_options, &option_index)) != -1) {
            
            switch (opt) {
                case 'f': { // --file
                    filePath = optarg;
                    hasFile = true;    // Marca que o arquivo foi especificado
                    break;
                }

                case 'd': { // --description
                    if (!hasFile) throw std::runtime_error("File not specified");
                    
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                    
                    std::cout << "Loading data from file: " << filePath << std::endl;
                    
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                    
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                    
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Retrieves the loaded data
                    auto loadedData = controller->getLoadedFileData();
                    if (!loadedData.has_value()) {
                        throw std::runtime_error("No data loaded to retrieve descriptions.");
                    }
                    
                    // Retrieves the descriptions
                    auto descriptions = controller->getSequenceDescriptions(loadedData);
                    std::cout << "Number of descriptions found: " << descriptions.size() << std::endl;
                    
                    // Displays the descriptions
                    for (const auto& desc : descriptions) {
                        std::cout << "ID: " << desc.sequenceId
                                  << ", Species: " << (desc.speciesName ? *desc.speciesName : "N/A")
                                  << ", Length: " << (desc.sequenceLength.has_value() ? std::to_string(desc.sequenceLength.value()) : "N/A")
                                  << ", Additional Info: " << (desc.additionalInfo ? *desc.additionalInfo : "N/A") << std::endl;
                    }
                    
                    break;
                }

                case 'n': { // --count-fasta
                    if (!hasFile) throw std::runtime_error("File not specified");
                    
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                    
                    std::cout << "Loading data from file: " << filePath << std::endl;
                    
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                    
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                    
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Counts the sequences in the multifasta
                    size_t totalSequences = controller->countSequencesInMultifasta();
                    std::cout << "Total number of sequences in the multifasta: " << totalSequences << std::endl;
                    
                    break;
                }
                
                case 'c': { // --complement
                    if (!hasFile) throw std::runtime_error("File not specified");
                    
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                    
                    std::cout << "Loading data from file: " << filePath << std::endl;
                
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                    
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                    
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Retrieves the loaded data and calculates the complementary sequences
                    auto loadedData = controller->getLoadedFileData();
                    if (!loadedData.has_value()) {
                        throw std::runtime_error("No data loaded to calculate complementary sequences.");
                    }
                    
                    controller->calculateComplementSequences(loadedData);
                    
                    std::cout << "Complementary sequences processed and saved successfully." << std::endl;
                    
                    break;
                }

                case 'r': { // --reverse-complement
                    if (!hasFile) throw std::runtime_error("File not specified");
                    
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                    
                    std::cout << "Loading data from file: " << filePath << std::endl;
                    
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                    
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                    
                    std::cout << "Data loaded into the controller: " 
                            << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Retrieves the loaded data and calculates the reverse-complement sequences
                    auto loadedData = controller->getLoadedFileData();
                    if (!loadedData.has_value()) {
                        throw std::runtime_error("No data loaded to calculate the reverse-complement sequences.");
                    }
                    
                    controller->calculateReverseComplementSequences(loadedData);
                    
                    std::cout << "Reverse-complement sequences processed and saved successfully." << std::endl;
                    
                    break;
                }

                case 'g': { // --gc-content
                    if (!hasFile) throw std::runtime_error("File not specified");
                    
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                    
                    std::cout << "Loading data from file: " << filePath << std::endl;
                    
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                    
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                    
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                    
                    // Retrieves the loaded data
                    auto loadedData = controller->getLoadedFileData();
                    if (!loadedData.has_value()) {
                        throw std::runtime_error("No data loaded to calculate GC content.");
                    }
                    
                    // Calculates the GC content
                    auto gcContent = controller->calculateGCContent(loadedData);
                    if (gcContent.has_value()) {
                        std::cout << "GC content calculated successfully:\n";
                        for (const auto& [sequenceId, gcValue] : *gcContent) {
                            std::cout << "ID: " << sequenceId << "\n GC Content: " << gcValue << "%" << std::endl;
                        }
                    } else {
                        std::cerr << "Error: Could not calculate GC content." << std::endl;
                    }
                    
                    break;
                }

                case 's': { // --start-stop
                    if (!hasFile) throw std::runtime_error("File not specified");
                
                    controller = new GeneticMaterialController(filePath);
                
                    // Parse the start and stop codons
                    std::string codonArgs(optarg);  // Example: "ATG,TAA,TAG,TGA"
                    std::stringstream ss(codonArgs);
                    std::string codon;
                
                    std::string startCodon;
                    std::vector<std::string> stopCodons;
                
                    if (std::getline(ss, codon, ',')) {
                        startCodon = codon;
                        while (std::getline(ss, codon, ',')) {
                            stopCodons.push_back(codon);
                        }
                    } else {
                        throw std::runtime_error("Invalid codons: provide at least 1 start and 1 stop codon.");
                    }
                
                    // Load and store the data
                    auto fileData = controller->loadFileContent(filePath);
                    controller->setLoadedFileData(fileData);
                
                    // Call the controller function
                    auto sequences = controller->getSequenceBetweenStartAndStop(
                        controller->getLoadedFileData(),
                        startCodon,
                        stopCodons
                    );
                
                    // Display the sequences found
                    for (const auto& seq : sequences) {
                        std::cout << "Sequence between codons: " << seq << std::endl;
                    }
                
                    break;
                }
                
                case 'x': { // --extract
                    if (!hasFile) throw std::runtime_error("File not specified");
                
                    controller = new GeneticMaterialController(filePath);
                
                    // Parse the parameters start, length, and sequenceId
                    std::string extractArgs(optarg);  // Example: "start=10,length=50,id=seq1"
                    size_t start = 0, length = 0;
                    std::optional<std::string> sequenceId = std::nullopt;
                
                    std::stringstream ss(extractArgs);
                    std::string param;
                    while (std::getline(ss, param, ',')) {
                        auto delimiterPos = param.find('=');
                        if (delimiterPos == std::string::npos) continue;
                
                        std::string key = param.substr(0, delimiterPos);
                        std::string value = param.substr(delimiterPos + 1);
                
                        if (key == "start") {
                            start = std::stoul(value);
                        } else if (key == "length") {
                            length = std::stoul(value);
                        } else if (key == "id") {
                            sequenceId = value;
                        }
                    }
                
                    if (length == 0) {
                        throw std::runtime_error("Parameter 'length' is required and must be greater than 0.");
                    }
                
                    // Load and store the data
                    auto fileData = controller->loadFileContent(filePath);
                    controller->setLoadedFileData(fileData);
                
                    // Call the controller function
                    auto result = controller->extractSubsequences(
                        controller->getLoadedFileData(),
                        start,
                        length,
                        sequenceId
                    );
                
                    // Display the extracted subsequences
                    if (result) {
                        auto [fasta, fastq] = *result;
                        if (fasta) {
                            for (const auto& file : *fasta) {
                                std::cout << "FASTA subsequence: " << *file.content << "\n";
                            }
                        }
                        if (fastq) {
                            for (const auto& file : *fastq) {
                                std::cout << "FASTQ subsequence: " << *file.content << "\n";
                            }
                        }
                    } else {
                        std::cout << "No subsequences found." << std::endl;
                    }
                
                    break;
                }
                
                
                case 'q': { // --average-quality
                    if (!hasFile) throw std::runtime_error("File not specified");
                    controller = new GeneticMaterialController(filePath);
                
                    std::string qualityStr = optarg;
                    double avg = controller->calculateAverageQuality(qualityStr);
                    std::cout << "Average quality: " << avg << std::endl;
                
                    break;
                }
                
                case 'F': { // --filter-fastq
                    if (!hasFile) throw std::runtime_error("File not specified");
                
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                
                    std::cout << "Loading data from file: " << filePath << std::endl;
                
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                        throw std::runtime_error("Error loading the file.");
                    }
                
                    // Sets the loaded data in the controller
                    controller->setLoadedFileData(fileData);
                
                    std::cout << "Data loaded into the controller: " 
                              << (controller->getLoadedFileData().has_value() ? "Yes" : "No") << std::endl;
                
                    // Retrieves the loaded data
                    auto loadedData = controller->getLoadedFileData();
                    if (!loadedData.has_value()) {
                        throw std::runtime_error("No data loaded to filter.");
                    }
                
                    // Accesses the FASTQ files (second element of the pair)
                    const auto& fastqFiles = loadedData->second;
                    if (fastqFiles.empty()) {
                        throw std::runtime_error("No valid FASTQ files loaded to filter.");
                    }
                
                    // Converts the FASTQ files to std::optional
                    std::optional<std::vector<GeneticMaterial::FileContent>> optionalFastqFiles = std::make_optional(fastqFiles);
                
                    // Calls the function to filter the FASTQ files
                    std::string qualityStr = optarg; // Minimum quality passed as an argument
                    auto filteredFiles = controller->filterFastqByQuality(optionalFastqFiles, qualityStr, "./");
                
                    if (filteredFiles.has_value()) {
                        std::cout << "Filtered files successfully saved in the current directory.\n";
                    } else {
                        std::cout << "No files met the quality criteria.\n";
                    }
                
                    break;
                }
                
                case 'a': { // --remove-polyAT
                    std::string value = optarg;
                    if (value == "true") {
                        config.removePolyAT = true;
                    } else if (value == "false") {
                        config.removePolyAT = false;
                    } else {
                        throw std::runtime_error("Invalid value for --remove-polyAT. Use 'true' or 'false'.");
                    }
                    break;
                }
                
                case 'b': { // --degenerate-bases
                    std::string value = optarg;
                    if (value == "true") {
                        config.removeDegenerate = true;
                    } else if (value == "false") {
                        config.removeDegenerate = false;
                    } else {
                        throw std::runtime_error("Invalid value for --degenerate-bases. Use 'true' or 'false'.");
                    }
                    break;
                }

                case 'B': { // --search-pattern
                    if (!hasFile) throw std::runtime_error("File not specified");
                
                    // Creates the controller with the file path
                    controller = new GeneticMaterialController(filePath);
                
                    std::cout << "Loading data from file: " << filePath << std::endl;
                
                    // Loads the file data
                    auto fileData = controller->loadFileContent(filePath);
                    if (!fileData.has_value()) {
                        throw std::runtime_error("Error loading the file.");
                    }
                
                    controller->setLoadedFileData(fileData);
                
                    // Builds the suffix array and LCP table
                    controller->buildIndex();
                
                    // Gets the pattern directly from optarg
                    std::string pattern = optarg;
                
                    // Performs the search
                    int position = controller->searchPattern(pattern);
                    if (position != -1) {
                        std::cout << "Pattern found at position: " << position << std::endl;
                    } else {
                        std::cout << "Pattern not found." << std::endl;
                    }
                
                    break;
                }
                
                case 'k': { // --prefix-to-remove
                    config.prefixToRemove = optarg;
                    break;
                }
                
                case 'm': { // --min-quality
                    config.minQuality = std::stod(optarg);
                    break;
                }


                case 'S': { // --suffix-table
                    if (!hasFile) throw std::runtime_error("File not specified");
                
                    controller = new GeneticMaterialController(filePath);
                
 
                    controller->buildAndSaveSuffixTable();
                
                    break;
                }
                
                case 'h':
                    showUsage();
                    exit(0);
                    
                case -1: // End of options
                    break;
                    
                default:
                    showUsage();
                    exit(1);
            }
        }

        if (!hasFile) {
            throw std::runtime_error("File path is required");
        }

        
        // Executa o comando apropriado
        if (processClean) {
            std::cout << "\nConfigured parameters:\n";
            std::cout << "File: " << filePath << "\n";
            std::cout << "Prefix to remove: " << (config.prefixToRemove ? *config.prefixToRemove : "None") << "\n";
            std::cout << "Remove poly-A/poly-T tails: " << (config.removePolyAT ? "Yes" : "No") << "\n";
            std::cout << "Minimum quality: " << config.minQuality << "\n";
            std::cout << "Remove degenerate bases: " << (config.removeDegenerate ? "Yes" : "No") << "\n";
            std::cout << "Process and clean sequences: " << (processClean ? "Yes" : "No") << "\n";
        
            if (!hasFile) throw std::runtime_error("File not specified");
        
            controller = new GeneticMaterialController(filePath);
        
            // Load the file data
            auto fileData = controller->loadFileContent(filePath);
            if (!fileData.has_value()) {
                std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
                throw std::runtime_error("Error loading the file.");
            }
        
            // Configura os dados carregados no controlador
            controller->setLoadedFileData(fileData);
        
            // Verifica o conteúdo dos dados carregados
            auto loadedData = controller->getLoadedFileData();
            if (!loadedData.has_value()) {
                throw std::runtime_error("No data loaded to process.");
            }
        
            // Combina os arquivos FASTA e FASTQ em um único vetor
            std::vector<GeneticMaterial::FileContent> combinedFiles;
            if (!loadedData->first.empty()) {
                combinedFiles.insert(combinedFiles.end(), loadedData->first.begin(), loadedData->first.end());
            }
            if (!loadedData->second.empty()) {
                combinedFiles.insert(combinedFiles.end(), loadedData->second.begin(), loadedData->second.end());
            }
        
            if (combinedFiles.empty()) {
                throw std::runtime_error("Nenhum dado válido carregado para processamento.");
            }
        
            // Processes and cleans the sequences
            auto cleanedFiles = controller->processAndCleanSequences(
                std::make_optional(combinedFiles), // Combine FASTA and FASTQ
                config.prefixToRemove,            // Prefix to remove
                config.removePolyAT,              // Remove poly-A/poly-T tails
                config.minQuality,                // Minimum quality
                config.removeDegenerate           // Remove degenerate characters
            );

            // Displays a summary of the processing
            if (cleanedFiles) {
                double totalQuality = 0.0;
                for (const auto& seq : *cleanedFiles) {
                    totalQuality += seq.avgQuality;
                }
                double avgQuality = totalQuality / cleanedFiles->size();

                std::cout << "\nProcessing completed successfully.\n";
                std::cout << "Total sequences processed: " << combinedFiles.size() << "\n";
                std::cout << "Total sequences cleaned: " << cleanedFiles->size() << "\n";
                std::cout << "Average quality of cleaned sequences: " << avgQuality << "\n";
            } else {
                std::cout << "\nNo sequences met the cleaning criteria.\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        showUsage();
        if (controller) delete controller;
        exit(1);
    }

    if (controller) delete controller;
}

