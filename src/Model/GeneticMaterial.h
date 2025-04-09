#ifndef GENETIC_MATERIAL_H  // Corrigido para corresponder ao nome do arquivo
#define GENETIC_MATERIAL_H

#include <string>
#include <vector>
#include <optional>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <variant>
#include <algorithm>
#include <unordered_map>



class GeneticMaterial {

    public:

        struct FileContent {
            std::string fileName;
            std::optional<std::string> fileType; // Pode ser "FASTA", "FASTQ" ou std::nullopt
            std::optional<std::string> content; // Conteúdo do arquivo, ou std::nullopt se não carregado
        };

            
        struct SequenceDescription {
            std::string sequenceId;
            std::optional<std::string> speciesName;
            std::optional<size_t> sequenceLength;
            std::optional<std::string> additionalInfo;
            std::string fileType;
        };

        // Estrutura para representar uma sequência processada
        struct ProcessedSequence {
            std::string header;
            std::string sequence;
            std::string quality;  // Para FASTQ
            double avgQuality;   
        };


        struct KmerMatch {
            int position;
            bool isReverseComplement;
        };

    private:
        std::string filePath;
        std::string fileContent;
        std::string fileExtension;
        FileContent fileData;
        std::vector<FileContent> m_loadedFiles;


        std::string sequence; // Sequência genética
        std::vector<int> suffixArray; // Tabela de sufixos
        std::vector<int> lcpArray; // Tabela de prefixos comuns (LCP)


    public:      
        explicit GeneticMaterial(const std::string& filePath);

        GeneticMaterial() = default;  // Construtor padrão

        std::string getFilePath() const; 

    
        std::pair<std::vector<FileContent>, std::vector<FileContent>> loadedData;

        // Método que retorna os dados carregados (FASTA e FASTQ)
        std::pair<std::vector<FileContent>, std::vector<FileContent>> getLoadedFileData() const;
        
        
       
        static std::pair<std::vector<FileContent>, std::vector<FileContent>> decompressZipFile(const std::string& filePath);
        static std::pair<std::vector<FileContent>, std::vector<FileContent>> openFile(const std::string& filePath);
        static std::optional<std::pair<std::optional<std::vector<FileContent>>, std::optional<std::vector<FileContent>>>> loadFileContent(const std::variant<std::string, std::vector<std::string>>& input);
        
        static bool verifySequence(const std::string& content, const std::string& fileType);
        static bool isValidNucleotide(char c);
        static bool isValidQualityChar(char c);
        static std::string normalizeNewlines(const std::string& input);
        
        size_t countSequencesInMultifasta(const std::vector<FileContent>& fastaFiles) const;

        std::vector<SequenceDescription> getSequenceDescriptions(const std::optional<std::pair<std::optional<std::vector<FileContent>>, std::optional<std::vector<FileContent>>>>& loadedFiles);

        static std::optional<std::vector<std::pair<std::string, double>>> calculateGCContent(const std::optional<std::pair<std::optional<std::vector<FileContent>>, std::optional<std::vector<FileContent>>>>& loadedFiles, std::optional<std::pair<size_t, size_t>> range = std::nullopt);
        
        static double calculateSequenceGC(const std::string& sequence, const std::optional<std::pair<size_t, size_t>>& range);
        
        static double calculateAverageQuality(const std::string& qualityStr);

        static std::string toString(const ProcessedSequence& sequence);
        
        static void removePrefix(std::string& sequence, std::string& quality, const std::string& prefix);
        
        static void trimPolyAT(std::string& sequence, std::string& quality);
        
        static bool containsDegenerateBases(const std::string& sequence);

        static std::string getComplementSequenceAndSave(const std::string& sequence, const std::string& outputFilePath);

        void calculateAndSaveComplementSequences(
            const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
            const std::string& outputDir = "");
        


        std::optional<std::vector<FileContent>> filterFastqByQuality(
            const std::optional<std::vector<FileContent>>& fastqFiles,
            const std::string& qualityStr,
            const std::string& outputDir);

            
        std::optional<std::vector<ProcessedSequence>> processAndCleanSequences(
            const std::optional<std::vector<FileContent>>& loadedFiles,
            const std::optional<std::string>& prefixToRemove,
            bool removePolyAT,
            double minQuality,
            bool removeDegenerate);

        static std::optional<std::pair<
            std::optional<std::vector<FileContent>>, 
            std::optional<std::vector<FileContent>>
            >> calculateComplementSequences(const std::optional<std::pair<
            std::optional<std::vector<FileContent>>, 
            std::optional<std::vector<FileContent>>
            >>& loadedFiles);
        
        static std::optional<std::pair<
            std::optional<std::vector<FileContent>>, 
            std::optional<std::vector<FileContent>>
        >> calculateReverseComplementSequences(const std::optional<std::pair<
            std::optional<std::vector<FileContent>>, 
            std::optional<std::vector<FileContent>>
        >>& loadedFiles);
        


        std::vector<std::string> getSequenceBetweenStartAndStop(
            const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
                                          std::optional<std::vector<GeneticMaterial::FileContent>>>>& fileData,
            const std::string& startCodon,
            const std::vector<std::string>& stopCodons);
        
        static std::optional<std::pair<
        std::optional<std::vector<GeneticMaterial::FileContent>>,
        std::optional<std::vector<GeneticMaterial::FileContent>>
        >> extractSubsequences(
        const std::optional<std::pair<
            std::optional<std::vector<GeneticMaterial::FileContent>>,
            std::optional<std::vector<GeneticMaterial::FileContent>>
        >>& loadedFiles,
        size_t start, 
        size_t length, 
        const std::optional<std::string>& sequenceId = std::nullopt);

        // Função auxiliar para obter o complemento de um nucleotídeo
        static char getComplement(char nucleotide);

        // Função para calcular o reverse complement de uma sequência
        static std::string getReverseComplement(const std::string& sequence);

        void buildAndSaveSuffixTable(); // Declaração do método
        // Métodos para manipulação de sequências
        std::vector<int> buildSuffixArray(const std::string& sequence);
        std::vector<int> buildLCPArray(const std::string& sequence, const std::vector<int>& suffixArray);
        int getSequenceLength() const;
        std::string getFactor(int i, int k) const;
        std::vector<int> findRepeatedFactors(int k) const;
        int searchPattern(const std::string& pattern);

        const std::vector<int>& getSuffixArray() const;

        const std::vector<int>& getLCPArray() const;

        const std::string& getSequence() const { return sequence; }
    
        
        std::string extractSequenceFromLoadedData(
            const std::optional<std::pair<
                std::vector<FileContent>,
                std::vector<FileContent>>>& loadedData) const;

        std::vector<int> findRepeatedFactors(int minLength);

        void setSequence(const std::string& seq);
        
        void buildIndex();
       
        void buildKmerIndex(int k);
        std::vector<KmerMatch> mapReadToGenome(const std::string& read, int k) const;
        
    private:
        static bool endsWith(const std::string& str, const std::string& suffix);
        static bool startsWith(const std::string& str, const std::string& prefix);
        static std::pair<std::vector<FileContent>, std::vector<FileContent>> decompressGzipFile(const std::string& filePath);
        static std::pair<std::vector<FileContent>, std::vector<FileContent>> descompressTarFile(const std::string& filePath);
        static std::string extractFromSequence(const std::string& sequence, size_t start, size_t length);




    std::string genomeSequence;
    std::unordered_map<std::string, std::vector<int>> kmerIndex;

    std::string getReverseComplementKmer(const std::string& kmer) const;
    void loadGenomeSequence();

    public:
    // Função para mapear reads no genoma
    std::vector<std::pair<std::string, std::vector<KmerMatch>>> mapReads(const std::vector<std::string>& reads, int k);

    // Função para obter o suporte de um k-mer
    int getKmerSupport(const std::string& kmer) const;

    // Função para diagnosticar variações em um read
    std::optional<std::string> diagnoseVariation(const std::string& read, const std::vector<KmerMatch>& matches, int k) const;
        
    void mapReadsAndReport(const std::vector<std::string>& reads, int k);
};

#endif  // GENETIC_MATERIAL_H
