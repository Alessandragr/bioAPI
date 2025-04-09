#ifndef GENETIC_MATERIAL_CONTROLLER_H
#define GENETIC_MATERIAL_CONTROLLER_H

#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <variant>

// Including the header of the GeneticMaterial class, as it is used here
#include "../Model/GeneticMaterial.h"

class GeneticMaterialController {
    private:
        GeneticMaterial geneticMaterial;  // Attribute to store the GeneticMaterial model
        std::optional<std::pair<
        std::vector<GeneticMaterial::FileContent>,
        std::vector<GeneticMaterial::FileContent>>> loadedFileData;


        GeneticMaterial& model; // Reference to the GeneticMaterial model
        std::vector<int> suffixArray; // Stores the suffix array
        std::vector<int> lcpArray;   // Stores the LCP (Longest Common Prefix) array


    public:

        // Constructor that initializes the controller with a file path
        explicit GeneticMaterialController(const std::string& filePath);

        // Loads and sets file data from the specified file path
        void loadAndSetFileData(const std::string& filePath);

        // Counts the number of sequences in a MULTIFASTA file
        size_t countSequencesInMultifasta() const;

        // Sets the loaded file data
        void setLoadedFileData(
            const std::optional<std::pair<
                std::optional<std::vector<GeneticMaterial::FileContent>>,
                std::optional<std::vector<GeneticMaterial::FileContent>>>>& data);

        // Retrieves the loaded file data
        std::optional<std::pair<
        std::vector<GeneticMaterial::FileContent>,
        std::vector<GeneticMaterial::FileContent>>>
        getLoadedFileData() const;

        // Retrieves the file path associated with the GeneticMaterial object
        std::string getFilePath() const;

        // Methods to interact with the GeneticMaterial object
        std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> openFile();
        std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> loadFileContent(const std::variant<std::string, std::vector<std::string>>& input);
        bool verifySequence(const std::string& content, const std::string& fileType);

        // Calculates the complement sequences for the loaded files
        std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
            std::optional<std::vector<GeneticMaterial::FileContent>>>>
             calculateComplementSequences(const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
                std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles);
        
        // Calculates the reverse complement sequences for the loaded files
        std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
        std::optional<std::vector<GeneticMaterial::FileContent>>>> 
        calculateReverseComplementSequences(const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
            std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles);

        // Saves the complement sequence to a specified output file
        void saveComplementSequence(const std::string& outputFilePath);

        // Extracts sequences between a start codon and stop codons
        std::vector<std::string> getSequenceBetweenStartAndStop(
            const std::optional<std::pair<
                std::optional<std::vector<GeneticMaterial::FileContent>>,
                std::optional<std::vector<GeneticMaterial::FileContent>>>>& fileData,
            const std::string& startCodon,
            const std::vector<std::string>& stopCodons);


        // Extracts sequences between a start codon and stop codons
        std::string getSequenceBetweenStartAndStop(const std::string& sequence, const std::string& startCodon, const std::vector<std::string>& stopCodons);
        
        // Retrieves sequence descriptions from the loaded files
        std::vector<GeneticMaterial::SequenceDescription> getSequenceDescriptions(const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles);
        
        std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> extractSubsequences(
            const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
            size_t start,
            size_t length,
            const std::optional<std::string>& sequenceId = std::nullopt
        );
        std::optional<std::vector<std::pair<std::string, double>>> calculateGCContent(
            const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
            std::optional<std::pair<size_t, size_t>> range = std::nullopt
        );
   
        
        std::string getReverseComplement(const std::string& sequence);

        // Static methods
        static double calculateAverageQuality(const std::string& qualityStr);
        
        static void removePrefix(std::string& sequence, 
            std::string& quality, 
            const std::string& prefix) {
            GeneticMaterial::removePrefix(sequence, quality, prefix);}
    
        static void trimPolyAT(std::string& sequence, 
            std::string& quality) {
            GeneticMaterial::trimPolyAT(sequence, quality);}


        std::optional<std::vector<GeneticMaterial::FileContent>> filterFastqByQuality(
            const std::optional<std::vector<GeneticMaterial::FileContent>>& fastqFiles,
            const std::string& qualityStr,
            const std::string& outputDir);

        

        std::optional<std::vector<GeneticMaterial::ProcessedSequence>> processAndCleanSequences(
            const std::optional<std::vector<GeneticMaterial::FileContent>>& loadedFiles,
            const std::optional<std::string>& prefixToRemove,
            bool removePolyAT,
            double minQuality,
            bool removeDegenerate);

        // Methods to build and access the index
        void buildIndex();
        void printSuffixArray() const;
        void printLCPArray() const;
        std::vector<int> findRepeatedFactors(int k) const;
        const std::vector<int>& getSuffixArray() const;

        // Auxiliary methods
        std::string getFactor(int i, int k) const;

        std::string getSequenceFromLoadedData(
            const std::optional<std::pair<
                std::vector<GeneticMaterial::FileContent>,
                std::vector<GeneticMaterial::FileContent>>>& loadedData) const;

        std::vector<int> buildSuffixArray(const std::string& sequence);
        std::vector<int> buildLCPArray(const std::string& sequence, const std::vector<int>& suffixArray);

        void buildAndSaveSuffixTable();

        int searchPattern(const std::string& pattern);

        std::vector<int> findRepeatedFactors(int minLength);

        void mapReads(const std::vector<std::string>& reads, int k);

        void mapReadsAndReport(const std::vector<std::string>& reads, int k);
};

#endif // GENETIC_MATERIAL_CONTROLLER_H
