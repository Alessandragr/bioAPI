#include <vector>
#include <optional>
#include <utility>
#include <zlib.h>
#include <archive.h>
#include <archive_entry.h>
#include <minizip/unzip.h>
#include <iostream>
#include <string>
#include <vector>

#include "../Model/GeneticMaterial.h"
#include "../View/MainView.h"
#include "GeneticMaterialController.h"






GeneticMaterialController::GeneticMaterialController(const std::string& filePath)
    : geneticMaterial(filePath), model(geneticMaterial) {}

std::string GeneticMaterialController::getFilePath() const {
    return geneticMaterial.getFilePath();
}

std::string GeneticMaterialController::getReverseComplement(const std::string& sequence) {
    return geneticMaterial.getReverseComplement(sequence);
}


std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> GeneticMaterialController::openFile() {
    return geneticMaterial.openFile(geneticMaterial.getFilePath());
}

std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> GeneticMaterialController::loadFileContent(const std::variant<std::string, std::vector<std::string>>& input) {
    return geneticMaterial.loadFileContent(input);
}

void GeneticMaterialController::loadAndSetFileData(const std::string& filePath) {
    auto fileData = loadFileContent(filePath);
    if (!fileData.has_value()) {
        std::cerr << "Error: Could not load the file content: " << filePath << std::endl;
        throw std::runtime_error("Error loading the file.");
    }

    setLoadedFileData(fileData);
}

bool GeneticMaterialController::verifySequence(const std::string& content, const std::string& fileType) {
    return geneticMaterial.verifySequence(content, fileType);
}

size_t GeneticMaterialController::countSequencesInMultifasta() const {
    if (loadedFileData.has_value()) { // Checks if loadedFileData contains a value
        const auto& fastaFiles = loadedFileData->first; // Accesses the 'first' member directly
        return geneticMaterial.countSequencesInMultifasta(fastaFiles);
    } else {
        throw std::runtime_error("No FASTA file loaded.");
    }
}

std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> 
GeneticMaterialController::calculateComplementSequences(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles) {
    return geneticMaterial.calculateComplementSequences(loadedFiles);
}

std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> GeneticMaterialController::calculateReverseComplementSequences(const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles) {
    return geneticMaterial.calculateReverseComplementSequences(loadedFiles);
}

std::vector<std::string> GeneticMaterialController::getSequenceBetweenStartAndStop(
    const std::optional<std::pair<
        std::optional<std::vector<GeneticMaterial::FileContent>>,
        std::optional<std::vector<GeneticMaterial::FileContent>>>>& fileData,
    const std::string& startCodon,
    const std::vector<std::string>& stopCodons) 
{
    return geneticMaterial.getSequenceBetweenStartAndStop(fileData, startCodon, stopCodons);
}

// Função para chamar o método getSequenceDescriptions da model GeneticMaterial
std::vector<GeneticMaterial::SequenceDescription> GeneticMaterialController::getSequenceDescriptions(
    const std::optional<std::pair<
        std::optional<std::vector<GeneticMaterial::FileContent>>,
        std::optional<std::vector<GeneticMaterial::FileContent>>
    >>& loadedFiles)
{
    return geneticMaterial.getSequenceDescriptions(loadedFiles);
}

// Função que chama o método extractSubsequences da model
std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>> 
GeneticMaterialController::extractSubsequences(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
    size_t start, 
    size_t length, 
    const std::optional<std::string>& sequenceId) {  // Remova o valor padrão aqui
    return geneticMaterial.extractSubsequences(loadedFiles, start, length, sequenceId);
}

// Função que chama o método calculateGCContent da model
std::optional<std::vector<std::pair<std::string, double>>> 
GeneticMaterialController::calculateGCContent(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, std::optional<std::vector<GeneticMaterial::FileContent>>>>& loadedFiles,
    std::optional<std::pair<size_t, size_t>> range) {  // Remova o valor padrão aqui
    return geneticMaterial.calculateGCContent(loadedFiles, range);
}

// Calcula a qualidade média de uma sequência FASTQ
double GeneticMaterialController::calculateAverageQuality(const std::string& qualityStr) {
    return GeneticMaterial::calculateAverageQuality(qualityStr);
}

// Remove um prefixo específico de uma sequência
void removePrefix(std::string& sequence, std::string& quality, const std::string& prefix);

// Remove caudas poly-A/poly-T de uma sequência
static void trimPolyAT(std::string& sequence, std::string& quality);

// Verifica se uma sequência contém caracteres degenerados
static bool containsDegenerateBases(const std::string& sequence);

std::optional<std::vector<GeneticMaterial::FileContent>> GeneticMaterialController::filterFastqByQuality(
    const std::optional<std::vector<GeneticMaterial::FileContent>>& fastqFiles,
    const std::string& qualityStr,
    const std::string& outputDir) {
    return geneticMaterial.filterFastqByQuality(fastqFiles, qualityStr, outputDir);
}

std::optional<std::pair<
    std::vector<GeneticMaterial::FileContent>,
    std::vector<GeneticMaterial::FileContent>>>
GeneticMaterialController::getLoadedFileData() const {
    return loadedFileData;
}

std::optional<std::vector<GeneticMaterial::ProcessedSequence>> GeneticMaterialController::processAndCleanSequences(
    const std::optional<std::vector<GeneticMaterial::FileContent>>& loadedFiles,
    const std::optional<std::string>& prefixToRemove,
    bool removePolyAT,
    double minQuality,
    bool removeDegenerate) {
   
    // Chama a função correspondente no modelo
    return geneticMaterial.processAndCleanSequences(
        loadedFiles,
        prefixToRemove,
        removePolyAT,
        minQuality,
        removeDegenerate
    );
}

void GeneticMaterialController::setLoadedFileData(
    const std::optional<std::pair<std::optional<std::vector<GeneticMaterial::FileContent>>, 
                                  std::optional<std::vector<GeneticMaterial::FileContent>>>>& fileData) {
        if (!fileData.has_value() || !fileData->first.has_value()) {
            throw std::runtime_error("No data loaded to configure.");
        }

    // Transformar fileData no tipo correto
    std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> transformedData;

    // Extrair e transformar o membro 'first'
    if (fileData->first.has_value()) {
        transformedData.first = fileData->first.value();
    }

    // Extrair e transformar o membro 'second'
    if (fileData->second.has_value()) {
        transformedData.second = fileData->second.value();
    }

    // Atribuir ao loadedFileData
    loadedFileData = transformedData;

    // Concatenar todas as linhas de sequência do arquivo FASTA
    std::string concatenatedSequence;
    for (const auto& line : transformedData.first) {
        if (line.content.has_value()) {
            concatenatedSequence += line.content.value();
        }
    }

    // Configurar a sequência concatenada na model
    geneticMaterial.setSequence(concatenatedSequence);
}



const std::vector<int>& GeneticMaterialController::getSuffixArray() const {
    return geneticMaterial.getSuffixArray();
}

std::vector<int> GeneticMaterialController::buildSuffixArray(const std::string& sequence) {
    return geneticMaterial.buildSuffixArray(sequence);
}

std::vector<int> GeneticMaterialController::buildLCPArray(const std::string& sequence, const std::vector<int>& suffixArray) {
    return geneticMaterial.buildLCPArray(sequence, suffixArray);
}

// Chama o método da model para construir as tabelas
void GeneticMaterialController::buildIndex() {
    geneticMaterial.buildIndex(); 
}

void GeneticMaterialController::printSuffixArray() const {
    const auto& suffixArray = geneticMaterial.getSuffixArray();
    const auto& sequence = geneticMaterial.getSequence();

    std::cout << "Suffix Table:\n";
    for (size_t i = 0; i < suffixArray.size(); ++i) {
        std::cout << i << ": " << sequence.substr(suffixArray[i]) << "\n";
    }
}

void GeneticMaterialController::printLCPArray() const {
    const auto& lcpArray = geneticMaterial.getLCPArray();

    std::cout << "Longest Common Prefix (LCP) Table:\n";
    for (size_t i = 1; i < lcpArray.size(); ++i) {
        std::cout << "LCP[" << i << "]: " << lcpArray[i] << "\n";
    }
}

std::vector<int> GeneticMaterialController::findRepeatedFactors(int k) const {
    return geneticMaterial.findRepeatedFactors(k);
}

std::string GeneticMaterialController::getFactor(int i, int k) const {
    return geneticMaterial.getFactor(i, k);
}

std::string GeneticMaterialController::getSequenceFromLoadedData(
    const std::optional<std::pair<
        std::vector<GeneticMaterial::FileContent>,
        std::vector<GeneticMaterial::FileContent>>>& loadedData) const {
    return geneticMaterial.extractSequenceFromLoadedData(loadedData);
}

void GeneticMaterialController::buildAndSaveSuffixTable() {
    geneticMaterial.buildAndSaveSuffixTable();
}

int GeneticMaterialController::searchPattern(const std::string& pattern) {
    return geneticMaterial.searchPattern(pattern);
}


std::vector<int> GeneticMaterialController::findRepeatedFactors(int minLength) {
    return geneticMaterial.findRepeatedFactors(minLength);
}


void GeneticMaterialController::mapReads(const std::vector<std::string>& reads, int k) {
    geneticMaterial.mapReadsAndReport(reads, k);
}


