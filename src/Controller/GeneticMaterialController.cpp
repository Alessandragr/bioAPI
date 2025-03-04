#include <iostream>
#include <string>
#include <vector>
#include "../Model/GeneticMaterial.cpp"

class GeneticMaterialController {
private:
    GeneticMaterial geneticMaterial; // Instância da model

public:
    // Construtor
    GeneticMaterialController(const std::string& filePath)
        : geneticMaterial(filePath) {}

    // Métodos para interagir com a model
    std::pair<std::vector<GeneticMaterial::FileContent>, std::vector<GeneticMaterial::FileContent>> openFile() {
        return geneticMaterial.openFile(geneticMaterial.getFilePath());
    }

    bool verifySequence(const std::string& content) {
        return geneticMaterial.verifySequence(content);
    }
};