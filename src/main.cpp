#include "Controller/GeneticMaterialController.cpp"
#include "View/MainView.cpp"

int main(int argc, char* argv[]) {
    std::string filePath;

    // Verifica se o caminho do arquivo foi passado como argumento
    if (argc < 2) {
        filePath = MainView::getFilePath(); // Solicita o caminho do arquivo ao usuário
    } else {
        filePath = argv[1]; // Usa o caminho passado como argumento
    }

    try {
        GeneticMaterialController controller(filePath);

        auto [fastaFiles, fastqFiles] = controller.openFile();

        // Verifica as sequências FASTA
        for (const auto& file : fastaFiles) {
            if (file.content.has_value()) {
                bool isValid = controller.verifySequence(file.content.value());
                MainView::displayMessage("O arquivo " + file.fileName + " contem uma sequencia FASTA " +
                                        (isValid ? "valida." : "invalida."));
            }
        }

        // Verifica as sequências FASTQ
        for (const auto& file : fastqFiles) {
            if (file.content.has_value()) {
                bool isValid = controller.verifySequence(file.content.value());
                MainView::displayMessage("O arquivo " + file.fileName + " contem uma sequencia FASTQ " +
                                        (isValid ? "valida." : "invalida."));
            }
        }

    } catch (const std::exception& e) {
        MainView::displayMessage("Erro: " + std::string(e.what()));
        return 1;
    }

    return 0;
}