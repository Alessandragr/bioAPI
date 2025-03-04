#include <iostream>
#include <string>

class MainView {
public:
    // Exibe mensagens para o usuário
    static void displayMessage(const std::string& message) {
        std::cout << message << std::endl;
    }

    // Solicita o caminho do arquivo ao usuário
    static std::string getFilePath() {
        std::string filePath;
        std::cout << "Digite o caminho do arquivo: ";
        std::getline(std::cin, filePath);
        return filePath;
    }
};