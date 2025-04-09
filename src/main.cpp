#include "View/MainView.h"  // Note: Inclua o .h, não o .cpp
#include "Controller/GeneticMaterialController.h"
#include <locale>
#include <clocale>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

void configureSystemLocale() {
    #ifdef _WIN32
    // Configuração específica para Windows
    SetConsoleOutputCP(65001);  // UTF-8
    std::system("chcp 65001 > nul");  // Altera a codificação do console
    #else
    // Configuração para Linux/Mac
    std::setlocale(LC_ALL, "en_US.UTF-8");
    std::locale::global(std::locale("en_US.UTF-8"));
    #endif
}

int main(int argc, char* argv[]) {
    configureSystemLocale();
    
    try {
        // Processa os comandos e executa a lógica principal
        MainView::processCommandLine(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}