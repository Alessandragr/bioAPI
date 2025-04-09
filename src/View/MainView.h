#ifndef MAIN_VIEW_H
#define MAIN_VIEW_H

#include <string>
#include <vector>

class MainView {
public:
    static void displayMessage(const std::string& message);
    static void showUsage();
    static void processCommandLine(int argc, char* argv[]);

private:
    #ifdef _WIN32
    static void initWindowsConsole();
    #endif
};

#endif // MAIN_VIEW_H