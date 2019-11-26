#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdio>

class GnuplotPipe {
public:
    inline GnuplotPipe(bool persist = true) {
        pipe = popen(persist ? "gnuplot -persist" : "gnuplot", "w");
        if (!pipe)
            std::cout << "Open gnuplot failed!" << std::endl;
    }
    inline virtual ~GnuplotPipe(){
        if (pipe) pclose(pipe);
    }

    void sendLine(const std::string& text, bool useBuffer = false){
        if (!pipe) return;
        if (useBuffer)
            buffer.push_back(text + "\n");
        else
            fputs((text + "\n").c_str(), pipe);
    }
    void sendEndOfData(unsigned repeatBuffer = 1){
        if (!pipe) return;
        for (unsigned i = 0; i < repeatBuffer; i++) {
            for (auto& line : buffer) fputs(line.c_str(), pipe);
            fputs("e\n", pipe);
        }
        fflush(pipe);
        buffer.clear();
    }
    void sendNewDataBlock(){
        sendLine("\n", !buffer.empty());
    }

    void writeBufferToFile(const std::string& fileName){
        std::ofstream fileOut(fileName);
        for (auto& line : buffer) fileOut << line;
        fileOut.close();
    }

private:
    GnuplotPipe(GnuplotPipe const&) = delete;
    void operator=(GnuplotPipe const&) = delete;

    FILE* pipe;
    std::vector<std::string> buffer;
};
