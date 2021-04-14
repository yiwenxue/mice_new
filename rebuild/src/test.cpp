/* #include <TPDF.h> */
#include <iostream>
#include <mathematics.h>
#include <utility.h>
#include <cstdio>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>

int main(int argc, char **argv)
{
    if (argc < 2){
        std::cout << "Usage: " << std::endl;
        std::cout << "      <" << argv[0] << "> [file]" << std::endl;
        exit(-1);
    }

    int lineNum = fileLines(std::string(argv[1]));
    std::cout << "line Num: " << lineNum << std::endl;
    FILE * readData = fopen(argv[1], "r");

    auto data = new double [lineNum];

    int i = 0; double temp = 0;
    while ((fscanf(readData, "%lf\n", &temp) == 1) && i < lineNum) {
        data[i++] = temp;
    }
    fclose(readData);
    std::cout << "Read data: " << i << std::endl;

    auto newdfa = new DFA(data, i, 2);
    std::cout << "dfa alpha: " << newdfa->index << std::endl;

    delete[] data;
    delete newdfa;
    return 0;
}
