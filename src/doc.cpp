#include <iostream>
#include <ostream>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc < 2){
        return -1;
    }
    FILE *config = fopen(argv[1],"r");
    char path[255];
    char mutant[10];
    char batch[10];
    int start, length;

    ofstream report("report.tex");
    {
        report << "\\documentclass{article}" << std::endl;
        report << "\\usepackage{graphicx}" << std::endl;
        report << "\\usepackage{subfigure}" << std::endl;
        report << "\\usepackage{rotating}" << std::endl;
        report << "\\usepackage{hyperref}" << endl;
        report << "\\hypersetup{colorlinks=true, linktoc=all, linkcolor=blue,urlcolor=red}" << endl;
        /* report << "\\usepackage[top=1cm, bottom=1.5cm, outer=1cm, inner=1cm, heightrounded, marginparwidth=1cm, marginparsep=1cm]{geometry}" << std::endl; */
        report << "\\usepackage{geometry}" << endl;
        report << "\\geometry{a4paper,scale=0.9}" << std::endl;
        report << "\n\n\\begin{document}" << std::endl;
        report << "\\tableofcontents" << endl;
        report << "\\newpage" << endl;
    }

    int mark;
    while(fscanf(config,"%s %s %s %d %d\n",path,mutant,batch,&start,&length) == 5){
        string str1(path);
        string str2(mutant);
        mark = str1.find_last_of("/");
        cout << mark <<"\n"<< endl;
        string str4 = str1.substr(mark+1);
        string name =  str4.substr(0,str4.find_first_of("."));
        report << "\\section{" + name +"\t\\,\\,\\,\\,\\,"+ (str2=="mutant"?"mutant":"wild") + "}\n" << endl;
        report << "Activity" << endl;
        report << "\\begin{figure}[htp]" << std::endl;
        report << "\\centering" << endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_heatmap" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_overview" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_powerspec" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_rhythm" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_std_rhythm" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Activity_dfa_rhythm" ".pdf}}" << std::endl;
        report << "\\end{figure}" << std::endl;
        report << "\\\\" << endl;
        report << "\\newpage " << endl;

        report << "Temperature" << endl;
        report << "\\begin{figure}[htp]" << std::endl;
        report << "\\centering" << endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Temperature_heatmap" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Temperature_overview" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Temperature_powerspec" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Temperature_rhythm" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Temperature_dfa_rhythm" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.42\\textwidth]{" + name +"_Sync" ".pdf}}" << std::endl;
        report << "\\end{figure}" << std::endl;
        report << "\\\\" << endl;
        report << "\\newpage " << endl;
    }

    {

        report << "\\begin{sidewaystable}" << endl;
        report << "%%Here will be the final table \n \\\\" << endl;
        report << "\\end{sidewaystable}" << endl;
        report << "\\end{document}" << std::endl;
        report.close();
        /* system(("pdflatex report.tex")); */
    }
    return 0;
}
