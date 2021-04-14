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
    char mouse[10];
    char batch[10];
    int start, length;

    /* ofstream report("report.tex"); */
    ofstream report("dfa_compare.tex");
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

    while(fscanf(config,"%s %s %s %s %d %d %*d %*d\n",path,mouse,mutant,batch,&start,&length) == 6){
        string str1(path);
        string str2(mutant);
        string name(mouse);

        report << "\\section{" + name + "\\,\\,\\,\\,\\," + "batch " << 
            batch << "\\,\\,\\," <<  (str2 == "mutant"?"mutant":"WT") << "}" << std::endl;
        /* name = "../outcome/" + name ; */
        name = "../" + name ;
        cout << name << endl;
        report << "\\begin{figure}[htp]" << std::endl;
        report << "\\centering" << endl;
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_Activity_heatmap" ".pdf}}" << std::endl; */
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_Temperature_heatmap" ".pdf}}" << std::endl; */
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_activity_DFA2" ".pdf}}" << std::endl; */
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_temperature_DFA2" ".pdf}}" << std::endl; */
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_period" ".pdf}}" << std::endl; */
        /* report << "\\subfigure{\\includegraphics[width=0.45\\textwidth]{" + name +"_DFA3" ".pdf}}" << std::endl; */
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_activity_DFA1" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_activity_DFA2" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_activity_DFA3" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_activity_DFA4" ".pdf}}" << std::endl;
        /* report << "\\end{figure}" << std::endl; */
        /* report << "\\newpage " << endl; */

        /* report << "\\begin{figure}[htp]" << std::endl; */
        /* report << "\\centering" << endl; */
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_temperature_DFA1" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_temperature_DFA2" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_temperature_DFA3" ".pdf}}" << std::endl;
        report << "\\subfigure{\\includegraphics[width=0.36\\textwidth]{" + name +"_temperature_DFA4" ".pdf}}" << std::endl;
        report << "\\end{figure}" << std::endl;
        report << "\\newpage " << endl;
    }

    {

        report << "\\end{document}" << std::endl;
        report.close();
        /* system(("pdflatex report.tex")); */
    }
    return 0;
}
