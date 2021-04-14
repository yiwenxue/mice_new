#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

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
    /* ofstream report("report.pdf"); */
    /* { */
    /*     report << "\\documentclass{article}" << std::endl; */
    /*     report << "\\usepackage{graphicx}" << std::endl; */
    /*     report << "\\usepackage{subfigure}" << std::endl; */
    /*     report << "\\usepackage[top=1cm, bottom=0cm, outer=1cm, inner=1cm, heightrounded, marginparwidth=1cm, marginparsep=1cm]{geometry}" << std::endl; */
    /*     report << "\\geometry{a4paper,scale=0.8}" << std::endl; */
    /*     report << "" << std::endl; */
    /* } */

    system("echo \"Batch,Mice,Type,Days,Act level,Act rsq,Act Std rsq,Act Dfa range,Act Dfa rsq,Temp level,Temp rsq,Temp Dfa range,Temp Dfa rsq,Temp Sync index\" > table.txt");
    while(fscanf(config,"%s %s %s %d %d\n",path,mutant,batch,&start,&length) == 5){
        string str1(path);
        string str2(mutant);
        string str3(batch);
        if(str3 == "batch6 ")
            system(("./main_batch6 "+str1+" "+str2+" "+str3+" "+to_string(start)+" "+to_string(length)).c_str());
        else 
            system(("./main "+str1+" "+str2+" "+str3+" "+to_string(start)+" "+to_string(length)).c_str());

        /* { */
        /*     report << "Activity:" << std::endl; */ 

        /* } */
    }

    /* { */
    /*     report << "\\end{document}" << std::endl; */
    /* } */
    fclose(config);
    /* report.close(); */
    return 0;
}
