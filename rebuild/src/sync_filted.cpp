#include <TGraph.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_long.h>
#include <iostream>
#include <iterator>
#include <ostream>
#include <random>
#include <string>
#include <vector>
#include "utility.h"

#include "sync_util.h"

#define MIN(x,y) (((x)>(y))?(y):(x))

using namespace std;

TimeSeq *readTimeseq(std::string infile);
void freeMeta();
int readMeta(std::string metaFile);
void Times_dump(std::string name, TimeSeq data);

using namespace std;

int main(int argc, char **argv)
{ 
    gErrorIgnoreLevel = kWarning;
    string select = string("");
    if (argc == 2 ){
        path_mask = string(argv[1]);
    } else if (argc == 3){
        path_mask = string(argv[1]);
        select = string(argv[2]);
    }
    cout << "path_mask: " + path_mask << endl;

    readMeta(path_mask + "./mice");
    filted = 0;
    phase_output = 1;

    int surro = 0;

    if (surro == 0)
    for (auto &i:metaData){

        if (select != "" && i.mice.name != select)
            continue;
        cout << "------------------------" << endl;
        cout << "  mice: " + i.mice.name << endl;
        cout << "  path: " + i.mice.path << endl;
        cout << "  n,m : " << i.n << "," << i.m << endl;

        /* analysis(&i); */
        analysis_ori(&i);

        cout << endl;
    }

    if (surro == 1)
    for (auto &i:metaData){
        for (auto &j:metaData){
            if (i.mice.name == j.mice.name)
                continue;
            if (i.mice.batch != j.mice.batch)
                continue;
            cout << "------------------------" << endl;
            cout << "  mice: " + i.mice.name <<"_"<< j.mice.name << endl;
            cout << "  path: " + i.mice.path << endl;
            cout << "  n,m : " << 1 << "," << 1 << endl;

            surrogate(&i, &j);
        }
    }

    printf("PSI table:\n");
    printf("%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "mouse", "rho", "lambda", "gamma", "W", "t_lag");
    for (auto i:table){
        std::cout << i << std::endl;
    }

    freeMeta();
    metaData.clear();
    return 0;
}
