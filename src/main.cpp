#include <cstdlib>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include "utility.h"

#define MIN(x,y) (((x)>(y))?(y):(x))

using namespace std;

int main(int argc, char *argv[])
{
    int seg_size = 90;
    std::string mutant(argv[2]);
    std::string batch(argv[3]);
    bool ifmutant = (mutant=="mutant")?(true):(false);
    std::string name(argv[1]);
    int start = atoi(argv[4]);
    int length = atoi(argv[5]);

    mice mice1(name,seg_size,start,length);

    mice1.print_details();

    mice1.plot_all();

    string new_mice;
    if(mice1.type == "Activity")
        new_mice = mice1.path + mice1.mice_name + "."
            + "Temperature.txt";
    else if(mice1.type == "Temperature")
        new_mice = mice1.path + mice1.mice_name + "."
            + "Activity.txt";

    mice mice2(new_mice,seg_size,start,length);
    mice2.print_details();
    mice2.plot_all();

    Sync new_sync(mice1.average,mice2.average,(48*360)/seg_size );
    new_sync.plot(mice1.mice_name, seg_size);
    std::cout << std::endl;
    mktable(&mice1,&mice2,new_sync.sync_index, batch, ifmutant);

    std::cout << "shutdown" << std::endl;
    return 0;
}
