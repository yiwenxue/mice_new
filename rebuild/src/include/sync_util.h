
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_statistics_double.h>
#include <time.h>

using namespace std;
struct meta {
    TMicemem mice;

    int start;
    int length;
    int n;
    int m;

    TimeSeq *x; 
    TimeSeq *y;
};


int filted = 1;
int phase_output = 0;
char buffer[255];
std::vector<struct meta> metaData;

std::string path_mask("./");
std::vector<std::string> table;

int readMeta(std::string metaFile){
    FILE *config = fopen(metaFile.c_str(), "r");

    char _name[255];
    char _path[255];
    char _mutant[10];
    int batch;
    int n,m;
    int start, length;

    TMicemem mice;
    struct meta readmeta;
    while(fscanf(config,"%s %s %s %d %d %d %d %d\n",_path,_name,_mutant,&batch,&start,&length,&n,&m) == 8){
        std::string path(_path);
        path = path_mask + path;
        std::string mutant(_mutant);
        std::string name(_name);

        readmeta.mice.name = name;
        readmeta.mice.path = path;
        readmeta.mice.batch = batch;
        readmeta.mice.ifmutant = (mutant == "mutant"? true:false);
        readmeta.start = start; 
        readmeta.length = length;
        readmeta.n = n;
        readmeta.m = m;

        readmeta.x = NULL;
        readmeta.y = NULL;

        metaData.push_back(readmeta);
    }
    return 0;
}

TimeSeq *readTimeseq(std::string infile){
    TimeSeq *indata = (TimeSeq *) malloc (sizeof(TimeSeq));
    indata->size = fileLines(infile);
    indata->time = (double *) malloc (sizeof(double) * indata->size);
    indata->data = (double *) malloc (sizeof(double) * indata->size);

    FILE *fp =  fopen(infile.c_str(), "r");

    double x, y;
    for (int i=0; i<indata->size; i++){
        fscanf(fp, "%lf\t%lf\n", &x, &y);
        indata->time[i] = x;
        indata->data[i] = y;
    }
    fclose(fp);
    return indata;
}


void freeMeta()
{
    for (auto const&i :metaData){
        if (i.x != NULL){
            free(i.x->time);
            free(i.x->data);
            free(i.x);
            free(i.y->time);
            free(i.y->data);
            free(i.y);
        }
    }
}

void Times_dump(std::string name, TimeSeq data){
    int size = data.size;
    FILE *fp = fopen(name.c_str(), "w");

    if (fp == NULL){
        perror("Err: On opening file.");
        exit(-1);
    }

    for (int i = 0; i < size; ++i) {
        fprintf(fp, "%lf\t%lf\n", data.time[i], data.data[i]);
    }

    fclose(fp);
}

TimeSeq *TimeSeq_MA(TimeSeq input, int winSize, int overLap){

    int steps = (input.size - winSize)/overLap;
    TimeSeq *output = (TimeSeq *)malloc(sizeof(TimeSeq));
    output->size= steps;
    output->time = new double[steps];
    output->data = new double[steps];

    for (int i = 0; i < steps; ++i) {
        output->time[i] = gsl_stats_mean(input.time+i*overLap, 1, winSize);
        output->data[i] = gsl_stats_mean(input.data+i*overLap, 1, winSize);
    }

    return output;
}

void surrogate(struct meta *member1, struct meta *member2){
    if (!filted){
        if (member1->x == NULL) 
            member1->x = readTimeseq(member1->mice.path+"/selected_data/"+member1->mice.name+".Activity.txt");
        if (member2->y == NULL) 
            member2->y = readTimeseq(member2->mice.path+"/selected_data/"+member2->mice.name+".Temperature.txt");
    } else {
        if (member1->x == NULL) 
            member1->x = readTimeseq(member1->mice.path+"/bandpass/"+member1->mice.name+"_act_conv.txt");
        if (member2->y == NULL) 
            member2->y = readTimeseq(member2->mice.path+"/bandpass/"+member2->mice.name+"_temp_conv.txt");
    }
    int winSize = 90;
    int winNum = 0;
    winNum = member1->x->size / winSize;
    TimeSeq seq1 = {
        .size = winNum,
        .data = (double *)malloc(sizeof(double) * winNum),
        .time = (double *)malloc(sizeof(double) * winNum),
    };
    for (int i = 0; i< winNum; i++){
        seq1.data[i] = gsl_stats_mean(member1->x->data+i*winSize, 1, winSize);
        seq1.time[i] = gsl_stats_mean(member1->x->time+i*winSize, 1, winSize);
    }

    winNum = member2->y->size / winSize;
    TimeSeq seq2 = {
        .size = winNum,
        .data = (double *)malloc(sizeof(double) * winNum),
        .time = (double *)malloc(sizeof(double) * winNum),
    };
    for (int i = 0; i< winNum; i++){
        seq2.data[i] = gsl_stats_mean(member2->y->data+i*winSize, 1, winSize);
        seq2.time[i] = gsl_stats_mean(member2->y->time+i*winSize, 1, winSize);
    }

    Sync theSync(seq1.data, seq2.data, min(seq1.size, seq2.size), 1, 1);

    /* theSync.phasePlot(member->mice.name); */
    /* data_overview(member->mice.name + "_dataOv.pdf", seq1, seq2, 0, 1000); */

    theSync.PSI_rho(0);
    theSync.PSI_gamma(0);
    theSync.PSI_lambda(0);
    double rho = theSync.rho_index;
    double gamma = theSync.gamma_index;
    double lambda = theSync.lambda_index;

    int day = 15;

    theSync.PSI_decay(member1->mice.name+"_"+member2->mice.name+"_suro", 4*24*day, 1, "rho");
    theSync.PSI_decay(member1->mice.name+"_"+member2->mice.name+"_suro", 4*24*day, 1, "lambda");
    theSync.PSI_decay(member1->mice.name+"_"+member2->mice.name+"_suro", 4*24*day, 1, "gamma");

    sprintf(buffer, "%10s\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf", (member1->mice.name+"_"+member2->mice.name).c_str(), rho, lambda, gamma, theSync.decay_radio, theSync.decay_max_time);
    table.push_back(string(buffer));

    free(seq1.data);
    free(seq1.time);
    free(seq2.data);
    free(seq2.time);
}



void analysis_ori(struct meta *member){
    if (!filted){
        member->x = readTimeseq(member->mice.path+"/selected_data/"+member->mice.name+".Activity.txt");
        member->y = readTimeseq(member->mice.path+"/selected_data/"+member->mice.name+".Temperature.txt");
    } else {
        member->x = readTimeseq(member->mice.path+"/bandpass/"+member->mice.name+"_act_conv.txt");
        member->y = readTimeseq(member->mice.path+"/bandpass/"+member->mice.name+"_temp_conv.txt");
    }

    int winSize = 1440;
    int overLap = 10;

    int start    = 36*24*0;
    int duration = 36 * 24 * 10;

    clock_t time = clock();
    TimeSeq *act = TimeSeq_MA(*(member->x), winSize, overLap);
    TimeSeq *temp = TimeSeq_MA(*(member->y), winSize, overLap);
    time = clock() - time;
    printf("Time used when moving average: %ldms\n", time/1000);
    /* start = act->size - duration; */

    int dayu = (int)(act->size /24.0/36.0);
    /* printf("Size: %lf day. Days used: %d\n", act->size/24.0/36.0, dayu); */

    Sync theSync(act->data+start, temp->data+start, dayu*24*36 /*act->size*/, member->n, member->m);
    if (phase_output){
        FILE *fp = fopen((member->mice.name+"_phase.txt").c_str(), "w");
        if (fp == NULL){
            perror("Opening phase data file");
            exit(-1);
        }
        for (int i = 0; i < theSync.size; ++i) {
            fprintf(fp, "%lf\t%lf\n", theSync.phasex[i], theSync.phasey[i]);
        }
        fclose(fp);}

    theSync.phasePlot(member->mice.name);
    data_overview(member->mice.name + "_dataOv.pdf", *act, *temp, start, start+duration);

    theSync.PSI_rho(0);
    theSync.PSI_gamma(0);
    theSync.PSI_lambda(0);
    double rho = theSync.rho_index;
    double gamma = theSync.gamma_index;
    double lambda = theSync.lambda_index;

    int day = 5;
    theSync.PSI_decay(member->mice.name, 36*24*day, 10, "rho");
    theSync.PSI_decay(member->mice.name, 36*24*day, 10, "lambda");
    theSync.PSI_decay(member->mice.name, 36*24*day, 10, "gamma");

    sprintf(buffer, "%10s\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf", member->mice.name.c_str(), rho, lambda, gamma, theSync.decay_radio, theSync.decay_max_time);
    table.push_back(string(buffer));

    free(act->data);
    free(act->time);
    free(temp->data);
    free(temp->time);
    free(act);
    free(temp);
}

void analysis(struct meta *member){
    if (!filted){
        member->x = readTimeseq(member->mice.path+"/selected_data/"+member->mice.name+".Activity.txt");
        member->y = readTimeseq(member->mice.path+"/selected_data/"+member->mice.name+".Temperature.txt");
    } else {
        member->x = readTimeseq(member->mice.path+"/bandpass/"+member->mice.name+"_act_conv.txt");
        member->y = readTimeseq(member->mice.path+"/bandpass/"+member->mice.name+"_temp_conv.txt");
    }

    int winSize = 90;
    int winNum = 0;
    winNum = member->x->size / winSize;
    TimeSeq seq1 = {
        .size = winNum,
        .data = (double *)malloc(sizeof(double) * winNum),
        .time = (double *)malloc(sizeof(double) * winNum),
    };
    for (int i = 0; i< winNum; i++){
        seq1.data[i] = gsl_stats_mean(member->x->data+i*winSize, 1, winSize);
        seq1.time[i] = gsl_stats_mean(member->x->time+i*winSize, 1, winSize);
    }

    winNum = member->y->size / winSize;
    TimeSeq seq2 = {
        .size = winNum,
        .data = (double *)malloc(sizeof(double) * winNum),
        .time = (double *)malloc(sizeof(double) * winNum),
    };
    for (int i = 0; i< winNum; i++){
        seq2.data[i] = gsl_stats_mean(member->y->data+i*winSize, 1, winSize);
        seq2.time[i] = gsl_stats_mean(member->y->time+i*winSize, 1, winSize);
    }

    /* Sync theSync(member->x->data, member->y->data,  60, member->n, member->m); */
    Sync theSync(seq1.data, seq2.data, min(seq1.size, seq2.size), member->n, member->m);

    /* Times_dump(member->mice.name+"_phase_act.txt", theSync.phasex); */
    /* Times_dump(member->mice.name+"_phase_temp.txt", theSync.phasex); */

    if (phase_output){
        FILE *fp = fopen((member->mice.name+"_phase.txt").c_str(), "w");
        if (fp == NULL){
            perror("Opening phase data file");
            exit(-1);
        }
        for (int i = 0; i < theSync.size; ++i) {
            fprintf(fp, "%lf\t%lf\n", theSync.phasex[i], theSync.phasey[i]);
        }
        fclose(fp);}

    theSync.phasePlot(member->mice.name);
    data_overview(member->mice.name + "_dataOv.pdf", seq1, seq2, 0, 360 * 24 );

    theSync.PSI_rho(0);
    theSync.PSI_gamma(0);
    theSync.PSI_lambda(0);
    double rho = theSync.rho_index;
    double gamma = theSync.gamma_index;
    double lambda = theSync.lambda_index;

    int day = 15;

    theSync.PSI_decay(member->mice.name, 4*24*day, 1, "rho");
    theSync.PSI_decay(member->mice.name, 4*24*day, 1, "lambda");
    theSync.PSI_decay(member->mice.name, 4*24*day, 1, "gamma");

    sprintf(buffer, "%10s\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf", member->mice.name.c_str(), rho, lambda, gamma, theSync.decay_radio, theSync.decay_max_time);
    table.push_back(string(buffer));

    free(seq1.data);
    free(seq1.time);
    free(seq2.data);
    free(seq2.time);
}
