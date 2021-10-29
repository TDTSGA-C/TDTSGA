#ifndef CSTCHANGE_CLASSDEFINE_H

#include "math.h"
#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>

#define CSTCHANGE_CLASSDEFINE_H
using namespace std;
//{file}
class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

//{task}
class Task {
public:
    double length;          //the length of task
    vector<int> ElgRsc;     //the set of resources which are eligible to  perform this task
    vector<int> parents;    //the set of parent tasks
    vector<int> children;   //the set of child tasks
    vector<vfile> IFile;    //the set of input files
    vector<vfile> OFile;    //the set of output files
};

//{resources}
class Resource {
public:
    vector<int> ElgTsk;     //the set of tasks which it is eligible to perform
    double pc, bw;          //processing capacity, bandwidth

    Resource(int id, double pc, double bw) {
        this->bw = bw;
        this->pc = pc;
    }
};

class chromosome {
public:
    vector<int> RscAlcLst;       //resources allocation, task-to-resource mapping, (Match String)
    vector<int> TskSchLst;       //task scheduling order (Scheduling String)
    vector<double> Code_RK;
    vector<list<int>> Code_TD;
    vector<double> EndTime;      //the finish time of task
    double FitnessValue;         //Fitness (makespan)
    //{sort chromosome and remove the same chromosome according to fitness value}
    bool operator<(const chromosome &otherChromosome)const {
        return this->FitnessValue + 1e-6 < otherChromosome.FitnessValue;
    }
//    //{sort chromosome according to fitness value and remove the same chromosome according to code }
//    bool operator<(const chromosome &otherChromosome)const {
//        int flag = -1;
//        for (int i = 0; i < this->TaskOrderList.size(); i++) {
//            if (this->VMAllocationList[i] != otherChromosome.VMAllocationList[i]) {
//                flag = i;
//                break;
//            }
//            if (this->TaskOrderList[i] != otherChromosome.TaskOrderList[i]) {
//                flag = i;
//                break;
//            }
//        }
//        if (flag == -1) {
//            return false;
//        } else {
//            if(fabs(this->FitnessValue-otherChromosome.FitnessValue)<1e-6) {
//                if(this->VMAllocationList[flag] != otherChromosome.VMAllocationList[flag]){
//                    return this->VMAllocationList[flag]<otherChromosome.VMAllocationList[flag];
//                }else {
//                    return this->TaskOrderList[flag]<otherChromosome.TaskOrderList[flag];
//                }
//            }else {
//                return this->FitnessValue < otherChromosome.FitnessValue;
//            }
//        }
//    }
};

class Paramet_TDTSGA {
public:
    int NumOfChormPerPop;
    double beta;
    double MutationStg1Rate;
    double MutationStg2Rate;
    double RunTimeRatioOfStg1;
};

class ComConst {
public:
    int NumOfTsk;             //the number of Tasks
    int NumOfRsc;             //the number of resources
};

class Orthogonal{
public:
    double NumOfChormPerPopRate;
    double beta;
    double Mtn1;
    double Mtn2;
    double RunTimeRatioOfStg1Rate;
};

extern vector<Task> Tasks;
extern vector<vector<int> > TskLstInLvl;
extern vector<vector<double> > ParChildTranFileSizeSum;
extern vector<int> LevelIdOfTask;
extern vector<Resource> Rscs;
extern vector<chromosome> population;
extern vector<vector<chromosome> > populations;
extern Paramet_TDTSGA Parameter_TDTSGA;
extern ComConst comConst;
extern vector<Orthogonal> orthogonal;
#endif //CSTCHANGE_CLASSDEFINE_H