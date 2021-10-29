#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "common.h"

using namespace std;

void GnrTskSchLst(chromosome& chrom , vector<int>& TskSchLst) {
    chromosome temch;
    IntChr(temch);
    temch = chrom;
    vector<int > upr(comConst.NumOfTsk,0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            if(!temch.Code_TD[j].empty() && upr[temch.Code_TD[j].front()] == 0 ){
                TskSchLst.push_back(temch.Code_TD[j].front());
                for (int l = 0; l < Tasks[temch.Code_TD[j].front()].children.size(); ++l) {
                    int childId = Tasks[temch.Code_TD[j].front()].children[l];
                    upr[childId] = upr[childId] - 1;
                }
                temch.Code_TD[j].erase(temch.Code_TD[j].begin());
                break;
            }
        }
    }
}

double DcdEvl(chromosome& ch, bool IsFrw) {
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        double StartTime = 0;
        if(IsFrw) {                              //forward-loading
            if (Tasks[TaskIndex].parents.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                    int ParentTask = Tasks[TaskIndex].parents[j];
                    int ParentRsc = ch.RscAlcLst[ParentTask];
                    double TransferTime = 0;
                    if(RscIndex != ParentRsc) {
                        TransferTime = ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                    }
                    double sum = ch.EndTime[ParentTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        } else {                                //backward-loading
            if (Tasks[TaskIndex].children.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                    int ChildTask = Tasks[TaskIndex].children[j];
                    int ChildRsc = ch.RscAlcLst[ChildTask];
                    double TransferTime = 0;
                    if(RscIndex != ChildRsc) {
                        TransferTime = ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    }
                    double sum = ch.EndTime[ChildTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        }
        set<double>::iterator pre  = ITL[RscIndex].begin();
        set<double>::iterator post = ITL[RscIndex].begin();
        ++post;
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        //{find an idle time-slot in ITL which can finish the task  at the earliest}
        while(post != ITL[RscIndex].end()) {
            if((*post - *pre) >= ExecutionTime && ReadyTime <= (*post)-ExecutionTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        ch.EndTime[TaskIndex] = StartTime + ExecutionTime;
        if (makespan < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(StartTime);
        } else {
            ITL[RscIndex].insert(StartTime);
        }

        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}

double Dcd(chromosome& chrom, bool IsFrw) {
    chromosome temch = chrom;
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0);
    if(IsFrw){
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].parents.size();
        }
    } else{
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i] = Tasks[i].children.size();
        }
    }

    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int taskIndex = -1;
        int RscIndex = -1;
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            if(!temch.Code_TD[j].empty() && upr[temch.Code_TD[j].front()] == 0){
                taskIndex = temch.Code_TD[j].front();
                RscIndex = j;
                chrom.RscAlcLst[taskIndex] = j;
                break;
            }
        }
        double ReadyTime = 0;
        double StartTime = 0;
        temch.Code_TD[RscIndex].erase(temch.Code_TD[RscIndex].begin());
        if(IsFrw) {                              //forward-loading
            for (int l = 0; l < Tasks[taskIndex].children.size(); ++l) {
                int childId = Tasks[taskIndex].children[l];
                upr[childId] = upr[childId] - 1;
            }
            if (Tasks[taskIndex].parents.size() != 0) {
                for (int j = 0; j < Tasks[taskIndex].parents.size(); ++j) {
                    int ParentTask = Tasks[taskIndex].parents[j];
                    int ParentRsc = chrom.RscAlcLst[ParentTask];
                    double TransferTime = 0;
                    if(RscIndex != ParentRsc) {
                        TransferTime = ParChildTranFileSizeSum[ParentTask][taskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                    }
                    double sum = chrom.EndTime[ParentTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        } else {                                //backward-loading
            for (int l = 0; l < Tasks[taskIndex].parents.size(); ++l) {
                int parentId = Tasks[taskIndex].parents[l];
                upr[parentId] = upr[parentId] - 1;
            }
            if (Tasks[taskIndex].children.size() != 0) {
                for (int j = 0; j < Tasks[taskIndex].children.size(); ++j) {
                    int ChildTask = Tasks[taskIndex].children[j];
                    int ChildRsc = chrom.RscAlcLst[ChildTask];
                    double TransferTime = 0;
                    if(RscIndex != ChildRsc) {
                        TransferTime = ParChildTranFileSizeSum[taskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    }
                    double sum = chrom.EndTime[ChildTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        }
        set<double>::iterator pre  = ITL[RscIndex].begin();
        set<double>::iterator post = ITL[RscIndex].begin();
        ++post;
        double ExecutionTime = Tasks[taskIndex].length / Rscs[RscIndex].pc;
        //{find an idle time-slot in ITL which can finish the task  at the earliest}
        while(post != ITL[RscIndex].end()) {
            if((*post - *pre) >= ExecutionTime && ReadyTime <= (*post)-ExecutionTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        chrom.EndTime[taskIndex] = StartTime + ExecutionTime;
        if (makespan < chrom.EndTime[taskIndex]) {
            makespan = chrom.EndTime[taskIndex];
        }
        //{update ITL}
        UpdateITL(ITL,RscIndex,StartTime,chrom.EndTime[taskIndex]);
    }
    chrom.FitnessValue = makespan;
    return chrom.FitnessValue;
}

//{calculate the cumulative probabilities for the population whose chromosome have been sorted}
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A , int& NumOfChormPerPop) {
    for (int n = 0; n < NumOfChormPerPop; ++n) {
        A[n] = pow(RtOfSltPrb,NumOfChormPerPop-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, NumOfChormPerPop) - 1);
    }
    for (int n = 1; n < NumOfChormPerPop; ++n){
        A[n] = A[n] + A[n - 1];
    }
}

//{select a chromosome using roulette wheel selection scheme}
int SltChr(vector<double>& A) {
    double lambda = RandomDouble(0, 1);
    for (int n = 0; n < A.size(); ++n)
        if (lambda <= A[n])
            return n;
}

void SlcTwoChrom(vector<double>& A,chromosome& Chrom1,chromosome& Chrom2 ){
    int ind1 = SltChr(A);
    int ind2 = SltChr(A);
    while (ind1 == ind2) {
        ind2 = SltChr(A);
    }
    if (population[ind1].FitnessValue + PrecisionValue < population[ind2].FitnessValue){
        Chrom1 = population[ind1];
        Chrom2 = population[ind2];
    } else {
        Chrom1 = population[ind2];
        Chrom2 = population[ind1];
    }
}

//{select two different chromosomes using the tournament method}
//it can only be used in the population where the chromosome have been sorted according fitness from good to bad
void SelectionTournament(int& parent_1, int& parent_2 , int& NumOfChormPerPop) {
    int P1 = rand() % NumOfChormPerPop;
    int P2 = rand() % NumOfChormPerPop;
    while (P1 == P2)
        P2 = rand() % NumOfChormPerPop;
    if (P1 < P2)
        parent_1 = P1;
    else
        parent_1 = P2;

    parent_2 = parent_1;
    while (parent_2 == parent_1) {
        P1 = rand() % NumOfChormPerPop;
        P2 = rand() % NumOfChormPerPop;
        while (P1 == P2)
            P2 = rand() % NumOfChormPerPop;
        if (P1 < P2)
            parent_2 = P1;
        else
            parent_2 = P2;
    }
}

void CrsMS_MP(chromosome& chrom1, chromosome& chrom2) {
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()]; //select a task randomly
        XY_SWAP(chrom1.RscAlcLst[TaskId], chrom2.RscAlcLst[TaskId], int);
    }
}

//{The SS mutation based TS:select a task randomly,then select a different position in its valid range to insert it}
void MtnSS_TS(chromosome& ch) {
    int pos = rand() % comConst.NumOfTsk;
    int TaskID = ch.TskSchLst[pos];
    int str = pos - 1, end = pos + 1;
    while (str > -1 && (find(Tasks[TaskID].parents.begin(), Tasks[TaskID].parents.end(), ch.TskSchLst[str]) ==
                        Tasks[TaskID].parents.end()))
        --str;
    ++str;
    while (end < comConst.NumOfTsk &&
           (find(Tasks[TaskID].children.begin(), Tasks[TaskID].children.end(), ch.TskSchLst[end]) ==
            Tasks[TaskID].children.end()))
        ++end;

    if (end - str <= 1) {
        return;
    }
    int InsertPoint = rand() % (end - str) + str;
    while (InsertPoint == pos) { // select a different position
        InsertPoint = rand() % (end - str) + str;
    }
    if (InsertPoint < pos) {
        for (int i = pos; i > InsertPoint; --i) {
            ch.TskSchLst[i] = ch.TskSchLst[i - 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    } else {
        for (int i = pos; i < InsertPoint; ++i) {
            ch.TskSchLst[i] = ch.TskSchLst[i + 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    }
}

//{The MS mutation based multiple point}
void MtnMS_MP(chromosome& ch) {
    int gamma = 1 + rand() % (int(comConst.NumOfTsk / 4) + 1);
    while (gamma--) {
        int i = rand() % comConst.NumOfTsk;
        int j = rand() % Tasks[i].ElgRsc.size();
        ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
    }
}

void MtnMS_SP(chromosome& ch) {
    int i = rand() % comConst.NumOfTsk;
    if (Tasks[i].ElgRsc.size() == 1) {
        return;
    }
    int j = rand() % Tasks[i].ElgRsc.size();
    while (ch.RscAlcLst[i] == Tasks[i].ElgRsc[j]) {
        j = rand() % Tasks[i].ElgRsc.size();
    }
    ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
}

void CrsSS_TS_R(chromosome& ch1, chromosome& ch2, int CrossPoint){
    vector<int> tem1 = ch1.TskSchLst;
    int delta = comConst.NumOfTsk-1;
    for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
        bool fd = false;
        for (int j = CrossPoint-1; j >= 0; --j) {
            if (ch2.TskSchLst[i] == ch1.TskSchLst[j]) {
                fd = true;
                break;
            }
        }
        if(!fd){
            ch1.TskSchLst[delta] = ch2.TskSchLst[i];
            delta--;
            if (delta < CrossPoint){
                break;
            }
        }
    }
    delta = comConst.NumOfTsk-1;
    for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
        bool fd = false;
        for (int j = CrossPoint-1; j >= 0; --j) {
            if (tem1[i] == ch2.TskSchLst[j]) {
                fd = true;
                break;
            }
        }
        if(!fd){
            ch2.TskSchLst[delta] = tem1[i];
            delta--;
            if (delta < CrossPoint){
                break;
            }
        }
    }
}

//{HGA: single point crossover }
void CrsMS_SP(chromosome& ch1, chromosome& ch2) {
    int CrossPoint = rand() % (comConst.NumOfTsk - 1) + 1;
    for (int i = 0; i < CrossPoint; ++i) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

//CGA crossover
void Crossover_CGA(chromosome& ch1, chromosome& ch2) {
    if (RandomDouble(0, 1) < Parameter_CGA.CrossoverRate) {
        CrsMS_SP(ch1, ch2);
    }
    if (RandomDouble(0, 1) < Parameter_CGA.CrossoverRate) {
        int CrossPoint = 1 + rand() % (comConst.NumOfTsk - 1);
        CrsSS_TS_R(ch1, ch2, CrossPoint);
    }
}

//CGA mutation
void Mutation_CGA(chromosome& ch) {
    if (RandomDouble(0, 1) < Parameter_CGA.MutationRate) {
        MtnMS_SP(ch);
    }
    if (RandomDouble(0, 1) < Parameter_CGA.MutationRate) {
        MtnSS_TS(ch);
    }
}

//{HGA:two point crossover, HGA}
void CrsMS_DP(chromosome& ch1, chromosome& ch2) {
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point1 == point2) {
        point2 = rand() % comConst.NumOfTsk;
    }
    if (point1 > point2) {
        XY_SWAP(point1, point2, int);
    }
    for (int i = point1;  i <= point2; ++i ) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

void GnrTskSchLst_HGA(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = ch.RscAlcLst[parent];
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
    }
}

//｛HGA crossover｝
void Crossover_HGA(chromosome& ch1, chromosome& ch2) {
    chromosome x1 = ch1;
    chromosome x2 = ch2;
    chromosome y1 = ch1;
    chromosome y2 = ch2;
    CrsMS_SP(x1, x2);
    GnrTskSchLst_HGA(x1);
    DcdEvl(x1, true);
    GnrTskSchLst_HGA(x2);
    DcdEvl(x2, true);

    CrsMS_DP(y1, y2);
    GnrTskSchLst_HGA(y1);
    DcdEvl(y1, true);
    GnrTskSchLst_HGA(y2);
    DcdEvl(y2, true);
    vector<chromosome> sub;
    sub.push_back(x1);
    sub.push_back(x2);
    sub.push_back(y1);
    sub.push_back(y2);
    sort(sub.begin(), sub.end(), SortPopOnFitValueByAscend);
    ch1 = sub[0];
    ch2 = sub[1];
}

void Mutation_HGA(chromosome& ch) {
    chromosome x = ch;
    chromosome y = ch;
    //{single point mutation on x}
    int point = rand() % comConst.NumOfTsk;
    x.RscAlcLst[point] = Tasks[point].ElgRsc[rand() % Tasks[point].ElgRsc.size()];
    //{double point mutation on y}
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point2 == point1) {
        point2 = rand() % comConst.NumOfTsk;
    }
    y.RscAlcLst[point1] = Tasks[point1].ElgRsc[rand() % Tasks[point1].ElgRsc.size()];
    y.RscAlcLst[point2] = Tasks[point2].ElgRsc[rand() % Tasks[point2].ElgRsc.size()];
    GnrTskSchLst_HGA(x);
    DcdEvl(x, true);
    GnrTskSchLst_HGA(y);
    DcdEvl(y, true);
    if ( y.FitnessValue + PrecisionValue < x.FitnessValue ) {
        ch = y;
    } else {
        ch = x;
    }
}

//load balance improvement for HGA
void RscLoadAdjust_HGA(vector<chromosome>& Pop) {
    vector<double> lb(Pop.size());
    #pragma omp parallel for
    for(int n =0 ;n < Pop.size(); ++n){
        //calculate the finish times of all task for each resource and find out the maximum
        vector<double> FT(comConst.NumOfRsc,0);
        for(int j = 0 ;j < Pop[n].RscAlcLst.size(); ++j){
            int IndexRsc = Pop[n].RscAlcLst[j];
            if(FT[IndexRsc] < Pop[n].EndTime[j]){
                FT[IndexRsc] = Pop[n].EndTime[j];
            }
        }
        //{find the minimum in FT}
        double min = 9999999;
        for(int j = 0; j < comConst.NumOfRsc; ++j){
            if(min>FT[j]){
                min = FT[j];
            }
        }
        lb[n] = Pop[n].FitnessValue - min;
    }
    //{sort lb}
    vector<int> IndexLb(Pop.size());
    IndexSort(IndexLb,lb);             //sorting chromosome by lb from small to large
    //{According to LB, select the last 50% from small to large to improve}
    #pragma omp parallel for
    for (int n = Pop.size() / 2; n < Pop.size(); ++n) {
        chromosome chrom = Pop[IndexLb[n]];
        vector<double> ld (comConst.NumOfRsc,0);
        vector<vector<int>> TSK(comConst.NumOfRsc);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            int RscIndex = chrom.RscAlcLst[i];
            ld[RscIndex] += (1.0 * Tasks[i].length) / Rscs[RscIndex].pc;
            TSK[RscIndex].push_back(i);
        }
        vector<int> Ind(comConst.NumOfRsc);
        IndexSort(Ind, ld);                                                                  //load sort
        int BigRsc = Ind[Ind.size() - 1];                                                            //obtain the maximum
        int RandTask = TSK[BigRsc][rand() % TSK[BigRsc].size()];                                     //select a task from the Rsc with the largest load
        chrom.RscAlcLst[RandTask] = Tasks[RandTask].ElgRsc[rand() % Tasks[RandTask].ElgRsc.size()];  //reallocation
        GnrTskSchLst_HGA(chrom);
        DcdEvl(chrom, true);
        if ( chrom.FitnessValue + PrecisionValue < Pop[IndexLb[n]].FitnessValue ) {
            Pop[IndexLb[n]] = chrom;
        }
    }
}

//{LWSGA: level swapping (exchange all tasks in the level) }
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2){
    int RandLevel = rand()%TskLstInLvl.size();
    //{Rsc swap}
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        int TaskIndex = TskLstInLvl[RandLevel][i];
        XY_SWAP(chrom1.RscAlcLst[TaskIndex], chrom2.RscAlcLst[TaskIndex], int);
    }
    //{find the start point of level }
    int pos = 0;
    for(int i = 0; i < RandLevel; ++i){
        pos += TskLstInLvl[i].size();
    }
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        XY_SWAP(chrom1.TskSchLst[pos], chrom2.TskSchLst[pos], int);
        ++pos;
    }
}

//{LWSGA: exchange two tasks in each level }
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2) {
    int p1 = 0, p2 = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()];
        //{Since the task with small level must be in front of the task with large level, the search can be started from the last recorded position}
        for (int j = p1; j < comConst.NumOfTsk; ++j) {
            if (chrom1.TskSchLst[j] == TaskId) {
                p1 = j;
                break;
            }
        }
        for (int j = p2; j < comConst.NumOfTsk; ++j) {
            if (chrom2.TskSchLst[j] == TaskId) {
                p2 = j;
                break;
            }
        }
        XY_SWAP(chrom1.TskSchLst[p1], chrom1.TskSchLst[p2], int);
        XY_SWAP(chrom2.TskSchLst[p1], chrom2.TskSchLst[p2], int);
    }
}

void Crossover_LWSGA(chromosome& ch1, chromosome& ch2) {
    int method = rand() % 3;
    if (method == 0) {
        Crs_Lvl(ch1, ch2);
    } else if (method == 1) {
        CrsSS_ExcTskInLvl(ch1, ch2);
    } else {
        CrsMS_MP(ch1, ch2);
    }
}

//(mutation: exchange two tasks in level)
void MtnSS_ExcTskInLvl(chromosome& chrom) {
    int RandLevel = rand() % TskLstInLvl.size();
    int t1 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int t2 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int p1 = -1, p2 = -1;
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t1) {
            p1 = i;
            break;
        }
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t2) {
            p2 = i;
            break;
        }
    }
    XY_SWAP(chrom.TskSchLst[p1], chrom.TskSchLst[p2], int);
}

//{select a level, rearrange these tasks in this level and reallocate the resources for these tasks in this level}
void Mtn_rebuild_level(chromosome& ch) {
    int SctLvl = rand() % TskLstInLvl.size();
    vector<int> TemTskLst = TskLstInLvl[SctLvl];
    random_shuffle(TemTskLst.begin(), TemTskLst.end()); // rearrange these tasks
    int StartIndex = 0;
    for(int i = 0; i < SctLvl; ++i) {
        StartIndex = StartIndex + TskLstInLvl[i].size();
    }
    for(int i = 0;i < TemTskLst.size(); ++i) {
        int index = StartIndex + i;
        int TaskId = TemTskLst[i];
        ch.TskSchLst[index] = TaskId;
        int RscIndex = Tasks[TaskId].ElgRsc[rand() % Tasks[TaskId].ElgRsc.size()];
        ch.RscAlcLst[TaskId] = RscIndex;
    }
}

void Mutation_LWSGA(chromosome& ch) {
    int method = rand() % 3;
    if (method == 0) {
        MtnSS_ExcTskInLvl(ch);
    } else if (method == 1) {
        MtnMS_MP(ch);
    } else {
        Mtn_rebuild_level(ch);
    }
}

void CrsLvl(chromosome& chrom1, chromosome& chrom2) {
    vector<list<int>> tem1 = chrom1.Code_TD;
    vector<list<int>> tem2 = chrom2.Code_TD;
    //层次值是从0开始的
    int crossoverSite = rand() % ( TskLstInLvl.size() - 1);
    vector<list<int>> TL1(comConst.NumOfRsc);
    vector<list<int>> TL2(comConst.NumOfRsc);
    list<int>::iterator iter;
    for(int j = 0; j < comConst.NumOfRsc; ++j) {
        for(iter = tem1[j].begin(); iter != tem1[j].end();){
            if(LevelIdOfTask[*iter] <= crossoverSite) {
                TL1[j].push_back(*iter);
                iter = tem1[j].erase(iter);
            }else{
                break;
            }
        }
    }

    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        for(iter = tem2[j].begin() ; iter != tem2[j].end() ;){
            if(LevelIdOfTask[*iter] <= crossoverSite) {
                TL2[j].push_back(*iter);
                iter = tem2[j].erase(iter);
            } else{
                break;
            }
        }
    }

    for(int j = 0; j < comConst.NumOfRsc; ++j) {
        for(iter = tem2[j].begin() ; iter != tem2[j].end() ; iter++){
            TL1[j].push_back(*iter);
        }
        for(iter = tem1[j].begin() ; iter != tem1[j].end() ; iter++){
            TL2[j].push_back(*iter);
        }
    }
    chrom1.Code_TD = TL1;
    chrom2.Code_TD = TL2;
}

void CrsTS(chromosome& chrom1, chromosome& chrom2) {
    vector<list<int> > OGL1 = chrom1.Code_TD;
    vector<list<int> > OGL2 = chrom2.Code_TD;
    vector<list<int> > NGL1(comConst.NumOfRsc);
    vector<list<int> > NGL2(comConst.NumOfRsc);

    vector<int > upr1(comConst.NumOfTsk,0);
    vector<int > upr2(comConst.NumOfTsk,0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr1[i] = Tasks[i].parents.size();
        upr2[i] = upr1[i];
    }

    int randNum = rand()%(comConst.NumOfTsk-1) + 1;
    while(randNum > 0) {
        --randNum;
        int task1 = -1;
        int Rsc1 = -1;
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            if (!OGL1[j].empty() && upr1[OGL1[j].front()] == 0) {
                task1 = OGL1[j].front();
                Rsc1 = j;
                OGL1[j].erase(OGL1[j].begin());
                break;
            }
        }
        int task2 = -1;
        int Rsc2 = -1;
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            if (!OGL2[j].empty() && upr2[OGL2[j].front()] == 0) {
                task2 = OGL2[j].front();
                Rsc2 = j;
                OGL2[j].erase(OGL2[j].begin());
                break;
            }
        }
        for (int k = 0; k < Tasks[task1].children.size(); ++k) {
            int childId = Tasks[task1].children[k];
            upr1[childId] = upr1[childId] - 1;
        }

        for (int k = 0; k < Tasks[task2].children.size(); ++k) {
            int childId = Tasks[task2].children[k];
            upr2[childId] = upr2[childId] - 1;
        }

        NGL1[Rsc1].push_back(task1);
        NGL2[Rsc2].push_back(task2);
        int flag = 0;

        for(int j = 0; j < comConst.NumOfRsc; ++j) {
            for(auto iter = chrom2.Code_TD[j].begin(); iter != chrom2.Code_TD[j].end();) {
                if(*iter == task1) {
                    chrom2.Code_TD[j].erase(iter);
                    flag = 1;
                    break;
                } else {
                    ++iter;
                }
            }
            if(flag==1){
                break;
            }
        }

        for(int j = 0; j < comConst.NumOfRsc; ++j) {
            for(auto iter = chrom1.Code_TD[j].begin(); iter != chrom1.Code_TD[j].end();) {
                if(*iter == task2) {
                    chrom1.Code_TD[j].erase(iter);
                    flag = 2;
                    break;
                } else {
                    ++iter;
                }
            }
            if(flag==2){
                break;
            }
        }
    }
    for(int j = 0; j < comConst.NumOfRsc; ++j) {
        for(auto iter = chrom2.Code_TD[j].begin() ; iter != chrom2.Code_TD[j].end() ; iter++){
            NGL1[j].push_back(*iter);
        }
    }
    for(int j = 0; j < comConst.NumOfRsc; ++j) {
        for(auto iter = chrom1.Code_TD[j].begin() ; iter != chrom1.Code_TD[j].end() ; iter++){
            NGL2[j].push_back(*iter);
        }
    }
    chrom1.Code_TD = NGL1;
    chrom2.Code_TD = NGL2;
}

void FindInsertRangeForLvlMtn(chromosome &Chrom, int &SlcNewRsc, int &SlcOldRsc, int &SlcTask, _List_iterator<int> &start,
                              _List_iterator<int> &end) {
    Chrom.Code_TD[SlcOldRsc].erase(find(Chrom.Code_TD[SlcOldRsc].begin(), Chrom.Code_TD[SlcOldRsc].end(), SlcTask));
    int RandTaskLevel = LevelIdOfTask[SlcTask];
    start = Chrom.Code_TD[SlcNewRsc].begin();
    end = Chrom.Code_TD[SlcNewRsc].end();
    for(auto iter = Chrom.Code_TD[SlcNewRsc].begin() ; iter != Chrom.Code_TD[SlcNewRsc].end() ; iter++){
        if(LevelIdOfTask[*iter] < RandTaskLevel) {
            start = iter;
            start++;
        } else if(LevelIdOfTask[*iter] > RandTaskLevel) {
            end = iter;
            break;
        }
    }
}

void FindInsertRangeForTplSrtMtn(chromosome &Chrom, int &SlcNewRsc, int &SlcOldRsc, int &SlcTask, _List_iterator<int> &start,
                                 _List_iterator<int> &end) {
    vector<int> TSL;
    GnrTskSchLst(Chrom, TSL);
    int RandPlace = -1;
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        if(SlcTask == TSL[i]){
            RandPlace = i;
            break;
        }
    }
    int temstr = RandPlace-1, temend = RandPlace+1;
    while (temstr >= 0 && find(Tasks[SlcTask].parents.begin(), Tasks[SlcTask].parents.end(),TSL[temstr]) == Tasks[SlcTask].parents.end()){
        temstr = temstr-1;
    }
    while (temend <= comConst.NumOfTsk-1 && find(Tasks[SlcTask].children.begin(), Tasks[SlcTask].children.end(),TSL[temend]) == Tasks[SlcTask].children.end()){
        temend = temend+1;
    }
    Chrom.Code_TD[SlcOldRsc].erase(find(Chrom.Code_TD[SlcOldRsc].begin(), Chrom.Code_TD[SlcOldRsc].end(), SlcTask));

    start = Chrom.Code_TD[SlcNewRsc].begin();
    list<int>::iterator iter = Chrom.Code_TD[SlcNewRsc].begin();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(iter == Chrom.Code_TD[SlcNewRsc].end()){
            break;
        }
        if (TSL[i] == *iter) {
            if (i <= temstr) {
                start = iter;
                start++;
                iter++;
            } else {
                break;
            }
        }
    }

    end = Chrom.Code_TD[SlcNewRsc].end();
    list<int>::iterator tend = end;
    list<int>::iterator PreBegin = Chrom.Code_TD[SlcNewRsc].begin();
    PreBegin--;
    tend--;
    for (int i = comConst.NumOfTsk-1; i >= 0; --i) {
        if(tend == PreBegin){
            break;
        }
        if (TSL[i] == *tend) {
            if (i >= temend){
                end = tend;
                tend--;
            } else {
                break;
            }
        }
    }
}

void Mtn(chromosome& chrom, double Pm, int stg) {
    double random = RandomDouble(0, 1);
    if (random < Pm) {
        vector<double> Id(comConst.NumOfRsc, 0);
        list<int>::iterator it;
        double sum = 0;
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            for (it = chrom.Code_TD[j].begin(); it != chrom.Code_TD[j].end(); it++) {
                Id[j] += Tasks[*it].length / Rscs[j].pc;
            }
            sum += Id[j];
        }
        vector<double> Prob(Id.size());
        for(int j = 0; j < Id.size(); ++j){
            Prob[j] = Id[j] / sum;
        }
        double randNum = rand() % 100 / 100.0;
        double ProbSum = 0;
        int SlcOldRsc = -1;
        for (int j = 0; j < Id.size(); ++j) {
            ProbSum += Prob[j];
            if (ProbSum > randNum) {
                SlcOldRsc = j;
                break;
            }
        }

        int randTaskIndex = rand() % chrom.Code_TD[SlcOldRsc].size();
        list<int>::iterator iterP = chrom.Code_TD[SlcOldRsc].begin();
        advance(iterP, randTaskIndex);
        int SlcTask = *iterP;
        int SlcNewRsc = -1;
        if (Tasks[SlcTask].ElgRsc.size() == 1) {
            SlcNewRsc = Tasks[SlcTask].ElgRsc[0];
        } else {
            Prob.resize(Tasks[SlcTask].ElgRsc.size());
            sum = 0;
            for (int j = 0; j < Tasks[SlcTask].ElgRsc.size(); ++j) {
                sum += Id[Tasks[SlcTask].ElgRsc[j]];
            }
            for (int j = 0; j < Tasks[SlcTask].ElgRsc.size(); ++j) {
                Prob[j] = Id[Tasks[SlcTask].ElgRsc[j]] / sum;
            }
            sum = 0;
            for (int j = 0; j < Prob.size(); ++j) {
                sum += 1 - Prob[j];
            }
            ProbSum = 0;
            randNum = rand() % 100 / 100.0;
            for (int j = 0; j < Prob.size(); ++j) {
                ProbSum += (1 - Prob[j]) / sum;
                if (ProbSum > randNum) {
                    SlcNewRsc = Tasks[SlcTask].ElgRsc[j];
                    break;
                }
            }
        }
        list<int>::iterator start;
        list<int>::iterator end;
        if (stg == 1) {
            FindInsertRangeForLvlMtn(chrom, SlcNewRsc, SlcOldRsc, SlcTask, start, end);
        } else {
            FindInsertRangeForTplSrtMtn(chrom, SlcNewRsc, SlcOldRsc, SlcTask, start, end);
        }
        int Num = 1;
        for (; start != end; --end) {
            Num++;
        }
        int a = (rand() % Num);
        advance(start, a);
        chrom.Code_TD[SlcNewRsc].insert(start, SlcTask);
    }
}