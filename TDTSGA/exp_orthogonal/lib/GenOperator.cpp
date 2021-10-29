#include <cstdlib>
#include <GenerateAChrom.h>
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

void CrsLvl(chromosome& chrom1, chromosome& chrom2) {
    vector<list<int>> tem1 = chrom1.Code_TD;
    vector<list<int>> tem2 = chrom2.Code_TD;
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
