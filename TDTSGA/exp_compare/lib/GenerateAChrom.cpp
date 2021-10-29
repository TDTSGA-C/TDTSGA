#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c) {
    for(int i =1 ;i < TskLstInLvl.size(); ++i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId = TskLstInLvl[i][j];
            for (int k = 0; k < Tasks[TaskId].parents.size(); ++k) {
                int tem = Tasks[TaskId].parents[k];
                double re = w[tem] + c[tem][TaskId] + RankList[tem];
                if (RankList[TaskId] < re) {
                    RankList[TaskId] = re;
                }
            }
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.Code_TD.resize(comConst.NumOfRsc);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

//{generate a task scheduling order by the levels of tasks from small to large}
//{Those haveing the same level are ranked arbitrarily among them}
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

double GnrMS_Evl(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
            double ReadyTime = 0;
            int v = Tasks[TaskIndex].ElgRsc[j];
            if(Tasks[TaskIndex].parents.size() != 0){
                for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                    int ParentIndex = Tasks[TaskIndex].parents[n];
                    int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                    double max = ch.EndTime[ParentIndex];
                    if(v != ParentRscIndex){
                        double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                    }
                    if (ReadyTime < max){
                        ReadyTime = max;
                    }
                }
            }
            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
            double StartTime = 0;
            double EndTime = 0;
            //{Find an idle time-slot as early as possible from ITL}
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                    StartTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            EndTime = StartTime + ExeTime;
            //{find/record the earliest finish time}
            if (EndTime < FinalEndTime) {
                FinalStartTime = StartTime;
                FinalEndTime = EndTime;
                RscIndex = v;
            }
        }
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        //{update ITL}
        if(ITL[RscIndex].find(FinalStartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(FinalStartTime);
        } else {
            ITL[RscIndex].insert(FinalStartTime);
        }
        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

void EFT(chromosome& ch,vector<set<double>>& ITL,int& TaskIndex ,int& RscIndex,double& FinalStartTime,double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int v = Tasks[TaskIndex].ElgRsc[j];
        if(Tasks[TaskIndex].parents.size() != 0){
            for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                int ParentIndex = Tasks[TaskIndex].parents[n];
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                double max = ch.EndTime[ParentIndex];
                if(v != ParentRscIndex){
                    double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                    max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                }
                if (ReadyTime < max){
                    ReadyTime = max;
                }
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
        double StartTime = 0;
        double EndTime = 0;
        //{Find an idle time-slot as early as possible from ITL}
        set<double>::iterator pre  = ITL[v].begin();
        set<double>::iterator post = ITL[v].begin();
        ++post;
        while(post != ITL[v].end()) {
            if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = v;
        }
    }
}

void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime){
    if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
        ITL[RscIndex].erase(StartTime);
    } else {
        ITL[RscIndex].insert(StartTime);
    }
    if(ITL[RscIndex].find(EndTime) != ITL[RscIndex].end()) {
        ITL[RscIndex].erase(EndTime);
    } else {
        ITL[RscIndex].insert(EndTime);
    }
}

void FindIdleTimeSlot(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITL[RscIndex].begin();
    set<double>::iterator post = ITL[RscIndex].begin();
    ++post;
    while(post != ITL[RscIndex].end()) {
        if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
            StartTime = XY_MAX(*pre, ReadyTime);
            break;
        } else {
            ++pre;
            ++pre;
            ++post;
            ++post;
        }
    }
}

chromosome GnrChr_HEFT_HGA(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

// {in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of the number of child tasks}
chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSort(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrMS_Evl(TemChrom);
   return TemChrom;
}

chromosome GnrChr_Lvl_Rnd(){
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<vector<int> > tem = TskLstInLvl;
    for (int k = 0; k < tem.size(); ++k) {
        random_shuffle(tem[k].begin(), (tem[k].end()-1));
        for (int i = 0; i < tem[k].size(); ++i) {
            int TaskId = tem[k][i];
            int RandNum = rand() % Tasks[TaskId].ElgRsc.size();
            TemChrom.Code_TD[Tasks[TaskId].ElgRsc[RandNum]].push_back(TaskId);
        }
    }
    Dcd(TemChrom, true);
    return TemChrom;
}

chromosome GnrChr_Tpl_Rnd(){
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<int > upr(comConst.NumOfTsk,0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0) {
            RTI.push_back(i);
        }
    }
    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        int RandNum = rand() % RTI.size();
        int taskIndex = RTI[RandNum];
        RTI.erase(find(RTI.begin(), RTI.end(), taskIndex));
        for (int l = 0; l < Tasks[taskIndex].children.size(); ++l) {
            int childId = Tasks[taskIndex].children[l];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                RTI.push_back(childId);
            }
        }
        int Rand = rand() % Tasks[taskIndex].ElgRsc.size();
        TemChrom.Code_TD[Tasks[taskIndex].ElgRsc[Rand]].push_back(taskIndex);
    }
    Dcd(TemChrom, true);
    return TemChrom;
}

chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int> ind(comConst.NumOfTsk);
    chromosome chrom;
    IntChr(chrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        int RscIndex = -1;
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        EFT(chrom,ITL,TaskIndex ,RscIndex,FinalStartTime,FinalEndTime);
        chrom.EndTime[TaskIndex] = FinalEndTime;
        chrom.RscAlcLst[TaskIndex] = RscIndex;
        chrom.Code_TD[RscIndex].push_back(TaskIndex);
        UpdateITL(ITL,RscIndex,FinalStartTime,chrom.EndTime[TaskIndex]);
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    chrom.FitnessValue = makespan;
    return chrom;
}

chromosome GnrChr_Lvl_EFT() {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    chromosome chrom;
    IntChr(chrom);
    vector<vector<int> > tem = TskLstInLvl;
    for (int k = 0; k < TskLstInLvl.size(); ++k) {
        random_shuffle(tem[k].begin(), tem[k].end());
        for (int i = 0; i < tem[k].size(); ++i) {
            int TaskIndex = tem[k][i];
            int RscIndex = -1;
            double FinalEndTime = 100000000000;
            double FinalStartTime = 0;
            EFT(chrom,ITL,TaskIndex ,RscIndex,FinalStartTime,FinalEndTime);
            chrom.EndTime[TaskIndex] = FinalEndTime;
            chrom.RscAlcLst[TaskIndex] = RscIndex;
            chrom.Code_TD[RscIndex].push_back(TaskIndex);
            UpdateITL(ITL,RscIndex,FinalStartTime,chrom.EndTime[TaskIndex]);
            makespan = XY_MAX(makespan, FinalEndTime);
        }
    }
    chrom.FitnessValue = makespan;
    return chrom;
}

chromosome GnrChr_Tpl_EFT() {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk,0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0) {
            RTI.push_back(i);
        }
    }
    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        int RandNum = rand() % RTI.size();
        int TaskIndex = RTI[RandNum];
        RTI.erase(find(RTI.begin(), RTI.end(), TaskIndex));
        for (int l = 0; l < Tasks[TaskIndex].children.size(); ++l) {
            int childId = Tasks[TaskIndex].children[l];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                RTI.push_back(childId);
            }
        }
        int RscIndex = -1;
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        EFT(chrom,ITL,TaskIndex ,RscIndex,FinalStartTime,FinalEndTime);
        chrom.EndTime[TaskIndex] = FinalEndTime;
        chrom.RscAlcLst[TaskIndex] = RscIndex;
        chrom.Code_TD[RscIndex].push_back(TaskIndex);
        //{update ITL}
        UpdateITL(ITL,RscIndex,FinalStartTime,chrom.EndTime[TaskIndex]);
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    chrom.FitnessValue = makespan;
    return chrom;
}

void SeletRsc_EFT(chromosome& ch,vector<set<double>>& ITL,int& TaskIndex ,int& RscIndex,double& FinalStartTime,double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int v = Tasks[TaskIndex].ElgRsc[j];
        if(Tasks[TaskIndex].parents.size() != 0){
            for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                int ParentIndex = Tasks[TaskIndex].parents[n];
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                double max = ch.EndTime[ParentIndex];
                if(v != ParentRscIndex){
                    double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                    max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                }
                if (ReadyTime < max){
                    ReadyTime = max;
                }
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
        double StartTime = 0;
        double EndTime = 0;
        //{Find an idle time-slot as early as possible from ITL}
        set<double>::iterator pre  = ITL[v].begin();
        set<double>::iterator post = ITL[v].begin();
        ++post;
        while(post != ITL[v].end()) {
            if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = v;
        }
    }
}

double IHEFT3_S(chromosome& ch) {
    list <int> TemTskSchLst;
    TemTskSchLst.assign(ch.TskSchLst.begin(),ch.TskSchLst.end());
    ch.RscAlcLst.resize(comConst.NumOfTsk,-1);
    ch.TskSchLst.resize(comConst.NumOfTsk,-1);
    vector<set<double> > ITL;
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
    }

    int IndexCount = 0;
    while (!TemTskSchLst.empty()){
        int CrnTask = TemTskSchLst.front();
        ch.TskSchLst[IndexCount] = CrnTask;
        IndexCount++;
        TemTskSchLst.erase(TemTskSchLst.begin());
        int FinalRscIdForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0;
        double FinalEndTimeOfCrnTask = 10000000000000;

        vector<int> NeedProcessChildTaskSet;
        for (int m = 0; m < Tasks[CrnTask].children.size(); ++m) {
            int childId = Tasks[CrnTask].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                NeedProcessChildTaskSet.push_back(childId);
            }
        }

        if(NeedProcessChildTaskSet.empty()){
            int FinalRscOfChildTask = -1;
            SeletRsc_EFT(ch,ITL,CrnTask, FinalRscOfChildTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscOfChildTask;
            ch.Code_TD[FinalRscOfChildTask].push_back(CrnTask);
            UpdateITL(ITL,FinalRscOfChildTask,FinalStartTimeOfCrnTask,ch.EndTime[CrnTask]);
            makespan = XY_MAX(makespan, FinalEndTimeOfCrnTask);
        } else{
            int FinalChildTask = -1 , FinalRscOfChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = 10000000000000;
            vector<set<double> > ITL_CrnTaskScheduled;
            for (int j = 0; j < Tasks[CrnTask].ElgRsc.size(); ++j) {
                double RscOfCrnTask = 0;
                int CrnTaskRsc = Tasks[CrnTask].ElgRsc[j];
                ch.RscAlcLst[CrnTask] = CrnTaskRsc;
                if(Tasks[CrnTask].parents.size() != 0){
                    for (int i2 = 0; i2 < Tasks[CrnTask].parents.size(); ++i2) {
                        int ParentTask = Tasks[CrnTask].parents[i2];
                        int RscOfParentTask = ch.RscAlcLst[ParentTask];
                        double max = ch.EndTime[ParentTask];
                        if(CrnTaskRsc != RscOfParentTask){
                            double TransferData = ParChildTranFileSizeSum[ParentTask][CrnTask];
                            max += TransferData / VALUE * 8 / (XY_MIN(Rscs[CrnTaskRsc].bw,Rscs[RscOfParentTask].bw));
                        }
                        if (RscOfCrnTask < max){
                            RscOfCrnTask = max;
                        }
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTask].length / Rscs[CrnTaskRsc].pc;
                double CrnTaskStartTime = 0, CrnTaskEndTime = 0;
                FindIdleTimeSlot(ITL,CrnTaskRsc,CrnTaskStartTime,ExeTimeOfCrnTask,RscOfCrnTask);
                CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL,CrnTaskRsc,CrnTaskStartTime,CrnTaskEndTime);
                ch.EndTime[CrnTask] = CrnTaskEndTime;
                for (int i3 = 0; i3 < NeedProcessChildTaskSet.size(); ++i3) {
                    int TemChildTask = NeedProcessChildTaskSet[i3];
                    for (int j1 = 0; j1 < Tasks[TemChildTask].ElgRsc.size(); ++j1) {
                        double TemChildReadyTime = 0;
                        int TemChildRsc = Tasks[TemChildTask].ElgRsc[j1];
                        for (int i4 = 0; i4 < Tasks[TemChildTask].parents.size(); ++i4) {
                            int TemParentTask = Tasks[TemChildTask].parents[i4];
                            int TemParRsc = ch.RscAlcLst[TemParentTask];
                            double max = ch.EndTime[TemParentTask];
                            if (TemChildRsc != TemParRsc) {
                                double TransferData = ParChildTranFileSizeSum[TemParentTask][TemChildTask];
                                max += TransferData / VALUE * 8 /
                                       (XY_MIN(Rscs[TemChildRsc].bw, Rscs[TemParRsc].bw));
                            }
                            if (TemChildReadyTime < max) {
                                TemChildReadyTime = max;
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildTask].length / Rscs[TemChildRsc].pc;
                        double TemChildStartTime = 0,TemChildEndTime = 0;
                        FindIdleTimeSlot(TemITL,TemChildRsc,TemChildStartTime,TemChildExeTime,TemChildReadyTime);
                        TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime) {
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalRscOfChildTask = TemChildRsc;
                            FinalChildTask = TemChildTask;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscIdForCrnTask = CrnTaskRsc;
                            FinalStartTimeOfCrnTask = CrnTaskStartTime;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscIdForCrnTask;
            ch.Code_TD[FinalRscIdForCrnTask].push_back(CrnTask);
            ch.TskSchLst[IndexCount] = FinalChildTask;
            ch.EndTime[FinalChildTask] = FinalEndTimeOfChildTask;
            ch.RscAlcLst[FinalChildTask] = FinalRscOfChildTask;
            ch.Code_TD[FinalRscOfChildTask].push_back(FinalChildTask);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL,FinalRscOfChildTask,FinalStartTimeOfChildTask,ch.EndTime[FinalChildTask]);
            makespan = XY_MAX(makespan, FinalEndTimeOfChildTask);
            TemTskSchLst.erase(find(TemTskSchLst.begin(),TemTskSchLst.end(),FinalChildTask));
            for (int m = 0; m < Tasks[FinalChildTask].children.size(); ++m) {
                int childId = Tasks[FinalChildTask].children[m];
                upr[childId] = upr[childId] - 1;
            }
            IndexCount++;
        }
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChrIHEFT3(vector<double> Rank_b,vector<double> Rank_t) {
    vector<int> Ind_Rankb(comConst.NumOfTsk);
    chromosome Chrom_Rankb;
    IntChr(Chrom_Rankb);
    IndexSort(Ind_Rankb, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        Chrom_Rankb.TskSchLst[i] = Ind_Rankb[comConst.NumOfTsk - i - 1];
    }
    vector<int> Ind_Rankt(comConst.NumOfTsk);
    chromosome Chrom_Rankt;
    IntChr(Chrom_Rankt);
    IndexSort(Ind_Rankt, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        Chrom_Rankt.TskSchLst[i] = Ind_Rankt[i];
    }
    IHEFT3_S(Chrom_Rankb);
    IHEFT3_S(Chrom_Rankt);
    if (Chrom_Rankb.FitnessValue + PrecisionValue < Chrom_Rankt.FitnessValue){
        return Chrom_Rankb;
    } else{
        return Chrom_Rankt;
    }
}