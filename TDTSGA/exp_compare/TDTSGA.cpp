//
// Created by yx on 21-2-25.
//

#include "TDTSGA.h"
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runTDTSGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration){
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);         //read model information  XmlFile, RscAlcFile
    ConfigParameter_TDTSGA();                //set the parameter values
    CalculateLevelList();                  //calculate the levels of tasks

    vector<double> Rank_t(comConst.NumOfTsk,0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double> > cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_t(Rank_t, ww, cc);
    Calculate_Rank_b(Rank_b, cc, ww);
    int stg = 1;
    set<chromosome> TemPop;
    int flag = 0;
    int num = 1;
    while (TemPop.size() < Parameter_TDTSGA.NumOfChormPerPop){
        chromosome TemChrom;
        if( num <= 2*Parameter_TDTSGA.NumOfChormPerPop ){
            TemChrom = GnrChr_Lvl_EFT();
        } else if (num <= 4*Parameter_TDTSGA.NumOfChormPerPop) {
            TemChrom = GnrChr_Lvl_Rnd();
        } else if(num <= 6*Parameter_TDTSGA.NumOfChormPerPop) {
            if (flag == 0) {
                TemChrom = GnrChrIHEFT3(Rank_b,Rank_t);
                flag = 1;
            } else {
                TemChrom = GnrChr_Tpl_EFT();
            }
        } else{
            TemChrom = GnrChr_Tpl_Rnd();
        }
        TemPop.insert(TemChrom);
        num++;
    }
    if(num > 4*Parameter_TDTSGA.NumOfChormPerPop){
        stg = 2;
    }
    population.assign(TemPop.begin(),TemPop.end());
    vector<double> A(Parameter_TDTSGA.NumOfChormPerPop);
    CalSlctProb_Rank(1+Parameter_TDTSGA.beta/Parameter_TDTSGA.NumOfChormPerPop, A ,Parameter_TDTSGA.NumOfChormPerPop);

    while ((double) (clock() - start) / CLOCKS_PER_SEC < SchTime){
        ++iteration;
        vector<chromosome> NewPopulation(Parameter_TDTSGA.NumOfChormPerPop);
        for (int n = 0; n < Parameter_TDTSGA.NumOfChormPerPop; n += 2) {
            chromosome chrom1 ;
            chromosome chrom2 ;
            SlcTwoChrom(A,chrom1,chrom2);
            if (stg == 1) {
                CrsLvl(chrom1, chrom2);
                NewPopulation[n] = chrom1;
                NewPopulation[n+1] = chrom2;
                Mtn(NewPopulation[n], Parameter_TDTSGA.MutationRate1, stg);
                Mtn(NewPopulation[n+1], Parameter_TDTSGA.MutationRate1, stg);
            } else {
                CrsTS( chrom1,chrom2);
                NewPopulation[n] = chrom1;
                NewPopulation[n+1] = chrom2;
                Mtn(NewPopulation[n], Parameter_TDTSGA.MutationRate2, stg);
                Mtn(NewPopulation[n+1], Parameter_TDTSGA.MutationRate2, stg);
            }
        }
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TDTSGA.NumOfChormPerPop; ++n) {
            Dcd(NewPopulation[n], true);
        }

        set<chromosome> NxtPop; // NxtSubPop is used for merging chromosomes, removing the same chromosomes and sorting chromosomes -xy6
        NxtPop.insert(population.begin(),population.end());
        NxtPop.insert(NewPopulation.begin(),NewPopulation.end());
        set<chromosome>::iterator iter = NxtPop.begin();
        advance(iter,Parameter_TDTSGA.NumOfChormPerPop);
        population.assign(NxtPop.begin(),iter);
        if((double) (clock() - start) / CLOCKS_PER_SEC > Parameter_TDTSGA.RunTimeRatioOfStg1*SchTime && stg == 1){
            stg = 2;
            chromosome chr = population[population.size()-1];
            population[population.size()-1] = GnrChrIHEFT3(Rank_b,Rank_t);
            set<chromosome> Tem;
            Tem.insert(population.begin(),population.end());
            if(Tem.size() != Parameter_TDTSGA.NumOfChormPerPop){
                population[population.size()-1] = chr;
            }
        }
    }
    ClearALL();
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    return population[0].FitnessValue;
}
