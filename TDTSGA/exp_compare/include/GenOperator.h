//
// Created by smallfish on 18-4-4.
//

#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H
void GnrTskSchLst(chromosome& chrom , vector<int>& TskSchLst);
double DcdEvl(chromosome& chrom, bool isForward);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SltChr(vector<double>& A);
void SlcTwoChrom(vector<double>& A,chromosome& Chrom1,chromosome& Chrom2 );
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void MtnSS_TS(chromosome& a);
void MtnMS_SP(chromosome& a);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);
double Dcd(chromosome& ch, bool IsFrw);
void CrsLvl(chromosome& chrom1, chromosome& chrom2);
void CrsTS(chromosome& chrom1, chromosome& chrom2);
void Mtn(chromosome& Chrom, double Pm, int stg);
#endif //CSTCHANGE_CROSSOVER_H
