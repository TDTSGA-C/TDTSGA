
#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

void GnrTskSchLst(chromosome& chrom , vector<int>& TskSchLst);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SltChr(vector<double>& A);
void SlcTwoChrom(vector<double>& A,chromosome& Chrom1,chromosome& Chrom2 );
double Dcd(chromosome& ch, bool IsFrw);
void CrsLvl(chromosome& chrom1, chromosome& chrom2);
void CrsTS(chromosome& chrom1, chromosome& chrom2);
void Mtn(chromosome& chrom, double Pm, int stg);

#endif //CSTCHANGE_CROSSOVER_H
