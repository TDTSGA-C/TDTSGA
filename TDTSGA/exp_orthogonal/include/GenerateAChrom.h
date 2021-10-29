

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_t(vector<double>& rankList, vector<double>& w, vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime);
chromosome GnrChr_Lvl_Rnd();
chromosome GnrChr_Tpl_Rnd();
chromosome GnrChr_HEFT(vector<double> Rank_b);
chromosome GnrChr_Lvl_EFT();
chromosome GnrChr_Tpl_EFT();
double IHEFT3_S(chromosome& ch);
chromosome GnrChrIHEFT3(vector<double> Rank_b,vector<double> Rank_t);
#endif //CSTCHANGE_GENERATEACHROM_H
