#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "TDTSGA.h"

using namespace std;

int main() {
    srand((int) time(0));
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] = 0.935;
    SchTime["Montage25_0.7"] = 1.553;
    SchTime["Montage25_1.0"] = 2.143;
    SchTime["Montage50_0.4"] = 2.999;
    SchTime["Montage50_0.7"] = 6.400;
    SchTime["Montage50_1.0"] = 9.421 ;
    SchTime["Montage100_0.4"] =10.060 ;
    SchTime["Montage100_0.7"] =17.142  ;
    SchTime["Montage100_1.0"] =24.437  ;

    SchTime["CyberShake30_0.4"] = 1.151 ;
    SchTime["CyberShake30_0.7"] = 2.279 ;
    SchTime["CyberShake30_1.0"] = 2.677 ;
    SchTime["CyberShake50_0.4"] = 3.021 ;
    SchTime["CyberShake50_0.7"] = 5.694 ;
    SchTime["CyberShake50_1.0"] = 7.347  ;
    SchTime["CyberShake100_0.4"] =9.627  ;
    SchTime["CyberShake100_0.7"] =15.455  ;
    SchTime["CyberShake100_1.0"] =21.186  ;

    SchTime["Epigenomics24_0.4"] = 0.880 ;
    SchTime["Epigenomics24_0.7"] = 1.202 ;
    SchTime["Epigenomics24_1.0"] = 2.028 ;
    SchTime["Epigenomics47_0.4"] = 3.217 ;
    SchTime["Epigenomics47_0.7"] = 5.158 ;
    SchTime["Epigenomics47_1.0"] = 7.269 ;
    SchTime["Epigenomics100_0.4"] =9.199 ;
    SchTime["Epigenomics100_0.7"] =20.431    ;
    SchTime["Epigenomics100_1.0"] =33.882    ;

    SchTime["Ligo30_0.4"] = 1.257 ;
    SchTime["Ligo30_0.7"] = 1.689 ;
    SchTime["Ligo30_1.0"] = 2.901 ;
    SchTime["Ligo50_0.4"] = 3.040 ;
    SchTime["Ligo50_0.7"] = 4.592 ;
    SchTime["Ligo50_1.0"] = 6.345  ;
    SchTime["Ligo100_0.4"] =9.867  ;
    SchTime["Ligo100_0.7"] =14.676   ;
    SchTime["Ligo100_1.0"] =20.346   ;

    string strLine;
    ifstream iFile("../exp.txt");
    if (!iFile) {
        cout << "filelist open failed!\n";
        exit(1);
    }
    while(getline(iFile,strLine)){
        istringstream is(strLine);
        Orthogonal TemOrthogonal;
        is >> TemOrthogonal.NumOfChormPerPopRate >> TemOrthogonal.Mtn1 >> TemOrthogonal.Mtn2 >> TemOrthogonal.beta
           >> TemOrthogonal.RunTimeRatioOfStg1Rate;
        orthogonal.push_back(TemOrthogonal);
    }
    iFile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string strLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, strLine);
        if (strLine.size() < 1) {
            cout << "Empty input file(fileList)" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(strLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        int index =0 ;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../print/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(5);
            for (int times = 0; times < 10; ++times) {
                double TDTSGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
                int TDTSGA_Iteration = 0;
                double TDTSGA_Result = runTDTSGA(XmlFile, RscAlcFile, TemOrthogonal, TDTSGA_SchTime, TDTSGA_Iteration );
                outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                        << TDTSGA_Result << " " << TDTSGA_SchTime << " " << TDTSGA_Iteration
                        << endl;
            }
            outfile.close();
        }
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    return 0;
}
