#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "CGA.h"
#include "LWSGA.h"
#include "HGA.h"
#include "IHEFT3.h"
#include "TDTSGA.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm}
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

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(0);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        double IHEFT3_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        double IHEFT3_Result = runIHEFT3(XmlFile, RscAlcFile,IHEFT3_SchTime);

        double HGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);

        double LWSGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);

        double CGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int CGA_Iteration = 0;
        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);

        double TDTSGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int TDTSGA_Iteration = 0;
        double TDTSGA_Result = runTDTSGA(XmlFile, RscAlcFile, TDTSGA_SchTime, TDTSGA_Iteration);



//        results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << IHEFT3_Result << " " << IHEFT3_SchTime << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << TDTSGA_Result  << " " << TDTSGA_SchTime  << " " << TDTSGA_Iteration  << " "
                << endl;
        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    return 0;
}
