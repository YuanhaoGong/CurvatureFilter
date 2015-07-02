#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <cmath>

using namespace cv;
using namespace std;

//default filter and iteration number
int ItNum = 10;
int Type = 2;

#include "DM.h"

int main(int argc, char** argv)
{
    
    DM DualMesh;
    if (argc!=4)
    {
       cout<<"usage: main filename filterType Iterations.\n For example: ./cf lena.bmp m 30, where m means Mean Curvature Filter. \n";
       return -1;
    }
    DualMesh.read(argv[1]);
    ItNum = atoi(argv[3]);

    char * filterType = argv[2];
    if (*filterType == 't') Type = 0;
    if (*filterType == 'm') Type = 1;
    if (*filterType == 'd') Type = 3;


    DualMesh.split();
    double mytime;
    DualMesh.Filter(Type, mytime, ItNum);
    cout<<"runtime is "<<mytime<<" milliseconds."<<endl;

    DualMesh.merge();
    DualMesh.write();

    DualMesh.read(argv[1]);
    DualMesh.FilterNoSplit(Type, mytime, ItNum);
    cout<<"runtime (noSplit) is "<<mytime<<" milliseconds."<<endl;
    DualMesh.write("CF_NoSplit_result.png");
    
    return 0;
}


