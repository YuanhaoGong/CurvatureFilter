/*=========================================================================
 *
 *                           Curvature Filter 
 *
 **************************************************************************  
 
            @phdthesis{gong:phd, 
             title={Spectrally regularized surfaces}, 
             author={Gong, Yuanhao}, 
             year={2015}, 
             school={ETH Zurich, Nr. 22616},
             note={http://dx.doi.org/10.3929/ethz-a-010438292}}

 *=========================================================================*/

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
float lambda = 1.0f;
float DataFitOrder = 1.0f;

#include "DM.h"

int main(int argc, char** argv)
{
    
    DM DualMesh;
    if ((argc!=4) && (argc!=6))
    {
       cout<<endl;
       cout<<" -------------------- Curvature Filter ------------------------- "<<endl;
       cout<<" Please cite Yuanhao's PhD thesis and related papers. Thank you! "<<endl;
       cout<<" --------------------------------------------------------------- \n\n";
       cout<<"usage: main imageName filterType Iterations.\n For example: ./cf lena.bmp m 30\n";
       cout<<"             or              "<<endl;
       cout<<"usage: main imageName filterType MaxItNum lambda DataFitOrder.\n For example: ./cf lena.bmp m 30 1.2 1.5\n";
       cout<<"************************************************\n";
       cout<<"Possible Filter Type: t (Total Variation) \n";
       cout<<"                      m (Mean Curvature) \n";
       cout<<"                      d (Difference Curvature) \n";
       cout<<"                      g (Gaussian Curvature) \n";
       cout<<"                      b (Bernstein Filter) \n";
       return -1;
    }

    DualMesh.read(argv[1]);

    char * filterType = argv[2];
    if (*filterType == 't') Type = 0;
    if (*filterType == 'm') Type = 1;
    if (*filterType == 'd') Type = 3;
    if (*filterType == 'b') Type = 4;

    ItNum = atoi(argv[3]);

    if (argc==6)
    {
    	lambda = (float)atof(argv[4]);
    	DataFitOrder = (float)atof(argv[5]);
    }

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
    
    if (argc==6)
    {
        //filter solver for the variational models
        DualMesh.read(argv[1]);
	DualMesh.Solver(Type, mytime, ItNum, lambda, DataFitOrder);
	cout<<"runtime is "<<mytime<<" milliseconds."<<endl;
	DualMesh.write("CF_Solver.png");
    }
    return 0;
}


