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
//If use these filters to solve a complex data fitting term, define the data fitting as the blackbox function
float BlackBox(int row, int col, Mat& U, Mat & img_orig, float & d)
{
    //this is an example of adaptive norm
    float diff = fabs(U.at<float>(row,col)+d - img_orig.at<float>(row,col));
    float order = 2 - (fabs(U.at<float>(row+1,col) - U.at<float>(row,col)) + fabs(U.at<float>(row,col+1) - U.at<float>(row,col)));
    return pow(diff, order);
}

int main(int argc, char** argv)
{
    
    DM DualMesh;
    if ((argc<4) || (argc>6))
    {
       cout<<endl;
       cout<<" -------------------- Curvature Filter ------------------------- "<<endl;
       cout<<" Please cite Yuanhao's PhD thesis and related papers. Thank you! "<<endl;
       cout<<" --------------------------------------------------------------- \n\n";
       cout<<"usage: main imageName filterType Iterations.\nFor example: ./cf lena.bmp m 30\n";
       cout<<"             or              "<<endl;
       cout<<"usage: main imageName filterType MaxItNum lambda DataFitOrder.\nFor example: ./cf lena.bmp m 30 1.2 1.5\n";
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
    switch(*filterType)
    {
        case 't':
            { Type = 0; break; }
        case 'm':
            { Type = 1; break; }
        case 'g':
            { Type = 2; break; }
        case 'd':
            { Type = 3; break; }
        case 'b':
            { Type = 4; break; }
        default:
            { cout<<"Filter Type is NOT correct."<<endl; return -1; }
    }

    ItNum = atoi(argv[3]);
    double mytime;
    
    //just smooth the image by the filter
    if (argc==4)
    {
      DualMesh.Filter(Type, mytime, ItNum);
      cout<<"runtime is "<<mytime<<" milliseconds."<<endl;
      DualMesh.write();

      DualMesh.read(argv[1]);
      DualMesh.FilterNoSplit(Type, mytime, ItNum);
      cout<<"runtime (noSplit) is "<<mytime<<" milliseconds."<<endl;
      DualMesh.write("CF_NoSplit_result.png");
    }

    //solve a variational model (data fitting term is blackbox)
    if (argc==5)
    {
      lambda = (float)atof(argv[4]);
     
      DualMesh.read(argv[1]);
      DualMesh.BlackBoxSolver(Type, mytime, ItNum, lambda, BlackBox);
      cout<<"runtime is "<<mytime<<" milliseconds."<<endl;
      DualMesh.write("CF_BlackBox.png");
    }
    
    //solve a variational model
    if (argc==6)
    {
      lambda = (float)atof(argv[4]);
      DataFitOrder = (float)atof(argv[5]);
      //filter solver for the variational models
      DualMesh.read(argv[1]);
      DualMesh.Solver(Type, mytime, ItNum, lambda, DataFitOrder);
      cout<<"runtime is "<<mytime<<" milliseconds."<<endl;
      DualMesh.write("CF_Solver.png");
    }

    return 0;
}


