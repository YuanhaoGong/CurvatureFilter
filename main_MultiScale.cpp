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
int ItNum = 10;
int Type = 1;
int levels = 3;

#include "DM.h"
//filter on one layer
void oneLayer(Mat& img, int Iterations, int Type);
//construct upper level
Mat upLevel(Mat& img);
//backward to down level
void downLevel(Mat& low, Mat& high);

int main(int argc, char** argv)
{
    
    if (argc!=5)
    {
       cout<<"usage: main filename levels filterType Iterations.\n For example: ./mcf lena.bmp 3 t 30 \n";
       return -1;
    }
    Mat img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);   // Read the file
    Mat imgF = Mat::zeros(img.rows, img.cols, CV_32FC1);
    img.convertTo(imgF, CV_32FC1);

    levels = atoi(argv[2]);
    ItNum = atoi(argv[4]);

    char * filterType = argv[3];
    if (*filterType == 't') Type = 0;
    if (*filterType == 'm') Type = 1;
    if (*filterType == 'd') Type = 3;
    if (*filterType == 'b') Type = 4;

    

    //construct multiscale
    vector<Mat> multiScale;
    Mat tmp = imgF;

    for (int i = 0; i < levels-1; ++i)
    {
        multiScale.push_back(tmp);
        tmp = upLevel(tmp);
    }
    multiScale.push_back(tmp);
    
    //filter on each layer and pull back
	bool debug = false;
	if (debug) namedWindow( "mcf", 1 );
    for (int i = 0; i < levels-1; ++i)
    {
        oneLayer(multiScale.at(levels-i-1), ItNum, Type);
        downLevel(multiScale.at(levels-i-1), multiScale.at(levels-i-2));
		if (debug)
		{
			imshow("mcf",multiScale.at(levels-i-2)/255);
            cvWaitKey();
		}
    }
    oneLayer(multiScale.at(0), ItNum, Type);

	//write the result
    Mat result = multiScale.at(0);
    Mat tmp2 = Mat::zeros(result.rows, result.cols, CV_8UC1);
    result.convertTo(tmp2, CV_8UC1);

    vector<int> params;
    params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    params.push_back(0);

    imwrite("mcf.png", tmp2, params);

    return 0;
}

void oneLayer(Mat& img, int Iterations, int Type)
{
    DM test;
    test.set(img);
    test.split();
    double mytime;
    test.Filter(Type, mytime, Iterations);
    test.merge();
    img = test.get()*255;
	cout<<"One Level running time: "<<mytime<<" millisecond"<<endl;
}

//construct upper level
Mat upLevel(Mat& img)
{
    float a, b, c, d;
    Mat up = Mat::zeros(img.rows/2, img.cols/2, CV_32FC1);
    float* p, * p_nex, *p_up;
    for (int i = 0; i < img.rows-1; ++i,++i)
    {
        p = img.ptr<float>(i);
        p_nex = img.ptr<float>(i+1);
        p_up = up.ptr<float>(i/2);
        for (int j = 0; j < img.cols-1; ++j, ++j)
        {
            a = p[j];
            b = p[j+1];
            c = p_nex[j];
            d = p_nex[j+1];

            p[j] = (a+b+c+d)/4;
            p[j+1] = (a-b+c-d)/4;
            p_nex[j] = (a+b-c-d)/4;
            p_nex[j+1] = (a-b-c+d)/4;
            p_up[j/2] = p[j];
        }
    }
    return up;
}
//backward to down level
void downLevel(Mat& low, Mat& high)
{
    float *p, *p_nex, *p_low;
    float x, y, z, t;
    for (int i = 0; i < high.rows-1; ++i,++i)
    {
        p = high.ptr<float>(i);
        p_nex = high.ptr<float>(i+1);
        p_low = low.ptr<float>(i/2);

        for (int j = 0; j < high.cols-1; ++j, ++j)
        {
            x = p_low[j/2]; 
            y = p[j+1];
            z = p_nex[j];
            t = p_nex[j+1];

            p[j] = x + y + z + t;
            p[j+1] = x + z - y - t;
            p_nex[j] = x + y - z - t;
            p_nex[j+1] = x - y + t - z;
        }
    }
}

