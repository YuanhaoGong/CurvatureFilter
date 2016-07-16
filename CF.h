/*=========================================================================
 *
 *                           Curvature Filter 
 *
 *                             Yuanhao Gong
 *                        gongyuanhao@gmail.com
 **************************************************************************  
 
@phdthesis{gong:phd, 
 title={Spectrally regularized surfaces}, 
 author={Gong, Yuanhao}, 
 year={2015}, 
 school={ETH Zurich, Nr. 22616},
 note={http://dx.doi.org/10.3929/ethz-a-010438292}}

 *=========================================================================*/

//Curvature Filter
class CF
{
public:
    /********************* basic IO *********************/
    //read one image from disk
    void read(const char* FileName);
    //set one image from memory
    void set(Mat& src);
    //get the filtered image 
    Mat get(){return imgF;};
    //write the result to disk
    void write();
    void write(const char* FileName);

    /********************* PSNR and Ns *********************/
    //PSNR
    double PSNR();
    double PSNR(const Mat& I1, const Mat& I2);
    //compute naturalness factor
    double Naturalness(){return Naturalness(imgF);}
    double Naturalness(const Mat & img);

    /********************* curvature and energy *********************/
    //compute TV: scheme=1, L1 norm; scheme=2, L2 norm
    void TV(const Mat & img, Mat & T, int scheme=1);
    //compute MC: 0, standard; 1, my scheme; 2, quadratic fitting; 3, my scheme; 4, isoline
    void MC(const Mat & img, Mat & MC, int scheme=0);
    //compute GC: 0, standard; 1, my scheme; 2, quadratic fitting
    void GC(const Mat & img, Mat & GC, int scheme=0);

    //compute energy for given TV, MC, or GC image by L1 norm
    double energy(const Mat& img, const int order=1);
    //compute data fitting energy between image and imgF
    double DataFitEnergy(Mat& tmp, double order);

    /**************************** curvature filters ******************************************
    Type=0, TV; Type=1, MC; Type=2, GC; Type=3, DC; Type=4, Bernstein;
    the stepsize parameter is in (0,1], stochastic step is in (0,stepsize)
    *********************************************************************************************/
    void Filter(const int Type, double & time, const int ItNum = 10, const float stepsize=1);
    void FilterNoSplit(const int Type, double & time, const int ItNum = 10, const float stepsize=1, const bool stochastic=false);

    /*********************************************************************************************
    ******************* generic solver for variational models *****************************
    solve variational models by the curvature filter, just guarantee to reduce the total energy
    and convergence, but not necessarily to a local minimum 
    *********************************************************************************************/
    //solve |U - I|^DataFitOrder + lambda * |curvature(U)|
    void Solver(const int Type, double & time, const int MaxItNum, 
                const float lambda = 2, const float DataFitOrder = 1, const float stepsize=1);
    //solve BlackBox(U,I) + lambda * |curvature(U)|
    void BlackBoxSolver(const int Type, double & time, const int MaxItNum, const float lambda, 
                            float (*BlackBox)(int row, int col, Mat& U, Mat & img_orig, float & d), const float stepsize=1);
    
    /*********************************************************************************************
    ***** the dm at given location, which is the gradient of regularization energy ****
    *********************************************************************************************/
    //Th result dm is only valid for the location type, NOT the full image
    //LocationType        start_row         start_col
    //   0 (BC),                  1,               1;
    //   1 (BT),                  2,               2;
    //   2 (WC),                  1,               2;
    //   3 (WT),                  2,               1;
    void DM(const int FilterType, const int LocationType, const Mat & img, Mat & dm);
    
    /*********************************************************************************************
    ****************************    Curvature Guided Filter   ************************************
    *********************************************************************************************/
    //compute the curvature from the guided image (scaled to the size of imgF)
    Mat GuideCurvature(const char * FileName, const int Type);
    //filter the image such that the result is close to the specified curvature
    void CurvatureGuidedFilter(const Mat & curv, const int Type, double & time, const int ItNum = 10, const float stepsize=1);

    /********************************************************************************************
    *********************************************************************************************
    *********************************  end of public functions   ********************************
    *********************************************************************************************
    *********************************************************************************************/
private:
    //padded original, actural float image, result
    Mat image, imgF, result;
    //four sets, be aware that position is fixed, see split() or Yuanhao Gong's PhD thesis
    Mat WC, WT, BC, BT;
    //image size
    int M, N, M_orig, N_orig, M_half, N_half;
    //six pointers
    float* p, *p_right, *p_down, *p_rd, *p_pre, *p_Corner;
    //pointer to the data
    const float* p_data;
private:
    //split imgF into four sets
    void split();
    //merge four sets back to imgF
    void merge();
    //naturalness evaluation
    double Naturalness_search(float* data, int N, int offset);
    //computing curvature by different schemes
    void MC_fit(const Mat & img, Mat & MC);
    void GC_new(const Mat & img, Mat & GC);
    void GC_fit(const Mat & img, Mat & GC);
    Mat GC_LUT_Init();
    void GC_LUT(const Mat & LUT, const Mat & img, Mat & GC);
    //fit coefficients for quad function
    void FiveCoefficient(const Mat & img, Mat & x2, Mat &y2, Mat & xy, Mat & x, Mat & y);
    
    //find the signed value with minimum abs value, dist contains FOUR floats
    inline float SignedMin(float * dist);
    inline float SignedMin_noSplit(float * dist);

    /*************************************** Split into 4 sets *********************************/
    //one is for BT and WC, two is for BC and WT
    inline void GC_one(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void GC_two(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);

    inline void MC_one(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void MC_two(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void TV_one(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void TV_two(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);

    inline void DC_one(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void DC_two(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void LS_one(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);
    inline void LS_two(float* p, const float* right, const float* down, const float *rd, const float* pre, const float* corner, const float& stepsize);

    /*************************************** Direct on imgF (no split) ***********************/
    inline float Scheme_GC(int i, const float * p_pre, const float * p, const float * p_nex, const float * p_guide = NULL);
    inline float Scheme_MC(int i, const float * p_pre, const float * p, const float * p_nex, const float * p_guide = NULL);
    inline float Scheme_TV(int i, const float * p_pre, const float * p, const float * p_nex, const float * p_guide = NULL);
    inline float Scheme_DC(int i, const float * p_pre, const float * p, const float * p_nex, const float * p_guide = NULL);
    inline float Scheme_LS(int i, const float * p_pre, const float * p, const float * p_nex, const float * p_guide = NULL);
};

/********************************************************************************************
*********************************************************************************************
*********************************  end of CF class   **************************************
*********************************************************************************************
*********************************************************************************************/

double CF::PSNR()
{
    return PSNR(image, imgF);
}

double CF::PSNR(const Mat& I1, const Mat& I2)
{
     Mat diff;
     absdiff(I1, I2, diff);           // |I1 - I2|
     diff = diff.mul(diff);           // |I1 - I2|^2

     Scalar d = sum(diff);            // sum elements per channel
     double sse = d.val[0];

     if( sse <= 1e-10)               // for small values return zero
         return 0;
     else
     {
         double mse =sse /(double)(I1.total());
         double psnr = - 10*log10(mse);
         return psnr;
     }
}

double CF::Naturalness(const Mat& imgF)
{
    //compute naturalness factor
    const float * p_row ;
    const float * pp_row;
    int indexX, indexY;
    int Offset = 256;
    int N = 512;
    double eps = 0.0001;
    double left(0), right(1), mid_left(0), mid_right(0);

    Mat GradCDF = Mat::zeros(2, N, CV_32FC1);
    float * Gradx = GradCDF.ptr<float>(0);
    float * Grady = GradCDF.ptr<float>(1);
    float f = 1.0f/((imgF.rows-1)*(imgF.cols-1));

    //not efficient but safe way
    for(int i = 0; i < imgF.rows - 1; i++)
    {
        p_row = imgF.ptr<float>(i);
        pp_row = imgF.ptr<float>(i+1);

        for(int j = 0; j < imgF.cols - 1; j++)
        {
            //scale back to 255
            indexX = Offset + int((p_row[j+1] - p_row[j])*255);
            indexY = Offset + int((pp_row[j] - p_row[j])*255);            
            Gradx[indexX] += f;
            Grady[indexY] += f;
        }
    }
    //convert Grad PDF into CDF
    for (int j = 1; j < N; ++j)
    {
        Gradx[j] += Gradx[j-1];
        Grady[j] += Grady[j-1];
    }
    //scale the data
    GradCDF -= 0.5f;
    GradCDF *= 3.14159f;

    //search parameter T for x component
    double Natural_x = Naturalness_search(Gradx, N, Offset);
    double Natural_y = Naturalness_search(Grady, N, Offset);

    //the final naturalness factor
    return (Natural_x+Natural_y)/(0.7508);
}

double CF::Naturalness_search(float* data, int N, int offset)
{
    //Ternary search
    float * p_d, *p_d2;
    double left(0), right(1), mid_left, mid_right;
    double eps(0.0001), tmp, tmp2, error_left, error_right;
    while(right-left>=eps)
    {
        mid_left=left+(right-left)/3;
        mid_right=right-(right-left)/3;

        error_left = 0; error_right = 0;
        p_d = data; p_d2 = data;
        for (int i=-offset+1; i<N-offset; ++i) 
        {
            tmp = atan(mid_left*(i)) - (*p_d++);
            tmp2 = atan(mid_right*(i)) - (*p_d2++);
            error_left += (tmp*tmp);
            error_right += (tmp2*tmp2);
        }

        if(error_left <= error_right)
            right=mid_right;
        else
            left=mid_left;
    }
    return (mid_left+mid_right)/2;
}

void CF::read(const char* FileName)
{
    //load the image and convert to float, pad to even size
    //the Dirichlet boundary is used.
    Mat tmp = imread(FileName, CV_LOAD_IMAGE_GRAYSCALE);   
    if(!tmp.data )                               
        {
            cout <<  "*********** Read Image Failed! ************" << std::endl ;
            return ;
        }

    Mat tmp2 = Mat::zeros(tmp.rows, tmp.cols, CV_32FC1);
    tmp.convertTo(tmp2, CV_32FC1);
    M_orig = tmp2.rows;
    N_orig = tmp2.cols;
    M = (int)ceil(M_orig/2.0)*2;
    N = (int)ceil(N_orig/2.0)*2;
    M_half = M/2;
    N_half = N/2;
    
    tmp2 /= 255.0f;

    image = Mat::zeros(M,N,CV_32FC1);
    tmp2.copyTo(image(Range(0,M_orig),Range(0,N_orig)));
    image.col(N-1) = image.col(N-2);
    image.row(M-1) = image.row(M-2);
    image.at<float>(M-1, N-1) = image.at<float>(M-2, N-2);

    //set the imgF
    imgF = Mat::zeros(M,N,CV_32FC1);
    image.copyTo(imgF);
}

void CF::set(Mat& file)
{
    //instead of loading image, it can be set directly from the memory.
    //the Dirichlet boundary is used.
    
    Mat tmp2 = Mat::zeros(file.rows, file.cols, CV_32FC1);
    file.convertTo(tmp2, CV_32FC1);
    M_orig = tmp2.rows;
    N_orig = tmp2.cols;
    M = (int)ceil(M_orig/2.0)*2;
    N = (int)ceil(N_orig/2.0)*2;
    M_half = M/2;
    N_half = N/2;
    
    tmp2 /= 255.0f;

    image = Mat::zeros(M,N,CV_32FC1);
    tmp2.copyTo(image(Range(0,M_orig),Range(0,N_orig)));
    image.col(N-1) = image.col(N-2);
    image.row(M-1) = image.row(M-2);
    image.at<float>(M-1, N-1) = image.at<float>(M-2, N-2);

    //set the imgF
    imgF = Mat::zeros(M,N,CV_32FC1);
    image.copyTo(imgF);
}
void CF::write()
{
    CF::write("CF.png");
}

void CF::write(const char* FileName)
{
    Mat tmp = Mat::zeros(M_orig,N_orig,CV_8UC1);
    Mat tmp2 = imgF*255.0f;
    tmp2(Range(0,M_orig),Range(0,N_orig)).convertTo(tmp, CV_8UC1);

    vector<int> params;
    params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    params.push_back(0);

    imwrite(FileName, tmp, params);
}

//compute Total Variation
void CF::TV(const Mat & imgF, Mat & T, int type)
{
    const float * p_row, * pn_row;
    float * p_t;
    switch(type)
    {
        case 1://default using L1 norm
            for(int i = 1; i < imgF.rows-1; i++)
            {
                p_row = imgF.ptr<float>(i);
                pn_row = imgF.ptr<float>(i+1);
                p_t = T.ptr<float>(i);
                for(int j = 1; j < imgF.cols-1; j++)
                {
                    p_t[j] = fabsf(p_row[j+1] - p_row[j]) + fabsf(pn_row[j] - p_row[j]);
                }   
            }
            break;
        case 2:
            float gx, gy;
            for(int i = 1; i < imgF.rows-1; i++)
            {
                p_row = imgF.ptr<float>(i);
                pn_row = imgF.ptr<float>(i+1);
                p_t = T.ptr<float>(i);
                for(int j = 1; j < imgF.cols-1; j++)
                {
                    gx = p_row[j+1] - p_row[j]; gy = pn_row[j] - p_row[j];
                    p_t[j] = sqrt(gx*gx + gy*gy);
                }   
            }
            break;
        default:
        cout<<"The type in TV is not correct."<<endl; 
        break;
    }
}

//compute Mean Curvature
void CF::MC(const Mat& imgF, Mat & MC, int type)
{
    //classical scheme is used by default
    const float * p_row, *pn_row, *pp_row;
    float *p_d;
    float Ix, Iy, Ixy, Ixx, Iyy, num, den, tmp;
    Mat kernel;
    switch(type)
    {
        case 0: //standard scheme
            for(int i = 1; i < imgF.rows-1; i++)
            {
                p_row = imgF.ptr<float>(i);
                pn_row = imgF.ptr<float>(i+1);
                pp_row = imgF.ptr<float>(i-1);
                p_d = MC.ptr<float>(i);
                
                for(int j = 1; j < imgF.cols-1; j++)
                {
                    Ix = (p_row[j+1] - p_row[j-1])/2;
                    Iy = (pn_row[j] - pp_row[j])/2;
                    Ixx = p_row[j+1] - 2*p_row[j] + p_row[j-1];
                    Iyy = pn_row[j] - 2*p_row[j] + pp_row[j];
                    Ixy = (pn_row[j-1] - pn_row[j+1]- pp_row[j-1] + pp_row[j+1])/4;
                    
                    num = (1+Ix*Ix)*Iyy - 2*Ix*Iy*Ixy + (1+Iy*Iy)*Ixx;
                    tmp = 1.0f + Ix*Ix + Iy*Iy;
                    den = sqrt(tmp)*tmp*2;
                    p_d[j] = num/den;
                }   
            }
            break;
        case 1: //separate kernel from Eq.6.12, add the center pixel later
            kernel = (Mat_<float>(1,3) << 0.25f, -1.25f, 0.25f); 
            sepFilter2D(imgF, MC, CV_32F, kernel, kernel,Point(-1,-1),0,BORDER_REPLICATE);
            MC *= -1;
            MC += (0.5625f*imgF);//add the center pixel
            break;
        case 2: //fit by a quadratic function
            MC_fit(imgF, MC); 
            break;
        case 3: //another linear kernel
            Laplacian(imgF, MC, MC.depth(), 1, 0.25);
            break;
        case 4: //isoline schemes, be aware the difference
            for(int i = 1; i < imgF.rows-1; i++)
            {
                p_row = imgF.ptr<float>(i);
                pn_row = imgF.ptr<float>(i+1);
                pp_row = imgF.ptr<float>(i-1);
                p_d = MC.ptr<float>(i);
                
                for(int j = 1; j < imgF.cols-1; j++)
                {
                    Ix = (p_row[j+1] - p_row[j-1])/2;
                    Iy = (pn_row[j] - pp_row[j])/2;
                    Ixx = p_row[j+1] - 2*p_row[j] + p_row[j-1];
                    Iyy = pn_row[j] - 2*p_row[j] + pp_row[j];
                    Ixy = (pn_row[j-1] - pn_row[j+1]- pp_row[j-1] + pp_row[j+1])/4;
                    
                    num = Ix*Ix*Iyy - 2*Ix*Iy*Ixy + Iy*Iy*Ixx;
                    den = Ix*Ix + Iy*Iy;
                    den = sqrt(den)*den + 0.000001f;
                    p_d[j] = num/den;
                }   
            }
            break;
        default:
        cout<<"The type in MC is not correct."<<endl; 
        break;
    }
}

//fit coefficients 
void CF::FiveCoefficient(const Mat & img, Mat & x2, Mat &y2, Mat & xy, Mat & x, Mat & y)
{
    Mat kernel_one = (Mat_<float>(1,3) << 1.0f, 1.0f, 1.0f);
    Mat kernel_one_h = (Mat_<float>(1,3) << 0.1666667f, -0.333333f, 0.1666667f);
    Mat kernel_three_h = (Mat_<float>(1,3) << -0.25f, 0.0f, 0.25f);
    Mat kernel_three_v = (Mat_<float>(1,3) << 1.0f, 0.0f, -1.0f);
    Mat kernel_four_h = - kernel_three_v;
    Mat kernel_four_v = (Mat_<float>(1,3)<<0.1666667f,0.1666667f,0.1666667f);

    sepFilter2D(img, x2, img.depth(), kernel_one_h, kernel_one);
    sepFilter2D(img, y2, img.depth(), kernel_one, kernel_one_h);
    sepFilter2D(img, xy, img.depth(), kernel_three_h, kernel_three_v);
    sepFilter2D(img, x, img.depth(), kernel_four_h, kernel_four_v);
    sepFilter2D(img, y, img.depth(), kernel_four_v, kernel_three_v);//reuse the same kernel
}

//compute Mean Curvature by fitting quad function
void CF::MC_fit(const Mat & img, Mat & MC)
{
    Mat x2, y2, xy, x, y;
    x2 = Mat::zeros(img.rows, img.cols, CV_32FC1);
    y2 = Mat::zeros(img.rows, img.cols, CV_32FC1);
    xy = Mat::zeros(img.rows, img.cols, CV_32FC1);
    x  = Mat::zeros(img.rows, img.cols, CV_32FC1);
    y  = Mat::zeros(img.rows, img.cols, CV_32FC1);

    FiveCoefficient(img, x2, y2, xy, x, y);
    float* p_x2, *p_y2, *p_xy, *p_x, *p_y, *p_d;
    float num, den;
    for (int i = 1; i < img.rows-1; ++i)
    {
        p_x2 = x2.ptr<float>(i);
        p_y2 = y2.ptr<float>(i);
        p_xy = xy.ptr<float>(i);
        p_x  = x.ptr<float>(i);
        p_y  = y.ptr<float>(i);
        p_d  = MC.ptr<float>(i);
        for (int j = 1; j < img.cols-1; ++j)
        {
            num = (1+p_x[j]*p_x[j])*p_y2[j] - p_x[j]*p_y[j]*p_xy[j] + (1+p_y[j]*p_y[j])*p_x2[j];
            den = 1+p_x[j]*p_x[j] + p_y[j]*p_y[j];
            den = sqrt(den)*den;//no multiply 2 here
            p_d[j] = num/den; 
        }
    }
}

//compute Gaussian Curvature by fitting quad function
void CF::GC_fit(const Mat & img, Mat & GC)
{
    Mat x2, y2, xy, x, y;
    x2 = Mat::zeros(img.rows, img.cols, CV_32FC1);
    y2 = Mat::zeros(img.rows, img.cols, CV_32FC1);
    xy = Mat::zeros(img.rows, img.cols, CV_32FC1);
    x  = Mat::zeros(img.rows, img.cols, CV_32FC1);
    y  = Mat::zeros(img.rows, img.cols, CV_32FC1);

    FiveCoefficient(img, x2, y2, xy, x, y);
    float* p_x2, *p_y2, *p_xy, *p_x, *p_y, *p_d;
    float num, den;
    for (int i = 1; i < img.rows-1; ++i)
    {
        p_x2 = x2.ptr<float>(i);
        p_y2 = y2.ptr<float>(i);
        p_xy = xy.ptr<float>(i);
        p_x  = x.ptr<float>(i);
        p_y  = y.ptr<float>(i);
        p_d  = GC.ptr<float>(i);
        for (int j = 1; j < img.cols-1; ++j)
        {
            num = p_x2[j] * p_y2[j] *4 - p_xy[j]*p_xy[j];
            den = 1+p_x[j]*p_x[j] + p_y[j]*p_y[j];
            den *= den;
            p_d[j] = num/den; 
        }
    }
}

//compute Gaussian curvature
void CF::GC(const Mat & imgF, Mat &GC, int type)
{
    //classical scheme is used by default
    const float * p_row, *pn_row, *pp_row;
    float *p_d;
    float Ix, Iy, Ixx, Iyy, Ixy, num, den;
    switch(type)
    {
        case 0: //classical schemes
            for(int i = 1; i < imgF.rows-1; i++)
            {
                p_row = imgF.ptr<float>(i);
                pn_row = imgF.ptr<float>(i+1);
                pp_row = imgF.ptr<float>(i-1);
                p_d = GC.ptr<float>(i);
                
                for(int j = 1; j < imgF.cols-1; j++)
                {
                    Ix = (p_row[j+1] - p_row[j-1])/2;
                    Iy = (pn_row[j] - pp_row[j])/2;
                    Ixx = p_row[j+1] - 2*p_row[j] + p_row[j-1];
                    Iyy = pn_row[j] -2*p_row[j] + pp_row[j];
                    Ixy = (pn_row[j-1] - pn_row[j+1]- pp_row[j-1] + pp_row[j+1])/4;

                    num = Ixx*Iyy - Ixy*Ixy;
                    den = (1.0f + Ix*Ix + Iy*Iy);
                    den *= den;
                    p_d[j] = num/den;
                }   
            }
            break;
        case 1: // the new scheme from my thesis
            GC_new(imgF, GC); 
            break;
        case 2: // fit the surface with a quadratic function
            GC_fit(imgF, GC); 
            break;
        default:
        cout<<"The type in GC is not correct."<<endl; 
        break;
    }
}

//new scheme, Eq.6.16 in my thesis
void CF::GC_new(const Mat & img, Mat & GC)
{
    Mat tmp = Mat::zeros(M,N,CV_32FC1);
    
    //six kernel from Eq.6.16
    Mat kernel_one = (Mat_<float>(3,3) << 
        0.002104f, 0.254187f, 0.002104f,
        0.254187f, -1.02516f, 0.254187f,
        0.002104f, 0.254187f, 0.002104f);
    Mat kernel_two = (Mat_<float>(3,3) << 
        0.25f, 0.0f, -0.25f,
        -0.5f, 0.0f, 0.5f,
        0.25f, 0.0f, -0.25f);
    Mat kernel_three = (Mat_<float>(3,3) << 
        -0.25f, 0.5f, -0.25f,
        0.0f, 0.0f, 0.0f,
        0.25f, -0.5f, 0.25f);
    Mat kernel_four = (Mat_<float>(3,3) << 
        -0.286419f, 0.229983f, -0.286419f,
        0.229983f, 0.225743f, 0.229983f,
        -0.286419f, 0.229983f, -0.286419f);
    Mat kernel_five = (Mat_<float>(3,3) << 
        -0.306186f, 0.0f, 0.306186f,
        0.0f, 0.0f, 0.0f,
        0.306186f, 0.0f, -0.306186f);
    Mat kernel_six = (Mat_<float>(3,3) << 
        0.0f, -0.176777f, 0.0f,
        0.176777f, 0.0f, 0.176777f,
        0.0f, -0.176777f, 0.0f);

    filter2D(img, GC, CV_32F, kernel_one);
    pow(GC, 2, GC);
    filter2D(img, tmp, CV_32F, kernel_two);
    pow(tmp, 2, tmp);
    GC -= tmp;
    filter2D(img, tmp, CV_32F, kernel_three);
    pow(tmp, 2, tmp);
    GC -= tmp;
    filter2D(img, tmp, CV_32F, kernel_four);
    pow(tmp, 2, tmp);
    GC -= tmp;
    filter2D(img, tmp, CV_32F, kernel_five);
    pow(tmp, 2, tmp);
    GC -= tmp;
    filter2D(img, tmp, CV_32F, kernel_six);
    pow(tmp, 2, tmp);
    GC -= tmp;
}

Mat CF::GC_LUT_Init()
{
    const int scale = 255*255;
    const int scale2 = 2*scale;
    int num;
    float den; 

    Mat LUT = Mat::zeros(512, 512, CV_32FC1);
    float * p_LUT;
    for (int i = -255; i < 1; ++i)
    {
        p_LUT = LUT.ptr<float>(i+255);
        for (int j = -255; j < 255; ++j)
        {
            num = i*j + scale;
            den = sqrt(float(i*i+scale))*sqrt(float(j*j+scale2));//avoid integer overflow
            p_LUT[j+255] = num/den;
            if (abs(p_LUT[j+255])>1)
            {
                cout<<i<<" "<<j<<" "<<p_LUT[j+255]<<endl;
            }
            p_LUT[j+255] = acos(p_LUT[j+255]);
            LUT.at<float>(255-i,j+255) = p_LUT[j+255];
        }
    }
    return LUT;
}

//approimate GC by LUT
void CF::GC_LUT(const Mat & LUT, const Mat & img, Mat & GC)
{
    const float * p_row, *pp_row, *pn_row, *p_LUT;
    float *g;
    float total;
    int row[4], col[4], offset;
    for (int i = 1; i < img.rows-1; ++i)
    {
        p_row = img.ptr<float>(i);
        pp_row = img.ptr<float>(i-1);
        pn_row = img.ptr<float>(i+1);
        g = GC.ptr<float>(i);
        for (int j = 1; j < img.cols-1; ++j)
        {
            total = 0;
            offset = int(255 - 255*p_row[j]);
            row[0] = int(pp_row[j]*255) + offset;
            col[0] = int(pp_row[j-1]*255) + offset;
            col[1] = int(pp_row[j+1]*255) + offset;
            row[1] = int(p_row[j-1]*255) + offset;
            row[2] = int(p_row[j+1]*255) + offset;
            row[3] = int(pn_row[j]*255) + offset;
            col[2] = int(pn_row[j-1]*255) + offset;
            col[3] = int(pn_row[j+1]*255) + offset;
            //top two
            p_LUT = LUT.ptr<float>(row[0]);
            total += (p_LUT[col[0]] + p_LUT[col[1]]);
            //left two
            p_LUT = LUT.ptr<float>(row[1]);
            total += (p_LUT[col[0]] + p_LUT[col[2]]);
            //right two
            p_LUT = LUT.ptr<float>(row[2]);
            total += (p_LUT[col[1]] + p_LUT[col[3]]);
            //bottom two
            p_LUT = LUT.ptr<float>(row[3]);
            total += (p_LUT[col[2]] + p_LUT[col[3]]);

            g[j] = 6.283185f - total;
        }
    }
}

//compute the curvature energy
double CF::energy(const Mat &img, const int order)
{
    switch(order)
    {
        case 1:
        {
            Scalar tmp = sum(cv::abs(img));
            return tmp(0);
        }
        case 2:
        {
            Scalar tmp = sum(cv::abs(img));
            return tmp(0);
        }
        default:
        cout <<"the order in curvature energy is not correct."<<endl;
        return 0;
    }
}

//compute the energy between image and imgF
double CF::DataFitEnergy(Mat & tmp, double order)
{
    tmp = abs(image - imgF);
    pow(tmp, order, tmp);
    Scalar tmp2 = sum(tmp);
    return tmp2(0);
}

//split the image into four sets
void CF::split()
{
    WC = Mat::zeros(M_half,N_half,CV_32FC1);
    WT = Mat::zeros(M_half,N_half,CV_32FC1);
    BC = Mat::zeros(M_half,N_half,CV_32FC1);
    BT = Mat::zeros(M_half,N_half,CV_32FC1);

    float *p;
    for (int i = 0; i < imgF.rows; ++i)
    {
        p=imgF.ptr<float>(i);
        for (int j = 0; j < imgF.cols; ++j)
        {
            if (i%2==0 && j%2==0) BT.at<float>(i/2,j/2) = p[j];
            if (i%2==0 && j%2==1) WT.at<float>(i/2,(j-1)/2) = p[j];
            if (i%2==1 && j%2==0) WC.at<float>((i-1)/2,j/2) = p[j];
            if (i%2==1 && j%2==1) BC.at<float>((i-1)/2,(j-1)/2) = p[j];
        }
    }
}

//merge the four sets into one image
void CF::merge()
{
    float *p;
    for (int i = 0; i < imgF.rows; ++i)
    {
        p=imgF.ptr<float>(i);
        for (int j = 0; j < imgF.cols; ++j)
        {
            if (i%2==0 && j%2==0) p[j] = BT.at<float>(i/2,j/2);
            if (i%2==0 && j%2==1) p[j] = WT.at<float>(i/2,(j-1)/2);
            if (i%2==1 && j%2==0) p[j] = WC.at<float>((i-1)/2,j/2);
            if (i%2==1 && j%2==1) p[j] = BC.at<float>((i-1)/2,(j-1)/2);
        }
    }
}

void CF::Filter(const int Type, double & time, const int ItNum, const float stepsize)
{
    clock_t Tstart, Tend;
    //split imgF into four sets
    split();

    void (CF::* Local_one)(float* p, const float* p_right, const float* p_down, const float *p_rd, 
                            const float* p_pre, const float* p_Corner, const float& stepsize);
    void (CF::* Local_two)(float* p, const float* p_right, const float* p_down, const float *p_rd, 
                            const float* p_pre, const float* p_Corner, const float& stepsize);

    switch(Type)
    {
        case 0:
        {
            Local_one = &CF::TV_one; Local_two = &CF::TV_two; 
            cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local_one = &CF::MC_one; Local_two = &CF::MC_two; 
            cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local_one = &CF::GC_one; Local_two = &CF::GC_two; 
            cout<<"GC Filter: "; break;
        }
        case 3:
        {
            Local_one = &CF::DC_one; Local_two = &CF::DC_two; 
            cout<<"DC Filter: "; break;
        }
        case 4:
        {
            Local_one = &CF::LS_one; Local_two = &CF::LS_two;
            cout<<"Bernstein Filter: "; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }

    Tstart = clock();
    for(int it=0;it<ItNum;++it)
    {
        //BC
        for (int i = 0; i < M_half-1; ++i)
        {
            p = BC.ptr<float>(i); p_right = WC.ptr<float>(i);
            p_down = WT.ptr<float>(i+1); p_rd = BT.ptr<float>(i+1); 
            p_pre = WT.ptr<float>(i); p_Corner = BT.ptr<float>(i);
            (this->*Local_two)(p, p_right, p_down, p_rd, p_pre, p_Corner, stepsize);
        }
        //BT
        for (int i = 1; i < M_half; ++i)
        {
            p = BT.ptr<float>(i); p_right = WT.ptr<float>(i);
            p_down = WC.ptr<float>(i); p_rd = BC.ptr<float>(i); 
            p_pre = WC.ptr<float>(i-1); p_Corner = BC.ptr<float>(i-1);
            (this->*Local_one)(p, p_right, p_down, p_rd, p_pre, p_Corner, stepsize);
        }
        //WC
        for (int i = 0; i < M_half-1; ++i)
        {
            p = WC.ptr<float>(i); p_right = BC.ptr<float>(i);
            p_down = BT.ptr<float>(i+1); p_rd = WT.ptr<float>(i+1); 
            p_pre = BT.ptr<float>(i); p_Corner = WT.ptr<float>(i);
            (this->*Local_one)(p, p_right, p_down, p_rd, p_pre, p_Corner, stepsize);
        }
        //WT
        for (int i = 1; i < M_half; ++i)
        {
            p = WT.ptr<float>(i); p_right = BT.ptr<float>(i);
            p_down = BC.ptr<float>(i); p_rd = WC.ptr<float>(i); 
            p_pre = BC.ptr<float>(i-1); p_Corner = WC.ptr<float>(i-1);
            (this->*Local_two)(p, p_right, p_down, p_rd, p_pre, p_Corner, stepsize);
        }
        
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);

    //merge four sets back to imgF
    merge();
}

//this nosplit is very useful for tasks like deconvolution, where the four sets need to be merged 
//every iteration if we use the split scheme.
void CF::FilterNoSplit(const int Type, double & time, const int ItNum, const float stepsize, const bool stochastic)
{
    clock_t Tstart, Tend;
    const float Max_rand_float = stepsize/RAND_MAX;
    float scaled_stepsize;

    float (CF::* Local)(int i, const float* p_pre, const float* p, const float* p_nex, const float * p_curv);

    switch(Type)
    {
        case 0:
        {
            Local = &CF::Scheme_TV; cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local = &CF::Scheme_MC; cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local = &CF::Scheme_GC; cout<<"GC Filter: "; break;
        }
        case 3:
        {
         Local = &CF::Scheme_DC; cout<<"DC Filter: "; break;
        }
        case 4:
        {
         Local = &CF::Scheme_LS; cout<<"Bernstein Filter: "; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }
    Tstart = clock();
    float d;
    for(int it=0;it<ItNum;++it)
    {
        //black circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,NULL);
                if (stochastic)
                {
                    scaled_stepsize = rand()*Max_rand_float;
                    p[j] += (scaled_stepsize*d);
                }
                else p[j] += (stepsize*d);
            }
        }

        //black triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,NULL);
                if (stochastic)
                {
                    scaled_stepsize = rand()*Max_rand_float;
                    p[j] += (scaled_stepsize*d);
                }
                else p[j] += (stepsize*d);
            }
        }

        //white circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,NULL);
                if (stochastic)
                {
                    scaled_stepsize = rand()*Max_rand_float;
                    p[j] += (scaled_stepsize*d);
                }
                else p[j] += (stepsize*d);
            }
        }

        //white triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,NULL);
                if (stochastic)
                {
                    scaled_stepsize = rand()*Max_rand_float;
                    p[j] += (scaled_stepsize*d);
                }
                else p[j] += (stepsize*d);
            }
        }
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
}

//compute the guide curvature from a given image
Mat CF::GuideCurvature(const char * FileName, const int Type)
{
    Mat tmp = imread(FileName, CV_8UC1);
    if (abs(tmp.rows - M)>2 || abs(tmp.cols - N)>2)
    {
        cout<<"warning: the guided image is resized."<<endl;
    }
    Mat tmp2 = Mat::zeros(M, N, CV_8UC1);
    resize(tmp, tmp2, tmp2.size());
    Mat tmp3 = Mat::zeros(M, N, CV_32FC1);
    tmp2.convertTo(tmp3, CV_32FC1);
    tmp3 /= 255.0f;
    Mat curv = Mat::zeros(M, N, CV_32FC1);
    switch(Type)
    {
        case 0:
        {
            TV(tmp3, curv); break;
        }
        case 1:
        {
            MC(tmp3, curv); break;
        }
        case 2:
        {
            GC(tmp3, curv); break;
        }
        default:
        {
            cout<<"The type in GuideCurvature is wrong. Do nothing."<<endl;
        }
    }
    return curv;
}

//curvature guided filter
void CF::CurvatureGuidedFilter(const Mat & curv, const int Type, double & time, const int ItNum, const float stepsize)
{
    clock_t Tstart, Tend;

    float (CF::* Local)(int i, const float* p_pre, const float* p, const float* p_nex, const float * p_curv);

    switch(Type)
    {
        case 0:
        {
            Local = &CF::Scheme_TV; cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local = &CF::Scheme_MC; cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local = &CF::Scheme_GC; cout<<"GC Filter: "; break;
        }
        case 3:
        {
         Local = &CF::Scheme_DC; cout<<"DC Filter: "; break;
        }
        case 4:
        {
         Local = &CF::Scheme_LS; cout<<"Bernstein Filter: "; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }
    Tstart = clock();
    float d;
    const float * p_curv;
    for(int it=0;it<ItNum;++it)
    {
        //black circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_curv = curv.ptr<float>(i);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,p_curv);
                p[j] += (stepsize*d);
            }
        }

        //black triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_curv = curv.ptr<float>(i);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,p_curv);
                p[j] += (stepsize*d);
            }
        }

        //white circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_curv = curv.ptr<float>(i);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,p_curv);
                p[j] += (stepsize*d);
            }
        }

        //white triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_curv = curv.ptr<float>(i);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = (this->*Local)(j,p_pre,p,p_down,p_curv);
                p[j] += (stepsize*d);
            }
        }
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
}

//generic filter solver for variational model |U - I|^DataFitOrder + lambda * Regularization
//the DataFitOrder can be fractional such as 1.5, which is not possible for other solvers.
void CF::Solver(const int Type, double & time, const int MaxItNum, const float lambda, const float DataFitOrder, const float stepsize)
{
    clock_t Tstart, Tend;
    float (CF::* Local)(int i, const float* p_pre, const float* p, const float* p_nex, const float * p_curv);
    void (CF::* curvature_compute)(const Mat& img, Mat& curv, int scheme);

    Mat curvature = Mat::zeros(M, N, CV_32FC1);
    Mat dataFit = Mat::zeros(M, N, CV_32FC1);
    std::vector<double> energyRecord_DataFit;
    std::vector<double> energyRecord_Curvature;

    std::cout<<"********************************************"<<endl;
    std::cout<<"*** Filter Solver for Variational Models ***\n    Lambda = "<<lambda<<" and DataFitOrder = "<<DataFitOrder<<endl;
    std::cout<<"********************************************"<<endl;

    switch(Type)
    {
        case 0:
        {
            Local = &CF::Scheme_TV; cout<<"TV Filter:(TVL1 by default) "; 
            curvature_compute = &CF::TV; break;
        }
        case 1:
        {
            Local = &CF::Scheme_MC; cout<<"MC Filter: "; 
            curvature_compute = &CF::MC; break;
        }
        case 2:
        {
            Local = &CF::Scheme_GC; cout<<"GC Filter: "; 
            curvature_compute = &CF::GC; break;
        }
        case 3:
        {
          Local = &CF::Scheme_DC; cout<<"DC Filter: "; 
            curvature_compute = &CF::GC; break;
        }
        case 4:
        {
          Local = &CF::Scheme_LS; cout<<"Bernstein Filter: "; 
            curvature_compute = &CF::MC; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }
    
    float d, energy_increase, dist_orig, dist_proj_orig;
    int count = 0;

    Tstart = clock();
    for(int it=0;it<MaxItNum;++it)
    {
        (this->*curvature_compute)(imgF, curvature,0);
        energyRecord_Curvature.push_back(lambda*energy(curvature));

        dataFit = imgF - image;
        energyRecord_DataFit.push_back(DataFitEnergy(dataFit,DataFitOrder));
        //if the energy starts to increase, stop the loop
        if(count>1 && (energyRecord_DataFit[it] + energyRecord_Curvature[it] > energyRecord_DataFit[it-1] + energyRecord_Curvature[it-1])) break;

        //black circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_data = image.ptr<float>(i);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = fabsf(p[j] - p_data[j]);
                dist_proj_orig = fabsf(p[j] + d - p_data[j]);
                energy_increase = powf(dist_proj_orig, DataFitOrder) - powf(dist_orig, DataFitOrder);
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //black triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_data = image.ptr<float>(i);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = fabsf(p[j] - p_data[j]);
                dist_proj_orig = fabsf(p[j] + d - p_data[j]);
                energy_increase = powf(dist_proj_orig, DataFitOrder) - powf(dist_orig, DataFitOrder);
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //white circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_data = image.ptr<float>(i);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = fabsf(p[j] - p_data[j]);
                dist_proj_orig = fabsf(p[j] + d - p_data[j]);
                energy_increase = powf(dist_proj_orig, DataFitOrder) - powf(dist_orig, DataFitOrder);
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //white triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            p_data = image.ptr<float>(i);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = fabsf(p[j] - p_data[j]);
                dist_proj_orig = fabsf(p[j] + d - p_data[j]);
                energy_increase = powf(dist_proj_orig, DataFitOrder) - powf(dist_orig, DataFitOrder);
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }
        count++;
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);

    cout<<"stop after "<<count<<" Iterations and ";

    (this->*curvature_compute)(imgF, curvature,0);
    dataFit = imgF - image;
    energyRecord_Curvature.push_back(lambda*energy(curvature));
    energyRecord_DataFit.push_back(DataFitEnergy(dataFit,DataFitOrder));

    //output the total energy profile
    ofstream energyProfile;
    energyProfile.open ("Energy.txt");
    energyProfile<<"### Iteration TotalEnergy DataFitEnergy RegularizationEnergy"<<endl;
    for (int i = 0; i <= count; ++i)
    {
    energyProfile<<i<<" "<<energyRecord_DataFit[i] + energyRecord_Curvature[i]<<" "<<energyRecord_DataFit[i]<<" "<<energyRecord_Curvature[i]<<endl;
    }
    energyProfile.close();
}


//solve BlackBox() + lambda * |curvature(U)|
 void CF::BlackBoxSolver(const int Type, double & time, const int MaxItNum, const float lambda, 
                        float (*BlackBox)(int row, int col, Mat& U, Mat & img_orig, float & d), const float stepsize)
 {
    clock_t Tstart, Tend;
    float (CF::* Local)(int i, const float* p_pre, const float* p, const float* p_nex, const float * p_curv);
    void (CF::* curvature_compute)(const Mat& img, Mat& curv, int scheme);

    Mat curvature = Mat::zeros(M, N, CV_32FC1);
    Mat dataFit = Mat::zeros(M, N, CV_32FC1);
    std::vector<double> energy_DataFit;
    std::vector<double> energy_Curvature;

    std::cout<<"********************************************"<<endl;
    std::cout<<"*** Filter Solver for Variational Models ***\n    Lambda = "<<lambda<<endl;
    std::cout<<"********************************************"<<endl;

    switch(Type)
    {
        case 0:
        {
            Local = &CF::Scheme_TV; cout<<"TV Filter:(TVL1 by default) "; 
            curvature_compute = &CF::TV; break;
        }
        case 1:
        {
            Local = &CF::Scheme_MC; cout<<"MC Filter: "; 
            curvature_compute = &CF::MC; break;
        }
        case 2:
        {
            Local = &CF::Scheme_GC; cout<<"GC Filter: "; 
            curvature_compute = &CF::GC; break;
        }
        case 3:
        {
          Local = &CF::Scheme_DC; cout<<"DC Filter: "; 
            curvature_compute = &CF::GC; break;
        }
        case 4:
        {
          Local = &CF::Scheme_LS; cout<<"Bernstein Filter: "; 
            curvature_compute = &CF::MC; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }
    
    float d, energy_increase, dist_orig, dist_proj_orig;
    int count = 0; float zero = 0.0f;

    Tstart = clock();
    for(int it=0;it<MaxItNum;++it)
    {
        (this->*curvature_compute)(imgF, curvature,0);
        energy_Curvature.push_back(lambda*energy(curvature));

        for (int i = 1; i < M-1; ++i)
            for (int j = 1; j < N-1; ++j)
            {
                dataFit.at<float>(i,j) = BlackBox(i,j, imgF, image,zero);
            }
        Scalar tmp_scalar = sum(dataFit);
        energy_DataFit.push_back(tmp_scalar(0));
        //if the energy starts to increase, stop the loop
        if(count>1 && (energy_DataFit[it] + energy_Curvature[it] > energy_DataFit[it-1] + energy_Curvature[it-1])) break;

        //black circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = BlackBox(i,j,imgF,image,zero);
                dist_proj_orig = BlackBox(i,j,imgF,image,d);
                energy_increase = dist_proj_orig - dist_orig;
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //black triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = BlackBox(i,j,imgF,image,zero);
                dist_proj_orig = BlackBox(i,j,imgF,image,d);
                energy_increase = dist_proj_orig - dist_orig;
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //white circle
        for (int i = 1; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 2; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = BlackBox(i,j,imgF,image,zero);
                dist_proj_orig = BlackBox(i,j,imgF,image,d);
                energy_increase = dist_proj_orig - dist_orig;
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }

        //white triangle
        for (int i = 2; i < M-1; ++i,++i)
        {
            p = imgF.ptr<float>(i);
            p_pre = imgF.ptr<float>(i-1);
            p_down = imgF.ptr<float>(i+1);
            for (int j = 1; j < N-1; ++j, ++j)
            {
                d = stepsize*(this->*Local)(j,p_pre,p,p_down,NULL);
                dist_orig = BlackBox(i,j,imgF,image,zero);
                dist_proj_orig = BlackBox(i,j,imgF,image,d);
                energy_increase = dist_proj_orig - dist_orig;
                if (energy_increase <= lambda*abs(d)) p[j] += d;
            }
        }
        count++;
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);

    cout<<"stop after "<<count<<" Iterations and ";

    (this->*curvature_compute)(imgF, curvature,0);
    energy_Curvature.push_back(lambda*energy(curvature));
    for (int i = 1; i < M-1; ++i)
        for (int j = 1; j < N-1; ++j)
        {
            dataFit.at<float>(i,j) = BlackBox(i,j, imgF, image,zero);
        }
    Scalar tmp_scalar = sum(dataFit);
    energy_DataFit.push_back(tmp_scalar(0));
    

    //output the total energy profile
    ofstream energyProfile;
    energyProfile.open ("Energy.txt");
    energyProfile<<"### Iteration TotalEnergy DataFitEnergy RegularizationEnergy"<<endl;
    for (int i = 0; i <= count; ++i)
    {
    energyProfile<<i<<" "<<energy_DataFit[i] + energy_Curvature[i]<<" "<<energy_DataFit[i]<<" "<<energy_Curvature[i]<<endl;
    }
    energyProfile.close();
 }

void CF::DM(const int FilterType, const int LocationType, const Mat & img, Mat & dm)
{
    const float *p, *p_pre, *p_down;
    float *p_d;
    float (CF::* Local)(int i, const float* p_pre, const float* p, const float* p_nex, const float * p_curv);
    switch(Type)
    {
        case 0:
        {
            Local = &CF::Scheme_TV; cout<<"TV dm: "; break;
        }
        case 1:
        {
            Local = &CF::Scheme_MC; cout<<"MC dm: "; break;
        }
        case 2:
        {
            Local = &CF::Scheme_GC; cout<<"GC dm: "; break;
        }
        case 4:
        {
            Local = &CF::Scheme_LS; cout<<"Bernstein dm: "; break;
        }
        default:
        {
            cout<<"The filter type is wrong. Do nothing."<<endl;
            return;
        }
    }
    int start_row, start_col;
    switch(LocationType)
    {
        case 0://black circle
        {
            start_row = 1; start_col = 1; break;
        }
        case 1://black triangle
        {
            start_row = 2; start_col = 2; break;
        }
        case 2://white circle
        {
            start_row = 1; start_col = 2; break;
        }
        case 3://white triangle
        {
            start_row = 2; start_col = 1; break;
        }
        default:
        {
            cout<<"The location type is wrong. Do nothing."<<endl;
            return;
        }
    }
    int M = img.rows;
    int N = img.cols; 
    
    for (int i = start_row; i < M-1; ++i,++i)
    {
        p = img.ptr<float>(i);
        p_pre = img.ptr<float>(i-1);
        p_down = img.ptr<float>(i+1);
        p_d = dm.ptr<float>(i);
        for (int j = start_col; j < N-1; ++j, ++j)
        {
            p_d[j] = (this->*Local)(j,p_pre,p,p_down,NULL);
        }
    }

    //statistically, following gives the gradient of curvature energy
    // if(Type==1) dm *= 2;
    // if(Type==2) dm *= 30; //GC filter is super efficient
 }

 //find the value with minimum abs value, 4 floats
inline float CF::SignedMin(float * dist)
{
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    if (fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
    if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
    return dist[0];
}

//find the value with minimum abs value, 4 floats
inline float CF::SignedMin_noSplit(float * dist)
{
    float absMin = fabsf(dist[0]);
    unsigned char index = 0;
    for (unsigned char i = 1; i < 4; ++i)
    {
        if (fabsf(dist[i])<absMin)
        {
            absMin = fabsf(dist[i]);
            index = i;
        }
    }
    return dist[index];
}

//*************************** Do NOT change anything! *****************************//
//************************* these filters are optimized ***************************//
//********** contact Yuanhao Gong if you need to change anything ******************//
//*************************** gongyuanhao@gmail.com *******************************//
inline void CF::GC_one(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float scaled_stepsize = stepsize/3;
    register float tmp, min_value, com_one, com_two, min_value2;
    for (int j = 1; j < N_half; ++j)
     {
        tmp = 2*p[j];
        dist[0] = p_pre[j]+p_down[j] - tmp;
        dist[1] = p_right[j-1]+p_right[j] - tmp;
        dist[2] = p_Corner[j-1]+p_rd[j] - tmp;
        dist[3] = p_Corner[j]+p_rd[j-1] - tmp;

        min_value = SignedMin(dist);

        tmp *= 1.5f;
        min_value *= 1.5f;
        com_one = p_pre[j] - tmp;
        com_two = p_down[j] - tmp;
        dist[0] = p_Corner[j-1] + p_right[j-1] + com_one;
        dist[1] = p_Corner[j] + p_right[j] + com_one;
        dist[2] = p_right[j-1] + p_rd[j-1] + com_two;
        dist[3] = p_right[j] + p_rd[j] +com_two;

        min_value2 = SignedMin(dist);

        if (fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void CF::GC_two(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float scaled_stepsize = stepsize/3;
    register float tmp, min_value, min_value2, com_one, com_two;
    for (int j = 0; j < N_half-1; ++j)
     {
        tmp = 2*p[j];
        dist[0] = p_pre[j]+p_down[j] - tmp;
        dist[1] = p_right[j]+p_right[j+1] - tmp;
        dist[2] = p_Corner[j]+p_rd[j+1] - tmp;
        dist[3] = p_Corner[j+1]+p_rd[j] - tmp;
        
        min_value = SignedMin(dist);

        tmp *= 1.5f;
        min_value *= 1.5f;
        com_one = p_pre[j] - tmp;
        com_two = p_down[j] - tmp;
        dist[0] = p_Corner[j] + p_right[j] + com_one;
        dist[1] = p_Corner[j+1] + p_right[j+1] + com_one;
        dist[2] = p_right[j] + p_rd[j] + com_two;
        dist[3] = p_right[j+1] + p_rd[j+1] + com_two;
        
        min_value2 = SignedMin(dist);
        if (fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void CF::MC_one(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8];
    register float scaled_stepsize = stepsize/8;
    register float min_value;
    for (int j = 1; j < N_half; ++j)
     {
         
        dist[4] = p[j]*8;
        dist[5] = (p_pre[j]+p_down[j])*2.5f - dist[4];
        dist[6] = (p_right[j-1]+p_right[j])*2.5f - dist[4];

        dist[0] = dist[5] + p_right[j]*5 -p_Corner[j]-p_rd[j];
        dist[1] = dist[5] + p_right[j-1]*5 -p_Corner[j-1]-p_rd[j-1];
        dist[2] = dist[6] + p_pre[j]*5 -p_Corner[j-1]-p_Corner[j];
        dist[3] = dist[6] + p_down[j]*5 -p_rd[j-1]-p_rd[j];

        min_value = SignedMin(dist);

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void CF::MC_two(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8];
    register float scaled_stepsize = stepsize/8;
    register float min_value;
    for (int j = 0; j < N_half-1; ++j)
    {
         
        dist[4] = p[j]*8;
        dist[5] = (p_pre[j]+p_down[j])*2.5f - dist[4];
        dist[6] = (p_right[j]+p_right[j+1])*2.5f - dist[4];

        dist[0] = dist[5] + p_right[j]*5 -p_Corner[j]-p_rd[j];
        dist[1] = dist[5] + p_right[j+1]*5 -p_Corner[j+1]-p_rd[j+1];
        dist[2] = dist[6] + p_pre[j]*5 -p_Corner[j]-p_Corner[j+1];
        dist[3] = dist[6] + p_down[j]*5 -p_rd[j]-p_rd[j+1];
        
        min_value = SignedMin(dist);
        
        p[j] += (scaled_stepsize*min_value);
    }
}

inline void CF::LS_one(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    //one is for BT and WC, two is for BC and WT
    register float dist[4];
    register float scaled_stepsize = stepsize/10;
    register float tmp;
    for (int j = 1; j < N_half; ++j)
    {
         
        tmp = p[j]*2;
        dist[0] = (p_pre[j]+p_down[j]) - tmp;
        dist[1] = (p_right[j-1]+p_right[j]) - tmp;
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        tmp *= 3.5f;
        dist[0] *= 3.3333f;

        dist[2] = 3*(p_Corner[j] + p_rd[j-1]) - tmp;
        dist[3] = 3*(p_Corner[j-1] + p_rd[j]) - tmp;

        dist[1] = p_right[j-1] + p_pre[j] - p_Corner[j-1] + dist[2];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[1] = p_right[j] + p_pre[j] - p_Corner[j] + dist[3];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[1] = p_right[j] + p_down[j] - p_rd[j] + dist[2];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[1] = p_right[j-1] + p_down[j]  - p_rd[j-1] + dist[3];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        p[j] += (scaled_stepsize*dist[0]);
    }
}

inline void CF::LS_two(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    //one is for BT and WC, two is for BC and WT
    register float dist[4];
    register float scaled_stepsize = stepsize/10;
    register float tmp;
    for (int j = 0; j < N_half-1; ++j)
    {
        tmp = p[j]*2;
        dist[0] = p_pre[j]+p_down[j] - tmp;
        dist[1] = p_right[j]+p_right[j+1] - tmp;
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        tmp *= 3.5f;
        dist[0] *= 3.3333f;

        dist[2] = 3*(p_Corner[j+1] + p_rd[j]) - tmp;
        dist[3] = 3*(p_Corner[j] + p_rd[j+1]) - tmp;

        dist[0] = p_right[j] + p_pre[j] - p_Corner[j] + dist[2];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[0] = p_right[j+1] + p_pre[j] - p_Corner[j+1] + dist[3];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[0] = p_right[j+1] + p_down[j] - p_rd[j+1] + dist[2];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        dist[0] = p_right[j] + p_down[j] - p_rd[j] + dist[3];
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        p[j] += (scaled_stepsize*dist[0]);
    }
}

inline void CF::TV_one(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8], dist2[4];
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five)
    register float scaled_stepsize = stepsize/6;
    register float scaledP, min_value, min_value2;
    for (int j = 1; j < N_half; ++j)
     {
        //temp var
        scaledP =  p[j]*5;
        dist[4] = p_pre[j]+p_down[j] - scaledP;
        dist[5] = p_right[j-1]+p_right[j] - scaledP;
        dist[6] = p_Corner[j-1]+p_Corner[j]+p_pre[j] - scaledP;
        dist[7] = p_rd[j-1]+p_rd[j]+p_down[j] - scaledP;

        dist[0] = dist[4] + p_Corner[j-1]+p_rd[j-1]+p_right[j-1];
        dist[1] = dist[4] + p_Corner[j]+p_rd[j]+p_right[j];
        dist[2] = dist[5] + p_Corner[j-1] + p_Corner[j] + p_pre[j];
        dist[3] = dist[5] + p_rd[j-1] + p_rd[j] + p_down[j];

        dist2[0] = dist[6] + p_right[j-1]+p_rd[j-1];
        dist2[1] = dist[6] + p_right[j]+p_rd[j];
        dist2[2] = dist[7] + p_right[j-1]+p_Corner[j-1];
        dist2[3] = dist[7] + p_right[j]+p_Corner[j];

        min_value = SignedMin(dist);
        min_value2 = SignedMin(dist2);
        if(fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void CF::TV_two(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8], tmp[4];
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five) 
    register float scaled_stepsize = stepsize/6;
    register float scaledP, min_value, min_value2;
    for (int j = 0; j < N_half-1; ++j)
     {
        scaledP = p[j]*5;
        tmp[0] = p_pre[j]+p_down[j] - scaledP;
        tmp[1] = p_right[j]+p_right[j+1] - scaledP;
        tmp[2] = p_Corner[j]+p_Corner[j+1]+p_pre[j] - scaledP;
        tmp[3] = p_rd[j]+p_rd[j+1]+p_down[j] - scaledP;

        dist[0] = tmp[0] + p_Corner[j]+p_rd[j]+p_right[j];
        dist[1] = tmp[0] + p_Corner[j+1]+p_rd[j+1]+p_right[j+1];
        dist[2] = tmp[1] + p_Corner[j] + p_Corner[j+1] + p_pre[j];
        dist[3] = tmp[1] + p_rd[j] + p_rd[j+1] + p_down[j];
        
        dist[4] = tmp[2] +p_right[j]+p_rd[j];
        dist[5] = tmp[2] +p_right[j+1]+p_rd[j+1];
        dist[6] = tmp[3] +p_right[j]+p_Corner[j];
        dist[7] = tmp[3] +p_right[j+1]+p_Corner[j+1];
        
        min_value = SignedMin(dist);
        min_value2 = SignedMin(dist+4);
        if(fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;
        p[j] += (scaled_stepsize*min_value);
     }
}


inline void CF::DC_one(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float weight = -0.225603f;
    
    for (int j = 1; j < N_half; ++j)
     {
        dist[1] = (p_pre[j]+p_down[j])/2 - p[j];
        dist[2] = (p_right[j-1]+p_right[j])/2 - p[j];
        
        dist[0] = dist[1] + (p_rd[j] + p_Corner[j] - 2*p_right[j])*weight;
        dist[3] = dist[1] + (p_rd[j-1] + p_Corner[j-1] - 2*p_right[j-1])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        dist[3] = dist[2] + (p_Corner[j-1] + p_Corner[j] - 2*p_pre[j])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        dist[3] = dist[2] + (p_rd[j-1] + p_rd[j] - 2*p_down[j])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        p[j] += (stepsize*dist[0]);
     }
}

inline void CF::DC_two(float* __restrict p, const float* __restrict p_right, const float* __restrict p_down, 
    const float * __restrict p_rd, const float* __restrict p_pre, const float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float weight = -0.225603f;

    for (int j = 0; j < N_half-1; ++j)
     {
        dist[1] = (p_pre[j]+p_down[j])/2 - p[j];
        dist[2] = (p_right[j]+p_right[j+1])/2 - p[j];
        
        dist[0] = dist[1] + (p_Corner[j+1] + p_rd[j+1] - 2*p_right[j+1])*weight;
        dist[3] = dist[1] + (p_Corner[j] + p_rd[j] - 2*p_right[j])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        dist[3] = dist[2] + (p_Corner[j] + p_Corner[j+1] - 2*p_pre[j])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        dist[3] = dist[2] + (p_rd[j] + p_rd[j+1] - 2*p_down[j])*weight;
        if(fabsf(dist[3]) < fabsf(dist[0])) dist[0] = dist[3];

        p[j] += (stepsize*dist[0]);
     }
}

/********************************************************************/
/********************** scheme at each pixel ************************/
/********************** only for noSplit case ***********************/
/********************************************************************/
inline float CF::Scheme_GC(int i, const float * __restrict p_pre, const float * __restrict p, 
                            const float * __restrict p_nex, const float * p_curv)
{
    register float dist[4];
    register float tmp, min_value, min_value2;
    tmp = -2*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) tmp = -2*(p[i] + p_curv[i]);
    dist[0] = p_pre[i] + p_nex[i] + tmp;
    dist[1] = p[i-1] + p[i+1] + tmp;
    dist[2] = p_pre[i-1] + p_nex[i+1] + tmp;
    dist[3]  = p_nex[i-1] + p_pre[i+1] + tmp;
    min_value = SignedMin_noSplit(dist);
    min_value *= 1.5f;
    
    tmp *= 1.5f;
    dist[0] = p_pre[i] + p_pre[i-1] + p[i-1] + tmp;
    dist[1] = p_pre[i] + p_pre[i+1] + p[i+1] + tmp;
    dist[2] = p_nex[i] + p_nex[i-1] + p[i-1] + tmp;
    dist[3] = p_nex[i] + p_nex[i+1] + p[i+1] + tmp;
    min_value2 = SignedMin_noSplit(dist);

    if(fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;
    min_value /= 3;
    
    return min_value;
}

inline float CF::Scheme_MC(int i, const float * __restrict p_pre, const float * __restrict p, 
                            const float * __restrict p_nex, const float * p_curv)
{
    //compute the movement according to half window
    //       a   b
    //       I   e
    //       c   d
    // return (2.5(a+c)+5*e)-b-d)/8.0;

    register float dist[4];
    register float tmp, com_one, com_two, min_value;
    tmp = 8*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) tmp = 8*(p[i] + p_curv[i]);

    com_one = (p_pre[i]+p_nex[i])*2.5f - tmp;
    com_two = (p[i-1]+p[i+1])*2.5f - tmp;

    dist[0] = com_one + 5.0f*p[i+1] - p_pre[i+1] - p_nex[i+1];
    dist[1] = com_one + 5.0f*p[i-1] - p_pre[i-1] - p_nex[i-1];
    dist[2] = com_two + 5.0f*p_nex[i] - p_nex[i-1] - p_nex[i+1];
    dist[3] = com_two + 5.0f*p_pre[i] - p_pre[i-1] - p_pre[i+1];

    min_value = SignedMin_noSplit(dist)/8;
    
    return min_value;
}

inline float CF::Scheme_LS(int i, const float * __restrict p_pre, const float * __restrict p, 
                            const float * __restrict p_nex, const float * p_curv)
{
    //compute the movement according to half window
    //   f   a   b            0 1/2 0               3/7 1/7 -1/7
    //       I   e               -1 0                    -1  1/7
    //       c   d              1/2 0                        3/7
    // or (include central pixel)
    //   f   a   b            0 1/3 0               3/10 1/10 -1/10
    //       I   e              -2/3 0                   -7/10  1/10
    //       c   d              1/3 0                        3/10

    register float dist[4];
    register float tmp, min_value, tmp_one, tmp_two;
    tmp = 2*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) tmp = 2*(p[i] + p_curv[i]);

    min_value = p_pre[i]+p_nex[i] - tmp;
    dist[0] = p[i-1] + p[i+1] - tmp;
    if(fabsf(dist[0])<fabsf(min_value)) min_value = dist[0];
    
    tmp *= 3.5f;
    min_value *= 3.3333f;

    tmp_one = 3*(p_nex[i-1] + p_pre[i+1]) - tmp;
    tmp_two = 3*(p_pre[i-1] + p_nex[i+1]) - tmp;

    dist[0] = p_pre[i] - p_pre[i-1] + p[i-1] + tmp_one;
    dist[1] = p_pre[i] - p_pre[i+1] + p[i+1] + tmp_two;
    dist[2] = p_nex[i] - p_nex[i+1] + p[i+1] + tmp_one;
    dist[3] = p_nex[i] - p_nex[i-1] + p[i-1] + tmp_two;
    
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    if(fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
    if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
    if(fabsf(dist[0])<fabsf(min_value)) min_value = dist[0];

    return min_value/10;//here 10 means including central pixel while 7 means exclusion
}

inline float CF::Scheme_TV(int i, const float * __restrict p_pre, const float * __restrict p, 
                            const float * __restrict p_nex, const float * p_curv)
{
    //       a   b
    //       I   e
    //       c   d
    // return (a+b+c+d+e)/5.0;
    register float dist[4], tmp[4];
     //old fashion, need 5*8 times plus or minus
    register float scaledP, min_value, min_value2;
    scaledP = 5*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) scaledP = 5*(p[i] + p_curv[i]);


    tmp[0] = p_pre[i-1]+p[i-1] + p_nex[i-1] - scaledP;
    tmp[1] = p_pre[i+1]+p[i+1] + p_nex[i+1] - scaledP;
    tmp[2] = p[i-1]+p[i+1] - scaledP;
    tmp[3] = p_pre[i]+p_nex[i];

    dist[0] = tmp[0] + tmp[3];
    dist[1] = tmp[1] + tmp[3];
    dist[2] = tmp[2] + p_pre[i-1] + p_pre[i] + p_pre[i+1];
    dist[3] = tmp[2] + p_nex[i-1] + p_nex[i] + p_nex[i+1];
    min_value = SignedMin_noSplit(dist);

    //diag
    dist[0] = tmp[0] + p_pre[i] + p_pre[i+1];
    dist[1] = tmp[0] + p_nex[i] + p_nex[i+1];
    dist[2] = tmp[1] + p_pre[i-1] + p_pre[i];
    dist[3] = tmp[1] + p_nex[i-1] + p_nex[i];
    min_value2 = SignedMin_noSplit(dist);

    if(fabsf(min_value2)<fabsf(min_value)) min_value = min_value2;
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five)

    return min_value/6;
}

inline float CF::Scheme_DC(int i, const float * __restrict p_pre, const float * __restrict p, 
                            const float * __restrict p_nex, const float * p_curv)
{
    float dist[2];
    float weight = -0.225603f;
    float scaledP = p[i];
    if (p_curv != NULL) scaledP = p[i] + p_curv[i];

    dist[0] = (p_pre[i] + p_nex[i])/2 + (p_pre[i+1] + p_nex[i+1] - 2*p[i+1])*weight - scaledP;
    dist[1] = (p_pre[i] + p_nex[i])/2 + (p_pre[i-1] + p_nex[i-1] - 2*p[i-1])*weight - scaledP;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = (p[i-1] + p[i+1])/2 + (p_pre[i-1] + p_pre[i+1] - 2*p_pre[i])*weight - scaledP;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = (p[i-1] + p[i+1])/2 + (p_nex[i-1] + p_nex[i+1] - 2*p_nex[i])*weight - scaledP;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    return dist[0];
}
