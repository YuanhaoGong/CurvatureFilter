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

//Dual Mesh sctructure
class DM
{
public:
    //read one image from disk
    void read(const char* FileName);
    //set one image 
    void set(Mat& file);
    //get the filtered image 
    Mat get(){return imgF;};
    //write the result to disk
    void write();
    void write(const char* FileName);
    //compute TV
    void TV(const Mat & img, Mat & T);
    //compute MC
    void MC(const Mat & img, Mat & MC);
    void MC_new(const Mat & img, Mat & MC);//new scheme, Eq.6.12 in my thesis
    void MC_fit(const Mat & img, Mat & MC);
    void MC_isoLine(const Mat & img, Mat & MC);
    //compute GC
    void GC(const Mat & img, Mat & GC);
    void GC_new(const Mat & img, Mat & GC);//new scheme, Eq.6.16 in my thesis
    void GC_fit(const Mat & img, Mat & GC);
    void GC_LUT_Init();
    void GC_LUT(const Mat & img, Mat & GC);
    //compute energy for given TV, MC, or GC image
    double energy(const Mat& img);
    //compute data fitting energy between image and imgF
    double DataFitEnergy(Mat& tmp, double order);
    //PSNR
    double PSNR();
    double PSNR(const Mat& I1, const Mat& I2);
    //compute naturalness factor
    double Naturalness(){return Naturalness(imgF);}
    double Naturalness(const Mat & img);
    /******************* filters *****************************/
    // Type=0, TV; Type=1, MC; Type=2, GC; (Type=3, DC, experimental);
    // the stepsize parameter is in (0,1]:smaller means more iterations, but reaches lower energy level; larger means less iterations, but converges at higher energy level
    //////////////////////////////////////////////////////////
    void Filter(const int Type, double & time, const int ItNum = 10, const float stepsize=1);//with split
    void FilterNoSplit(const int Type, double & time, const int ItNum = 10, const float stepsize=1);//direct on imgF 
    /******************* Curvature Guided Filter *****************************/
    //compute the curvature from the guided image (scaled to the size of imgF)
    Mat GuideCurvature(const char * FileName, const int Type);
    //filter the image such that the result is close to the specified curvature
    void CurvatureGuidedFilter(const Mat & curv, const int Type, double & time, const int ItNum = 10, const float stepsize=1);
    /******************* generic solver for variational models *****************************/
    //solve |U - I|^DataFitOrder + lambda * |curvature(U)|
    void Solver(const int Type, double & time, const int MaxItNum, const float lambda = 2, const float DataFitOrder = 1, const float stepsize=1);
    //solve BlackBox(U,I) + lambda * |curvature(U)|
    void BlackBoxSolver(const int Type, double & time, const int MaxItNum, const float lambda, float (*BlackBox)(int row, int col, Mat& U, Mat & img_orig, float & d), const float stepsize=1);

private:
    //padded original, tmp, result
    Mat image, imgF, result;
    //four sets, be aware that position is fixed, see split() or Yuanhao Gong's PhD thesis
    Mat WC, WT, BC, BT;
    //image size
    int M, N, M_orig, N_orig, M_half, N_half;
    //six pointers
    float* p, *p_right, *p_down, *p_rd, *p_pre, *p_Corner;
    //pointer to the data
    const float* p_data;
    //Look Up Table for fast computing GC
    Mat LUT;
private:
    //split imgF into four sets
    void split();
    //merge four sets back to imgF
    void merge();
    //naturalness evaluation
    double Naturalness_search(float* data, int N, int offset);
    //fit coefficients 
    void FiveCoefficient(const Mat & img, Mat & x2, Mat &y2, Mat & xy, Mat & x, Mat & y);
    /*************************************** Split into 4 sets *********************************/
    //one is for BT and WC, two is for BC and WT
    inline void GC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void GC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);

    inline void MC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void MC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void TV_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void TV_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);

    inline void DC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void DC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void LS_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    inline void LS_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);

    /*************************************** Direct on imgF (no split) ***********************/
    inline float Scheme_GC(int i, float * p_pre, float * p, float * p_nex, const float * p_curv = NULL);
    inline float Scheme_MC(int i, float * p_pre, float * p, float * p_nex, const float * p_curv = NULL);
    inline float Scheme_TV(int i, float * p_pre, float * p, float * p_nex, const float * p_curv = NULL);
    inline float Scheme_DC(int i, float * p_pre, float * p, float * p_nex, const float * p_curv = NULL);
    inline float Scheme_LS(int i, float * p_pre, float * p, float * p_nex, const float * p_curv = NULL);
};

double DM::PSNR()
{
    return PSNR(image, imgF);
}

double DM::PSNR(const Mat& I1, const Mat& I2)
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
         double  mse =sse /(double)(I1.total());
         double psnr = - 10*log10(mse);
         return psnr;
     }
}

double DM::Naturalness(const Mat& imgF)
{
    //compute naturalness factor
    const float * p_row ;
    const float * pp_row;
    short int indexX, indexY;
    short int Offset = 256;
    short int N = 512;
    double eps = 0.0001;
    double left(0), right(1), mid_left(0), mid_right(0);

    Mat GradCDF = Mat::zeros(2, N, CV_32FC1);
    float * Gradx = GradCDF.ptr<float>(0);
    float * Grady = GradCDF.ptr<float>(1);
    double f = 1.0/((imgF.rows-1)*(imgF.cols-1));

    //not efficient but safe way
    for(int i = 0; i < imgF.rows - 1; i++)
    {
        p_row = imgF.ptr<float>(i);
        pp_row = imgF.ptr<float>(i+1);

        for(int j = 0; j < imgF.cols - 1; j++)
        {
            //scale back to 255
            indexX = Offset + (p_row[j+1] - p_row[j])*255;
            indexY = Offset + (pp_row[j] - p_row[j])*255;            
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

double DM::Naturalness_search(float* data, int N, int offset)
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

void DM::read(const char* FileName)
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

void DM::set(Mat& file)
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
void DM::write()
{
    DM::write("CF_result.png");
}

void DM::write(const char* FileName)
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
void DM::TV(const Mat & imgF, Mat & T)
{
    const float * p_row, * pn_row;
    float * p_t;
    if (true) //the switch between TVL1 and TVL2
    {
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
    }else //TVL2
    {
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
    }
}

//compute Mean Curvature
void DM::MC(const Mat& imgF, Mat & MC)
{
    //classical scheme is used
    const float * p_row, *pn_row, *pp_row;
    float *p_d;
    float Ix, Iy, Ixy, Ixx, Iyy, num, den, tmp;
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
}

//compute Mean Curvature, Eq.6.12 in my thesis
void DM::MC_new(const Mat& imgF, Mat & MC)
{
    //separate kernel from Eq.6.12, add the center pixel later
    Mat kernel = (Mat_<float>(1,3) << 0.25f, -1.25f, 0.25f); 
    sepFilter2D(imgF, MC, CV_32F, kernel, kernel,Point(-1,-1),0,BORDER_REPLICATE);
    MC *= -1;
    MC += (0.5625f*imgF);//the center pixel
}

//fit coefficients 
void DM::FiveCoefficient(const Mat & img, Mat & x2, Mat &y2, Mat & xy, Mat & x, Mat & y)
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
void DM::MC_fit(const Mat & img, Mat & MC)
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
            den = sqrt(den)*den;
            p_d[j] = num/den; 
        }
    }
}

void DM::MC_isoLine(const Mat & img, Mat & MC)
{
    //compute the MC by iso lines scheme
    const float * p_row, *pn_row, *pp_row;
    float *p_d;
    float Ix, Iy, Ixy, Ixx, Iyy, num, den, tmp;
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
            den = Ix*Ix + Iy*Iy + 0.000001f;
            p_d[j] = num/den;
        }   
    }
}

//compute Gaussian Curvature by fitting quad function
void DM::GC_fit(const Mat & img, Mat & GC)
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
void DM::GC(const Mat & imgF, Mat &GC)
{
    //classical scheme is used
    const float * p_row, *pn_row, *pp_row;
    float *p_d;
    float Ix, Iy, Ixx, Iyy, Ixy, num, den;
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
}

//new scheme, Eq.6.16 in my thesis
void DM::GC_new(const Mat & img, Mat & GC)
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
    /*
    //
    Mat kernel_one_sep = (Mat_<float>(1,3)<<0.0458709f, 5.54135f,0.0458709f); //center pixel -31.7317
    Mat kernel_two_sep_h = (Mat_<float>(1,3)<<1.0f, 0.0f, -1.0f);
    Mat kernel_two_sep_v = (Mat_<float>(1,3)<<0.25f, -0.5f, 0.25f);
    Mat kernel_three_sep_h = kernel_two_sep_v;
    Mat kernel_three_sep_v = (Mat_<float>(1,3)<<-1.0f, 0.0f, 1.0f);
    Mat kernel_four_sep_h = (Mat_<float>(1,3)<<0.535181f, -0.429729f, 0.535181f);//0.410410
    Mat kernel_four_sep_v = - kernel_four_sep_h;
    Mat kernel_five_sep_h = (Mat_<float>(1,3)<<-0.553341f,0.0f,0.553341f);
    Mat kernel_five_sep_v = - kernel_five_sep_h;
    
    sepFilter2D(img, GC, img.depth(), kernel_one_sep, kernel_one_sep);
    GC -= (31.7317f*img);
    pow(GC, 2, GC);
    sepFilter2D(img, tmp, img.depth(), kernel_two_sep_h, kernel_two_sep_v);
    pow(tmp, 2, tmp);
    GC -= tmp;
    sepFilter2D(img, tmp, img.depth(), kernel_three_sep_h, kernel_three_sep_v);
    pow(tmp, 2, tmp);
    GC -= tmp;
    sepFilter2D(img, tmp, img.depth(), kernel_four_sep_h, kernel_four_sep_v);
    tmp += 0.41041f*img;
    pow(tmp, 2, tmp);
    GC -= tmp;
    sepFilter2D(img, tmp, img.depth(), kernel_five_sep_h, kernel_five_sep_v);
    pow(tmp, 2, tmp);
    GC -= tmp;
    */
    
}

void DM::GC_LUT_Init()
{
    const int scale = 255*255;
    const int scale2 = 2*scale;
    int num;
    float den; 

    LUT = Mat::zeros(512, 512, CV_32FC1);
    float * p_LUT;
    for (int i = -255; i < 1; ++i)
    {
        p_LUT = LUT.ptr<float>(i+255);
        for (int j = -255; j < 255; ++j)
        {
            num = i*j + scale;
            den = sqrt(i*i+scale)*sqrt(j*j+scale2);//avoid integer overflow
            p_LUT[j+255] = num/den;
            if (abs(p_LUT[j+255])>1)
            {
                cout<<i<<" "<<j<<" "<<p_LUT[j+255]<<endl;
            }
            p_LUT[j+255] = acos(p_LUT[j+255]);
            LUT.at<float>(255-i,j+255) = p_LUT[j+255];
        }
    }
}

//approimate GC by LUT
void DM::GC_LUT(const Mat & img, Mat & GC)
{
    const float * p_row, *pp_row, *pn_row;
    float *g, *p_LUT;
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
double DM::energy(const Mat &img)
{
    Scalar tmp = sum(cv::abs(img));
    return tmp(0);
}

//compute the energy between image and imgF
double DM::DataFitEnergy(Mat & tmp, double order)
{
    tmp = abs(image - imgF);
    pow(tmp, order, tmp);
    Scalar tmp2 = sum(tmp);
    return tmp2(0);
}

//split the image into four sets
void DM::split()
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
void DM::merge()
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

void DM::Filter(const int Type, double & time, const int ItNum, const float stepsize)
{
    clock_t Tstart, Tend;
    //split imgF into four sets
    split();

    void (DM::* Local_one)(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);
    void (DM::* Local_two)(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner, const float& stepsize);

    switch(Type)
    {
        case 0:
        {
            Local_one = &DM::TV_one; Local_two = &DM::TV_two; 
            cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local_one = &DM::MC_one; Local_two = &DM::MC_two; 
            cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local_one = &DM::GC_one; Local_two = &DM::GC_two; 
            cout<<"GC Filter: "; break;
        }
        case 3:
        {
            Local_one = &DM::DC_one; Local_two = &DM::DC_two; 
            cout<<"DC Filter: "; break;
        }
        case 4:
        {
            Local_one = &DM::LS_one; Local_two = &DM::LS_two;
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
void DM::FilterNoSplit(const int Type, double & time, const int ItNum, const float stepsize)
{
    clock_t Tstart, Tend;

    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex, const float * p_curv);

    switch(Type)
    {
        case 0:
        {
            Local = &DM::Scheme_TV; cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local = &DM::Scheme_MC; cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local = &DM::Scheme_GC; cout<<"GC Filter: "; break;
        }
        case 3:
        {
         Local = &DM::Scheme_DC; cout<<"DC Filter: "; break;
        }
        case 4:
        {
         Local = &DM::Scheme_LS; cout<<"Bernstein Filter: "; break;
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
                p[j] += (stepsize*d);
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
                p[j] += (stepsize*d);
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
                p[j] += (stepsize*d);
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
                p[j] += (stepsize*d);
            }
        }
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
}

//compute the guide curvature from a given image
Mat DM::GuideCurvature(const char * FileName, const int Type)
{
    Mat tmp = imread(FileName, CV_8UC1);
    if (abs(tmp.rows - M)>2 || abs(tmp.cols - N)>2)
    {
        cout<<"warning: the guided image is rescaled."<<endl;
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
void DM::CurvatureGuidedFilter(const Mat & curv, const int Type, double & time, const int ItNum, const float stepsize)
{
    clock_t Tstart, Tend;

    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex, const float * p_curv);

    switch(Type)
    {
        case 0:
        {
            Local = &DM::Scheme_TV; cout<<"TV Filter: "; break;
        }
        case 1:
        {
            Local = &DM::Scheme_MC; cout<<"MC Filter: "; break;
        }
        case 2:
        {
            Local = &DM::Scheme_GC; cout<<"GC Filter: "; break;
        }
        case 3:
        {
         Local = &DM::Scheme_DC; cout<<"DC Filter: "; break;
        }
        case 4:
        {
         Local = &DM::Scheme_LS; cout<<"Bernstein Filter: "; break;
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
void DM::Solver(const int Type, double & time, const int MaxItNum, const float lambda, const float DataFitOrder, const float stepsize)
{
    clock_t Tstart, Tend;
    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex, const float * p_curv);
    void (DM::* curvature_compute)(const Mat& img, Mat& curv);

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
            Local = &DM::Scheme_TV; cout<<"TV Filter:(TVL1 by default) "; 
            curvature_compute = &DM::TV; break;
        }
        case 1:
        {
            Local = &DM::Scheme_MC; cout<<"MC Filter: "; 
            curvature_compute = &DM::MC; break;
        }
        case 2:
        {
            Local = &DM::Scheme_GC; cout<<"GC Filter: "; 
            curvature_compute = &DM::GC; break;
        }
        case 3:
        {
          Local = &DM::Scheme_DC; cout<<"DC Filter: "; 
            curvature_compute = &DM::GC; break;
        }
        case 4:
        {
          Local = &DM::Scheme_LS; cout<<"Bernstein Filter: "; 
            curvature_compute = &DM::MC; break;
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
        (this->*curvature_compute)(imgF, curvature);
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

    (this->*curvature_compute)(imgF, curvature);
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
 void DM::BlackBoxSolver(const int Type, double & time, const int MaxItNum, const float lambda, float (*BlackBox)(int row, int col, Mat& U, Mat & img_orig, float & d), const float stepsize)
 {
    clock_t Tstart, Tend;
    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex, const float * p_curv);
    void (DM::* curvature_compute)(const Mat& img, Mat& curv);

    Mat curvature = Mat::zeros(M, N, CV_32FC1);
    Mat dataFit = Mat::zeros(M, N, CV_32FC1);
    std::vector<double> energyRecord_DataFit;
    std::vector<double> energyRecord_Curvature;

    std::cout<<"********************************************"<<endl;
    std::cout<<"*** Filter Solver for Variational Models ***\n    Lambda = "<<lambda<<endl;
    std::cout<<"********************************************"<<endl;

    switch(Type)
    {
        case 0:
        {
            Local = &DM::Scheme_TV; cout<<"TV Filter:(TVL1 by default) "; 
            curvature_compute = &DM::TV; break;
        }
        case 1:
        {
            Local = &DM::Scheme_MC; cout<<"MC Filter: "; 
            curvature_compute = &DM::MC; break;
        }
        case 2:
        {
            Local = &DM::Scheme_GC; cout<<"GC Filter: "; 
            curvature_compute = &DM::GC; break;
        }
        case 3:
        {
          Local = &DM::Scheme_DC; cout<<"DC Filter: "; 
            curvature_compute = &DM::GC; break;
        }
        case 4:
        {
          Local = &DM::Scheme_LS; cout<<"Bernstein Filter: "; 
            curvature_compute = &DM::MC; break;
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
        (this->*curvature_compute)(imgF, curvature);
        energyRecord_Curvature.push_back(lambda*energy(curvature));

        for (int i = 1; i < M-1; ++i)
            for (int j = 1; j < N-1; ++j)
            {
                dataFit.at<float>(i,j) = BlackBox(i,j, imgF, image,zero);
            }
        Scalar tmp_scalar = sum(dataFit);
        energyRecord_DataFit.push_back(tmp_scalar(0));
        //if the energy starts to increase, stop the loop
        if(count>1 && (energyRecord_DataFit[it] + energyRecord_Curvature[it] > energyRecord_DataFit[it-1] + energyRecord_Curvature[it-1])) break;

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

    (this->*curvature_compute)(imgF, curvature);
    energyRecord_Curvature.push_back(lambda*energy(curvature));
    for (int i = 1; i < M-1; ++i)
        for (int j = 1; j < N-1; ++j)
        {
            dataFit.at<float>(i,j) = BlackBox(i,j, imgF, image,zero);
        }
    Scalar tmp_scalar = sum(dataFit);
    energyRecord_DataFit.push_back(tmp_scalar(0));
    

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

//*************************** Do NOT change anything! *****************************//
//************************* these filters are optimized ***************************//
//********** contact Yuanhao Gong if you need to change anything ******************//
//*************************** gongyuanhao@gmail.com *******************************//
inline void DM::GC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float scaled_stepsize = stepsize/3;
    register float tmp, min_value, com_one, com_two;
    for (int j = 1; j < N_half; ++j)
     {
        tmp = 2*p[j];
        dist[0] = p_pre[j]+p_down[j] - tmp;
        dist[1] = p_right[j-1]+p_right[j] - tmp;
        dist[2] = p_Corner[j-1]+p_rd[j] - tmp;
        dist[3] = p_Corner[j]+p_rd[j-1] - tmp;

        if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if (fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        min_value = dist[0];

        tmp *= 1.5f;
        min_value *= 1.5f;
        com_one = p_pre[j] - tmp;
        com_two = p_down[j] - tmp;
        dist[0] = p_Corner[j-1] + p_right[j-1] + com_one;
        dist[1] = p_Corner[j] + p_right[j] + com_one;
        dist[2] = p_right[j-1] + p_rd[j-1] + com_two;
        dist[3] = p_right[j] + p_rd[j] +com_two;

        if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if (fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        if (fabsf(dist[0])<fabsf(min_value)) min_value = dist[0];

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void DM::GC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[4];
    register float scaled_stepsize = stepsize/3;
    register float tmp, min_value, com_one, com_two;
    for (int j = 0; j < N_half-1; ++j)
     {
        tmp = 2*p[j];
        dist[0] = p_pre[j]+p_down[j] - tmp;
        dist[1] = p_right[j]+p_right[j+1] - tmp;
        dist[2] = p_Corner[j]+p_rd[j+1] - tmp;
        dist[3] = p_Corner[j+1]+p_rd[j] - tmp;
        
        if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if (fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        min_value = dist[0];

        tmp *= 1.5f;
        min_value *= 1.5f;
        com_one = p_pre[j] - tmp;
        com_two = p_down[j] - tmp;
        dist[0] = p_Corner[j] + p_right[j] + com_one;
        dist[1] = p_Corner[j+1] + p_right[j+1] + com_one;
        dist[2] = p_right[j] + p_rd[j] + com_two;
        dist[3] = p_right[j+1] + p_rd[j+1] + com_two;
        
        if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if (fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        if (fabsf(dist[0])<fabsf(min_value)) min_value = dist[0];

        p[j] += (scaled_stepsize*min_value);
     }
}

inline void DM::MC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8];
    register float scaled_stepsize = stepsize/8;
    for (int j = 1; j < N_half; ++j)
     {
         
        dist[4] = p[j]*8;
        dist[5] = (p_pre[j]+p_down[j])*2.5f - dist[4];
        dist[6] = (p_right[j-1]+p_right[j])*2.5f - dist[4];

        dist[0] = dist[5] + p_right[j]*5 -p_Corner[j]-p_rd[j];
        dist[1] = dist[5] + p_right[j-1]*5 -p_Corner[j-1]-p_rd[j-1];
        dist[2] = dist[6] + p_pre[j]*5 -p_Corner[j-1]-p_Corner[j];
        dist[3] = dist[6] + p_down[j]*5 -p_rd[j-1]-p_rd[j];

        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if(fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];

        p[j] += (scaled_stepsize*dist[0]);
     }
}

inline void DM::MC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8];
    register float scaled_stepsize = stepsize/8;
    for (int j = 0; j < N_half-1; ++j)
    {
         
        dist[4] = p[j]*8;
        dist[5] = (p_pre[j]+p_down[j])*2.5f - dist[4];
        dist[6] = (p_right[j]+p_right[j+1])*2.5f - dist[4];


        dist[0] = dist[5] + p_right[j]*5 -p_Corner[j]-p_rd[j];
        dist[1] = dist[5] + p_right[j+1]*5 -p_Corner[j+1]-p_rd[j+1];
        dist[2] = dist[6] + p_pre[j]*5 -p_Corner[j]-p_Corner[j+1];
        dist[3] = dist[6] + p_down[j]*5 -p_rd[j]-p_rd[j+1];
        
        if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
        if(fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
        if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        
        p[j] += (scaled_stepsize*dist[0]);
    }
}

inline void DM::LS_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
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

inline void DM::LS_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
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

inline void DM::TV_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8], tmp[4];
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five)
    register float scaled_stepsize = stepsize/6;
    register float scaledP, absMin;
    int index;
    for (int j = 1; j < N_half; ++j)
     {
        //temp var
        scaledP =  p[j]*5;
        tmp[0] = p_pre[j]+p_down[j] - scaledP;
        tmp[1] = p_right[j-1]+p_right[j] - scaledP;
        tmp[2] = p_Corner[j-1]+p_Corner[j]+p_pre[j] - scaledP;
        tmp[3] = p_rd[j-1]+p_rd[j]+p_down[j] - scaledP;

        dist[0] = tmp[0] + p_Corner[j-1]+p_rd[j-1]+p_right[j-1];
        dist[1] = tmp[0] + p_Corner[j]+p_rd[j]+p_right[j];
        dist[2] = tmp[1] + p_Corner[j-1] + p_Corner[j] + p_pre[j];
        dist[3] = tmp[1] + p_rd[j-1] + p_rd[j] + p_down[j];

        dist[4] = tmp[2] + p_right[j-1]+p_rd[j-1];
        dist[5] = tmp[2] + p_right[j]+p_rd[j];
        dist[6] = tmp[3] + p_right[j-1]+p_Corner[j-1];
        dist[7] = tmp[3] + p_right[j]+p_Corner[j];

        index = 0; absMin = fabsf(dist[0]);
        for (int i = 1; i < 8; ++i)
        {
            if(fabsf(dist[i])<absMin) {absMin = fabsf(dist[i]); index = i;}
        }

        p[j] += (scaled_stepsize*dist[index]);
     }
}

inline void DM::TV_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
{
    register float dist[8], tmp[4];
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five) 
    register float scaled_stepsize = stepsize/6;
    register float scaledP, absMin;
    int index;
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
        
        index = 0; absMin = fabsf(dist[0]);
        for (int i = 1; i < 8; ++i)
        {
            if(fabsf(dist[i])<absMin) {absMin = fabsf(dist[i]); index = i;}
        }
        
        p[j] += (scaled_stepsize*dist[index]);
     }
}


inline void DM::DC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
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

inline void DM::DC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner, const float& stepsize)
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
inline float DM::Scheme_GC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex, const float * p_curv)
{
    register float dist[4];
    register float tmp, min_value;
    tmp = -2*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) tmp = -2*(p[i] + p_curv[i]);
    min_value = p_pre[i] + p_nex[i] + tmp;
    dist[1] = p[i-1] + p[i+1] + tmp;
    dist[2] = p_pre[i-1] + p_nex[i+1] + tmp;
    dist[3]  = p_nex[i-1] + p_pre[i+1] + tmp;
    for (int j = 1; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }
    
    tmp *= 1.5f;
    min_value *= 1.5f;
    dist[0] = p_pre[i] + p_pre[i-1] + p[i-1] + tmp;
    dist[1] = p_pre[i] + p_pre[i+1] + p[i+1] + tmp;
    dist[2] = p_nex[i] + p_nex[i-1] + p[i-1] + tmp;
    dist[3] = p_nex[i] + p_nex[i+1] + p[i+1] + tmp;
    for (int j = 0; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }

    min_value /= 3;
    
    return min_value;
}

inline float DM::Scheme_MC(int i, float* p_pre, float* p, float* p_nex, const float * p_curv)
{
    //compute the movement according to half window
    //       a   b
    //       I   e
    //       c   d
    // return (2.5(a+c)+5*e)-b-d)/8.0;
    register float dist[4];
    register float tmp, com_one, com_two;
    tmp = 8*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) tmp = 8*(p[i] + p_curv[i]);

    com_one = 2.5f*(p_pre[i]+p_nex[i]) - tmp;
    com_two = 2.5f*(p[i-1]+p[i+1]) - tmp;

    dist[0] = com_one + 5*p[i+1] - p_pre[i+1] - p_nex[i+1];
    dist[1] = com_one + 5*p[i-1] - p_pre[i-1] - p_nex[i-1];
    dist[2] = com_two - p_nex[i-1] + 5*p_nex[i] - p_nex[i+1];
    dist[3] = com_two - p_pre[i-1] + 5*p_pre[i] - p_pre[i+1];

    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    if(fabsf(dist[3])<fabsf(dist[2])) dist[2] = dist[3];
    if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];

    dist[0] /= 8;
    
    return dist[0];
}

inline float DM::Scheme_LS(int i, float* p_pre, float* p, float* p_nex, const float * p_curv)
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

inline float DM::Scheme_TV(int i, float* p_pre, float* p, float* p_nex, const float * p_curv)
{
    //       a   b
    //       I   e
    //       c   d
    // return (a+b+c+d+e)/5.0;
    register float dist[4], tmp[4];
     //old fashion, need 5*8 times plus or minus
    register float scaledP, min_value;
    scaledP = 5*p[i];
    //specify the curvature if provided
    if (p_curv != NULL) scaledP = 5*(p[i] + p_curv[i]);


    tmp[0] = p_pre[i-1]+p[i-1] + p_nex[i-1] - scaledP;
    tmp[1] = p_pre[i+1]+p[i+1] + p_nex[i+1] - scaledP;
    tmp[2] = p[i-1]+p[i+1] - scaledP;
    tmp[3] = p_pre[i]+p_nex[i];

    min_value = tmp[0] + tmp[3];
    dist[1] = tmp[1] + tmp[3];
    dist[2] = tmp[2] + p_pre[i-1] + p_pre[i] + p_pre[i+1];
    dist[3] = tmp[2] + p_nex[i-1] + p_nex[i] + p_nex[i+1];
    for (int j = 1; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }

    //diag
    dist[0] = tmp[0] + p_pre[i] + p_pre[i+1];
    dist[1] = tmp[0] + p_nex[i] + p_nex[i+1];
    dist[2] = tmp[1] + p_pre[i-1] + p_pre[i];
    dist[3] = tmp[1] + p_nex[i-1] + p_nex[i];
    for (int j = 0; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }

    min_value/=6;
    // if use 6, the scheme includes central pixel;
    // if use 5, the scheme does not include central pixel (my PhD thesis uses five)

    return min_value;

/*
    tmp[0] = p_pre[i-1]+p_pre[i] + p_pre[i+1] - scaledP;
    tmp[1] = p_nex[i-1]+p_nex[i] + p_nex[i+1] - scaledP;
    tmp[2] = p[i-1]+p[i+1];
    tmp[3] = p_pre[i]+p_nex[i] - scaledP;

    min_value = tmp[3]+ p_pre[i+1]+p[i+1] + p_nex[i+1];
    dist[1] = tmp[3]+ p[i-1] + p_pre[i-1] + p_nex[i-1];
    dist[2] = tmp[2]+tmp[1];
    dist[3] = tmp[2]+tmp[0];
    for (int j = 1; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }

    //diag
    dist[0] = tmp[0] + p[i-1] + p_nex[i-1];
    dist[1] = tmp[0] + p[i+1] + p_nex[i+1];
    dist[2] = tmp[1] + p[i-1] + p_pre[i-1];
    dist[3] = tmp[1] + p[i+1] + p_pre[i+1];
    for (int j = 0; j < 4; ++j)
    {
        if(fabsf(dist[j])<fabsf(min_value)) min_value = dist[j];
    }
    
    /*
    //only need 8 + 8*3 = 32 times plut or minus, but slower than above code on my MacBook
    //this leads to a new idea that using box filter for TV filter to further reduce computation
    float all = p_pre[i-1] + p_pre[i] + p_pre[i+1] + p[i-1] + p[i+1] + p_nex[i-1] + p_nex[i] + p_nex[i+1] - 5*p[i];
    dist[0] = all - p_pre[i+1] - p[i+1] - p_nex[i+1];
    dist[1] = all - p_pre[i-1] - p[i-1] - p_nex[i-1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = all - p_pre[i-1] - p_pre[i] - p_pre[i+1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = all - p_nex[i-1] - p_nex[i] - p_nex[i+1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    //diag
    dist[1] = all - p_pre[i-1] - p_pre[i] - p[i-1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = all - p_pre[i] - p_pre[i+1] - p[i+1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = all - p_nex[i] - p_nex[i+1] - p[i+1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = all - p_nex[i-1] - p_nex[i] - p[i-1];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    return dist[0]/5;
    */
    
}

inline float DM::Scheme_DC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex, const float * p_curv)
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
