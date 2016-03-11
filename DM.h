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
    //compute GC
    void GC(const Mat & img, Mat & GC);
    void GC_new(const Mat & img, Mat & MC);//new scheme, Eq.6.16 in my thesis
    //compute energy for given TV, MC, or GC image
    double energy(const Mat& img);
    //compute data fitting energy between image and imgF
    double DataFitEnergy(Mat& tmp, double order);
    //PSNR
    double PSNR();
    double PSNR(const Mat& I1, const Mat& I2);
    /******************* filters *****************************/
    // Type=0, TV; Type=1, MC; Type=2, GC; (Type=3, DC, experimental);
    // the stepsize parameter is in (0,1]:smaller means more iterations, but reaches lower energy level; larger means less iterations, but converges at higher energy level
    //////////////////////////////////////////////////////////
    void Filter(const int Type, double & time, const int ItNum = 10, const float stepsize=1);//with split
    void FilterNoSplit(const int Type, double & time, const int ItNum = 10, const float stepsize=1);//direct on imgF 
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
    inline float Scheme_GC(int i, float * p_pre, float * p, float * p_nex);
    inline float Scheme_MC(int i, float * p_pre, float * p, float * p_nex);
    inline float Scheme_TV(int i, float * p_pre, float * p, float * p_nex);
    inline float Scheme_DC(int i, float * p_pre, float * p, float * p_nex);
    inline float Scheme_LS(int i, float * p_pre, float * p, float * p_nex);
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
            den = sqrt(tmp)*tmp/2;
            p_d[j] = num/den;
        }   
    }
}

//compute Mean Curvature, Eq.6.12 in my thesis
void DM::MC_new(const Mat& imgF, Mat & MC)
{
    //new scheme used
    const float * p_row, *pn_row, *pp_row;
    float * p_d;
    float one, five;
    for(int i = 1; i < imgF.rows-1; i++)
    {
        p_row = imgF.ptr<float>(i);
        pn_row = imgF.ptr<float>(i+1);
        pp_row = imgF.ptr<float>(i-1);
        p_d = MC.ptr<float>(i);
        
        for(int j = 1; j < imgF.cols-1; j++)
        {
            one = pp_row[j-1] + pp_row[j+1] + pn_row[j-1] + pn_row[j+1];
            five = pp_row[j] + p[j-1] + p[j+1] + pn_row[j];

            p_d[j] = 0.3125f*five - 0.0625f*one - p_row[j];
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
void DM::GC_new(const Mat & img, Mat & MC)
{
    const float * p_row, *pn_row, *pp_row;
    float * p_d;
    float six[6];
    float sum_ud, sum_lr, sum_diag_one, sum_diag_two;
    float diff_ud, diff_lr, diff_diag_one, diff_diag_two;
    float tmp, tmp2;
    for(int i = 1; i < imgF.rows-1; ++i)
    {
        p_row = imgF.ptr<float>(i);
        pn_row = imgF.ptr<float>(i+1);
        pp_row = imgF.ptr<float>(i-1);
        p_d = MC.ptr<float>(i);
        
        for(int j = 1; j < imgF.cols-1; ++j)
        {
            sum_lr = p_row[j-1] + p_row[j+1];
            diff_lr= p_row[j-1] - p_row[j+1];
            sum_ud = pp_row[j] + pn_row[j];
            diff_ud= pp_row[j] - pn_row[j];
            sum_diag_one = pp_row[j-1] + pn_row[j+1];
            diff_diag_one= pp_row[j-1] - pn_row[j+1];
            sum_diag_two = pp_row[j+1] + pn_row[j-1];
            diff_diag_two= pp_row[j+1] - pn_row[j-1];

            tmp = sum_lr + sum_ud;
            tmp2= sum_diag_one + sum_diag_two;

            six[0] = tmp*0.254187f + tmp2*0.00210414f - 1.02516f*p[j];
            six[1] = tmp*0.229983f - tmp2*0.286419f - 0.225743f*p[j]; 
            six[2] = (sum_diag_two - sum_diag_one)*0.306186f;
            six[3] = (sum_lr - sum_ud)*0.176777f;
            six[4] = (diff_diag_one - diff_diag_two)/4 - diff_lr/2;
            six[5] = diff_ud/2 - (diff_diag_one+diff_diag_two)/4;

            six[0]*=six[0];
            six[1]*=six[1];
            six[2]*=six[2];
            six[3]*=six[3];
            six[4]*=six[4];
            six[5]*=six[5];

            six[4] += six[5];
            six[3] += six[4];
            six[2] += six[3];
            six[1] += six[2];
            six[0] -= six[1];
            p_d[j] = six[0];
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

    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex);

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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
                d = stepsize*(this->*Local)(j,p_pre,p,p_down);
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
inline float DM::Scheme_GC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex)
{
    register float dist[4];
    register float tmp, min_value;
    tmp = -2*p[i];
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

inline float DM::Scheme_MC(int i, float* p_pre, float* p, float* p_nex)
{
    //compute the movement according to half window
    //       a   b
    //       I   e
    //       c   d
    // return (2.5(a+c)+5*e)-b-d)/8.0;
    register float dist[4];
    register float tmp = 8*p[i];
    register float com_one, com_two;
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

inline float DM::Scheme_LS(int i, float* p_pre, float* p, float* p_nex)
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
    register float tmp = 2*p[i], min_value;
    register float tmp_one, tmp_two;

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

inline float DM::Scheme_TV(int i, float* p_pre, float* p, float* p_nex)
{
    //       a   b
    //       I   e
    //       c   d
    // return (a+b+c+d+e)/5.0;
    register float dist[4], tmp[4];
     //old fashion, need 5*8 times plus or minus
    register float scaledP = 5*p[i], min_value;
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

inline float DM::Scheme_DC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex)
{
    float dist[2];
    float weight = -0.225603f;

    dist[0] = (p_pre[i] + p_nex[i])/2 + (p_pre[i+1] + p_nex[i+1] - 2*p[i+1])*weight - p[i];
    dist[1] = (p_pre[i] + p_nex[i])/2 + (p_pre[i-1] + p_nex[i-1] - 2*p[i-1])*weight - p[i];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = (p[i-1] + p[i+1])/2 + (p_pre[i-1] + p_pre[i+1] - 2*p_pre[i])*weight - p[i];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = (p[i-1] + p[i+1])/2 + (p_nex[i-1] + p_nex[i+1] - 2*p_nex[i])*weight - p[i];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    return dist[0];
}
