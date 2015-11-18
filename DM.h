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
    //split imgF into four sets
    void split();
    //merge four sets back to imgF
    void merge();
    //compute TV
    void TV(Mat & img, Mat & T);
    //compute M
    void MC(Mat & img, Mat & MC);
    //compute MC
    void GC(Mat & img, Mat & GC);
    //compute energy for given TV, MC, or GC image
    double energy(Mat& img);
    //compute data fitting energy between image and imgF
    double DataFitEnergy(Mat& tmp, double order);
    //PSNR
    double PSNR();
    double PSNR(const Mat& I1, const Mat& I2);
    /******************* filters *****************************/
    // Type=0, TV; Type=1, MC; Type=2, GC; (Type=3, DC, experimental);
    void Filter(int Type, double & time, int ItNum = 10);//with split
    void FilterNoSplit(int Type, double & time, int ItNum = 10);//direct on imgF 
    /******************* generic solver for variational models *****************************/
    void Solver(int Type, double & time, int MaxItNum,  float lambda = 2, float DataFitOrder = 1);

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
	/*************************************** Split into 4 sets *********************************/
	//one is for BT and WC, two is for BC and WT
	inline void GC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	inline void GC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	
	inline void MC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	inline void MC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	
	inline void TV_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	inline void TV_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);

    inline void DC_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
    inline void DC_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	
	inline void LS_one(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	inline void LS_two(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
	
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
     //in our case, I1 and I2 are float type, in [0, 1].
     double v_min, v_max;
     minMaxLoc(I1, &v_min, &v_max);
     if (v_max>1.01)
     {
         cout<<"Input in PSNR is not correct.\n";
         return 0;
     }
     minMaxLoc(I2, &v_min, &v_max);
     if (v_max>1.01)
     {
         cout<<"Input in PSNR is not correct.\n";
         return 0;
     }

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
         double psnr = -10.0*log10(mse);
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
    tmp2 /=255;
    M_orig = tmp2.rows;
    N_orig = tmp2.cols;
    M = (int)ceil(M_orig/2.0)*2;
    N = (int)ceil(N_orig/2.0)*2;
    M_half = M/2;
    N_half = N/2;


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
    tmp2 /=255;
    M_orig = tmp2.rows;
    N_orig = tmp2.cols;
    M = (int)ceil(M_orig/2.0)*2;
    N = (int)ceil(N_orig/2.0)*2;
    M_half = M/2;
    N_half = N/2;


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
    Mat tmp2 = imgF*255;
    tmp2(Range(0,M_orig),Range(0,N_orig)).convertTo(tmp, CV_8UC1);

    vector<int> params;
    params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    params.push_back(0);

    imwrite(FileName, tmp, params);
}

//compute Total Variation
void DM::TV(Mat & imgF, Mat & T)
{
    float * p_row, * pn_row, * p_t;
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
void DM::MC(Mat& imgF, Mat & MC)
{
	//classical scheme is used
    float * p_row, *pn_row, *pp_row, *p_d;
    float Ix, Iy, Ixx, Iyy, num, den;
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
            Iyy = pn_row[j] -2*p_row[j] + pp_row[j];
            
            num = Ixx + Iyy;
            den = (1.0f + Ix*Ix + Iy*Iy);
            p_d[j] = num/den;
        }   
    }
}

//compute Gaussian curvature
void DM::GC(Mat & imgF, Mat &GC)
{
	//classical scheme is used
    float * p_row, *pn_row, *pp_row, *p_d;
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
            den = powf(den, 1.5f);
            p_d[j] = num/den;
        }   
    }
}

//compute the curvature energy
double DM::energy(Mat &img)
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

void DM::Filter(int Type, double & time, int ItNum )
{
    clock_t Tstart, Tend;

    void (DM::* Local_one)(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);
    void (DM::* Local_two)(float* p, float* p_right, float* p_down, float *p_rd, float* p_pre, float* p_Corner);

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
                (this->*Local_two)(p, p_right, p_down, p_rd, p_pre, p_Corner);
        }
        //BT
        for (int i = 1; i < M_half; ++i)
        {
                p = BT.ptr<float>(i); p_right = WT.ptr<float>(i);
                p_down = WC.ptr<float>(i); p_rd = BC.ptr<float>(i); 
                p_pre = WC.ptr<float>(i-1); p_Corner = BC.ptr<float>(i-1);
                (this->*Local_one)(p, p_right, p_down, p_rd, p_pre, p_Corner);
        }
        //WC
        for (int i = 0; i < M_half-1; ++i)
        {
                p = WC.ptr<float>(i); p_right = BC.ptr<float>(i);
                p_down = BT.ptr<float>(i+1); p_rd = WT.ptr<float>(i+1); 
                p_pre = BT.ptr<float>(i); p_Corner = WT.ptr<float>(i);
                (this->*Local_one)(p, p_right, p_down, p_rd, p_pre, p_Corner);
        }
    	//WT
        for (int i = 1; i < M_half; ++i)
        {
                p = WT.ptr<float>(i); p_right = BT.ptr<float>(i);
                p_down = BC.ptr<float>(i); p_rd = WC.ptr<float>(i); 
                p_pre = BC.ptr<float>(i-1); p_Corner = WC.ptr<float>(i-1);
                (this->*Local_two)(p, p_right, p_down, p_rd, p_pre, p_Corner);
        }
        
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
}

//this nosplit is very useful for tasks like deconvolution, where the four sets need to be merged 
//every iteration if we use the split scheme.
void DM::FilterNoSplit(int Type, double & time, int ItNum )
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
                p[j] += d;
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
                p[j] += d;
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
                p[j] += d;
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
                p[j] += d;
            }
        }
    }
    Tend = clock() - Tstart;   
    time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
}

//generic filter solver for variational model |U - I|^DataFitOrder + lambda * Regularization
//the DataFitOrder can be fractional such as 1.5, which is not possible for other solvers.
void DM::Solver(int Type, double & time, int MaxItNum, float lambda, float DataFitOrder)
{
    clock_t Tstart, Tend;
    float (DM::* Local)(int i, float* p_pre, float* p, float* p_nex);
    void (DM::* curvature_compute)(Mat& img, Mat& curv);

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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
                d = (this->*Local)(j,p_pre,p,p_down);
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
inline void DM::GC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[10];
	for (int j = 1; j < N_half; ++j)
     {
		dist[0] = 2*p[j];
     	dist[1] = (p_pre[j]+p_down[j]) - dist[0];
		dist[2] = (p_right[j-1]+p_right[j]) - dist[0];
		dist[3] = (p_Corner[j-1]+p_rd[j]) - dist[0];
		dist[4] = (p_Corner[j]+p_rd[j-1]) - dist[0];
		

		dist[5] = 3*p[j];
		dist[6] = (p_Corner[j-1] + p_pre[j] + p_right[j-1]) - dist[5];
		dist[7] = (p_Corner[j] + p_pre[j] + p_right[j]) - dist[5];
		dist[8] = (p_right[j-1] + p_rd[j-1] + p_down[j]) - dist[5];
		dist[9] = (p_right[j] + p_down[j] + p_rd[j]) - dist[5];
        
		for (int i = 2; i < 5; ++i)
		{
		    if (fabsf(dist[i])<fabsf(dist[1])) dist[1] = dist[i];
		    if (fabsf(dist[i+5])<fabsf(dist[6])) dist[6] = dist[i+5];
		}
		dist[1] /= 2;
		dist[6] /= 3;
		if (fabsf(dist[6])<fabsf(dist[1])) dist[1] = dist[6];

		p[j] += dist[1];
     }
}

inline void DM::GC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[10];
	for (int j = 0; j < N_half-1; ++j)
     {
		dist[0] = 2*p[j];
     	dist[1] = (p_pre[j]+p_down[j]) - dist[0];
		dist[2] = (p_right[j]+p_right[j+1]) - dist[0];
		dist[3] = (p_Corner[j]+p_rd[j+1]) - dist[0];
		dist[4] = (p_Corner[j+1]+p_rd[j]) - dist[0];
		

		dist[5] = 3*p[j];
		dist[6] = (p_Corner[j] + p_pre[j] + p_right[j]) - dist[5];
		dist[7] = (p_Corner[j+1] + p_pre[j] + p_right[j+1]) - dist[5];
		dist[8] = (p_right[j] + p_rd[j] + p_down[j]) - dist[5];
		dist[9] = (p_right[j+1] + p_down[j] + p_rd[j+1]) - dist[5];
		
		for (int i = 2; i < 5; ++i)
		{
		    if (fabsf(dist[i])<fabsf(dist[1])) dist[1] = dist[i];
		    if (fabsf(dist[i+5])<fabsf(dist[6])) dist[6] = dist[i+5];
		}
		dist[1] /= 2;
		dist[6] /= 3;
		if (fabsf(dist[6])<fabsf(dist[1])) dist[1] = dist[6];

		p[j] += dist[1];
     }
}

inline void DM::MC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[8];
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
          if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
          if(fabsf(dist[3])<fabsf(dist[0])) dist[0] = dist[3];

		dist[0] /= 8;

		p[j] += dist[0];
     }
}

inline void DM::MC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[8];
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
        if(fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
        if(fabsf(dist[3])<fabsf(dist[0])) dist[0] = dist[3];
        
		dist[0] /= 8;
		
		p[j] += dist[0];
     }
}

inline void DM::LS_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	//one is for BT and WC, two is for BC and WT
	register float dist[5];
	for (int j = 1; j < N_half; ++j)
     {
     	
		dist[4] = p[j]*2;
		dist[0] = (p_pre[j]+p_down[j]) - dist[4];
		dist[1] = (p_right[j-1]+p_right[j]) - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
		dist[1] = (p_Corner[j-1]+p_rd[j]) - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
		dist[1] = (p_Corner[j]+p_rd[j-1]) - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

		dist[4] *= 3.5;
		dist[0] *= 3.5;

		dist[2] = 3*(p_Corner[j] + p_rd[j-1]) - dist[4];
		dist[3] = 3*(p_Corner[j-1] + p_rd[j]) - dist[4];


		dist[1] = p_right[j-1] + p_pre[j] - p_Corner[j-1] + dist[2];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

		dist[1] = p_right[j] + p_pre[j] - p_Corner[j] + dist[3];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

		dist[1] = p_right[j] + p_down[j] - p_rd[j] + dist[2];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

		dist[1] = p_right[j-1] + p_down[j]  - p_rd[j-1] + dist[3];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];


		dist[0] /= 7;

		p[j] += dist[0];
     }
}

inline void DM::LS_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	//one is for BT and WC, two is for BC and WT
	register float dist[5];
	for (int j = 0; j < N_half-1; ++j)
     {
     	
		dist[4] = p[j]*2;
		dist[0] = p_pre[j]+p_down[j] - dist[4];
		dist[1] = p_right[j]+p_right[j+1] - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
		dist[1] = p_Corner[j]+p_rd[j+1] - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
		dist[1] = p_Corner[j+1]+p_rd[j] - dist[4];
		if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

		dist[4] *= 3.5;
		dist[0] *= 3.5;

		dist[2] = 3*(p_Corner[j+1] + p_rd[j]) - dist[4];
		dist[3] = 3*(p_Corner[j] + p_rd[j+1]) - dist[4];

     	dist[0] = p_right[j] + p_pre[j] - p_Corner[j] + dist[2];
     	if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

     	dist[0] = p_right[j+1] + p_pre[j] - p_Corner[j+1] + dist[3];
     	if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

     	dist[0] = p_right[j+1] + p_down[j] - p_rd[j+1] + dist[2];
     	if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

     	dist[0] = p_right[j] + p_down[j] - p_rd[j] + dist[3];
     	if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

        
		dist[0] /= 7;
		
		p[j] += dist[0];
     }
}

inline void DM::TV_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[13]; // first eight are distances, the rest are temp
	for (int j = 1; j < N_half; ++j)
     {
		//temp var
		dist[8] =  p[j]*5;
		dist[9] = p_pre[j]+p_down[j] - dist[8];
		dist[10] = p_right[j-1]+p_right[j] - dist[8];
		dist[11] = p_Corner[j-1]+p_Corner[j]+p_pre[j] - dist[8];
		dist[12] = p_rd[j-1]+p_rd[j]+p_down[j] - dist[8];

     	dist[0] = dist[9] + p_Corner[j-1]+p_rd[j-1]+p_right[j-1];
		dist[1] = dist[9] + p_Corner[j]+p_rd[j]+p_right[j];
		dist[2] = dist[10] + p_Corner[j-1] + p_Corner[j] + p_pre[j];
		dist[3] = dist[10] + p_rd[j-1] + p_rd[j] + p_down[j];

		dist[4] = dist[11] + p_right[j-1]+p_rd[j-1];
		dist[5] = dist[11] + p_right[j]+p_rd[j];
		dist[6] = dist[12] + p_right[j-1]+p_Corner[j-1];
		dist[7] = dist[12] + p_right[j]+p_Corner[j];


        for (int i = 0; i < 8; i+=2)
		{
		    if (fabsf(dist[i+1])<fabsf(dist[i])) dist[i] = dist[i+1];
		}
		if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
		if (fabsf(dist[6])<fabsf(dist[4])) dist[4] = dist[6];
		if (fabsf(dist[4])<fabsf(dist[0])) dist[0] = dist[4];

		dist[0] /= 5;
		p[j] += dist[0];
     }
}

inline void DM::TV_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float* __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
{
	register float dist[13]; // first eight are distances, the rest are temp
	for (int j = 0; j < N_half-1; ++j)
     {
		dist[8] = p[j]*5;
		dist[9] = p_pre[j]+p_down[j] - dist[8];
		dist[10] = p_right[j]+p_right[j+1] - dist[8];
		dist[11] = p_Corner[j]+p_Corner[j+1]+p_pre[j] - dist[8];
		dist[12] = p_rd[j]+p_rd[j+1]+p_down[j] - dist[8];

     	dist[0] = dist[9] + p_Corner[j]+p_rd[j]+p_right[j];
		dist[1] = dist[9] + p_Corner[j+1]+p_rd[j+1]+p_right[j+1];
		dist[2] = dist[10] + p_Corner[j] + p_Corner[j+1] + p_pre[j];
		dist[3] = dist[10] + p_rd[j] + p_rd[j+1] + p_down[j];
        
		dist[4] = dist[11] +p_right[j]+p_rd[j];
		dist[5] = dist[11] +p_right[j+1]+p_rd[j+1];
		dist[6] = dist[12] +p_right[j]+p_Corner[j];
		dist[7] = dist[12] +p_right[j+1]+p_Corner[j+1];

        for (int i = 0; i < 8; i+=2)
		{
		    if (fabsf(dist[i+1])<fabsf(dist[i])) dist[i] = dist[i+1];
		}
		if (fabsf(dist[2])<fabsf(dist[0])) dist[0] = dist[2];
		if (fabsf(dist[6])<fabsf(dist[4])) dist[4] = dist[6];
		if (fabsf(dist[4])<fabsf(dist[0])) dist[0] = dist[4];
        

		dist[0] /= 5;
		p[j] += dist[0];
     }
}


inline void DM::DC_one(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
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

        p[j] += dist[0];
     }
}

inline void DM::DC_two(float* __restrict p, float* __restrict p_right, float* __restrict p_down, float * __restrict p_rd, float* __restrict p_pre, float* __restrict p_Corner)
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


        p[j] += dist[0];
     }
}

/********************************************************************/
/********************** scheme at each pixel ************************/
/********************** only for noSplit case ***********************/
/********************************************************************/
inline float DM::Scheme_GC(int i, float * __restrict p_pre, float * __restrict p, float * __restrict p_nex)
{
	float dist[4];
    dist[2] = 2*p[i];
    dist[0] = (p_pre[i] + p_nex[i]) - dist[2];
    dist[1] = (p[i-1] + p[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    
    dist[1] = (p_pre[i-1] + p_nex[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1]  = (p_nex[i-1] + p_pre[i+1]) - dist[2];
    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[0] /= 2;
    
    
    
    dist[3] = 3*p[i];
    dist[1] = (p_pre[i] + p_pre[i-1] + p[i-1])- dist[3];
    
    dist[2] = (p_pre[i] + p_pre[i+1] + p[i+1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[2] = (p_nex[i] + p_nex[i-1] + p[i-1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[2] = (p_nex[i] + p_nex[i+1] + p[i+1])- dist[3];
    if (fabsf(dist[2])<fabsf(dist[1])) dist[1] = dist[2];
    dist[1] /= 3;

    if (fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    
    return dist[0];
}

inline float DM::Scheme_MC(int i, float* p_pre, float* p, float* p_nex)
{
    //compute the movement according to half window
    //       a   b
    //       I   e
    //       c   d
    // return (2.5(a+c)+5*e)-b-d)/8.0;
	float dist[2];
    float tmp = 8*p[i];

    dist[0] = 2.5f*(p_pre[i]+p_nex[i]) + 5*p[i+1] - p_pre[i+1] - p_nex[i+1] - tmp;
    dist[1] = 2.5f*(p_pre[i]+p_nex[i]) + 5*p[i-1] - p_pre[i-1] - p_nex[i-1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = 2.5f*(p[i-1]+p[i+1]) + 5*p_nex[i] - p_nex[i-1] - p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = 2.5f*(p[i-1]+p[i+1]) + 5*p_pre[i]- p_pre[i-1] - p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[0] /= 8.0f;
    
    return dist[0];
}

inline float DM::Scheme_LS(int i, float* p_pre, float* p, float* p_nex)
{
    //compute the movement according to half window
    //   f   a   b            0 1/2 0               3/7 1/7 -1/7
    //       I   e               -1 0                    -1  1/7
    //       c   d              1/2 0                        3/7
    // 
    float dist[4];
    float tmp = 2*p[i];

    dist[0] = p_pre[i]+p_nex[i] - tmp;
    dist[1] = p[i-1] + p[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_pre[i-1] + p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1] + p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    

    
    tmp *= 3.5;
    dist[0] *= 3.5;

    dist[2] = 3*(p_pre[i+1] + p_nex[i-1]) - tmp;
    dist[3] = 3*(p_pre[i-1] + p_nex[i+1]) - tmp;

    dist[1] = p_pre[i] - p_pre[i-1] +  p[i-1] + dist[2];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_pre[i] - p_pre[i+1] + p[i+1] + dist[3];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_nex[i] - p_nex[i-1] + p[i-1] + dist[3];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p_nex[i]  - p_nex[i+1] + p[i+1] + dist[2];
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[0] /= 7;

    return dist[0];
}

inline float DM::Scheme_TV(int i, float* p_pre, float* p, float* p_nex)
{
    //       a   b
    //       I   e
    //       c   d
    // return (a+b+c+d+e)/5.0;
	float dist[2];
    float tmp = 5*p[i];

    dist[0] = p_pre[i]+ p_pre[i+1]+p[i+1]+p_nex[i] + p_nex[i+1] - tmp;

    dist[1] = p_pre[i]+p[i-1]+p_nex[i] + p_pre[i-1] + p_nex[i-1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    dist[1] = p[i-1]+p[i+1]+p_nex[i] + p_nex[i-1] + p_nex[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p[i-1]+p[i+1]+p_pre[i] + p_pre[i-1] + p_pre[i+1] - tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];

    //diag
    dist[1] = p_pre[i-1]+p_pre[i] + p_pre[i+1] + p[i-1] + p_nex[i-1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_pre[i-1]+p_pre[i] + p_pre[i+1] + p[i+1] + p_nex[i+1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1]+p_nex[i] + p_nex[i+1] + p[i-1] + p_pre[i-1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];
    dist[1] = p_nex[i-1]+p_nex[i] + p_nex[i+1] + p[i+1] + p_pre[i+1]- tmp;
    if(fabsf(dist[1])<fabsf(dist[0])) dist[0] = dist[1];


    dist[0]/=5.0;
    return dist[0];
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
