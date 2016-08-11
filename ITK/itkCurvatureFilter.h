/*=========================================================================
 *
 *   Curvature Filter in ITK, please cite following work in you paper!
 *
 **************************************************************************  
 
            @phdthesis{gong:phd,
             title={Spectrally regularized surfaces}, 
             author={Gong, Yuanhao}, 
             year={2015}, 
             school={ETH Zurich, Nr. 22616},
             note={http://dx.doi.org/10.3929/ethz-a-010438292}}

 *=========================================================================*/
//2D curvature filter
void CurvatureFilter(itk::Image< float, 2 >::Pointer image, int FilterType, int IterationNumber=10, float stepsize=1);


//estimate the d_m by different curvature filters for a 2D image
float GC(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn);
float MC(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn);
float TV(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn);
float (*projection)(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn);

void CurvatureFilter(itk::Image< float, 2 >::Pointer image, int FilterType, int IterationNumber, float stepsize)
{

    if (FilterType == 0) {projection = TV; std::cout<<"TV filter: with "<<IterationNumber<<" Iterations";}
    if (FilterType == 1) {projection = MC; std::cout<<"MC filter: with "<<IterationNumber<<" Iterations";}
    if (FilterType == 2) {projection = GC; std::cout<<"GC filter: with "<<IterationNumber<<" Iterations";}
    float d, v; 
    //image size and neighbor offset
    itk::Image< float, 2 >::SizeType size = image->GetLargestPossibleRegion().GetSize();
    itk::Image< float, 2 >::IndexType index;
    itk::Image< float, 2 >::OffsetType off_prev, off_rigUp, off_right, off_rigDn, off_next, off_leftDn, off_left, off_leftUp;

    off_prev[0] = 0; off_prev[1] = -1;
    off_rigUp[0] = 1; off_rigUp[1] = -1;
    off_right[0] = 1; off_right[1] = 0;
    off_rigDn[0] = 1; off_rigDn[1] = 1;
    off_next[0] = 0; off_next[1] = 1;
    off_leftDn[0] = -1; off_leftDn[1] = 1;
    off_left[0] = -1; off_left[1] = 0;
    off_leftUp[0] = -1; off_leftUp[1] = -1;

    clock_t Tstart, Tend;
    Tstart = clock();

    for(int it = 0; it < IterationNumber; ++it)
    {
      //domain decomposition, four sets
      for (int start_r = 1; start_r < 3; start_r++)
        for (int start_c = 1; start_c <3; start_c++)
          for (int i = start_r; i < size[1]-1; ++i, ++i)//row
          {
              index[1]=i;
              for (int j = start_c; j < size[0]-1; ++j, ++j)//col
              {
                  index[0] = j;
                  v = image->GetPixel(index);

                  d = (*projection)(image->GetPixel(index + off_leftUp), image->GetPixel(index + off_prev), 
                      image->GetPixel(index + off_rigUp), image->GetPixel(index + off_left), v,
                      image->GetPixel(index + off_right), image->GetPixel(index + off_leftDn), 
                      image->GetPixel(index + off_next), image->GetPixel(index + off_rigDn));

                  image->SetPixel(index, v + stepsize*d);
              }
          }
    }
    Tend = clock() - Tstart;   
    double time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
    std::cout<<" running time is "<<time<<" milliseconds"<<std::endl;
}

// *************************** functions ********************************************
//estimate d_m (minimal projection signed distance)
float GC(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn)
{
    register float TM = 2*cur;
    float d_m = pre + next - TM;
    float tmp = left + right - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = lefUp + rigDn - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = lefDn + rigUp - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;

    //hybrid
    d_m *= 1.5f;
    TM *= 1.5f;
    tmp = left + lefUp + pre - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = pre + rigUp + right - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = right + rigDn + next - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;  
    tmp = next + lefDn + left - TM;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;

    return d_m/3;
}
float MC(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn)
{
    float TM = 2*cur;
    float var1 = (pre + next)/2 - TM;
    float var2 = (left + right)/2 - TM;

    float d_m = var1 + right;
    float tmp = var1 + left;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var2 + pre;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var2 + next;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;

    return d_m/2;
}
float TV(float & lefUp, float & pre, float & rigUp, float & left, float & cur, float & right, float & lefDn, float & next, float & rigDn)
{
    register float TM = 5*cur;
    float var1 = pre + lefUp + left - TM;
    float var2 = next + rigDn + right - TM;

    float d_m = var1 + next + lefDn;
    float tmp = var2 + pre + rigUp;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var1 + rigUp + right;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var2 + left + lefDn;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;

    //diag
    var1 = lefUp + pre + rigUp - TM;
    var2 = lefDn + next + rigDn - TM;

    tmp = var1 + left + lefDn;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var1 + right + rigDn;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var2 + left + lefUp;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;
    tmp = var2 + right + rigUp;
    if (fabs(tmp) < fabs(d_m)) d_m = tmp;

    return d_m/5;
}

