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
//estimate the d_m by different curvature filters
float GC(float & cur, float & pre, float & rigUp, float & right, float & rigDn, float & next, float & lefDn, float & left, float & lefUp);
float MC(float & cur, float & pre, float & rigUp, float & right, float & rigDn, float & next, float & lefDn, float & left, float & lefUp);
float TV(float & cur, float & pre, float & rigUp, float & right, float & rigDn, float & next, float & lefDn, float & left, float & lefUp);
float (*projection)(float & cur, float & pre, float & rigUp, float & right, float & rigDn, float & next, float & lefDn, float & left, float & lefUp);
void CurvatureFilter(itk::Image< float, 2 >::Pointer image, int FilterType, int IterationNumber)
{

  if (FilterType == 0) {projection = TV; std::cout<<"TV filter: with "<<IterationNumber<<" Iterations";}
  if (FilterType == 1) {projection = MC; std::cout<<"MC filter: with "<<IterationNumber<<" Iterations";}
  if (FilterType == 2) {projection = GC; std::cout<<"GC filter: with "<<IterationNumber<<" Iterations";}
  float d, cur; 
  //image size and neighbor index
  itk::Image< float, 2 >::SizeType size = image->GetLargestPossibleRegion().GetSize();
  itk::Image< float, 2 >::IndexType index_cur, index_prev, index_rigUp, index_right, index_rigDn, index_next, index_leftDn, index_left, index_leftUp;

  clock_t Tstart, Tend;
  Tstart = clock();
  
  for(int it = 0; it < IterationNumber; ++it)
  {
	  //domain decomposition, four sets
	  for (int i = 1; i < size[0]; ++i, ++i)
  		for (int j = 1; j < size[1]; ++j, ++j)
  		{
  			index_leftUp[0] = j-1, index_leftUp[1] = i-1, index_prev[0] = j, index_prev[1] = i-1,  index_rigUp[0]=j+1, index_rigUp[1] = i-1;
  			index_left[0]=j-1, index_left[1]=i, index_cur[0] = j, index_cur[1] = i, index_right[0]=j+1, index_right[1]=i;
  			index_leftDn[0]=j-1, index_leftDn[1]=i+1, index_next[0]=j, index_next[1] = i+1, index_rigDn[0]=j+1, index_rigDn[1] = i+1;

			cur = image->GetPixel(index_cur);
  			d = (*projection)(cur, image->GetPixel(index_prev), image->GetPixel(index_rigUp), image->GetPixel(index_right), image->GetPixel(index_rigDn), image->GetPixel(index_next), image->GetPixel(index_leftDn), image->GetPixel(index_left), image->GetPixel(index_leftUp));
			image->SetPixel(index_cur, cur + d);
  		}
	  for (int i = 2; i < size[0]; ++i, ++i)
  		for (int j = 2; j < size[1]; ++j, ++j)
  		{
  			index_leftUp[0] = j-1, index_leftUp[1] = i-1, index_prev[0] = j, index_prev[1] = i-1,  index_rigUp[0]=j+1, index_rigUp[1] = i-1;
  			index_left[0]=j-1, index_left[1]=i, index_cur[0] = j, index_cur[1] = i, index_right[0]=j+1, index_right[1]=i;
  			index_leftDn[0]=j-1, index_leftDn[1]=i+1, index_next[0]=j, index_next[1] = i+1, index_rigDn[0]=j+1, index_rigDn[1] = i+1;

			cur = image->GetPixel(index_cur);
  			d = (*projection)(cur, image->GetPixel(index_prev), image->GetPixel(index_rigUp), image->GetPixel(index_right), image->GetPixel(index_rigDn), image->GetPixel(index_next), image->GetPixel(index_leftDn), image->GetPixel(index_left), image->GetPixel(index_leftUp));
			image->SetPixel(index_cur, cur + d);
  		}
	  for (int i = 1; i < size[0]; ++i, ++i)
  		for (int j = 2; j < size[1]; ++j, ++j)
  		{
  			index_leftUp[0] = j-1, index_leftUp[1] = i-1, index_prev[0] = j, index_prev[1] = i-1,  index_rigUp[0]=j+1, index_rigUp[1] = i-1;
  			index_left[0]=j-1, index_left[1]=i, index_cur[0] = j, index_cur[1] = i, index_right[0]=j+1, index_right[1]=i;
  			index_leftDn[0]=j-1, index_leftDn[1]=i+1, index_next[0]=j, index_next[1] = i+1, index_rigDn[0]=j+1, index_rigDn[1] = i+1;

			cur = image->GetPixel(index_cur);
  			d = (*projection)(cur, image->GetPixel(index_prev), image->GetPixel(index_rigUp), image->GetPixel(index_right), image->GetPixel(index_rigDn), image->GetPixel(index_next), image->GetPixel(index_leftDn), image->GetPixel(index_left), image->GetPixel(index_leftUp));
			image->SetPixel(index_cur, cur + d);
  		}
	  for (int i = 2; i < size[0]; ++i, ++i)
  		for (int j = 1; j < size[1]; ++j, ++j)
  		{
  			index_leftUp[0] = j-1, index_leftUp[1] = i-1, index_prev[0] = j, index_prev[1] = i-1,  index_rigUp[0]=j+1, index_rigUp[1] = i-1;
  			index_left[0]=j-1, index_left[1]=i, index_cur[0] = j, index_cur[1] = i, index_right[0]=j+1, index_right[1]=i;
  			index_leftDn[0]=j-1, index_leftDn[1]=i+1, index_next[0]=j, index_next[1] = i+1, index_rigDn[0]=j+1, index_rigDn[1] = i+1;

			cur = image->GetPixel(index_cur);
  			d = (*projection)(cur, image->GetPixel(index_prev), image->GetPixel(index_rigUp), image->GetPixel(index_right), image->GetPixel(index_rigDn), image->GetPixel(index_next), image->GetPixel(index_leftDn), image->GetPixel(index_left), image->GetPixel(index_leftUp));
			image->SetPixel(index_cur, cur + d);
  		}
  }

  Tend = clock() - Tstart;   
  double time = double(Tend)/(CLOCKS_PER_SEC/1000.0);
  std::cout<<" running time is "<<time<<" milliseconds"<<std::endl;
}

// *************************** functions ********************************************
//estimate d_m (minimal projection signed distance)
float GC(float & cur, float & prev, float & rigUp, float & right, float & rigDn, float & down, float & lefDn, float & left, float& lefUp)
{
	float d_m = (prev + down)/2 - cur;
	float tmp = (left, right)/2 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (lefUp, rigDn)/2 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (lefDn, rigUp)/2 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	//hybrid
	tmp = (left + lefUp + prev)/3 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (prev + rigUp + right)/3 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (right + rigDn + down)/3 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;  
	tmp = (down + lefDn + left)/3 - cur;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	return d_m;
}
float MC(float & cur, float & prev, float & rigUp, float & right, float & rigDn, float & down, float & lefDn, float & left, float& lefUp)
{
	float TM = 8*cur;

	float d_m = (prev + down)*2.5f - rigUp - rigDn + 5*right - TM;
	float tmp = (prev + down)*2.5f - lefUp - lefDn + 5*left - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (left + right)*2.5f - lefUp - rigUp + 5*prev - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = (left + right)*2.5f - lefDn - rigDn + 5*down - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	return d_m/8;
}
float TV(float & cur, float & prev, float & rigUp, float & right, float & rigDn, float & down, float & lefDn, float & left, float& lefUp)
{
	float TM = 5*cur;

	float d_m = prev + down + lefUp + left + lefDn - TM;
	float tmp = prev + down + rigUp + rigDn + right - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = left + lefUp + prev + rigUp + right - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = left + lefDn + down + rigDn + right - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	//diag
	tmp = lefUp + prev + rigUp + left + lefDn - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = lefUp + prev + rigUp + right + rigDn - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = lefDn + down + rigDn + left + lefUp - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = lefDn + down + rigDn + right + rigUp - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	return d_m/5;
}

