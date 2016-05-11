/*=========================================================================
 *
 *   Curvature Filter in ITK, please cite following work in you paper!
 *
 **************************************************************************  
 
            @phdthesis{gong:phd, 
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
	register float TM = 2*cur;
	float d_m = prev + down - TM;
	float tmp = left + right - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = lefUp + rigDn - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = lefDn + rigUp - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	//hybrid
	d_m *= 1.5f;
	TM *= 1.5f;
	tmp = left + lefUp + prev - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = prev + rigUp + right - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = right + rigDn + down - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;  
	tmp = down + lefDn + left - TM;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	return d_m/3;
}
float MC(float & cur, float & prev, float & rigUp, float & right, float & rigDn, float & down, float & lefDn, float & left, float& lefUp)
{
	register float TM = 8*cur;
	float var1 = (prev + down)*2.5f - TM;
	float var2 = (left + right)*2.5f - TM;

	float d_m = var1 - rigUp - rigDn + 5*right;
	float tmp = var1 - lefUp - lefDn + 5*left;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = var2 - lefUp - rigUp + 5*prev;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = var2 - lefDn - rigDn + 5*down;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	return d_m/8;
}
float TV(float & cur, float & prev, float & rigUp, float & right, float & rigDn, float & down, float & lefDn, float & left, float& lefUp)
{
	register float TM = 5*cur;
	float var1 = prev + lefUp + left - TM;
	float var2 = down + rigDn + right - TM;

	float d_m = var1 + down + lefDn;
	float tmp = var2 + prev + rigUp;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = var1 + rigUp + right;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;
	tmp = var2 + left + lefDn;
	if (fabs(tmp) < fabs(d_m)) d_m = tmp;

	//diag
	var1 = lefUp + prev + rigUp - TM;
	var2 = lefDn + down + rigDn - TM;

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

