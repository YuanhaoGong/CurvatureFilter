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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include <time.h>
#include "itkCurvatureFilter.h"

int FType = 2;

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  outputImageFile filterType numberOfIterations" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int numberOfIterations = atoi( argv[4] );
  const char *   filterType = argv[3];

  if (*filterType == 't') {FType = 0; }
  if (*filterType == 'm') {FType = 1; }
  if (*filterType == 'g') {FType = 2; }

  //  read image
  typedef float                             PixelType;
  typedef itk::Image< PixelType, 2 >        ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }
  ImageType::Pointer image = reader->GetOutput();

  /*********** curvature filter ******************/
  CurvatureFilter(image, FType, numberOfIterations);
  /*********** curvature filter ******************/

  // save the result
  typedef unsigned char                          WritePixelType;
  typedef itk::Image< WritePixelType, 2 >        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType > WriterType;

  typedef itk::RescaleIntensityImageFilter<
               ImageType, WriteImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput(image);

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(rescaler->GetOutput());
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  return EXIT_SUCCESS;
}

