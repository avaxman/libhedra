// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_SCALAR2RGB_H
#define HEDRA_SCALAR2RGB_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 


namespace hedra
{
  // computes a cool-warm set of RGB [0.0,1.0] colors. The presented values are normalized to [minValue, maxValue] and values beyond this range are truncated.
  
  //Inputs and outputs are then straightforward
  
  //Assumption: minValue < maxValue
  IGL_INLINE bool scalar2RGB(const Eigen::VectorXd& scalar,
                             const double minValue,
                             const double maxValue,
                             Eigen::MatrixXd& C)
  {
    using namespace Eigen;
    
    MatrixXd ColorMap(32, 3);
    ColorMap<<85,72,193,
    94,87,205,
    103,101,217,
    112,115,227,
    121,129,236,
    130,142,243,
    139,154,248,
    149,166,252,
    158,177,255,
    167,186,255,
    176,195,254,
    185,203,252,
    194,209,248,
    202,214,242,
    210,218,234,
    217,220,226,
    225,219,215,
    231,214,203,
    236,207,190,
    240,200,178,
    242,191,165,
    243,181,152,
    243,170,139,
    241,157,126,
    237,144,114,
    232,130,101,
    226,115,90,
    219,99,78,
    210,81,68,
    200,62,58,
    189,40,48,
    177,1,39;
    
    ColorMap/=255.0;
    
    C.conservativeResize(scalar.size(),3);
    for (int i=0;i<scalar.size();i++){
      double CurrValue=31.0*(scalar(i)-minValue)/(maxValue-minValue);
      /*std::cout<<"CurrValue: "<<CurrValue<<std::endl;
      std::cout<<"scalar(i): "<<scalar(i)<<std::endl;
      std::cout<<"minValue: "<<minValue<<std::endl;
      std::cout<<"maxValue: "<<maxValue<<std::endl;*/
      if (CurrValue<=0.0) C.row(i)<<ColorMap.row(0);
      else
        if (CurrValue>=31.0) C.row(i)<<ColorMap.row(31);
        else{
          int Entry=(int)floor(CurrValue);
          double Residue=CurrValue-(double)Entry;
          C.row(i)=ColorMap.row(Entry)*(1.0-Residue)+ColorMap.row(Entry+1)*Residue;
        }
      
    }
    return true;
  }
}




#endif


