// Copyright (C) 2008 Wim Meeussen <meeussen at willowgarage com>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//  
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  


#ifndef __NON_LINEAR_SYSTEM_CONDITIONAL_GAUSSIAN_ODO__
#define __NON_LINEAR_SYSTEM_CONDITIONAL_GAUSSIAN_ODO__

#include <bfl/pdf/analyticconditionalgaussian_additivenoise.h>

/*
bfl用到的几种数据类型
#define MyMatrix          MatrixWrapper::Matrix
#define MyColumnVector    MatrixWrapper::ColumnVector
#define MyRowVector       MatrixWrapper::RowVector
#define MySymmetricMatrix MatrixWrapper::SymmetricMatrix


#### MatrixWrapper::ColumnVector ####
	  ColumnVector();
	  ColumnVector(int nrows);

#### MatrixWrapper::RowVector ####
	  RowVector();
	  RowVector(int ncols);
	  
#### MatrixWrapper::Matrix ####
      Matrix();
      Matrix(int m, int n);

#### MatrixWrapper::SymmetricMatrix ####
      SymmetricMatrix();
      SymmetricMatrix(int n);
	  
包含了矩阵的基本操作，例如inverse，transpose,determinant,矩阵的运算
以及矩阵的求解SVD/cholesky_semidefinite（乔里斯基半正定分解）等
*/



namespace BFL
{
  /// Non Linear Conditional Gaussian
  /**
     - \f$ \mu = Matrix[1] . ConditionalArguments[0] +
     Matrix[2]. ConditionalArguments[1]  + ... + Noise.\mu \f$
     - Covariance is independent of the ConditionalArguments, and is
     the covariance of the Noise pdf
  */
  class NonLinearAnalyticConditionalGaussianOdo : public AnalyticConditionalGaussianAdditiveNoise
  {
    public:
      /// Constructor
      /** @pre:  Every Matrix should have the same amount of rows!
	  This is currently not checked.  The same goes for the number
	  of columns, which should be equal to the number of rows of
	  the corresponding conditional argument!
	  @param additiveNoise Pdf representing the additive Gaussian uncertainty
      */
		
	 //Additive :because it is added to any noise that might be intrinsic to the information system.
      NonLinearAnalyticConditionalGaussianOdo( const Gaussian& additiveNoise);

      /// Destructor
      virtual ~NonLinearAnalyticConditionalGaussianOdo();

      // redefine virtual functions
	  //Get the expected value E[x] of the pdf.
      virtual MatrixWrapper::ColumnVector    ExpectedValueGet() const;
	  
	  //returns derivative from function to n-th conditional variable
      virtual MatrixWrapper::Matrix          dfGet(unsigned int i)       const;

    private:
	  //mutable 只能用来修饰类的数据成员；而被 mutable 修饰的数据成员，可以在 const成员函数中修改
      mutable MatrixWrapper::Matrix df;
    };

} // End namespace BFL
 
#endif //  
