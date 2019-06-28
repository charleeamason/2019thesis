// PolyOblModelNHQS.cpp
//
// Neutron and Hybrid Quark Star oblateness model
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#include "PolyOblModelNHQS.h"
#include <iostream>
PolyOblModelNHQS::PolyOblModelNHQS( const double& Rspot_nounits, const double& Req_nounits, const double& zeta, const double& eps )
  : PolyOblModelBase(Rspot_nounits, Req_nounits, zeta, eps) { }

double PolyOblModelNHQS::a0() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  //return double(-0.18*eps + 0.23*zeta*eps - 0.05*eps*eps);
   return double( 0.5*a2());

}

double PolyOblModelNHQS::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  //return double(-0.39*eps + 0.29*zeta*eps + 0.13*eps*eps);
  //std::cout<<" eps = "<<eps<<" zeta = "<<zeta<<std::endl;
  return double(2.0/3.0 * eps*(-0.788 + 1.030*zeta));

}

double PolyOblModelNHQS::a4() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  //  return double(0.04*eps - 0.15*zeta*eps + 0.07*eps*eps );

  //return double( -8.0/3.0*(-0.18*eps + 0.23*zeta*eps - 0.05*eps*eps - 0.5*(-0.39*eps + 0.29*zeta*eps + 0.13*eps*eps)));
    return double(0.0);


}
