#ifndef CUSTOM_SITESET_H
#define CUSTOM_SITESET_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

SiteSet
custom_spin_boson(const int N , const int max_occ);

SiteSet
custom_spin_boson_doubling(const int N , const int max_occ);
#endif
