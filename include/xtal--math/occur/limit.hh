#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

struct  LIMIT;

template <unsigned int N=0, typename ...As>
using   limit_t = inferred_t<unsigned int, LIMIT, bond::word<N>, As...>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
