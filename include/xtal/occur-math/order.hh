#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

struct   ORDER;

template <unsigned int N=0, typename ...As>
using    order_t = inferred_t<unsigned int, ORDER, bond::word<N>, As...>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
