#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

XTAL_TYP LIMIT;

template <unsigned int N=0, typename ...As>
XTAL_USE limit_t = inferred_t<unsigned int, LIMIT, bond::word<N>, As...>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
