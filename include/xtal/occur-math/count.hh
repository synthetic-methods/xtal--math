#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

XTAL_TYP COUNT;

template <unsigned int N=0, typename ...As>
XTAL_USE count_t = inferred_t<unsigned int, COUNT, bond::word<N>, As...>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
