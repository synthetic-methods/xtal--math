#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

XTAL_TYP SHAPE;

template <unsigned int N=0, typename ...As>
XTAL_USE shape_t = inferred_t<unsigned int, SHAPE, bond::word<N>, As...>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
