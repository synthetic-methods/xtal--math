#pragma once
#include "./any.hh"

#include "./logarithm.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Produces the base-two logarithm and its inverse. \

///\todo\
Either define `subdivision` to apply `1/N_ind` prescaling, \
or supply an `attach`ed parameter to do so...


template <int M_ism=2, int M_car=0>
struct  octave;


////////////////////////////////////////////////////////////////////////////////

template <int M_car>
struct octave<+2, M_car>
:	process::lift<void
	,	dilate<[] XTAL_1FN_(to) (one*logarithm_f<1>(bond::fit<>::alpha_f(2)))>
	,	logarithm<+1, M_car>
	>
{
};
template <int M_car>
struct octave<-2, M_car>
:	process::lift<void
	,	logarithm<-1, M_car>
	,	dilate<[] XTAL_1FN_(to) (one/logarithm_f<1>(bond::fit<>::alpha_f(2)))>
	>
{
};

template <int M_car>
struct octave<+1, M_car>
:	process::lift<void
	,	dilate<bond::operate{[] XTAL_1FN_(to)
			(one*bond::fit<>::patio_f(2))}>
	,	imagine<-1>
	,	logarithm<+1, M_car>
	>
{
};
template <int M_car>
struct octave<-1, M_car>
:	process::lift<void
	,	logarithm<-1, M_car>
	,	imagine<+1>
	,	dilate<bond::operate{[] XTAL_1FN_(to)
			(one/bond::fit<>::patio_f(2))}>
	>
{
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=2, int M_car=0>
using octave_t = process::confined_t<octave<M_ism, M_car>>;

template <int M_ism=2, int M_car=0, int ...Ns>
XTAL_DEF_(let)
octave_f = [] XTAL_1FN_(call) (octave_t<M_ism, M_car>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
