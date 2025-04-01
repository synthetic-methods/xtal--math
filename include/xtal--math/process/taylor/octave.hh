#pragma once
#include "./any.hh"

#include "./logarithm.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Produces the base-two logarithm and its inverse.
\todo    Either define `subdivision` to apply `1/N_ind` prescaling, or supply an `attach`ed parameter to do so.
\todo    Allow for recentering w.r.t. `A4`, e.g. `440*Exp2[-69/12]*Exp2[#/12]&`,
         or `440*0.1858136117191751604527105712350021e-1L*octave_f<-2>(n/m_div)`.
*/
template <int M_ism=2, int M_div=1>
struct  octave;


////////////////////////////////////////////////////////////////////////////////

template <int M_div>
struct octave<+2, M_div>
:	process::lift<void
	,	dilate<[] XTAL_1FN_(to) (one*_detail::base_f<M_div>(2.L))>
	,	logarithm<+1, 1>
	>
{
};
template <int M_div>
struct octave<-2, M_div>
:	process::lift<void
	,	logarithm<-1, 1>
	,	dilate<[] XTAL_1FN_(to) (one/_detail::base_f<M_div>(2.L))>
	>
{
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=2, int M_div=1>
using octave_t = process::confined_t<octave<M_ism, M_div>>;

template <int M_ism=2, int M_div=1, int ...Ns>
XTAL_DEF_(let)
octave_f = [] XTAL_1FN_(call) (octave_t<M_ism, M_div>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
