#pragma once
#include "./any.hh"

#include "../dilating.hh"
#include "../square.hh"
#include "../root.hh"


XTAL_ENV_(push)
namespace xtal::process::math::gudermannian
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0> XTAL_TYP tangent {static_assert(M_ism);};
template <int M_ism=1, int M_car=0> XTAL_USE tangent_t = process::confined_t<tangent<M_ism, M_car>>;
template <int M_ism=1, int M_car=0>
XTAL_FN2 tangent_f(auto &&o)
XTAL_0EX
{
	return tangent_t<M_ism, M_car>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
///\
Defines a class of square-root-based-approximations of Gudermannian-related functions, \
indexed by `M_ism`: \
\
	-2 -> InverseGudermannian[#*Pi]/Pi&\
	 2 ->        Gudermannian[#*Pi]/Pi&\
	-1 ->              ArcTan[#*Pi]/Pi&\
	 1 ->                 Tan[#*Pi]/Pi&\

///\note\
The indexing is intended to reflect that used by the trigonometric functions, \
where `+/- 1` and `+/- 2` respectively indicate circular and hyperbolic evaluations. \

template <int M_ism>
struct tangent<M_ism, -0>
{
	using subkind = bond::compose<discarding<1>, tangent<M_ism, -1>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&u)
		XTAL_0EX
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
			XTAL_0IF (N_lim <  0) {
				using _op = bond::operate<decltype(u)>;
				auto const up =   _op::patio_1;
				auto const dn = 1/_op::patio_1;
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return atan(sinh(XTAL_REF_(u)*up))*dn;}// `Gudermannian`
				XTAL_0IF (M_ism ==  1) {return      (tan(XTAL_REF_(u)*up))*dn;}
				XTAL_0IF (M_ism == -1) {return     (atan(XTAL_REF_(u)*up))*dn;}
				XTAL_0IF (M_ism == -2) {return asinh(tan(XTAL_REF_(u)*up))*dn;}// `InverseGudermannian`
			}
		}

	};
};
template <int M_ism>
struct tangent<M_ism, -1>
:	bond::compose<discarding<2>, tangent<M_ism, -2>>
{
};
template <int M_ism>
struct tangent<M_ism, -2>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&w)
		XTAL_0EX
		{
			if constexpr (N_lim < 0) {
				auto const u = root_f<2>(XTAL_REF_(w));
				return tangent<M_ism, -0>::template function<N_lim>(u)/(u);
			}
			else {
				using _op = bond::operate<decltype(w)>;
				auto const _1 = _op::alpha_1*1;
				auto const _2 = _op::alpha_1*2;
				auto const _4 = _op::alpha_1*4;
				auto const _8 = _op::alpha_1*8;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return root_f<2>(horner::term_f<  >(_1,     w))/horner::term_f<  >(_1, _2, XTAL_REF_(w));}
				XTAL_0IF (M_ism ==  1) {return root_f<2>(horner::term_f<-1>(_1, _2, w))/horner::term_f<-1>(_1, _4, XTAL_REF_(w));}
				XTAL_0IF (M_ism == -1) {auto const m = horner::term_f<  >(_1, _8, XTAL_REF_(w)); return _1/root_f<2>(_op::haplo_1*(root_f<2>(m) + m));}
				XTAL_0IF (M_ism == -2) {auto const m = horner::term_f<-1>(_1, _4, XTAL_REF_(w)); return _1/root_f<2>(_op::haplo_1*(root_f<2>(m) + m));}
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
