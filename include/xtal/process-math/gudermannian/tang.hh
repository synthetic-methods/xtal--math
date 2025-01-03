#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::gudermannian
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines a class of Gudermannian-related function/approximations indexed by `M_ism`: \
	 1 ->                 Tan[#*Pi]/Pi -> # * Sqrt[(1 - 2 #^2)]/(1 - 4 #^2)\
	 2 ->        Gudermannian[#*Pi]/Pi -> # * Sqrt[(1 +   #^2)]/(1 + 2 #^2)\
	-1 ->              ArcTan[#*Pi]/Pi -> # * Sqrt[2]/Sqrt[(1 + 8 #^2) + Sqrt[(1 + 8 #^2)]]\
	-2 -> InverseGudermannian[#*Pi]/Pi -> # * Sqrt[2]/Sqrt[(1 - 4 #^2) + Sqrt[(1 - 4 #^2)]]\

///\note\
The (co)domain is normalized around `+/- 1/2`, with derivative `1` at `0`. \

///\example\
	using   Tanh = process::confined_t<dilated<2>, tang< 2>>;\
	using ArTanh = process::confined_t<dilated<2>, tang<-2>>;\

template <int M_ism=1, int M_car=0, typename ...As>
	requires in_q<M_ism, 1, 2, -1, -2> and in_q<M_car, -0, -1, -2>
struct   tang
:	process::lift<tang<M_ism, M_car>, bond::compose<As...>>
{
};
template <int M_ism=1, typename ...As>
using    tang_t = process::confined_t<tang<M_ism, bond::seek_constant_n<As..., constant_t<0>>, As...>>;

template <int M_ism=1, typename ...As>
XTAL_DEF_(short)
XTAL_LET tang_f(auto &&o)
noexcept -> decltype(auto)
{
	return tang_t<M_ism, As...>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct tang<M_ism, -0>
{
	using superkind = bond::compose<discarded<1>, tang<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> decltype(auto)
		{
			static_assert(N_lim <= 0);

			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
			XTAL_0IF (N_lim <  0) {
				using _op = bond::operate<decltype(u)>;
				auto constexpr up =     _op::patio_1;
				auto constexpr dn = one/_op::patio_1;
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return atan(sinh(XTAL_REF_(u)*up))*dn;}// `Gudermannian`
				XTAL_0IF (M_ism ==  1) {return       tan(XTAL_REF_(u)*up) *dn;}
				XTAL_0IF (M_ism == -1) {return      atan(XTAL_REF_(u)*up) *dn;}
				XTAL_0IF (M_ism == -2) {return asinh(tan(XTAL_REF_(u)*up))*dn;}// `InverseGudermannian`
			//	Alternative formulations...
			//	`M_ism == -1`: `function[2]@Log[x + Sqrt[x^2 + 1]]`
			//	`M_ism == -2`: `ArSinh[Tan[x]] = Log[(1 + Sin[up x])/Cos[up x]]*dn`

			}
		}

	};
};
template <int M_ism>
struct tang<M_ism, -1>
:	bond::compose<discarded<2>, tang<M_ism, -2>>
{
};
template <int M_ism>
struct tang<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> decltype(auto)
		{
			static_assert(N_lim <= 0);

			if constexpr (N_lim < 0) {
				auto const u = root_f<2>(XTAL_REF_(w));
				return tang<M_ism, -0>::template function<N_lim>(u)/(u);
			}
			else {
				using _op = bond::operate<decltype(w)>;
				auto const _1 = _op::diplo_f(0);
				auto const _2 = _op::diplo_f(1);
				auto const _4 = _op::diplo_f(2);
				auto const _8 = _op::diplo_f(3);
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return root_f<2>(term_f(one,  _1, w))/term_f(one, _2, XTAL_REF_(w));}
				XTAL_0IF (M_ism ==  1) {return root_f<2>(term_f(one, -_2, w))/term_f(one,-_4, XTAL_REF_(w));}
				XTAL_0IF (M_ism == -1) {auto const m = term_f(one, _8, XTAL_REF_(w)); return root_f<-2>(_op::haplo_1*(root_f<2>(m) + m));}
				XTAL_0IF (M_ism == -2) {auto const m = term_f(one,-_4, XTAL_REF_(w)); return root_f<-2>(_op::haplo_1*(root_f<2>(m) + m));}
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
