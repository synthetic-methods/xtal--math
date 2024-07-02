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
	using   Tanh = process::confined_t<dilating<1>, sigmoid< 2>>;\
	using ArTanh = process::confined_t<dilating<1>, sigmoid<-2>>;\

template <int M_ism=1, int M_car=0, typename ...As>
	requires inclusive_q<M_ism, 1, 2, -1, -2> and inclusive_q<M_car, -0, -1, -2>
XTAL_TYP sigmoid
:	process::lift<sigmoid<M_ism, M_car>, bond::compose<As...>>
{
};
template <int M_ism=1, typename ...As>
XTAL_USE sigmoid_t = process::confined_t<sigmoid<M_ism, bond::seek_constant_n<As..., nominal_t<0>>, As...>>;

template <int M_ism=1, typename ...As>
XTAL_DEF_(return,inline)
XTAL_RET sigmoid_f(auto &&o)
XTAL_0EX
{
	return sigmoid_t<M_ism, As...>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct sigmoid<M_ism, -0>
{
	using subkind = bond::compose<discarding<1, +1>, sigmoid<M_ism, -1>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&u)
		XTAL_0EX -> decltype(auto)
		{
			static_assert(N_lim <= 0);

			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
			XTAL_0IF (N_lim <  0) {
				using _op = bond::operate<decltype(u)>;
				auto constexpr _1 =    _op::alpha_1;
				auto constexpr up =    _op::patio_1;
				auto constexpr dn = _1/_op::patio_1;
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
struct sigmoid<M_ism, -1>
:	bond::compose<discarding<1, +2>, sigmoid<M_ism, -2>>
{
};
template <int M_ism>
struct sigmoid<M_ism, -2>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&w)
		XTAL_0EX -> decltype(auto)
		{
			static_assert(N_lim <= 0);

			if constexpr (N_lim < 0) {
				auto const u = root_f<2>(XTAL_REF_(w));
				return sigmoid<M_ism, -0>::template function<N_lim>(u)/(u);
			}
			else {
				using _op = bond::operate<decltype(w)>;
				auto const _1 = _op::diplo_f(0);
				auto const _2 = _op::diplo_f(1);
				auto const _4 = _op::diplo_f(2);
				auto const _8 = _op::diplo_f(3);
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
