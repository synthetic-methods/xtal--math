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
	using   Tanh = process::confined_t<dilated<1>, argument< 2>>;\
	using ArTanh = process::confined_t<dilated<1>, argument<-2>>;\

template <int M_ism=1, int M_car=0, typename ...As>
	requires in_n<M_ism, 1, 2, -1, -2> and in_n<M_car, -0, -1, -2>
XTAL_TYP argument
:	process::lift<argument<M_ism, M_car>, bond::compose<As...>>
{
};
template <int M_ism=1, typename ...As>
XTAL_USE argument_t = process::confined_t<argument<M_ism, bond::seek_constant_n<As..., nominal_t<0>>, As...>>;

template <int M_ism=1, typename ...As>
XTAL_DEF_(return,inline)
XTAL_LET argument_f(auto &&o)
XTAL_0EX -> decltype(auto)
{
	return argument_t<M_ism, As...>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct argument<M_ism, -0>
{
	using subkind = bond::compose<discarded<1, +1>, argument<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_number_q auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			using horner::term_f;
			using _op = bond::operate<decltype(o)>;
			using U_aphex = typename _op::aphex_type;
			using U_alpha = typename _op::alpha_type;
			using W_alpha = algebra::scalar_t<U_alpha[2]>;

			auto const &[o_re, o_im] = invalued_f(o);

			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
			XTAL_0IF (N_lim <  0) {
				auto o_re_abs = o_re, o_re_sgn = _op::design_f(o_re_abs);
				auto o_im_abs = o_im, o_im_sgn = _op::design_f(o_im_abs);
				W_alpha mode(o_im_abs < o_re_abs);

				W_alpha up{o_im, o_re_abs}; up *= mode;
				W_alpha dn{o_re_abs,-o_im}; dn *= mode;

				U_alpha u = up.sum()/dn.sum();

				U_alpha constexpr pi2 = 1.5707963267948966192313216916397514420985846996875529104874722961;
				U_alpha constexpr sq2 = 1.4142135623730950488016887242096980785696718753769480731766797379;
				U_alpha constexpr co_ = 0.8735141470922506811011384624884910162505195635622737479437808099;
				
				auto const x = term_f(_op::alpha_1, co_, u, u);
				auto const y = sq2*u*root_f<-2>(x + root_f<2>(x)) + mode[1]*pi2*o_re_sgn*o_im_sgn;
			//	auto const y = sq2*u*root_f<-2>(x + root_f<2>(x)) + n*pi2*_op::assigned_f(o_re);

			}
		}

	};
};
template <int M_ism>
struct argument<M_ism, -1>
:	bond::compose<discarded<1, +2>, argument<M_ism, -2>>
{
};
template <int M_ism>
struct argument<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&w)
		XTAL_0EX -> decltype(auto)
		{
			static_assert(N_lim <= 0);

			if constexpr (N_lim < 0) {
				auto const u = root_f<2>(XTAL_REF_(w));
				return argument<M_ism, -0>::template function<N_lim>(u)/(u);
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
