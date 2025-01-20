#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Implements the function `#^m &`. \

template <int M_exp=1, auto ...Ms>
struct   power;

template <int M_exp=1, auto ...Ms>
using    power_t = process::confined_t<power<M_exp, Ms...>>;

template <int M_exp=1, int ...Ns>
XTAL_DEF_(short)
XTAL_LET power_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return power_t<M_exp>::template static_method<Ns...>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////

template <int M_exp, auto ...Ms> requires un_n<M_exp>
struct   power<M_exp, Ms...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	public:

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &o)
		noexcept -> auto
		{
			return XTAL_ALL_(o) {one};
		}

	};
};
template <int M_exp, auto ...Ms> requires above_p<0, M_exp>
struct   power<M_exp, Ms...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	public:

		template <int ...Ns> requires in_n<0, bond::fixture<>::template expound_f<2>(M_exp)>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fix = bond::fixture<decltype(o)>;
			if constexpr (M_exp == 2) {
				if constexpr (complex_variable_q<decltype(o)>) {
					auto const &[x, y] = destruct_f(o);
					auto const      u  = _xtd::accumulator(x*x, y*y, _fix::alpha_f(-1));
					auto const      v  = two*x*y;
					return {u, v};
				}
				else {
					return o*o;
				}
			}
			else {
				auto oo = objective_f(o);
				#pragma unroll
				for (int i{M_exp}; i >>= 1; oo *= power_f<1>(XTAL_MOV_(oo))); return oo;
			}
		}

		template <int ...Ns> requires in_n<0, bond::fixture<>::template expound_f<3>(M_exp)>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fix = bond::fixture<decltype(o)>;
			if constexpr (M_exp == 3) {
				if constexpr (complex_variable_q<decltype(o)>) {
					auto const &[x, y] = destruct_f(o);
					auto const      xx =  x*x;
					auto const     _yy = -y*y;
					auto const      u  =  x*_xtd::accumulator( xx, _yy, _fix::alpha_f(3));
					auto const      v  =  y*_xtd::accumulator(_yy,  xx, _fix::alpha_f(3));
					return {u, v};
				}
				else {
					return o*o*o;
				}
			}
			else {
				auto oo = objective_f(o);
				#pragma unroll
				for (int i{M_exp}; i  /= 3; oo *= power_f<2>(XTAL_MOV_(oo))); return oo;
			}
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto const &o)
		noexcept -> auto
		requires un_n<0, bond::fixture<int>::template expound_f<2>(M_exp)>
		and      un_n<0, bond::fixture<int>::template expound_f<3>(M_exp)>
		{
			using O = XTAL_ALL_(o);
			XTAL_IF0
			XTAL_0IF (arrange::collate_q<O>) {
				return O::template map_f<XTAL_FUN_(static_method<Ns...>)>(o);
			}
			XTAL_0IF (M_exp   == 0) {return                             O{1};}
			XTAL_0IF (M_exp   == 1) {return                             o   ;}
			XTAL_0IF (M_exp%3 == 0) {return power_f<M_exp/3>(power_f<3>(o)) ;}
			XTAL_0IF (M_exp%2 == 0) {return power_f<M_exp/2>(power_f<2>(o)) ;}
			XTAL_0IF_(else)         {return        o*power_f<M_exp - 1>(o)  ;}
		}

	};
};
template <int M_exp, auto ...Ms> requires below_p<0, M_exp>
struct   power<M_exp, Ms...>
{
	using superkind = power<-M_exp, Ms...>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

	public:

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&o)
		noexcept -> auto
		{
			return S_::template static_method<Ns...>(one/XTAL_REF_(o));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
