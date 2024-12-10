#pragma once
#include "./any.hh"

#include "./root.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Implements the function `#^m^Sgn@m &`. \


////////////////////////////////////////////////////////////////////////////////

template <int M_exp=1, auto ...Ms>
struct   monomial;

template <auto ...Ms>
struct   monomial<0, Ms...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			return XTAL_ALL_(o){1};
		}

	};
};
template <int M_exp, auto ...Ms> requires (0 < M_exp) and in_n<bond::operating::bit_count_f(M_exp), 1>
struct   monomial<M_exp, Ms...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto o)
		noexcept -> auto
		{
			#pragma unroll
			for (int i{M_exp}; i >>= 1; o *= o); return o;
		}

	};
};
template <int M_exp, auto ...Ms> requires (2 < M_exp) and un_n<bond::operating::bit_count_f(M_exp), 1>
struct   monomial<M_exp, Ms...>
{
	XTAL_SET M_par = M_exp&1;
	XTAL_SET M_per = M_exp >> 1;

	using superkind = monomial<M_per, Ms...>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &o)
		noexcept -> auto
		{
			XTAL_LET oo_f = [] (auto const &u) XTAL_0FN_(u*u);
			XTAL_IF0
			XTAL_0IF (1 == M_par) {return o*oo_f(S_::template function<Ns...>(o));}
			XTAL_0IF (0 == M_par) {return   oo_f(S_::template function<Ns...>(o));}
		}

	};
};
template <int M_exp, auto ...Ms> requires (M_exp < 0) and complete_q<root<M_exp, Ms...>>
struct   monomial<M_exp, Ms...>
{
	using superkind = root<-M_exp, Ms...>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &o)
		noexcept -> auto
		{
			return S_::template function<Ns...>(o);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    monomial_t = process::confined_t<monomial<Ms...>>;

template <int ...Ms>
XTAL_DEF_(short)
XTAL_LET monomial_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return monomial_t<Ms...>::function(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
