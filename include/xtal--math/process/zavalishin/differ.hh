#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Differentiates the incoming signal. \
If the modulator and/or phasor are supplied, \
the result is renormalized w.r.t. the chain-rule of differentiation. \

///\note\
The `method` uses local arena-like allocation to manage state. \
Unless `As...` includes, `stowed` the maximum `sizeof(input) + sizeof(modulator) <= 64`. \

///\todo\
Enable logarithmic differentiation by including the original signal with the normalization factor? \

template <typename ...As> struct   differ;
template <typename ...As> using    differ_t = process::confined_t<differ<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, int N_pole>
struct differ<U_pole[N_pole]>
{
//	using sample_type = occur::sample_t<>;
	using  order_type = occur::inferred_t<struct ORDER, unsigned int, bond::seek_s<N_pole + 1>>;

	using superkind = bond::compose<bond::tag<differ_t>
	,	provision::stowed<U_pole[N_pole << 1]>
//	,	typename sample_type::template   attach<>
	,	typename  order_type::template dispatch<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <int N_ord=1> requires un_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto &&...oo)
		noexcept -> auto
		{
			return u;
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u)
		noexcept -> auto
		{
			auto [u_] = S_::stow(u);
			return (u - u_);
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, atom::math::phason_q auto const &t_)
		noexcept -> auto
		{
			auto [u_] = S_::stow(u);
			return (u - u_)*root_f<-1>(t_(1));
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto const &v, atom::math::phason_q auto const &t_)
		noexcept -> auto
		{
			auto [u_, v_] = S_::stow(u, v);
			return (u - u_)*root_f<-1>(t_(1) + v - v_);
		}

	};
};
template <>
struct differ<>
:	differ<typename bond::fixture<>::aphex_type[2]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
