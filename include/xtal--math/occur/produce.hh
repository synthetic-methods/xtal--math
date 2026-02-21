#pragma once
#include "./any.hh"

#include "../atom/dot.hh"
#include "./indent.hh"



XTAL_ENV_(push)
namespace xtal::occur::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides in-situ dot-production for simple messages.
*/
template <         typename ..._s>	struct  produce;
template <class S, typename ..._s>	using   produce_s = bond::compose_s<S, produce<_s...>>;
template <         typename ..._s>	concept produce_q = bond::tag_in_p<produce, _s...> or
	indent_q<_s...> and bond::tag_in_p<produce, typename _s::tail_type...>;


////////////////////////////////////////////////////////////////////////////////

template <typename ..._s>
struct produce
{

public:
	template <class S         > using innertype = XTAL_ALL_(XTAL_ANY_(S).product());
	template <class S, class T> using innerkind = _detail::navigate<T, innertype<S>>;
	using superkind = bond::compose<bond::tag<produce>, _s...>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind, innerkind<S, subtype<S>>>
	{
		//\
		static_assert(atom::math::dot_q<S> and S::rank() == 2);
		static_assert(atom::math::dot_q<S>);
		using S_ = bond::compose_s<S, superkind, innerkind<S, subtype<S>>>;
		using T_ = typename S_::self_type;
		using U_ = typename S_::tail_type;
	//	using V_ = typename U_::value_type;
	//	using W_ = produce_s<U_, _s...>;

		XTAL_DEF_(set) valve = [] (auto &&o)
			XTAL_0FN_(to) (XTAL_REF_(o).product());

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS

		XTAL_FN0_(go) (XTAL_DEF_(return,inline,get) basis, [] (auto &&o)
		XTAL_0FN_(to) (static_cast<S_&&>(XTAL_REF_(o)).tail().template element<0>()))

		XTAL_FN0_(go) (XTAL_DEF_(return,inline,get) point, [] (auto &&o)
		XTAL_0FN_(to) (static_cast<S_&&>(XTAL_REF_(o)).tail().template element<1>()))

		XTAL_FN0_(go) (XTAL_DEF_(return,inline,let)              value, valve)
		XTAL_FN0_(go) (XTAL_DEF_(return,inline,let)               head, valve)
		XTAL_FN0_(go) (XTAL_DEF_(return,inline,implicit) operator auto, valve)

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
