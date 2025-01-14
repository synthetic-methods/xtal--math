#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./dots.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct signum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &o)
		noexcept -> auto
		requires un_n<anyplex_variable_q<decltype(o)>>
		{
			return (0 < o) - (o < 0) + (o == 0);
		}
		/**/
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(ordinal_variable_q auto o)
		noexcept -> auto
		{
			using _op    = bond::operate<decltype(o)>;
			using _alpha = typename _op::alpha_type;
			o  -= 1;
			o >>= _op::positive.depth;
			//\
			if constexpr (_op::IEC&559) {
			if constexpr (false) {
				o &= _op::sign.mask;
				o |= _op::unit.mask;
				return _xtd::bit_cast<_alpha>(o);
			}
			else {
				o |= 1;
				return    static_cast<_alpha>(o);
			}
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(cardinal_variable_q auto o)
		noexcept -> auto
		{
			using _op    = bond::operate<decltype(o)>;
			using _alpha = typename _op::alpha_type;
			using _delta = typename _op::delta_type;
			//\
			if constexpr (_op::IEC&559) {
			if constexpr (false) {
				o <<= _op::positive.depth;
				o  ^= _op::sign.mask|_op::unit.mask;
				return _xtd::bit_cast<_alpha>(o);
			}
			else {
				return function(_xtd::bit_cast<_delta>(o&one));
			}
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return _xtd::copysign(_op::alpha_1, o);
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return o*dot_f<-2>(o);
		}

		template <int N_side=0> requires in_n<N_side, 1, 0, -1>
		XTAL_DEF_(long)
		XTAL_SET edit(ordinal_variable_q auto &u)
		noexcept -> auto
		{
			auto const v = bond::bit_sign_f(u);
			XTAL_IF0
			XTAL_0IF (N_side == -1) {
				u &= v;
				u += v;
				u ^= v;
			}
			XTAL_0IF (N_side ==  0) {
				u += v;
				u ^= v;
			}
			XTAL_0IF (N_side ==  1) {
				u &=~v;
			}
			return function(v);
		}
		template <int N_side=0> requires in_n<N_side, 1, 0, -1>
		XTAL_DEF_(long)
		XTAL_SET edit(cardinal_variable_q auto &u)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u)>;
			return edit<N_side>(reinterpret_cast<typename _op::delta_type &>(u));
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(real_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto const o_sgn = function<Ns...>(o); o *= o_sgn; return o_sgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET edit(complex_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			auto [u, v] = dots_f<2>(o);
			auto const o_sgn = o*v;
			auto const o_mgn = XTAL_ALL_(o){u};
			o =    o_mgn;
			return o_sgn;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using    signum_t = process::confined_t<signum<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET signum_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return signum_t<Ms...>::function(XTAL_REF_(oo)...);
	return signum_t<>::template function<Ns...>(XTAL_REF_(oo)...);
}
template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET signum_e(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return signum_t<Ms...>::function(XTAL_REF_(oo)...);
	return signum_t<>::template     edit<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
