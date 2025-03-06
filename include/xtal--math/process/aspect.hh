#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./dots.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Resolves the aspect indicated by the supplied type.

Providing related utility w.r.t. e.g. `signed`/`unsigned` i.o.w. sign/magnitude.
*/
template <typename ...As>
struct aspect;


////////////////////////////////////////////////////////////////////////////////

template <>
struct aspect<signed>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &o)
		noexcept -> auto
		requires un_n<numeric_variable_q<decltype(o)>>
		{
			return (0 < o) - (o < 0) + (o == 0);
		}
		/**/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(ordinal_variable_q auto o)
		noexcept -> auto
		{
			using _fit    = bond::fit<decltype(o)>;
			using _alpha = typename _fit::alpha_type;
			o  -= 1;
			o >>= _fit::positive.depth;
			//\
			if constexpr (_fit::IEC&559) {
			if constexpr (false) {
				o &= _fit::sign.mask;
				o |= _fit::unit.mask;
				return _xtd::bit_cast<_alpha>(o);
			}
			else {
				o |= 1;
				return    static_cast<_alpha>(o);
			}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(cardinal_variable_q auto o)
		noexcept -> auto
		{
			using _fit    = bond::fit<decltype(o)>;
			using _alpha = typename _fit::alpha_type;
			using _delta = typename _fit::delta_type;
			//\
			if constexpr (_fit::IEC&559) {
			if constexpr (false) {
				o <<= _fit::positive.depth;
				o  ^= _fit::sign.mask|_fit::unit.mask;
				return _xtd::bit_cast<_alpha>(o);
			}
			else {
				return method_f(_xtd::bit_cast<_delta>(o&one));
			}
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(real_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fit = bond::fit<decltype(o)>;
			return _xtd::copysign(_fit::alpha_1, o);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(complex_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fit = bond::fit<decltype(o)>;
			return o*root_f<-2>(dot_f(o));
		}

		template <int N_side=0> requires in_n<N_side, 1, 0, -1>
		XTAL_DEF_(return,inline,set)
		edit_f(ordinal_variable_q auto &u)
		noexcept -> auto
		{
			auto const v = bond::math::bit_sign_f(u);
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
			return method_f(v);
		}
		template <int N_side=0> requires in_n<N_side, 1, 0, -1>
		XTAL_DEF_(return,inline,set)
		edit_f(cardinal_variable_q auto &u)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(u)>;
			return edit_f<N_side>(reinterpret_cast<typename _fit::delta_type &>(u));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(real_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto const o_sgn = method_f<Ns...>(o); o *= o_sgn; return o_sgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(complex_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fit = bond::fit<decltype(o)>;
			auto [u, v] = dots_f<2>(o);
			auto const o_sgn = o*v;
			auto const o_mgn = XTAL_ALL_(o){u};
			o =    o_mgn;
			return o_sgn;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <>
struct aspect<unsigned>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(_std::unsigned_integral auto const &o)
		noexcept -> auto
		{
			return o;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(_std::  signed_integral auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto const v = o >> _fit::sign.shift;
			return (o^v) - v;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(real_variable_q auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			return _xtd::copysign(o, _fit::alpha_1);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(complex_variable_q auto const &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			return root_f<2>(dot_f(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(real_variable_q auto &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto const o_sgn = _xtd::copysign(_fit::alpha_1, o);
			auto const o_mgn = o*o_sgn;
			o =    o_sgn;
			return o_mgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(complex_variable_q auto &o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto [u, v] = dots_f<2>(o); o *= v; return u;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
using aspect_t = process::confined_t<aspect<As...>>;

template <typename ...As>
XTAL_DEF_(let)
aspect_f = [] XTAL_1FN_(call) (aspect_t<As...>::method_f);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
