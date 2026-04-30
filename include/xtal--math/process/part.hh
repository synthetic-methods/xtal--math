#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./root.hh"
#include "./roots.hh"
#include "../atom/phason.hh"

XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Resolves the value indicated by the supplied type,
         e.g. `signed`/`unsigned` resolves the sign/magnitude, respectively.

\note    When applied to floating-point values,
         the `part<signed>` excludes zero.

\todo    Provide for smoothing via `(Sqrt[#1^2 + #2^2]&)`.
\todo    Include projections e.g. `real` and `imag`.

*/
template <typename ...As>
struct part;


////////////////////////////////////////////////////////////////////////////////

template <>
struct part<signed>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(integral_constant_q auto const &o)
		const noexcept -> auto
		{
		//	return (0 < o) - (o < 0) + (o == 0);
			return constant_t<method(XTAL_ALL_(o)::value)>{};
		}
		/**/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(integral_variable_q auto o)
		const noexcept -> auto
		{
			using U       = XTAL_ALL_(o);
			using U_fit   = bond::fit<U>;
			using U_alpha = typename U_fit::alpha_type;
			using U_delta = typename U_fit::delta_type;
			using U_sigma = typename U_fit::sigma_type;
			XTAL_IF0
			XTAL_0IF ( ordinal_variable_q<U>) {
				o  -= 1;
				o >>= U_fit::positive.depth;
				//\
				if constexpr (U_fit::IEC&559) {
				if constexpr (false) {
					o &= U_fit::sign.mask;
					o |= U_fit::unit.mask;
					return _xtd::bit_cast<U_alpha>(o);
				}
				else {
					o |= 1;
					return    static_cast<U_alpha>(o);
				}
			}
			XTAL_0IF (cardinal_variable_q<U>) {
				//\
				if constexpr (U_fit::IEC&559) {
				if constexpr (false) {
					o <<= U_fit::positive.depth;
					o  ^= U_fit::sign.mask|U_fit::unit.mask;
					return _xtd::bit_cast<U_alpha>(o);
				}
				else {
					return method(_xtd::bit_cast<U_delta>(o&one));
				}
			}
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(real_variable_q auto const &o)
		const noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			return _xtd::copysign(U_fit::alpha_1, o);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(complex_variable_q auto const &o)
		const noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			return o*root_f<-2>(dot_f(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::groupoid_q auto &&o)
		const noexcept -> decltype(auto)
		{
			using  U = XTAL_ALL_(o);
			return U::template zip_from<[] XTAL_1FN_(call) (subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(_std::in_place_t, auto &&...oo)
		const noexcept -> decltype(auto)
		{
			return edit(XTAL_REF_(oo)...);
		}

		template <int N_side=0> requires in_v<N_side, 1, 0, -1>
		XTAL_DEF_(return,inline,let)
		edit(ordinal_variable_q auto &u)
		const noexcept -> auto
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
			return method(v);
		}
		template <int N_side=0> requires in_v<N_side, 1, 0, -1>
		XTAL_DEF_(return,inline,let)
		edit(cardinal_variable_q auto &u)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(u)>;
			return edit<N_side>(reinterpret_cast<typename U_fit::delta_type &>(u));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		edit(real_variable_q auto &o)
		const noexcept -> XTAL_ALL_(o)
		{
			auto const o_sgn = method<Ns...>(o); o *= o_sgn; return o_sgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		edit(complex_variable_q auto &o)
		const noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			auto [u, v] = roots_f<2>(dot_f(o));
			auto const o_sgn = o*v;
			auto const o_mgn = XTAL_ALL_(o){u};
			o =    o_mgn;
			return o_sgn;
		}

	};
};
template <>
struct part<unsigned>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(cardinal_variable_q auto const &o)
		const noexcept -> auto
		{
			return o;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method( ordinal_variable_q auto const &o)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto const v = o >> U_fit::sign.shift;
			return (o^v) - v;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(real_variable_q auto const &o)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			return _xtd::copysign(o, U_fit::alpha_1);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(complex_variable_q auto const &o)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			return root_f<2>(dot_f(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::groupoid_q auto &&o)
		const noexcept -> decltype(auto)
		{
			using  U = XTAL_ALL_(o);
			return U::template zip_from<[] XTAL_1FN_(call) (subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		edit(real_variable_q auto &o)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto const o_sgn = _xtd::copysign(U_fit::alpha_1, o);
			auto const o_mgn = o*o_sgn;
			o =    o_sgn;
			return o_mgn;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		edit(complex_variable_q auto &o)
		const noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto [u, v] = roots_f<2>(dot_f(o)); o *= v; return u;
		}

	};
};
template <>
struct part<_xtd::real<>>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		const noexcept -> auto
		{
			return _std::real(XTAL_REF_(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::math::phason_q auto &&t_)
		const noexcept -> auto
		{
			return method<Ns...>(t_(0));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
using part_t = process::confined_t<part<As...>>;

template <typename ...As>
XTAL_DEF_(let)
part_f = [] XTAL_1FN_(call) (part_t<As...>{}.method);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
