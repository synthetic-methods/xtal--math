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
		XTAL_DEF_(return,inline,set)
		method(integral_constant_q auto const &o)
		noexcept -> auto
		{
		//	return (0 < o) - (o < 0) + (o == 0);
			return constant_t<method(XTAL_ALL_(o)::value)>{};
		}
		/**/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(integral_variable_q auto o)
		noexcept -> auto
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
					return xtd::bit_cast<U_alpha>(o);
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
					return xtd::bit_cast<U_alpha>(o);
				}
				else {
					return method(xtd::bit_cast<U_delta>(o&one));
				}
			}
		}
		/***/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(real_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			return xtd::copysign(U_fit::alpha_1, o);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(complex_variable_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			return o*root_f<-2>(dot_f(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(atom::quantify_q auto &&o)
		noexcept -> decltype(auto)
		{
			using  U = XTAL_ALL_(o);
			return U::template zip_from<[] XTAL_1FN_(call)
				(subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}

		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(ordinal_variable_q auto &u, integral_constant_q auto N_dir)
		noexcept -> auto
		{
		//	static_assert(in_v<N_dir, 1, 0, -1>);
			auto const v = bond::math::bit_sign_f(u);
			XTAL_IF0
			XTAL_0IF (N_dir == -1) {
				u &= v;
				u += v;
				u ^= v;
			}
			XTAL_0IF (N_dir ==  0) {
				u += v;
				u ^= v;
			}
			XTAL_0IF (N_dir ==  1) {
				u &=~v;
			}
			return method(v);
		}
		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(cardinal_variable_q auto &u, integral_constant_q auto N_dir)
		noexcept -> auto
		{
			using U = XTAL_ALL_(u);
			using V = xtd::signed_case_t<U>;
			return method<std::in_place>(reinterpret_cast<V &>(u), N_dir);
		}
		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(integral_variable_q auto &u)
		noexcept -> auto
		{
			return method<std::in_place>(u, constant_t<0>{});
		}
		//\
		template <exvariable_q auto N_dir>
		template <auto N_dir> requires requires {exvariable_f(N_dir);}// Alt. for LLVM 16
		XTAL_DEF_(return,inline,set)
		method(integral_variable_q auto &u)
		noexcept -> auto
		{
			return method<std::in_place>(u, exvariable_f(N_dir));
		}

		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(real_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto const o_sgn = method(o); o *= o_sgn; return o_sgn;
		}
		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(complex_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
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
		XTAL_DEF_(return,inline,set)
		method(cardinal_variable_q auto const &o)
		noexcept -> auto
		{
			return o;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method( ordinal_variable_q auto const &o)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto const v = o >> U_fit::sign.shift;
			return (o^v) - v;
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(real_variable_q auto const &o)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			return xtd::copysign(o, U_fit::alpha_1);
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(complex_variable_q auto const &o)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			return root_f<2>(dot_f(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(atom::quantify_q auto &&o)
		noexcept -> decltype(auto)
		{
			using  U = XTAL_ALL_(o);
			return U::template zip_from<[] XTAL_1FN_(call)
				(subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}

		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(real_variable_q auto &o)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto const o_sgn = xtd::copysign(U_fit::alpha_1, o);
			auto const o_mgn = o*o_sgn;
			o =    o_sgn;
			return o_mgn;
		}
		template <std::in_place_t>
		XTAL_DEF_(return,inline,set)
		method(complex_variable_q auto &o)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(o)>;
			auto [u, v] = roots_f<2>(dot_f(o)); o *= v; return u;
		}

	};
};
template <>
struct part<xtd::real<>>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> auto
		{
			return std::real(XTAL_REF_(o));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(atom::math::phason_q auto &&t_)
		noexcept -> auto
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
part_f = [] XTAL_1FN_(call) (part_t<As...>::method);

template <typename ...As>
XTAL_DEF_(let)
part_e = [] XTAL_1FN_(call) (part_t<As...>::template method<std::in_place>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
