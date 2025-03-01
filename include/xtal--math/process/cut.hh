#pragma once
#include "./any.hh"
#include "./nearing.hh"
#include "./magnum.hh"
#include "./signum.hh"
#include "./dots.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


namespace _detail
{///////////////////////////////////////////////////////////////////////////////

template <auto N>
auto constexpr flank_n = (int) signum_f(signum_f(N())*nearing_f(magnum_f(N())));


}///////////////////////////////////////////////////////////////////////////////

template <auto M_stop_, int M_side=_detail::flank_n<M_stop_>> requires in_n<M_side, 1, -1>
struct  cut;

template <auto M_stop_, int M_side=_detail::flank_n<M_stop_>> requires in_n<M_side, 1, -1>
using   cut_t = process::confined_t<cut<M_stop_, M_side>>;

template <auto M_stop_, int M_side=_detail::flank_n<M_stop_>, auto ...Ns>
XTAL_DEF_(return,inline,let)
cut_f(auto &&o)
noexcept -> decltype(auto)
{
	return cut_t<M_stop_, M_side>::template method_f<Ns...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <auto M_stop_, int M_side> requires in_n<M_side, 1, -1>
struct cut
{
	static auto constexpr M_stop = M_stop_();

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,set)
		edit_f(bond::pack_q auto &o)
		noexcept -> auto
		{
		//	TODO: Accommodate returning materialized `atom::block` from `span`s...
			return [&]<auto ...I> (bond::seek_t<I...>)
				XTAL_0FN_(to) (bond::pack_f(edit_f(get<I>(o))...))
			(bond::antiseek_s<bond::pack_size<decltype(o)>{}()>{});
		}
		template <auto ...Ns>
		XTAL_DEF_(return,set)
		edit_f(real_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _xtd::bit_cast;

			using _fit = bond::fit<XTAL_ALL_(o)>;
			using U_delta = typename _fit::delta_type;
			using U_sigma = typename _fit::sigma_type;
			using U_alpha = typename _fit::alpha_type;

			auto constexpr o_stop = static_cast<U_alpha>(magnum_f(M_stop));
			auto constexpr N_side = static_cast<U_delta>(M_side);
			//\
			if (_std::is_constant_evaluated() or not _fit::IEC) {
			if (_std::is_constant_evaluated() or not _fit::IEC or XTAL_ENV_(GNUC)) {
				U_alpha const s = signum_t<>::edit_f(o);
				XTAL_IF0
				XTAL_0IF (M_side <= 0) {o = _fit::minimum_f(XTAL_MOV_(o), o_stop);}
				XTAL_0IF (M_side == 1) {o = _fit::maximum_f(XTAL_MOV_(o), o_stop);}
				U_alpha const q = o == o_stop;
				o *= s; return q*s;
			}
			else {
				U_sigma constexpr K_stop = bit_cast<U_sigma>(o_stop);
				U_sigma constexpr K_side = N_side >> _fit::sign.shift;
				auto constexpr sign = _fit::sign;
				auto constexpr unit = _fit::unit;

				auto   &t = reinterpret_cast<U_sigma &>(o);
				U_sigma s = t&sign.mask;
				U_sigma r = t^s;
				s  |=   unit.mask;
				r  -= K_stop;
				U_delta q = bit_cast<U_delta>(r) >> sign.shift;
				q  ^= K_side;
				r  &= q;
				s  &= q;
				t  -= r;
				return bit_cast<U_alpha>(s);
			}
		}
		/**/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(complex_variable_q auto &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto &[x, y] = destruct_f(o);
			return {edit_f<Ns...>(x), edit_f<Ns...>(y)};
		}
		/*/
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(complex_variable_q auto &o)
		noexcept -> XTAL_ALL_(abs(o))
		{
			using _fit = bond::fit<decltype(o)>;

			auto &[x, y] = destruct_f(o);
			auto  [w, m] = dots_f<2>(o);
			auto r = edit_f<Ns...>(w);
			w *= m;
			x *= w;
			y *= w;
			return r;
		}
		/***/

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto o)
		noexcept -> auto
		{
			(void) edit_f(o); return o;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
