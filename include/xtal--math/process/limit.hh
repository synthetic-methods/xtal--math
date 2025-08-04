#pragma once
#include "./any.hh"
#include "./nearing.hh"
#include "./component.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto M_lim>
XTAL_TYP_(new) limit;

template <auto M_lim>
XTAL_TYP_(let) limit_t = process::confined_t<limit<M_lim>>;

template <auto M_lim, auto ...Ns>
XTAL_DEF_(let) limit_f = [] XTAL_1FN_(call) (limit_t<M_lim>::template method_f<Ns...>);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <applicable_q auto M_app>
struct limit<M_app>
{
	static auto constexpr L_abs = component_f<unsigned>(M_app());
	static auto constexpr L_sgn = component_f<  signed>(M_app());
	static auto constexpr L_dir = (int) component_f<signed>(nearing_f(L_abs)*(L_sgn));

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
		requires un_n<complex_variable_q<decltype(o)>>
		{
		//	TODO: Accommodate returning materialized `atom::block` from `span`s...
			return [&]<auto ...I> (bond::seek_t<I...>)
				XTAL_0FN_(to) (bond::pack_f(edit_f(get<I>(o))...))
			(bond::seek_reverse_s<bond::pack_size<decltype(o)>{}()>{});
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

			auto constexpr o_stop = static_cast<U_alpha>(L_abs);
			auto constexpr N_side = static_cast<U_delta>(L_dir);
			//\
			if (_std::is_constant_evaluated() or not _fit::IEC) {
			if (_std::is_constant_evaluated() or not _fit::IEC or XTAL_ENV_(GNUC)) {
				U_alpha const s = component_t<signed>::edit_f(o);
				XTAL_IF0
				XTAL_0IF (L_dir <= 0) {o = _fit::minimum_f(XTAL_MOV_(o), o_stop);}
				XTAL_0IF (L_dir == 1) {o = _fit::maximum_f(XTAL_MOV_(o), o_stop);}
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
			auto  [w, m] = roots_f<2>(dot_f(o));
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
template <int M_dir>
struct limit<M_dir>
{
	static constexpr auto N_dir = component_f<unsigned>(M_dir);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&u)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(u);
			using U_fit = bond::fit<U>;
			auto constexpr N_min = [] XTAL_1FN_(to) (+U_fit::minilon_f(N_dir));
			auto constexpr N_max = [] XTAL_1FN_(to) (-U_fit::maxilon_f(N_dir));
			XTAL_IF0
			XTAL_0IF (0 < M_dir) {return limit_t<N_min>::template method_f<Ns...>(XTAL_REF_(u));}
			XTAL_0IF (M_dir < 0) {return limit_t<N_max>::template method_f<Ns...>(XTAL_REF_(u));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		edit_f(auto &u)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(u);
			using U_fit = bond::fit<U>;
			auto constexpr N_min = [] XTAL_1FN_(to) (+U_fit::minilon_f(N_dir));
			auto constexpr N_max = [] XTAL_1FN_(to) (-U_fit::maxilon_f(N_dir));
			XTAL_IF0
			XTAL_0IF (0 < M_dir) {return limit_t<N_min>::template   edit_f<Ns...>(XTAL_REF_(u));}
			XTAL_0IF (M_dir < 0) {return limit_t<N_max>::template   edit_f<Ns...>(XTAL_REF_(u));}
		}

	};
};
template <fixed_shaped_q<null_type[2]> auto M_range>
struct limit<M_range>
{
	static auto constexpr M_dn = get<0>(M_range);
	static auto constexpr M_up = get<1>(M_range);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Implement `edit`?
	//	TODO: Implement smooth version (together with smooth `component<unsigned>`) (via `Ns...`)?

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&u)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(u);
			using U_fit = bond::fit<U>;

			U constexpr m_up(M_up);
			U constexpr m_dn(M_dn);
			U constexpr m_vs(M_dn + M_up);
			U const   w{m_vs + component_f<unsigned>(u - m_dn) - component_f<unsigned>(u - m_up)};
			if constexpr (integral_variable_q<U>) {
				return XTAL_MOV_(w) >> 1;
			}
			else {
				return XTAL_MOV_(w)*U_fit::haplo_1;
			}
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
