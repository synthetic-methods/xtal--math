#pragma once
#include "./any.hh"
#include "./nearing.hh"
#include "./decompose.hh"



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
	static auto constexpr M_val = M_app();
	static auto constexpr M_abs = decompose_f<unsigned>(M_val);
	static auto constexpr M_sgn = decompose_f<  signed>(M_val);
	static auto constexpr M_dir = (int) decompose_f<signed>(nearing_f(M_abs)*(M_sgn));
	static auto constexpr M_opp = [] XTAL_1FN_(to) (-one/M_val);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

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

			auto constexpr o_stop = static_cast<U_alpha>(M_abs);
			auto constexpr N_side = static_cast<U_delta>(M_dir);
			//\
			if (_std::is_constant_evaluated() or not _fit::IEC) {
			if (_std::is_constant_evaluated() or not _fit::IEC or XTAL_ENV_(GNUC)) {
				U_alpha const s = decompose_t<signed>::edit_f(o);
				XTAL_IF0
				XTAL_0IF (M_dir <= 0) {o = _fit::minimum_f(XTAL_MOV_(o), o_stop);}
				XTAL_0IF (M_dir == 1) {o = _fit::maximum_f(XTAL_MOV_(o), o_stop);}
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
		XTAL_DEF_(return,set)
		edit_f(bond::pack_q auto &w_)
		noexcept -> auto
		requires un_n<complex_variable_q<decltype(w_)>>
		{
			XTAL_TYP_(let) W_ = XTAL_ALL_(w_);
			auto constexpr N_ = bond::pack_size_n<W_>;
			XTAL_IF0
			XTAL_0IF (0 == N_) {
				return W_{};
			}
			XTAL_0IF (1 == N_) {
				return W_{edit_f<Ns...>(get<0>(w_))};
			}
			XTAL_0IF (2 == N_ and atom::couple_q<W_>) {
			//	auto &[w0, w1] = w_;
				auto &w0 = get<0>(w_); using W0 = objective_t<decltype(w0)>;
				auto &w1 = get<1>(w_); using W1 = objective_t<decltype(w1)>;
				XTAL_IF0
				//\
				XTAL_0IF (wniplex_q<W_>) {
				XTAL_0IF (complex_variable_q<W0> and atom::couple_q<W1>) {
					auto const u0 =                          edit_f<Ns...>(get<0>(w1));
					auto const u1 = limit_t<M_opp>::template edit_f<Ns...>(get<1>(w1));
					return W_{W0{one, zero}, W1{XTAL_MOV_(u0), XTAL_MOV_(u1)}};
				}
				XTAL_0IF (atom::couple_q<W0> and complex_variable_q<W1>) {
					auto const u0 =                          edit_f<Ns...>(get<0>(w0));
					auto const u1 = limit_t<M_opp>::template edit_f<Ns...>(get<1>(w0));
					return W_{W0{XTAL_MOV_(u0), XTAL_MOV_(u1)}, W1{one, zero}};
				}
				//\
				XTAL_0IF (same_q<W0, W1>) {
				XTAL_0IF_(else) {
					return W_{edit_f<Ns...>(w0), edit_f<Ns...>(w1)};
				}
			}
			XTAL_0IF_(else) {
			//	TODO: Accommodate returning materialized `atom::block` from `span`s...
				return [&]<auto ...I> (bond::seek_t<I...>)
					XTAL_0FN_(to) (W_(edit_f<Ns...>(get<I>(w_))...))
				(bond::seek_reverse_s<N_>{});
			}
		}

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
	static constexpr auto N_dir = decompose_f<unsigned>(M_dir);

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
	//	TODO: Implement smooth version (together with smooth `decompose<unsigned>`) (via `Ns...`)?

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
			U const   w{m_vs + decompose_f<unsigned>(u - m_dn) - decompose_f<unsigned>(u - m_up)};
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
