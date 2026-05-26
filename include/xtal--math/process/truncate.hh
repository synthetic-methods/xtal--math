#pragma once
#include "./any.hh"
#include "./nearing.hh"
#include "./part.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto M_lim>
XTAL_TYP_(new) truncate;

template <auto M_lim>
XTAL_TYP_(let) truncate_t = process::confined_t<truncate<M_lim>>;

template <auto M_lim, auto ...Ns>
XTAL_DEF_(let) truncate_f = [] XTAL_1FN_(call) (truncate_t<M_lim>::template method<Ns...>);

template <auto M_lim>
XTAL_DEF_(let) truncate_e = [] XTAL_1FN_(call) (truncate_t<M_lim>::template method<std::in_place>);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <applicable_q auto M_app>
struct truncate<M_app>
{
	static auto constexpr M_val = M_app();
	static auto constexpr M_abs = part_f<unsigned>(M_val);
	static auto constexpr M_sgn = part_f<  signed>(M_val);
	static auto constexpr M_dir = (int) part_f<signed>(nearing_f(M_abs)*(M_sgn));
	static auto constexpr M_opp = [] XTAL_1FN_(to) (-one/M_val);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <std::in_place_t>
		XTAL_DEF_(inline,set)
		method(real_variable_q auto &o)
		noexcept -> decltype(o)
		{
			using U_fit   = bond::fit<XTAL_ALL_(o)>;
			using U_delta = typename U_fit::delta_type;
			using U_sigma = typename U_fit::sigma_type;
			using U_alpha = typename U_fit::alpha_type;

			auto constexpr o_stop =   static_cast<U_alpha>(M_abs);
			auto constexpr N_side =   static_cast<U_delta>(M_dir);
			auto constexpr K_stop = xtd::bit_cast<U_sigma>(o_stop);
			auto constexpr K_side = xtd::bit_cast<U_sigma>(-N_side >> U_fit::sign.shift);
			//\
			if (std::is_constant_evaluated() or not (U_fit::IEC&559)) {
			if (std::is_constant_evaluated() or not (U_fit::IEC&559) or XTAL_ENV_(GNUC)) {
				U_alpha const s = part_e<signed>(o);
				XTAL_IF0
				XTAL_0IF (M_dir <= 0) {o = U_fit::minimum_f(XTAL_MOV_(o), o_stop);}
				XTAL_0IF (M_dir == 1) {o = U_fit::maximum_f(XTAL_MOV_(o), o_stop);}
				U_alpha const q = o == o_stop;
				o *= s;
			}
			else {
				auto  &t = reinterpret_cast<U_sigma &>(o), r = t&U_fit::sign.mask|K_stop;
				r  -=  t;  r &= K_side^bond::math::bit_sign_f(r);
				t  +=  r;
			}
			return o;
		}
		template <std::in_place_t>
		XTAL_DEF_(inline,set)
		method(complex_variable_q auto &o)
		noexcept -> decltype(o)
		{
			auto &[x, y] = destruct_f(o);
			(void) method<std::in_place>(x);
		//	(void) method<std::in_place>(y);
			return o;
		}
		template <std::in_place_t>
		XTAL_DEF_(inline,set)
		method(bond::pack_q auto &o)
		noexcept -> decltype(o)
		requires un_v<complex_variable_q<decltype(o)>>
		{
			XTAL_TYP_(let) W_ = XTAL_ALL_(o);
			auto constexpr N_ = (int) bond::pack_size_v<W_>;
			XTAL_IF0
			XTAL_0IF (0 == N_) {
			}
			XTAL_0IF (1 == N_) {
				(void) method<std::in_place>(get<0>(o));
			}
			XTAL_0IF (2 == N_ and atom::quantity_multiplies_q<W_>) {
				auto &x = get<0>(o); using X = XTAL_ALL_(x);
				auto &y = get<1>(o); using Y = XTAL_ALL_(y);
				XTAL_IF0
			//	XTAL_0IF (atom::math::pade::uniplex_q<W_>) {
				XTAL_0IF (complex_variable_q<X> and atom::couple_q<Y>) {(void) method<std::in_place>(y);}
				XTAL_0IF (atom::couple_q<X> and complex_variable_q<Y>) {(void) method<std::in_place>(x);}
				XTAL_0IF (     same_q<X, Y>) {
					(void)                              method<std::in_place>(x);
					(void) truncate_t<M_opp>::template method<std::in_place>(y);
				}
				XTAL_0IF (different_q<X, Y>) {
					(void) method<std::in_place>(x);
					(void) method<std::in_place>(y);
				}
			}
			XTAL_0IF_(else) {
			//	TODO: Accommodate returning materialized `atom::bucket` from `span`s...
				[&]<auto ...I> (bond::seek_in_t<I...>)
					XTAL_0FN_(do) (W_(method<std::in_place>(get<I>(o))...))
				(bond::seek_to_t<-N_>{});
			}
			return o;
		}

		template <invariable_q auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto o)
		noexcept -> auto
		{
			(void) method<std::in_place>(o); return o;
		}

	};
};
template <int M_dir>
struct truncate<M_dir>
{
	static constexpr auto N_dir = part_f<unsigned>(M_dir);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <invariable_q auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(o);
			using U_fit = bond::fit<U>;
			auto constexpr N_min = [] XTAL_1FN_(to) (+U_fit::minilon_f(N_dir));
			auto constexpr N_max = [] XTAL_1FN_(to) (-U_fit::maxilon_f(N_dir));
			XTAL_IF0
			XTAL_0IF (0 < M_dir) {return truncate_t<N_min>::template method<Ns...>(XTAL_REF_(o));}
			XTAL_0IF (M_dir < 0) {return truncate_t<N_max>::template method<Ns...>(XTAL_REF_(o));}
		}
		template <std::in_place_t>
		XTAL_DEF_(inline,set)
		method(auto &o)
		noexcept -> decltype(o)
		{
			using U     = XTAL_ALL_(o);
			using U_fit = bond::fit<U>;
			auto constexpr N_min = [] XTAL_1FN_(to) (+U_fit::minilon_f(N_dir));
			auto constexpr N_max = [] XTAL_1FN_(to) (-U_fit::maxilon_f(N_dir));
			XTAL_IF0
			XTAL_0IF (0 < M_dir) {truncate_t<N_min>::template method<std::in_place>(XTAL_REF_(o));}
			XTAL_0IF (M_dir < 0) {truncate_t<N_max>::template method<std::in_place>(XTAL_REF_(o));}
			return o;
		}

	};
};
template <fixed_shaped_q<null_type[2]> auto M_range>
struct truncate<M_range>
{
	static auto constexpr M_dn = get<0>(M_range);
	static auto constexpr M_up = get<1>(M_range);

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Implement `method<std::in_place>`?
	//	TODO: Implement smooth version (together with smooth `part<unsigned>`) (via `Ns...`)?

		template <invariable_q auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> auto
		{
			using U     = XTAL_ALL_(o);
			using U_fit = bond::fit<U>;

			U constexpr m_up(M_up);
			U constexpr m_dn(M_dn);
			U constexpr m_vs(M_dn + M_up);
			U const   w{m_vs + part_f<unsigned>(o - m_dn) - part_f<unsigned>(o - m_up)};
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
