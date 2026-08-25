#pragma once
#include "./any.hh"


#include "../../process/roots.hh"
#include "../../process/pade/unity.hh"
#include "../../process/taylor/octarithm.hh"

XTAL_ENV_(push)
namespace xtal::atom::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Represents a complex reciprocal-pair via the sign and reciprocal magnitudes.
*/
template <class   ..._s>	XTAL_TYP_(new) uniplex;
template <class   ..._s>	XTAL_TYP_(let) uniplex_t = typename uniplex<_s...>::type;
template <class   ...Ts>	XTAL_TYP_(ask) uniplex_q = bond::tag_inner_p<uniplex_t, Ts...>;

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

XTAL_VAL_(return,inline,let)
uniplex_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return atom::_detail::factory<uniplex_t>::make(XTAL_REF_(oo)...);
}
XTAL_VAL_(return,inline,let)
uniplex_f(decltype(std::in_place), simplex_variable_q auto &&o)
noexcept -> decltype(auto)
{
	using  W = uniplex_t<XTAL_ALL_(o)>;
	return W{{one, zero}, process::math::roots_f<1>(XTAL_REF_(o))};
}
XTAL_VAL_(return,inline,let)
uniplex_f(decltype(std::in_place), complex_variable_q auto &&o)
noexcept -> decltype(auto)
{
	using  W = uniplex_t<XTAL_ALL_(o)>;
	auto   const vs = process::math::roots_f<2>(process::math::dot_f(o));
	auto   const dn = get<1>(vs);
	return W{XTAL_REF_(o)*XTAL_MOV_(dn), XTAL_MOV_(vs)};
}
XTAL_VAL_(return,inline,let)
uniplex_f(complex_variable_q auto &&o, decltype(std::in_place))
noexcept -> decltype(auto)
{
	using  W = uniplex_t<XTAL_ALL_(o)>;
	return W{XTAL_REF_(o), {one, one}};
}


}///////////////////////////////////////////////////////////////////////////////

XTAL_VAL_(let) uniplex_f = [] XTAL_1FN_(call) (_detail::uniplex_f);


////////////////////////////////////////////////////////////////////////////////

template <class A>
struct uniplex<A>
:	uniplex<typename fixed<A>::value_type>
{
};
template <extra_scalar_q A> requires simplex_field_q<A>
struct uniplex<A>
{
	using simplex_type = A;
	using complex_type = std::complex <simplex_type   >;
	using couplex_type =      couple_t<complex_type[2]>;
	using  duplex_type =      couple_t<simplex_type[2]>;
//	using couplex_type = XTAL_ALL_(process::math::roots_f<1>(XTAL_ANY_(complex_type)));
//	using  duplex_type = XTAL_ALL_(process::math::roots_f<1>(XTAL_ANY_(simplex_type)));


private:
	XTAL_VAL_(set) sig_f = process::math::  pade::    unity_f<+1>;
	XTAL_VAL_(set) mag_f = process::math::taylor::octarithm_f<-1>;

	template <class T>
	using endotype = typename quantity_multiplies<complex_type, duplex_type>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<uniplex_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;
		using I_ = typename S_::difference_type;

	public:// TYPE
		using signum_type = complex_type;
		using magnum_type =  duplex_type;

	public:// CONSTRUCT
	//	using S_::S_;
		XTAL_VAL_(delete) (homotype, noexcept=default)
		XTAL_VAL_(create) (homotype, noexcept=default)
		XTAL_VAL_(move)   (homotype, noexcept=default)
		XTAL_VAL_(copy)   (homotype, noexcept=default)
		XTAL_VAL_(induce) (homotype, noexcept,homotype)
		XTAL_VAL_(reduce) (homotype, noexcept,S_)

		XTAL_VAL_(inline)
		XTAL_VAL_(new,explicit)
		homotype(signum_type &&o, magnum_type &&q)
		noexcept
		:	S_{XTAL_MOV_(o), XTAL_MOV_(q)}
		{
		}
		XTAL_VAL_(inline)
		XTAL_VAL_(new,explicit)
		homotype(complex_variable_q auto &&x)
		noexcept
		:	homotype{
				signum_type{sig_f(x.real())}
			,	magnum_type{mag_f(x.imag()), std::in_place}
			}
		{
		}
		XTAL_VAL_(inline)
		XTAL_VAL_(new,explicit)
		homotype(simplex_variable_q auto &&x)
		noexcept
		:	homotype{
				signum_type{sig_f(x)}
			,	magnum_type{one, one}
			}
		{
		}
		XTAL_VAL_(inline)
		XTAL_VAL_(new,explicit)
		homotype(decltype(one))
		noexcept
		:	homotype{
				signum_type{one     }
			,	magnum_type{one, one}
			}
		{
		}

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_VAL_(return,inline,set)
		magnum_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return xtd::qualify_cast<T>(XTAL_REF_(o)).template element<1>();
		}
		XTAL_VAL_(return,inline,set)
		signum_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return xtd::qualify_cast<T>(XTAL_REF_(o)).template element<0>();
		}

		template <int N_pow>
		XTAL_VAL_(return,inline,set)
		magnum_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (N_pow ==  0) {return simplex_type{one};}
			XTAL_0IF (N_pow ==  1) {return  get<0>(magnum_f(XTAL_REF_(o)));}
			XTAL_0IF (N_pow == -1) {return  get<1>(magnum_f(XTAL_REF_(o)));}
		}
		template <int N_pow>
		XTAL_VAL_(return,inline,set)
		signum_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (N_pow ==  0) {return complex_type{one};}
			XTAL_0IF (N_pow ==  1) {return        (signum_f(XTAL_REF_(o)));}
			XTAL_0IF (N_pow == -1) {return    conj(signum_f(XTAL_REF_(o)));}
		}

		XTAL_FN1_(go) (template <auto ...Ns> XTAL_VAL_(return,inline,get) magnum, magnum_f<Ns...>)
		XTAL_FN1_(go) (template <auto ...Ns> XTAL_VAL_(return,inline,get) signum, signum_f<Ns...>)


	public:// RECONSTRUCT

		template <int N_dir=0>
		XTAL_VAL_(return,inline,set)
		resolution_f(auto &&that)
		noexcept -> auto
		{
			auto &&[o, q_] = xtd::qualify_cast<T>(XTAL_REF_(that));
			auto const q_up = get<0>(q_);
			auto const q_dn = get<1>(q_);
			XTAL_IF0
			XTAL_0IF (N_dir ==  0) {return couplex_type {{q_up*o.real(), q_up*o.imag()}, {q_dn*o.real(),-q_dn*o.imag()}};}
			XTAL_0IF (N_dir ==  1) {return complex_type  {q_up*o.real(), q_up*o.imag()}                                 ;}
			XTAL_0IF (N_dir == -1) {return complex_type                                  {q_dn*o.real(),-q_dn*o.imag()} ;}
		}


		template <int N_dir=0>
		XTAL_VAL_(return,inline,set)
		reflection_f(auto &&that)
		noexcept -> auto
		{
			return reflection_f<N_dir>(XTAL_REF_(that), zero);
		}
		template <int N_dir=0> requires un_v<N_dir, 0>
		XTAL_VAL_(return,inline,set)
		reflection_f(auto &&that, auto const &plus)
		noexcept -> complex_type
		{
			using process::math::term_f;
			auto &&[o, q_] = xtd::qualify_cast<T>(XTAL_REF_(that));
			XTAL_IF0
			XTAL_0IF (same_q<decltype(zero), decltype(plus)>) {return {
				o.real()*q_.template sum<+N_dir>()
			,	o.imag()*q_.template sum<-N_dir>()
			};}
			XTAL_0IF (simplex_field_q<       decltype(plus)>) {return {
				term_f(plus       , o.real(), q_.template sum<+N_dir>())
			,	term_f(zero       , o.imag(), q_.template sum<-N_dir>())
			};}
			XTAL_0IF (complex_field_q<       decltype(plus)>) {return {
				term_f(plus.real(), o.real(), q_.template sum<+N_dir>())
			,	term_f(plus.imag(), o.imag(), q_.template sum<-N_dir>())
			};}
		}
		template <int N_dir=0> requires in_v<N_dir, 0>
		XTAL_VAL_(return,inline,set)
		reflection_f(auto &&that, auto const &plus=zero)
		noexcept -> couplex_type
		{
			using process::math::term_f;
			auto &&[o, q_] = xtd::qualify_cast<T>(XTAL_REF_(that));
			auto const q_up = q_.template sum<+1>();
			auto const q_dn = q_.template sum<-1>();
			XTAL_IF0
			XTAL_0IF (same_q<decltype(zero), decltype(plus)>) {return {
				{q_up*o.real(), q_dn*o.imag()},
				{q_dn*o.real(), q_up*o.imag()}
			};}
			XTAL_0IF (simplex_field_q<       decltype(plus)>) {return {
				{term_f(plus       , o.real(), q_up), term_f(zero       , o.imag(), q_dn)},
				{term_f(plus       , o.real(), q_dn), term_f(zero       , o.imag(), q_up)}
			};}
			XTAL_0IF (complex_field_q<       decltype(plus)>) {return {
				{term_f(plus.real(), o.real(), q_up), term_f(plus.imag(), o.imag(), q_dn)},
				{term_f(plus.real(), o.real(), q_dn), term_f(plus.imag(), o.imag(), q_up)}
			};}
		}
		template <int N_dir=0>
		XTAL_VAL_(return,inline,set)
		reflection_f(auto &&that, auto &&plus, auto &&...times_)
		noexcept -> auto
		requires (1 <= sizeof...(times_))
		{
			auto &&[o, q_] = xtd::qualify_cast<T>(XTAL_REF_(that));
			return reflection_f<N_dir>(T{XTAL_REF_(o), XTAL_REF_(q_)*(times_ *...* one)}, XTAL_REF_(plus));
		}

		XTAL_FN1_(go) (template <int N_dir=0> XTAL_VAL_(return,inline,get) resolution, resolution_f<N_dir>)
		XTAL_FN1_(go) (template <int N_dir=0> XTAL_VAL_(return,inline,get) reflection, reflection_f<N_dir>)
		XTAL_FN1_(go) (template <int N_dir=1> XTAL_VAL_(return,inline,get)        sum, reflection_f<N_dir>)

	public:// OPERATE

		XTAL_VAL_(return,inline,let)
		flipped(simplex_type const w) const
		noexcept -> auto
		{
			auto const &o  = signum();
			auto const &q_ = magnum();
			return S_::form_f(complex_type{o.real(), w*o.imag()}, q_.flipped(w));
		}
		XTAL_VAL_(return,inline,let)
		flipped() const
		noexcept -> auto
		{
			auto const &o  = signum();
			auto const &q_ = magnum();
			return S_::form_f(conj(o), q_.flipped());
		}
		XTAL_VAL_(return,inline,let)
		operator ~ () const
		noexcept -> auto
		{
			return flipped();
		}

		using S_::operator*=;
		using S_::operator/=;
		XTAL_VAL_(mutate,inline,get) operator *=(complex_variable_q          auto &&t) noexcept {auto &s = signum(); s *=      XTAL_REF_(t) ; return self();}
		XTAL_VAL_(mutate,inline,get) operator /=(complex_variable_q          auto &&t) noexcept {auto &s = signum(); s *= conj(XTAL_REF_(t)); return self();}
		XTAL_VAL_(mutate,inline,get) operator +=(                   homotype      &&t) noexcept {auto &s =   self(); s *=      XTAL_MOV_(t) ; return self();}
		XTAL_VAL_(mutate,inline,get) operator -=(                   homotype      &&t) noexcept {auto &s =   self(); s /=      XTAL_MOV_(t) ; return self();}
		XTAL_VAL_(mutate,inline,get) operator +=(                   homotype const &t) noexcept {auto &s =   self(); s *=      XTAL_REF_(t) ; return self();}
		XTAL_VAL_(mutate,inline,get) operator -=(                   homotype const &t) noexcept {auto &s =   self(); s /=      XTAL_REF_(t) ; return self();}
		XTAL_VAL_(return,inline,met) operator + (homotype const &s, homotype const &t) noexcept {return s * t;}
		XTAL_VAL_(return,inline,met) operator - (homotype const &s, homotype const &t) noexcept {return s / t;}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
