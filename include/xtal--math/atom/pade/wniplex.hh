#pragma once
#include "./any.hh"



#include "../../process/pade/unity.hh"
#include "../../process/taylor/octarithm.hh"

XTAL_ENV_(push)
namespace xtal::atom::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  wniplex;
template <class   ..._s>	using   wniplex_t = typename wniplex<_s...>::type;
template <class   ...Ts>	concept wniplex_q = bond::tag_in_p<wniplex_t, Ts...>;

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(return,inline,let)
wniplex_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return atom::_detail::factory<wniplex_t>::make(XTAL_REF_(oo)...);
}
XTAL_DEF_(return,inline,let)
wniplex_f(decltype(_std::in_place), simplex_variable_q auto &&o)
noexcept -> decltype(auto)
{
	using  W = wniplex_t<XTAL_ALL_(o)>;
	return W{{one, zero}, process::math::roots_f<1>(XTAL_REF_(o))};
}
XTAL_DEF_(return,inline,let)
wniplex_f(decltype(_std::in_place), complex_variable_q auto &&o)
noexcept -> decltype(auto)
{
	using  W = wniplex_t<XTAL_ALL_(o)>;
	auto   const vs = process::math::roots_f<2>(process::math::dot_f(o));
	auto   const dn = get<1>(vs);
	return W{XTAL_REF_(o)*XTAL_MOV_(dn), XTAL_MOV_(vs)};
}
XTAL_DEF_(return,inline,let)
wniplex_f(complex_variable_q auto &&o, decltype(_std::in_place))
noexcept -> decltype(auto)
{
	using  W = wniplex_t<XTAL_ALL_(o)>;
	return W{XTAL_REF_(o), {one, one}};
}


}///////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(let) wniplex_f = [] XTAL_1FN_(call) (_detail::wniplex_f);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides a reciprocal complex coupling,
optimized for multiplication and scalar reflection.
*/
template <class A>
struct wniplex<A>
:	wniplex<typename fixed<A>::value_type>
{
};
template <scalar_q A> requires simplex_field_q<A>
struct wniplex<A>
{
	using simplex_type =  A;
	using complex_type = _std::complex<simplex_type>;
	using couplex_type =      couple_t<complex_type[2]>;
	using  duplex_type =      couple_t<simplex_type[2]>;

private:
	XTAL_DEF_(set) sig_f = process::math::  pade::    unity_f<+1>;
	XTAL_DEF_(set) mag_f = process::math::taylor::octarithm_f<-1>;

	template <class T>
	//\
	using endotype = typename group_multiplication<complex_type, duplex_type>::template homotype<T>;
	using endotype = typename couple<complex_type, duplex_type>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<wniplex_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;
		using I_ = typename S_::difference_type;

	public:// TYPE
		using signum_type = complex_type;
		using magnum_type =  duplex_type;

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		signum(), S_::template element<0>())

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		magnum(), S_::template element<1>())

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		signum(constant_q auto &&n), signum<XTAL_ALL_(n){}>())

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		magnum(constant_q auto &&n), magnum<XTAL_ALL_(n){}>())

		XTAL_FX4_(do) (template <int N_pow>
		XTAL_DEF_(return,inline,get)
		signum(),
		{
			XTAL_IF0
			XTAL_0IF (N_pow ==  0) {return complex_type{one};}
			XTAL_0IF (N_pow ==  1) {return        (signum());}
			XTAL_0IF (N_pow == -1) {return    conj(signum());}
		})
		XTAL_FX4_(do) (template <int N_pow>
		XTAL_DEF_(return,inline,get)
		magnum(),
		{
			XTAL_IF0
			XTAL_0IF (N_pow ==  0) {return   simplex_type{one};}
			XTAL_0IF (N_pow ==  1) {return  get<0>(magnum());}
			XTAL_0IF (N_pow == -1) {return  get<1>(magnum());}
		})

	public:// CONSTRUCT
	//	using S_::S_;
		XTAL_NEW_(else) (homotype, noexcept:S_)

		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(complex_type o, duplex_type w)
		noexcept
		:	S_{XTAL_MOV_(o), XTAL_MOV_(w)}
		{
		}
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(complex_variable_q auto &&x)
		noexcept
		:	homotype{
				signum_type{sig_f(x.real())}
			,	magnum_type{mag_f(x.imag()), _std::in_place}
			}
		{
		}
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(simplex_variable_q auto &&x)
		noexcept
		:	homotype{
				signum_type{sig_f(x)}
			,	magnum_type{one, one}
			}
		{
		}
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(decltype(one))
		noexcept
		:	homotype{
				signum_type{one     }
			,	magnum_type{one, one}
			}
		{
		}

	public:// DECONSTRUCT


		template <int N_dir=0>
		XTAL_DEF_(return,inline,let)
		resolution() const
		noexcept -> complex_type
		{
			auto const &[o, q_] = self();
			auto const q_up = get<0>(q_);
			auto const q_dn = get<1>(q_);
			XTAL_IF0
			XTAL_0IF (N_dir ==  1) {return {q_up*o.real(), q_up*o.imag()};}
			XTAL_0IF (N_dir == -1) {return {q_dn*o.real(),-q_dn*o.imag()};}
		}
		XTAL_DEF_(return,inline,let)
		resolution(constant_q auto const n) const
		noexcept -> decltype(auto)
		{
			return resolution<XTAL_ALL_(n){}>();
		}

		XTAL_DEF_(return,inline,let)
		resolution() const
		noexcept -> couplex_type
		{
			auto const &[o, q_] = self();
			auto const q_up = get<0>(q_);
			auto const q_dn = get<1>(q_);
			return {
				{q_up*o.real(), q_up*o.imag()},
				{q_dn*o.real(),-q_dn*o.imag()}
			};
		}


		template <int N_dir>
		XTAL_DEF_(return,inline,let)
		reflection(complex_type const x) const
		noexcept -> complex_type
		{
			//\
			using _xtd::plus_multiplies_f;
			using process::math::term_f;
			auto const &[o, q_] = self();
			return {
				term_f(x.real(), o.real(), q_.template sum<+N_dir>())
			,	term_f(x.imag(), o.imag(), q_.template sum<-N_dir>())
			};
		}
		template <int N_dir>
		XTAL_DEF_(return,inline,let)
		reflection() const
		noexcept -> complex_type
		{
			auto const &[o, q_] = self();
			return {
				o.real()*q_.template sum<+N_dir>()
			,	o.imag()*q_.template sum<-N_dir>()
			};
		}
		XTAL_DEF_(return,inline,let)
		reflection(constant_q auto const n) const
		noexcept -> decltype(auto)
		{
			return reflection<XTAL_ALL_(n){}>();
		}

		XTAL_DEF_(return,inline,let)
		reflection(complex_type const x) const
		noexcept -> couplex_type
		{
			//\
			using _xtd::plus_multiplies_f;
			using process::math::term_f;
			auto const &[o, q_] = self();
			auto const q_up = q_.template sum<+1>();
			auto const q_dn = q_.template sum<-1>();
			return {
				{term_f(x.real(), o.real(), q_up), term_f(x.imag(), o.imag(), q_dn)},
				{term_f(x.real(), o.real(), q_dn), term_f(x.imag(), o.imag(), q_up)}
			};
		}
		XTAL_DEF_(return,inline,let)
		reflection() const
		noexcept -> couplex_type
		{
			auto const &[o, q_]  = self();
			auto const      q_up = q_.template sum<+1>();
			auto const      q_dn = q_.template sum<-1>();
			return {
				{q_up*o.real(), q_dn*o.imag()},
				{q_dn*o.real(), q_up*o.imag()}
			};
		}


	public:// OPERATE

		template <int N_side=1>
		XTAL_DEF_(return,inline,let)
		sum(auto &&w=simplex_type{}) const
		noexcept -> auto
		{
			return XTAL_REF_(w) + reflection<N_side>();
		}

		XTAL_DEF_(return,inline,let)
		flipped(simplex_type const w) const
		noexcept -> auto
		{
			auto const &o  = signum();
			auto const &q_ = magnum();
			return S_::form(complex_type{o.real(), w*o.imag()}, q_.flipped(w));
		}
		XTAL_DEF_(return,inline,let)
		flipped() const
		noexcept -> auto
		{
			auto const &o  = signum();
			auto const &q_ = magnum();
			return S_::form(conj(o), q_.flipped());
		}
		XTAL_DEF_(return,inline,let)
		operator ~ () const
		noexcept -> auto
		{
			return flipped();
		}

		using S_::operator*=;
		using S_::operator/=;
		XTAL_DEF_(mutate,inline,get) operator *=(complex_variable_q          auto &&t) noexcept {auto &s = signum(); s *=      XTAL_REF_(t) ; return self();}
		XTAL_DEF_(mutate,inline,get) operator /=(complex_variable_q          auto &&t) noexcept {auto &s = signum(); s *= conj(XTAL_REF_(t)); return self();}
		XTAL_DEF_(mutate,inline,get) operator +=(                   homotype      &&t) noexcept {auto &s =   self(); s *=      XTAL_MOV_(t) ; return self();}
		XTAL_DEF_(mutate,inline,get) operator -=(                   homotype      &&t) noexcept {auto &s =   self(); s /=      XTAL_MOV_(t) ; return self();}
		XTAL_DEF_(mutate,inline,get) operator +=(                   homotype const &t) noexcept {auto &s =   self(); s *=      XTAL_REF_(t) ; return self();}
		XTAL_DEF_(mutate,inline,get) operator -=(                   homotype const &t) noexcept {auto &s =   self(); s /=      XTAL_REF_(t) ; return self();}
		XTAL_DEF_(return,inline,met) operator + (homotype const &s, homotype const &t) noexcept {return s * t;}
		XTAL_DEF_(return,inline,met) operator - (homotype const &s, homotype const &t) noexcept {return s / t;}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
