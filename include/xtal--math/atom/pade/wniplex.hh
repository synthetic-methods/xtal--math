#pragma once
#include "./any.hh"

#include "../../process/term.hh"
#include "../../process/pade/unity.hh"
#include "../../process/taylor/octarithm.hh"


XTAL_ENV_(push)
namespace xtal::atom::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  wniplex;
template <class   ..._s>	using   wniplex_t = typename wniplex<_s...>::type;
template <class   ...Ts>	concept wniplex_q = bond::tag_in_p<wniplex_t, Ts...>;

XTAL_DEF_(let) wniplex_f = [] XTAL_1FN_(call) (_detail::factory<wniplex_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides a reciprocal complex coupling,
optimized for multiplication and scalar reflection.
*/
template <class A>
struct wniplex<A>
:	wniplex<fixed_valued_u<A>>
{
};
template <scalar_q A> requires simplex_field_q<A>
struct wniplex<A>
{
private:
	using value_type = A;

	using complex_type = _std::complex<A>;
	using couplex_type = couple_t<complex_type[2]>;
	using duplex_type  = couple_t<  value_type[2]>;

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

	protected:
		template <int N_dir>
		using resolve_t = _std::conditional_t<N_dir == 0, couplex_type, complex_type>;

	public:// TYPE
	//	using typename S_::value_type;

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		signum(), S_::template element<0>())

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		magnum(), S_::template element<1>())

		XTAL_FX4_(to) (template <int N> requires ((N&1) == 0) XTAL_DEF_(return,inline,get)
		signum(),     (S_::template element<0>()))

		XTAL_FX4_(to) (template <int N> requires ((N&1) == 1) XTAL_DEF_(return,inline,get)
		signum(), conj(S_::template element<0>()))

		XTAL_FX4_(to) (template <int N> XTAL_DEF_(return,inline,get)
		magnum(), get<N&1>(S_::template element<1>()))

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
				complex_type{process::math::  pade::    unity_f<(+1)>(x.real())}
			,	 duplex_type{process::math::taylor::octarithm_f<(-1)>(x.imag()), _std::in_place}
			}
		{
		}
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(simplex_variable_q auto &&x)
		noexcept
		:	homotype{
				complex_type{process::math::  pade::    unity_f<(+1)>(x)}
			,	 duplex_type{one, one}
			}
		{
		}
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(decltype(one))
		noexcept
		:	homotype{
				complex_type{one}
			,	 duplex_type{one, one}
			}
		{
		}

	public:// CONSTRUCT

		template <int N_dir=0>
		XTAL_DEF_(return,inline,let)
		resolution() const
		noexcept -> resolve_t<N_dir>
		{
			auto const &[o, q_] = self();
			auto const   o_re = o.real(), q_up = get<0>(q_);
			auto const   o_im = o.imag(), q_dn = get<1>(q_);
			XTAL_IF0
			XTAL_0IF (N_dir ==  0) {
				return {{o_re*q_up,  o_im*q_up}, {o_re*q_dn, -o_im*q_dn}};
			}
			XTAL_0IF (N_dir ==  1) {
				return  {o_re*q_up,  o_im*q_up}                          ;
			}
			XTAL_0IF (N_dir == -1) {
				return                           {o_re*q_dn, -o_im*q_dn} ;
			}
		}
		XTAL_DEF_(return,inline,let)
		resolution(constant_q auto const n) const
		noexcept -> decltype(auto)
		{
			return resolution<XTAL_ALL_(n){}>();
		}

		template <int N_dir=0>
		XTAL_DEF_(return,inline,let)
		reflection(complex_type const x) const
		noexcept -> resolve_t<N_dir>
		{
			//\
			using _xtd::accumulator;
			using process::math::term_f;
			auto const &[o, q_] = self();
			auto const   q_up = q_.template sum<+1>();
			auto const   q_dn = q_.template sum<-1>();
			auto const &[o_re, o_im] = destruct_f(o);
			auto const &[x_re, x_im] = destruct_f(x);
			XTAL_IF0
			XTAL_0IF (N_dir ==  0) {
				return {
					{term_f(x_re, o_re, q_up), term_f(x_im, o_im, q_dn)},
					{term_f(x_re, o_re, q_dn), term_f(x_im, o_im, q_up)}
				};
			}
			XTAL_0IF (N_dir ==  1) {
				return {term_f(x_re, o_re, q_up), term_f(x_im, o_im, q_dn)};
			}
			XTAL_0IF (N_dir == -1) {
				return {term_f(x_re, o_re, q_dn), term_f(x_im, o_im, q_up)};
			}
		}
		template <int N_dir=0>
		XTAL_DEF_(return,inline,let)
		reflection() const
		noexcept -> resolve_t<N_dir>
		{
			//\
			using _xtd::accumulator;
			using process::math::term_f;
			auto const &[o, q_] = self();
			auto const   q_up = q_.template sum<+1>();
			auto const   q_dn = q_.template sum<-1>();
			auto const &[o_re, o_im] = destruct_f(o);
			XTAL_IF0
			XTAL_0IF (N_dir ==  0) {
				return {
					{o_re*q_up, o_im*q_dn},
					{o_re*q_dn, o_im*q_up}
				};
			}
			XTAL_0IF (N_dir ==  1) {
				return {o_re*q_up, o_im*q_dn};
			}
			XTAL_0IF (N_dir == -1) {
				return {o_re*q_dn, o_im*q_up};
			}
		}
		XTAL_DEF_(return,inline,let)
		reflection(constant_q auto const n) const
		noexcept -> decltype(auto)
		{
			return reflection<XTAL_ALL_(n){}>();
		}

		template <int N_side=1>
		XTAL_DEF_(return,inline,let)
		sum(auto &&w=value_type{}) const
		noexcept -> auto
		{
			return XTAL_REF_(w) + reflection<N_side>();
		}

		XTAL_DEF_(return,inline,let)
		flipped(value_type const w) const
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

		XTAL_DEF_(mutate,inline,get) operator +=(                   homotype      &&t) noexcept {auto &s = self(); s *= XTAL_MOV_(t); return s;}
		XTAL_DEF_(mutate,inline,get) operator -=(                   homotype      &&t) noexcept {auto &s = self(); s /= XTAL_MOV_(t); return s;}
		XTAL_DEF_(mutate,inline,get) operator +=(                   homotype const &t) noexcept {auto &s = self(); s *= XTAL_REF_(t); return s;}
		XTAL_DEF_(mutate,inline,get) operator -=(                   homotype const &t) noexcept {auto &s = self(); s /= XTAL_REF_(t); return s;}
		XTAL_DEF_(return,inline,met) operator + (homotype const &s, homotype const &t) noexcept {return s * t;}
		XTAL_DEF_(return,inline,met) operator - (homotype const &s, homotype const &t) noexcept {return s / t;}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
