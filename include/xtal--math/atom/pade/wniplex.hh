#pragma once
#include "./any.hh"

#include "../../process/pade/unity.hh"
#include "../../process/taylor/octarithm.hh"
#include "../../process/limit.hh"
#include "../../process/term.hh"

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
\brief   Extends `couple` with Dirichlet characterization and modulo access.
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
		using typename S_::value_type;

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let)
		unicycle(), S_::template element<0>())

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,let)
		magnitudes(), S_::template element<1>())

	public:// CONSTRUCT
		using S_::S_;

		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(complex_variable_q auto &&x)
		noexcept
		:	S_{
				process::math::  pade::    unity_f<(+1)>(x.real())
			,	process::math::taylor::octarithm_f<(-1)>(x.imag())
			}
		{
		}

	public:// CONSTRUCT

		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		limit()
		noexcept -> auto &
		{
			auto &s = self();
			(void) process::math::limit_t<-4>::edit_f(bond::pack_item_f<1, 1>(s));
			return s;
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		limited() const
		noexcept -> auto
		{
			return twin().template limit<N_lim>();
		}
		template <int N_sgn=0>
		XTAL_DEF_(return,inline,let)
		resolved() const
		noexcept -> _std::conditional_t<N_sgn == 0, couplex_type, complex_type>
		{
			auto const &[o, q_] = self();
			auto const   o_re = o.real(), q_up = get<0>(q_);
			auto const   o_im = o.imag(), q_dn = get<1>(q_);
			XTAL_IF0
			XTAL_0IF (N_sgn ==  0) {
				return {{o_re*q_up, o_im*q_up}, {o_re*q_up, o_im*q_up}};
			}
			XTAL_0IF (N_sgn ==  1) {
				return {o_re*q_up, o_im*q_up};
			}
			XTAL_0IF (N_sgn == -1) {
				return {o_re*q_up, o_im*q_up};
			}
		}
		template <int N_sgn=0>
		XTAL_DEF_(return,inline,let)
		reflexed() const
		noexcept -> _std::conditional_t<N_sgn == 0, couplex_type, complex_type>
		{
			//\
			using _xtd::accumulator;
			using process::math::term_f;
			auto const &[o, q_] = self();
			auto const   q_up = q_.template sum<+1>();
			auto const   q_dn = q_.template sum<-1>();
			auto const &[o_re, o_im] = destruct_f(o);
			XTAL_IF0
			XTAL_0IF (N_sgn ==  0) {
				return {
					{o_re*q_up, o_im*q_dn},
					{o_re*q_dn, o_im*q_up}
				};
			}
			XTAL_0IF (N_sgn ==  1) {
				return {o_re*q_up, o_im*q_dn};
			}
			XTAL_0IF (N_sgn == -1) {
				return {o_re*q_dn, o_im*q_up};
			}
		}
		template <int N_sgn=0>
		XTAL_DEF_(return,inline,let)
		reflexed(complex_type const x) const
		noexcept -> _std::conditional_t<N_sgn == 0, couplex_type, complex_type>
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
			XTAL_0IF (N_sgn ==  0) {
				return {
					{term_f(x_re, o_re, q_up), term_f(x_im, o_im, q_dn)},
					{term_f(x_re, o_re, q_dn), term_f(x_im, o_im, q_up)}
				};
			}
			XTAL_0IF (N_sgn ==  1) {
				return {term_f(x_re, o_re, q_up), term_f(x_im, o_im, q_dn)};
			}
			XTAL_0IF (N_sgn == -1) {
				return {term_f(x_re, o_re, q_dn), term_f(x_im, o_im, q_up)};
			}
		}
		template <int N_side=1>
		XTAL_DEF_(return,inline,let)
		sum(auto &&w=value_type{}) const
		noexcept -> auto
		{
			return XTAL_REF_(w) + reflexed<N_side>();
		}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
