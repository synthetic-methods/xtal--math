#pragma once
#include "./any.hh"

#include "../../process/pade/unity.hh"
#include "../../process/taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::atom::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  uniplex;
template <class   ..._s>	using   uniplex_t = typename uniplex<_s...>::type;
template <class   ...Ts>	concept uniplex_q = bond::tag_in_p<uniplex_t, Ts...>;

XTAL_DEF_(let) uniplex_f = [] XTAL_1FN_(call) (_detail::factory<uniplex_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `couple` with Dirichlet characterization and modulo access.
*/
template <class A>
struct uniplex<A>
:	uniplex<fixed_valued_u<A>>
{
};
template <scalar_q A> requires simplex_field_q<A>
struct uniplex<A>
{
private:
	using value_type = A;

	using complex_type = _std::complex<A>;
	using couplex_type = couple_t<complex_type[2]>;
	using duplex_type  = couple_t<  value_type[2]>;

	XTAL_DEF_(return,inline,set)  cis_2pi (auto &&o) {return process::math::  pade::    unity_t< 1   >::template method_f< 4>(XTAL_REF_(o));};
	XTAL_DEF_(return,inline,set) _exp_2pi (auto &&o) {return process::math::taylor::logarithm_t<-1, 1>::template method_f< 2>(XTAL_REF_(o) *bond::fit<decltype(o)>::template patio_f<+1>(-2));};

	template <class T>
	using endotype = typename couple<complex_type, duplex_type>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<uniplex_t>>;

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

	public:// CONSTRUCT
		using S_::S_;

		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(complex_variable_q auto &&x)
		noexcept
		:	S_{cis_2pi(x.real()), _exp_2pi(x.imag())}
		{
		}

	public:// CONSTRUCT

		template <int N=0>
		XTAL_DEF_(return,inline,let)
		resolved()
		noexcept -> _std::conditional_t<N == 0, couplex_type, complex_type>
		{
			auto const &[o, q_] = self();
			auto const   o_re = o.real(), q_up = get<0>(q_);
			auto const   o_im = o.imag(), q_dn = get<1>(q_);
			XTAL_IF0
			XTAL_0IF (N ==  0) {
				return {{o_re*q_up, o_im*q_up}, {o_re*q_up, o_im*q_up}};
			}
			XTAL_0IF (N ==  1) {
				return {o_re*q_up, o_im*q_up};
			}
			XTAL_0IF (N == -1) {
				return {o_re*q_up, o_im*q_up};
			}
		}
		template <int N=0>
		XTAL_DEF_(return,inline,let)
		desolved()
		noexcept -> _std::conditional_t<N == 0, couplex_type, complex_type>
		{
			auto const &[o, q_] = self();
			auto const   o_re = o.real(), q_up = q_.template sum<+1>();
			auto const   o_im = o.imag(), q_dn = q_.template sum<-1>();
			XTAL_IF0
			XTAL_0IF (N ==  0) {
				return {{o_re*q_up, o_im*q_dn}, {o_re*q_dn, o_im*q_up}};
			}
			XTAL_0IF (N ==  1) {
				return {o_re*q_up,  o_im*q_dn};
			}
			XTAL_0IF (N == -1) {
				return {o_re*q_dn,  o_im*q_up};
			}
		}
		template <int N_side=1>
		XTAL_DEF_(return,inline,let)
		sum(auto &&w=value_type{})
		noexcept -> auto
		{
			return XTAL_REF_(w) + desolved<N_side>();
		}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
