#pragma once
#include "./any.hh"

#include "../bond/bit.hh"
#include "../process/term.hh"
#include "../process/square.hh"
#include "../process/nearest.hh"

XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  loop;
template <class   ..._s>	using   loop_t = typename loop<_s...>::type;
template <class   ...Ts>	concept loop_q = bond::tag_inner_p<loop_t, Ts...>;

XTAL_DEF_(let) loop_f = [] XTAL_1FN_(call) (_detail::factory<loop_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `group_addition` with Dirichlet characterization and modulo access.
*/
template <scalar_q ..._s> requires same_q<_s...>
struct loop<_s ...>
:	loop<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_q A>
struct loop<A>
{
private:
	using _fit = bond::fit<A>;
	
	template <class T>
	using endotype = typename group_addition<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<loop_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// TYPE
		using typename S_::value_type;
		using typename S_::scale_type;
		using typename S_::index_type;

	//\
	public:// ACCESS
	protected:// ACCESS
		typename S_::index_type pushdex{};

	public:// ACCESS
	//	using S_::element_f;
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_DEF_(set) mask = size_constant_t<size - 1>{};

		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		element_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return T ::element_f(XTAL_REF_(o), I);
		}
		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		element_f(auto &&o, index_type i)
		noexcept -> decltype(auto)
		{
			i += qualify_f<homotype>(XTAL_REF_(o)).pushdex;
			return S_::template element_f<I%size + size>(XTAL_REF_(o), XTAL_MOV_(i));
		}

		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		coelement_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return S_::template coelement_f<I>(XTAL_REF_(o));
		}
		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		coelement_f(auto &&o, integral_q auto &&i)
		noexcept -> decltype(auto)
		{
			return S_::template coelement_f<I>(XTAL_REF_(o), XTAL_REF_(i));
		}
		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		coelement_f(auto &&o, real_q auto k)
		noexcept -> decltype(auto)
		{
			using process::math::   term_f;
			using process::math:: square_f;
			using process::math::nearest_f;
			using K = XTAL_ALL_(k);
			auto constexpr two     =     K(2);
			auto constexpr four    =     K(4);
			auto constexpr quarter = one/K(4);

		//	TODO: Customize interpolation either within `atom` or `process`.
		//	TODO: Customize out-of-bounds handling.

			auto const k_dn = nearest_f<-1>(k), t_dn = k_dn - k;
			auto const k_up =       k_dn + one, t_up = k - k_up;
			auto const i_dn = static_cast<index_type>(k_dn);
			auto const i_up = static_cast<index_type>(k_up);
			auto const u_dn = element_f<I>(o, i_dn), w_dn = element_f<I>(o, i_dn - 1);
			auto const u_up = element_f<I>(o, i_up), w_up = element_f<I>(o, i_up + 1);
			auto const v_dn = term_f(u_up, term_f(u_up, u_dn - w_up, quarter), term_f(two, two, t_dn))*square_f(t_dn);
			auto const v_up = term_f(u_dn, term_f(u_dn, u_up - w_dn, quarter), term_f(two, two, t_up))*square_f(t_up);
			return v_dn + v_up;
		}

		XTAL_DEF_(let)
		push(auto &&e)
		noexcept -> auto
		{
			--pushdex;
			XTAL_IF0
			XTAL_0IF (1 == _std::popcount(size())) {
				pushdex &= mask;
			}
			XTAL_0IF (2 <= _std::popcount(size())) {
				pushdex += size&bond::math::bit_sign_f(pushdex);
			}
			return _std::exchange(S_::element(), XTAL_REF_(e));
		}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
