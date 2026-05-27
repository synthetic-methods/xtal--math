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
\brief   Extends `quantity_plus` with Dirichlet characterization and modulo access.
*/
template <scalar_array_q ..._s> requires same_q<_s...>
struct loop<_s ...>
:	loop<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_array_q A>
struct loop<A>
{
private:
	using U_fit = bond::fit<A>;
	
	template <class T>
	using endotype = typename quantity_plus<A>::template homotype<T>;

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
		typename S_::index_type m_index{};
		typename S_::index_type m_limit{};

	public:// ACCESS
	//	using S_::element_f;
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

	protected:
		template <index_type I=0>
		XTAL_DEF_(inline,set)
		wrap_e(index_type &i)
		noexcept -> void
		{
			i += I;
			XTAL_IF0
			XTAL_0IF (1 == std::popcount(size())) {
				i &= mask;
			}
			XTAL_0IF (2 <= std::popcount(size())) {
				i %= size;
				i += size;
				i %= size;
			}
		}
		template <index_type I=0>
		XTAL_DEF_(return,inline,set)
		wrap_f(index_type  i)
		noexcept -> decltype(auto)
		{
			(void) wrap_e<I>(i); return i;
		}

	public:
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
			i += qualify_f<homotype>(o).m_index;
			return S_::element_f(XTAL_REF_(o), wrap_f(XTAL_MOV_(i)));
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
		coelement_f(auto &&o, real_q auto h)
		noexcept -> decltype(auto)
		{
			using process::math::   term_f;
			using process::math:: square_f;
			using process::math::nearest_f;
			using H = XTAL_ALL_(h);
			auto constexpr two     =     H(2);
			auto constexpr four    =     H(4);
			auto constexpr quarter = one/H(4);

		//	TODO: Customize interpolation either within `atom` or `process`.
		//	TODO: Customize out-of-bounds handling.

			auto const k_dn = nearest_f<-1>(h), t_dn = k_dn - h;
			auto const k_up =       k_dn + one, t_up = h - k_up;
			auto const i_dn = static_cast<index_type>(k_dn);
			auto const i_up = static_cast<index_type>(k_up);
			auto const u_dn = element_f<I>(o, i_dn), w_dn = element_f<I>(o, i_dn - 1);
			auto const u_up = element_f<I>(o, i_up), w_up = element_f<I>(o, i_up + 1);
			auto const v_dn = term_f(u_up, term_f(u_up, u_dn - w_up, quarter), term_f(two, two, t_dn))*square_f(t_dn);
			auto const v_up = term_f(u_dn, term_f(u_dn, u_up - w_dn, quarter), term_f(two, two, t_up))*square_f(t_up);
			return v_dn + v_up;
		}

		XTAL_DEF_(inline,let)
		peek(atom::quantify_q auto const &h_)
		noexcept -> auto
		{
			using H_ = objective_t<decltype(h_)>;
			return [&]<auto ...I> (bond::seek_in_t<I...>)
			XTAL_0FN_(to) (H_{peek(get<I>(h_))...})
				(bond::seek_to_t<H_::size()>{});
		}
		XTAL_DEF_(inline,let)
		peek(real_q auto h)
		noexcept -> auto
		{
			return S_::coelement(h);
		}
		XTAL_DEF_(inline,let)
		push(auto &&e)
		noexcept -> auto
		{
			--m_index; (void) wrap_e(m_index);
			return std::exchange(S_::element(), XTAL_REF_(e));
		}

	};
	using type = bond::derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
