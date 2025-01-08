#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math::hadamard
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct   series;
template <class   ..._s>	using    series_t = typename series<_s...>::type;
template <class   ...Ts>	concept  series_q = bond::tag_p<series_t, Ts...>;
template <class  V=void>
XTAL_DEF_(short)
XTAL_LET series_f(auto &&...oo)
noexcept -> auto
{
	return _detail::initialize<series_t>::template via<V>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////
///\
Extends `atom::order` with point-wise multiplication. \

template <vector_q A>
struct series<A>
{
	using A_op = bond::operate<A>;
	using A_sigma = typename A_op::sigma_type;
	using A_alpha = typename A_op::alpha_type;
	using A_aphex = typename A_op::aphex_type;
	
	template <class T>
	using endotype = typename atom::order<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<series_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using  S_ = holotype<T>;
	
	protected:
		using          S_::N_data;
		using typename S_::U_data;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::self;
		using S_::twin;

	public:// OPERATE
		using S_::operator*=;
		using S_::operator/=;
		using S_::operator%=;

		XTAL_DEF_(short)  XTAL_LET operator * (auto const &t)              const noexcept -> auto   {return S_::twin() *=   t ;}
		XTAL_DEF_(short)  XTAL_LET operator / (auto const &t)              const noexcept -> auto   {return S_::twin() /=   t ;}
		XTAL_DEF_(short)  XTAL_LET operator % (auto const &t)              const noexcept -> auto   {return S_::twin() %=   t ;}

		XTAL_DEF_(inline) XTAL_LET operator *=(_std::initializer_list<U_data> t) noexcept -> auto & {return S_::self() *= T(t);}
		XTAL_DEF_(inline) XTAL_LET operator /=(_std::initializer_list<U_data> t) noexcept -> auto & {return S_::self() /= T(t);}
		XTAL_DEF_(inline) XTAL_LET operator %=(_std::initializer_list<U_data> t) noexcept -> auto & {return S_::self() %= T(t);}

		XTAL_DEF_(inline)
		XTAL_LET operator *=(array_q<N_data> auto const &t)
		noexcept -> T &
		{
			return S_::template pointwise<[] (auto &u, auto const &v)
				XTAL_0FN {u *= v;}>(XTAL_REF_(t));
		}
		XTAL_DEF_(inline)
		XTAL_LET operator /=(array_q<N_data> auto const &t)
		noexcept -> T &
		{
			return S_::template pointwise<[] (auto &u, auto const &v)
				XTAL_0FN {u /= v;}>(XTAL_REF_(t));
		}
		XTAL_DEF_(inline)
		XTAL_LET operator %=(array_q<N_data> auto const &t)
		noexcept -> T &
		{
			return S_::template pointwise<[] (auto &u, auto const &v)
				XTAL_0FN {u %= v;}>(XTAL_REF_(t));
		}

	};
	using type = bond::isotype<homotype>;

};
static_assert(atomic_q<series_t<float[2]>>);

static_assert(fungible_q<_std::array<float, 2>,
	XTAL_ALL_(XTAL_ANY_(series_t<float(&)[2]>)*XTAL_ANY_(series_t<float(&)[2]>))>
);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
