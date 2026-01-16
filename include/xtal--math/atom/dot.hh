#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `group`-addition with `operator*` defined by the scalar product.

Indended to act as a coefficient of a similar type where a scalar result is required.

\todo    Either define `std::complex` construction/operation,
or create a similar complex sentinel with multiplication/projection.

\todo    Specialize `plus_multiplies_f` or `fma`?
*/

template <class ..._s>	struct  dot;
template <class ..._s>	using   dot_t = typename dot<_s...>::type;
template <class ..._s>	concept dot_q = bond::tag_infixed_p<dot_t, _s...>;

XTAL_DEF_(let) dot_f = [] XTAL_1FN_(call) (_detail::factory<dot_t>::make);


////////////////////////////////////////////////////////////////////////////////

template <scalar_q ..._s> requires same_q<_s...>
struct dot<_s ...>
:	dot<common_t<_s...>[sizeof...(_s)]>
{
};
template <matrix_q A>
struct dot<A>
:	dot<dot_t<typename fixed<A>::value_type>[fixed<A>::extent()]>
{
};
template <class ..._s>
struct dot
{
private:
	template <class T>
	//\
	using endotype = typename group_addition<_s...>::template homotype<T>;
	using endotype = typename group<wrap_s<_s, _std::plus>...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<dot_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

		using typename S_::revalue_type;
		using typename S_::  value_type;
		using typename S_::  scale_type;

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

	public:// OPERATE

		XTAL_DEF_(return,inline,let)
		product() const
		noexcept -> auto
		{
			static_assert(dot_q<value_type> and size == 2);
			auto &s = self();
			return get<0>(s) * get<1>(s);
		}

		XTAL_DEF_(return,inline,let)
		operator() () const noexcept {product();}

		XTAL_DEF_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_unless) (t.size())
		{
			return S_::operator*(t);
		}
		XTAL_DEF_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_if) (t.size()) and XTAL_TRY_(to_if)     (t.capacity())
		{
			auto &s = self();
			auto  u = revalue_type{};
			auto  n = bond::math::bit_extremum_f<-1>(t.size(), size());
			for (int i{}; i < n; ++i) {
				u = _xtd::plus_multiplies_f(XTAL_MOV_(u), s[i], t[i]);
			}
			return u;
		}
		XTAL_DEF_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_if) (t.size()) and XTAL_TRY_(to_unless) (t.capacity())
		{
			auto &s = self();
			revalue_type u{};
			bond::seek_until_f<size>([&]<constant_q I> (I) XTAL_0FN {
				u = _xtd::plus_multiplies_f(XTAL_MOV_(u), got<I{}>(s), got<I{}>(t));
			});
			return u;
		}
		template <class U>
		XTAL_DEF_(return,inline,met)
		operator * (U const &u, homotype const &s)
		noexcept
		requires un_v<_std::derived_from<U, homotype>>
		{
			return s.operator*(u);
		}

	};
	using type = bond::derive_t<homotype>;

};
template <scalar_q U>
struct dot<U>
:	dot<U[2]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
