#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `quantity`-addition with `operator*` defined by the scalar product.

Indended to act as a coefficient of a similar type where a scalar result is required.

\todo    Either define `std::complex` construction/operation,
or create a similar complex sentinel with multiplication/projection.

\todo    Specialize `plus_multiplies` or `fma`?
*/

template <class ..._s>	struct  dot;
template <class ..._s>	using   dot_t = typename dot<_s...>::type;
template <class ..._s>	concept dot_q = bond::classify_fixed_tag_p<dot_t, _s...>;

XTAL_VAL_(let) dot_f = [] XTAL_1FN_(call) (_detail::factory<dot_t>::make);


////////////////////////////////////////////////////////////////////////////////

template <class ...Us> requires un_v<bond::devise_condensed_p<     Us...>>
struct dot<Us...> :                  bond::devise_condensed_s<dot, Us...>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and un_q<xtd::plus<>, Us...>
struct dot<Us...> : dot<Us..., xtd::plus<>>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and in_q<xtd::plus<>, Us...>
and un_v<true, extra_matrix_q<Us>...>
struct dot<Us...>
{
private:
	template <class T>
	using endotype = typename quantity<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<dot_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_    = holotype<T>;
		using A_    = typename S_:: archetype;
		using I_    = typename S_::index_type;
		using U_    = typename S_::value_type;
		using V_    = typename S_::scale_type;
		using V_fit = bond::fit<V_>;
		static_assert(fixed_q<A_>);

	public:// CONSTRUCT
		using S_::S_;

		XTAL_VAL_(new,implicit)
		homotype(std::initializer_list<U_> xs)
		noexcept
		requires complete_q<U_>
		:	S_(xs)
		{
		}

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

	public:// OPERATE
		using S_::product;

		XTAL_VAL_(return,inline,let)
		product() const
		noexcept -> auto
		requires dot_q<U_> and in_v<size, 2>
		{
			auto &s = self();
			return get<0>(s) * get<1>(s);
		}

		XTAL_VAL_(return,inline,let)
		operator() () const noexcept {product();}

		XTAL_VAL_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_unless) (t.size())
		{
			return S_::operator*(t);
		}
		XTAL_VAL_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_if) (t.size()) and XTAL_TRY_(to_if)     (t.capacity())
		{
			auto &s = self();
			using U = XTAL_ALL_(s[0]*t[0]);
			auto  u = U{};
			auto  n = bond::math::bit_extremum_f<-1>(t.size(), size());
			for (int i{}; i < n; ++i) {
				u = xtd::plus_multiplies{} (XTAL_MOV_(u), s[i], t[i]);
			}
			return u;
		}
		XTAL_VAL_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		requires XTAL_TRY_(to_if) (t.size()) and XTAL_TRY_(to_unless) (t.capacity())
		{
			auto &s = self();
			using U = XTAL_ALL_(get<0>(s)*get<0>(t));
			U u{};
			bond::seek_to_e<size>([&]<constant_q I> (I) XTAL_0FN -> void {
				u = xtd::plus_multiplies{} (XTAL_MOV_(u), got<I{}>(s), got<I{}>(t));
			});
			return u;
		}
		template <class U>
		XTAL_VAL_(return,inline,met)
		operator * (U const &u, homotype const &s)
		noexcept
		requires un_v<std::derived_from<U, homotype>>
		{
			return s.operator*(u);
		}

	};
	using type = bond::derive_t<homotype>;

};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
