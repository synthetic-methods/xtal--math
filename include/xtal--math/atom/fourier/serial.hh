#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math::fourier
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class        ..._s>	struct  serial;
template <class        ..._s>	using   serial_t = typename serial<_s...>::type;
template <class        ...Ts>	concept serial_q = bond::classify_tag_p<serial_t, Ts...>;
template <int N, class ...Ts>	concept serial_p = serial_q<Ts...> and (...and (N == Ts::size()));

XTAL_VAL_(let) serial_f = [] XTAL_1FN_(call) (_detail::factory<serial_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `differential` with multiplication via linear convolution.
*/
template <class ...Us> requires un_v<bond::devise_condensed_p<        Us...>>
struct serial<Us...> :               bond::devise_condensed_s<serial, Us...>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and un_q<xtd::plus<>, Us...>
struct serial<Us...> : serial<Us..., xtd::plus<>>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and in_q<xtd::plus<>, Us...>
struct serial<Us...>
{
	template <class T>
	using endotype = typename differential<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<serial_t>>;

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

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;

	public:// OPERATE
		using S_::operator*=;

		XTAL_VAL_(return,inline,let)  operator  * (auto const &                      t) const noexcept -> auto   {return twin() *=   t ;}
		XTAL_VAL_(inline,let)         operator  *=(std::initializer_list<U_> t)       noexcept -> auto & {return self() *= T(t);}

		/*!
		\brief   Multiplication by linear convolution, truncated by `size`.
		*/
		XTAL_VAL_(let)
		operator *=(T const &t)
		noexcept -> T &
		{
			auto &s = self();

			if constexpr (V_fit::alignment < size) {
				for (auto i = size(); ~--i;) {element(i) *= get<0>(t);
				for (auto j =      i; j-- ;) {element(i) += t.element(j)*element(i - j);}}
			}
			else {
				int constexpr N{size};
				bond::seek_to_e<-N, 0>([&, this] (auto I) XTAL_0FN -> void {get<I>(s) *= get<0>(t);
				bond::seek_to_e<-I, 1>([&, this] (auto J) XTAL_0FN -> void {get<I>(s) += get<J>(t)*get<I - J>(s);});});
			}
			return s;
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
