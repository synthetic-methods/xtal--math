#pragma once
#include "./any.hh"

#include "../../process/pade/unity.hh"
#include "../../process/proot.hh"



XTAL_ENV_(push)
namespace xtal::atom::math::dirichlet
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  symbol;
template <class   ..._s>	using   symbol_t = typename symbol<_s...>::type;
template <class   ...Ts>	concept symbol_q = bond::classify_tag_p<symbol_t, Ts...>;

XTAL_VAL_(let) symbol_f = [] XTAL_1FN_(call) (_detail::factory<symbol_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `couple` with Dirichlet characterization and modulo access.
*/
template <class ...Us> requires un_v<bond::devise_condensed_p<        Us...>>
struct symbol<Us...> :               bond::devise_condensed_s<symbol, Us...>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and un_q<xtd::multiplies<>, Us...>
struct symbol<Us...> : symbol<Us..., xtd::multiplies<>>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and in_q<xtd::multiplies<>, Us...>
struct symbol<Us...>
{
private:
	template <class T>
	using endotype = typename couple<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<symbol_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// TYPE
		using typename S_::value_type;
		using typename S_::index_type;

	public:// ACCESS
		using S_::element;
		using S_::mask;
		using S_::size;
		using S_::self;
		using S_::twin;

		XTAL_VAL_(return,inline,set)
		deindex_f(index_type i)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (1 == std::popcount(size())) {
				i &= mask;
			}
			XTAL_0IF (2 <= std::popcount(size())) {
				i %= size;
				i += size;
				i %= size;
			}
			return i;
		}

	public:// CONSTRUCT
		using S_::S_;

		/*!
		\brief   Dirichlet character generation.
		*/
		template <int N_subscript=1> requires ((bool) (1&size))
		XTAL_VAL_(let)
		characterize()
		noexcept -> T &
		{
			using namespace process::math;
			using U_fit = bond::fit<value_type>;
			extent_type constexpr N = size;
			extent_type constexpr M = size - 1;
			extent_type constexpr K = M >> 1U;
			extent_type           k = N_subscript;
			element(0) = {};

			if constexpr (integral_variable_q<value_type>) {
				bond::seek_to_e<K>([&, this] (auto i) XTAL_0FN -> void {
					auto const o = k%N;
					element(    o) =  i;
					element(N - o) =  i - K;
					k *= process::math::proot_f<>(size);
				});
				element(1) = 0;
			}
			else {
				value_type w =  1;
				value_type u = -1;
				if constexpr (complex_field_q<value_type>) {
					u = process::math::pade::unity_f<1>(U_fit::ratio_f(1, 2*K));
				}
				bond::seek_to_e<K>([&, this] (auto i) XTAL_0FN -> void  {
					auto const o = k%N;
					element(    o) =  w;
					element(N - o) = -w;
					w *= u;
					k *= process::math::proot_f<>(size);
				});
			}
			return self();
		}
		template <int N_subscript=1>
		XTAL_VAL_(let)
		subcharacterize()
		noexcept -> T &
		{
			using namespace process::math;
			using U_fit = bond::fit<value_type>;
			extent_type constexpr N = size*2 + 1;
			extent_type constexpr M = size*2 + 0;
			extent_type constexpr K = size;
			extent_type           k = size;

			if constexpr (integral_variable_q<value_type>) {
				bond::seek_to_e<K>([&, this] (auto i) XTAL_0FN -> void  {
					auto const o = k%N;
					if (K < o) {
						element(M - o) = (1 + i) - K;
					}
					else {
						element(o - 1) = (1 + i);
					}
					k *= K;
				});
				element(0) = 0;
			}
			else {
				value_type w, u;
				if constexpr (complex_field_q<value_type>) {
					u = process::math::pade::unity_f<1>(U_fit::ratio_f(1, 2*K));
				}
				else {
					u = 1;
				}
				w = u;
				bond::seek_to_e<K>([&, this] (auto i) XTAL_0FN -> void  {
					auto const o = k%N;
					if (K < o) {
						element(M - o) = -w;
					}
					else {
						element(o - 1) =  w;
					}
					w *= u;
					k *= K;
				});
			}
			return self();
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
