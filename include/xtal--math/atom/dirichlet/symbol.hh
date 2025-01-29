#pragma once
#include "./any.hh"

#include "../../process/pade/unity.hh"




XTAL_ENV_(push)
namespace xtal::atom::math::dirichlet
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct   symbol;
template <class   ..._s>	using    symbol_t = typename symbol<_s...>::type;
template <class   ...Ts>	concept  symbol_q = bond::tag_p<symbol_t, Ts...>;


XTAL_FX0_(alias) (template <auto f=_std::identity{}>
XTAL_DEF_(return,inline,let)
symbol_f(auto &&...oo),
	_detail::factory<symbol_t>::
		template make<f>(XTAL_REF_(oo)...))


////////////////////////////////////////////////////////////////////////////////
///\
Extends `couple` with Dirichlet characterization and modulo access. \

template <scalar_q ..._s> requires common_q<_s...>
struct symbol<_s ...>
:	symbol<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_q A>
struct symbol<A>
{
	using _fix = bond::fixture<A>;
	
	template <class T>
	using endotype = typename couple<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<symbol_t>>;

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

		XTAL_DEF_(return,inline,let) element(I_ i) const &&noexcept -> decltype(auto) {return XTAL_MOV_(S_::operator[](process::math::modulo_f<size>(i)));}
		XTAL_DEF_(return,inline,let) element(I_ i) const  &noexcept -> decltype(auto) {return           S_::operator[](process::math::modulo_f<size>(i)) ;}

		XTAL_DEF_(return,inline,let) element(I_ i)       &&noexcept -> decltype(auto) {return XTAL_MOV_(S_::operator[](process::math::modulo_f<size>(i)));}
		XTAL_DEF_(return,inline,let) element(I_ i)        &noexcept -> decltype(auto) {return           S_::operator[](process::math::modulo_f<size>(i)) ;}


	public:// CONSTRUCT
		using S_::S_;

		///\
		Dirichlet character generation. \

		template <int N_subscript=1> requires ((bool) (1&size))
		XTAL_DEF_(let)
		characterize()
		noexcept -> T &
		{
			using namespace process::math;
			extent_type constexpr N = size;
			extent_type constexpr M = size - 1;
			extent_type constexpr K = M >> 1U;
			extent_type           k = N_subscript;
			element(0) = {};

			if constexpr (integral_variable_q<value_type>) {
				bond::seek_forward_f<K>([&, this] (auto i) XTAL_0FN {
					auto const o = k%N;
					element(    o) =  i;
					element(N - o) =  i - K;
					k *= K;
				});
				element(1) = 0;
			}
			else {
				value_type w =  1;
				value_type u = -1;
				if constexpr (complex_field_q<value_type>) {
					u = pade::unity_t<1>::template static_method<6>(_fix::ratio_f(1, 2*K));
				}
				bond::seek_forward_f<K>([&, this] (auto i) XTAL_0FN {
					auto const o = k%N;
					element(    o) =  w;
					element(N - o) = -w;
					w *= u;
					k *= K;
				});
			}
			return self();
		}
		template <int N_subscript=1>
		XTAL_DEF_(let)
		subcharacterize()
		noexcept -> T &
		{
			using namespace process::math;
			extent_type constexpr N = size*2 + 1;
			extent_type constexpr M = size*2 + 0;
			extent_type constexpr K = size;
			extent_type           k = size;

			if constexpr (integral_variable_q<value_type>) {
				bond::seek_forward_f<K>([&, this] (auto i) XTAL_0FN {
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
					u = pade::unity_t<1>::template static_method<6>(_fix::ratio_f(1, 2*K));
				}
				else {
					u = 1;
				}
				w = u;
				bond::seek_forward_f<K>([&, this] (auto i) XTAL_0FN {
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
	using type = derive_t<homotype>;

};

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
