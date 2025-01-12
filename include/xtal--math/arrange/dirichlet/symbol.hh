#pragma once
#include "./any.hh"

#include "../../process/pade/unity.hh"




XTAL_ENV_(push)
namespace xtal::arrange::math::dirichlet
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct   symbol;
template <class   ..._s>	using    symbol_t = typename symbol<_s...>::type;
template <class   ...Ts>	concept  symbol_q = bond::tag_p<symbol_t, Ts...>;
template <class  V=void>
XTAL_DEF_(short)
XTAL_LET symbol_f(auto &&...oo)
noexcept -> auto
{
	return _detail::initialize<symbol_t>::template via<V>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////
///\
Extends `couple` with Dirichlet characterization and modulo access. \

template <vector_q A>
struct symbol<A>
{
	using _op = bond::operate<A>;
	
	template <class T>
	using endotype = typename couple<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<symbol_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		friend T;
		using  S_ = holotype<T>;
		using  I_ = typename S_::difference_type;

	protected:
		using          S_::N_data;
		using typename S_::U_data;

		static constexpr auto modulo = [] (I_ i) XTAL_0FN_(((i%N_data) + N_data)%N_data);

	public:// ACCESS
		using S_::element;
		using S_::self;
		using S_::twin;

	public:// OPERATE

		XTAL_DEF_(short) XTAL_LET element(I_ i) const &&noexcept -> decltype(auto) {return XTAL_MOV_(S_::operator[](modulo(i)));}
		XTAL_DEF_(short) XTAL_LET element(I_ i) const  &noexcept -> decltype(auto) {return           S_::operator[](modulo(i)) ;}

		XTAL_DEF_(short) XTAL_LET element(I_ i)       &&noexcept -> decltype(auto) {return XTAL_MOV_(S_::operator[](modulo(i)));}
		XTAL_DEF_(short) XTAL_LET element(I_ i)        &noexcept -> decltype(auto) {return           S_::operator[](modulo(i)) ;}


	public:// CONSTRUCT
		using S_::S_;

		///\
		Dirichlet character generation. \

		template <int N_subscript=1> requires ((bool) (1&N_data))
		XTAL_LET characterize()
		noexcept -> T &
		{
			using namespace process::math;
			extent_type constexpr N = N_data;
			extent_type constexpr M = N_data - 1;
			extent_type constexpr K = M >> 1U;
			extent_type           k = N_subscript;
			element(0) = {};

			if constexpr (integral_variable_q<U_data>) {
				bond::seek_forward_f<K>([&, this] (auto i) XTAL_0FN {
					auto const o = k%N;
					element(    o) =  i;
					element(N - o) =  i - K;
					k *= K;
				});
				element(1) = 0;
			}
			else {
				U_data w =  1;
				U_data u = -1;
				if constexpr (complex_field_q<U_data>) {
					u = pade::unity_t<1>::template function<6>(_op::ratio_f(1, 2*K));
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
		XTAL_LET subcharacterize()
		noexcept -> T &
		{
			using namespace process::math;
			auto constexpr N = N_data*2 + 1;
			auto constexpr M = N_data*2 + 0;
			auto constexpr K = N_data;
			auto           k = N_data;

			if constexpr (integral_variable_q<U_data>) {
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
				U_data w, u;
				if constexpr (complex_field_q<U_data>) {
					u = pade::unity_t<1>::template function<6>(_op::ratio_f(1, 2*K));
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
