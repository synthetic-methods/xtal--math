#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::arrange::math::fourier
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class        ..._s>	struct   serial;
template <class        ..._s>	using    serial_t = typename serial<_s...>::type;
template <class        ...Ts>	concept  serial_q = bond::tag_p<serial_t, Ts...>;
template <int N, class ...Ts>	concept  serial_p = serial_q<Ts...> and (...and (N == Ts::size()));
template <class  V=void>
XTAL_DEF_(return,inline,let)
serial_f(auto &&...oo)
noexcept -> auto
{
	return _detail::initialize<serial_t>::template via<V>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////
///\
Extends `grade` with multiplication via linear convolution. \

template <vector_q A>
struct serial<A>
{
	using _fix = bond::fixture<A>;
	
	template <class T>
	using endotype = typename grade<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<serial_t>>;

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
		using S_::element;
		using S_::self;
		using S_::twin;
	
	public:// OPERATE
		using S_::operator*=;

		XTAL_DEF_(return,inline,let)  operator  * (auto           const &t) const noexcept -> auto   {return twin() *=   t ;}
		XTAL_DEF_(inline,let) operator  *=(initializer_s<U_data> t)       noexcept -> auto & {return self() *= T(t);}

		///\
		Multiplication by linear convolution, truncated by `N_data`. \

		XTAL_DEF_(let)
		operator *=(T const &t)
		noexcept -> T &
		{
			auto &s = self();
			
			if constexpr (_fix::alignment::value < N_data) {
				for (auto i = N_data; ~--i;) {element(i) *= get<0>(t);
				for (auto j = i; j-- ;) {element(i) += t.element(j)*element(i - j);}}
			}
			else {
				bond::seek_backward_f<N_data, 0>([&, this] (auto I) XTAL_0FN {get<I>(s) *= get<0>(t);
				bond::seek_backward_f<     I, 1>([&, this] (auto J) XTAL_0FN {get<I>(s) += get<J>(t)*get<I - J>(s);});});
			}
			return s;
		}

	};
	using type = derive_t<homotype>;

};
static_assert(atomic_q<serial_t<float[2]>>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
