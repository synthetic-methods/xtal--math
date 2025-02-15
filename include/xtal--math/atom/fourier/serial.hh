#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math::fourier
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class        ..._s>	struct   serial;
template <class        ..._s>	using    serial_t = typename serial<_s...>::type;
template <class        ...Ts>	concept  serial_q = bond::tag_p<serial_t, Ts...>;
template <int N, class ...Ts>	concept  serial_p = serial_q<Ts...> and (...and (N == Ts::size()));

XTAL_FX0_(to) (template <auto f=_std::identity{}>
XTAL_DEF_(return,inline,let)
serial_f(auto &&...oo),
	_detail::factory<serial_t>::
		template make<f>(XTAL_REF_(oo)...))


////////////////////////////////////////////////////////////////////////////////
///\
Extends `grade` with multiplication via linear convolution. \

template <scalar_q ..._s> requires common_q<_s...>
struct serial<_s ...>
:	serial<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_q A>
struct serial<A>
{
	using _fit = bond::fit<A>;
	
	template <class T>
	using endotype = typename grade<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<serial_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// TYPE
		using typename S_::value_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::element;
		using S_::size;
		using S_::self;
		using S_::twin;
	
	public:// OPERATE
		using S_::operator*=;

		XTAL_DEF_(return,inline,let)  operator  * (auto const &                       t) const noexcept -> auto   {return twin() *=   t ;}
		XTAL_DEF_(inline,let)         operator  *=(_std::initializer_list<value_type> t)       noexcept -> auto & {return self() *= T(t);}

		///\
		Multiplication by linear convolution, truncated by `size`. \

		XTAL_DEF_(let)
		operator *=(T const &t)
		noexcept -> T &
		{
			auto &s = self();
			
			if constexpr (typename _fit::alignment{}() < size()) {
				for (auto i = size(); ~--i;) {element(i) *= get<0>(t);
				for (auto j =      i; j-- ;) {element(i) += t.element(j)*element(i - j);}}
			}
			else {
				bond::seek_backward_f<size(), 0>([&, this] (auto I) XTAL_0FN {get<I>(s) *= get<0>(t);
				bond::seek_backward_f<     I, 1>([&, this] (auto J) XTAL_0FN {get<I>(s) += get<J>(t)*get<I - J>(s);});});
			}
			return s;
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
