#pragma once
#include "./any.hh"

#include "./term.hh"
#include "./square.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Produces the dot-product of the given arguments. \
For single arguments, the dot-square or `norm` is produced. \

///\note\
The `M_alt` parameter determines the sign of the even terms. \

///\note\
For `M_alt=-1` and `complex` products, \
the `real` component w.r.t. multiplication is returned. \

///\todo\
Integrate or redefine with `ìmagine`. \

template <int M_alt=1> requires in_n<M_alt, 0, 1,-1>
struct  dot;

template <int M_alt=1> requires in_n<M_alt, 0, 1,-1>
using   dot_t = process::confined_t<dot<M_alt>>;

template <int M_alt=1>
XTAL_DEF_(return,inline,let)
dot_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return dot_t<M_alt>::method_f(XTAL_REF_(oo)...);
}

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

template <class ...Ts>
using multiplied_t = common_t<return_t<_std::multiplies<void>, Ts, Ts>...>;

template <class ...Ts>
using multiplied_u = multiplied_t<fluid_u<Ts>...>;


}///////////////////////////////////////////////////////////////////////////////

template <int M_alt> requires in_n<M_alt, 0, 1,-1>
struct dot
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x)
		noexcept -> auto
		requires un_n<fixed_shaped_q<decltype(x)>>
		{
			using _std::norm;// In case...
		//	static_assert(M_alt == 1);
			XTAL_IF1_(to) (norm(XTAL_REF_(x)))
			XTAL_0IF_(to) (square_f(abs(x)))
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x)
		noexcept -> auto
		requires in_n<fixed_shaped_q<decltype(x)>>
		{
			using _std::norm;// In case...
			using X = XTAL_ALL_(x);
			auto constexpr N = fixed_shaped<X>::extent();

			XTAL_IF0
			XTAL_0IF (complex_field_q<decltype(x)> and 1 == M_alt) {
				return norm(XTAL_REF_(x));
			}
			XTAL_0IF (1 == N) {
				auto const &[x0] = destruct_f(XTAL_REF_(x));
				return x0;
			}
			XTAL_0IF (2 == N) {
				auto const &[x0, x1] = destruct_f(XTAL_REF_(x));
				return term_f<M_alt, 2>(square_f(x0), x1);
			}
			XTAL_0IF (3 <= N) {
				_detail::multiplied_u<X> w{0};
				bond::seek_out_f<N>([&]<constant_q I> (I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 < M_alt) {w = term_f<               2>(XTAL_MOV_(w), get<I{}>(x));}
					XTAL_0IF (M_alt < 0) {w = term_f<cosign_v<I{}>, 2>(XTAL_MOV_(w), get<I{}>(x));}
				});
				return w;
			}
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x, auto &&y)
		noexcept -> auto
		requires in_n<fixed_shaped_q<decltype(x), decltype(y)>>
		{
			using X = XTAL_ALL_(x);
			using Y = XTAL_ALL_(y);
			auto constexpr N = fixed_shaped<X>::extent();

			XTAL_IF0
			XTAL_0IF (1 == N) {
				auto const &[x0] = destruct_f(XTAL_REF_(x));
				auto const &[y0] = destruct_f(XTAL_REF_(y));
				return x0*y0;
			}
			XTAL_0IF (2 == N) {
				auto const &[x0, x1] = destruct_f(XTAL_REF_(x));
				auto const &[y0, y1] = destruct_f(XTAL_REF_(y));
				return term_f<M_alt>(x0*y0, x1,y1);
			}
			XTAL_0IF (3 <= N) {
				_detail::multiplied_u<X, Y> w{0};
				bond::seek_out_f<N>([&]<constant_q I> (I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 < M_alt) {w = term_f<               1>(XTAL_MOV_(w), get<I{}>(x), get<I{}>(y));}
					XTAL_0IF (M_alt < 0) {w = term_f<cosign_v<I{}>, 1>(XTAL_MOV_(w), get<I{}>(x), get<I{}>(y));}
				});
				return w;
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
