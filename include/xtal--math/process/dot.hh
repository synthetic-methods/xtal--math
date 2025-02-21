#pragma once
#include "./any.hh"

#include "./term.hh"




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


////////////////////////////////////////////////////////////////////////////////

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
		{
			if constexpr (complex_field_q<decltype(x)> and 1 == M_alt) {
				return norm(XTAL_REF_(x));
			}
			else {
			//	static_assert(fixed_shaped_q<decltype(x)>);
			//	auto const &[x0, x1] = destruct_f(XTAL_REF_(x));
			//	return term_f<M_alt>(x0*x0, x1,x1);
				return method_f(x, x);
			}
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x, auto &&y)
		noexcept -> auto
		{
			using X = XTAL_ALL_(x);
			using Y = XTAL_ALL_(y);
			static_assert(fixed_shaped_q<X, Y>);
			auto constexpr N = fixed_shaped<X>::extent();
			if constexpr (N == 2) {
				auto const &[x0, x1] = destruct_f(XTAL_REF_(x));
				auto const &[y0, y1] = destruct_f(XTAL_REF_(y));
				return term_f<M_alt>(x0*y0, x1,y1);
			}
			using value_type = fixed_valued_u<X, Y>;
			using scale_type = absolve_u<value_type>;
			value_type u{0};
			
			bond::seek_out_f<N>([&]<constant_q I> (I) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 < M_alt) {u = term_f(                    XTAL_MOV_(u), get<I{}>(x), get<I{}>(y));}
				XTAL_0IF (M_alt < 0) {u = term_f<-sign_v<I{}&1, -1>>(XTAL_MOV_(u), get<I{}>(x), get<I{}>(y));}
			});
			
			return u;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
