#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::lambert
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_car=0> struct   sine {static_assert(M_ism);};
template <int M_ism=1, int N_car=0> using    sine_t = process::confined_t<sine<M_ism, N_car>>;

template <int M_ism>
struct sine<M_ism, -0>
{
	using superkind = bond::compose<discarded<1>, sine<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;
	//	using S_::function;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			return function<N_lim>(XTAL_REF_(o), _op::alpha_1);
		}
		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o, auto &&a)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;

			XTAL_LET N = below_m<(1<<2), (unsigned) N_lim>;
			XTAL_IF0
			XTAL_0IF (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
			//	XTAL_0IF (M_ism == -1) {return asin (XTAL_REF_(o));}
			//	XTAL_0IF (M_ism ==  1) {return  sin (XTAL_REF_(o));}
				XTAL_0IF (M_ism == -2) {return asinh(XTAL_REF_(o));}
				XTAL_0IF (M_ism ==  2) {return  sinh(XTAL_REF_(o));}
			}
			XTAL_0IF (M_ism == -2) {
				auto const u = o + root_f<2>(term_f(_op::alpha_1, o, o));
				return roots_f<(1 << N)>(u).template sum<-1>()*_op::diplo_f(N - 1);
			}
			XTAL_0IF (M_ism == +2) {
				return S_::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism,-2>
struct sine<M_ism, -1>
{
	using supertype = sine_t<M_ism, 0>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			return function<N_lim>(XTAL_REF_(o), _op::alpha_1);
		}
		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o, auto &&a)
		noexcept -> decltype(auto)
		{
			return supertype::template function<N_lim>(XTAL_REF_(o), XTAL_REF_(a))*root_f<-1, 4>(o);
		}

	};
};
template <int M_ism> requires in_n<M_ism, 2>
struct sine<M_ism, -1>
{
	using superkind = bond::compose<discarded<2>, sine<M_ism, -2>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;
	//	using S_::function;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o, auto &&a)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			XTAL_IF0
			XTAL_0IF (N_lim < 0 or 3 < N_lim) {
				auto const t = term_f(o, a, a);
				auto const s = o + _op::alpha_1;

				auto const [u, n] = roots_f<1, 4>(XTAL_REF_(o));
				return root_f<2>(t)*root_f<-2>(s)*sinh(u)*n;
			}
			XTAL_0IF_(else) {
				return S_::template function<N_lim>(XTAL_REF_(o), XTAL_REF_(a));
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			return function<N_lim>(XTAL_REF_(o), _op::alpha_1);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 2>
struct sine<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o, auto &&a)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;

			XTAL_LET N = below_m<(1<<2), (unsigned) N_lim>;
			auto const w = XTAL_REF_(o)*_op::haplo_f(N << 1);
			auto const t = term_f(w, a, a);
			auto const s = w + _op::alpha_1;

			XTAL_IF0
			XTAL_0IF (N == 0) {
				return a;
			}
			XTAL_0IF (N == 1) {
				return root_f<2>(t);
			}
			XTAL_0IF (N == 2) {
				return root_f<2>(t)*
					term_f(_op::alpha_1, _op::diplo_f(1), w);
			}
			XTAL_0IF (N == 3) {
				return root_f<2>(t)*
					term_f(_op::alpha_1, _op::diplo_f(1),   w)*
					term_f(_op::alpha_1, _op::diplo_f(3), t*w);
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			return function<N_lim>(XTAL_REF_(o), _op::alpha_1);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
