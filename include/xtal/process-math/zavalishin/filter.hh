#pragma once
#include "./any.hh"

#include "./sigmoid.hh"
#include "../root.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> struct   filter;
template <typename ...As> using    filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct filter
:	process::link<filter<>, As...>
{
};
template <>
struct filter<>
{
	using  alpha_type = typename bond::operating::alpha_type;
	using  aphex_type = typename bond::operating::aphex_type;

	using   zoom_type = occur::inferred_t<struct   ZOOM, float>;
	using  order_type = occur::inferred_t<struct  ORDER, bond::seek_t<0, 1, 2, 3>>;
	using  limit_type = occur::inferred_t<struct  LIMIT, bond::seek_t<0, 1>>;
	using select_type = occur::inferred_t<struct SELECT, bond::seek_t<0, 1>>;

	XTAL_SET N_order = order_type::cardinality() - 1;
	XTAL_SET N_cache = N_order*2;

//	TODO: Implement `maximum`? \

//	NOTE: Assuming a base-type of `std::complex`, \
	the minimum storage required is `(2)*(2*8) == 32` bytes-per-pole. \

	//\
	using subject = provision::context<void
	using subject = bond::compose<void
	,	typename   zoom_type::template   attach<>
	,	typename  limit_type::template dispatch<>
	,	typename  order_type::template dispatch<>
	,	typename select_type::template dispatch<>
	>;
	using superkind = bond::compose<void
	,	provision::cached<aphex_type[N_cache]>
	,	provision::example<>
	,	subject
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		XTAL_DEF_(short)
		XTAL_LET infuse(auto &&o)
		noexcept -> signed
		{
			XTAL_IF0
			XTAL_0IF (is_q<decltype(o), order_type>) {
				S_::cache(0);
			}
			return S_::infuse(XTAL_REF_(o));
		}
		template <int N_lim=0, int N_ord=0, int N_kit=0>
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&x
		,	real_number_q auto &&o
		,	real_number_q auto &&a
		)
		noexcept -> XTAL_ALL_(x)
		{
			using X  = XTAL_ALL_(x); using _op = bond::operate<X>;
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return XTAL_REF_(x);
			}
			XTAL_0IF (0 == N_kit) {
				auto constexpr N  = N_ord + 1;
				auto constexpr N_ = N_ord + 0;

				using A  = algebra::sector_t<typename _op::alpha_type[N ]>;
				using A_ = algebra::sector_t<typename _op::alpha_type[N_]>;
				using W  = algebra::sector_t<X[N ]>;
				using W_ = algebra::sector_t<X[N_]>;

				auto constexpr _1 = _op::alpha_1;
				auto constexpr _2 = _op::alpha_2;

				XTAL_IF0
				XTAL_0IF (1 == N_ord) {
					A_ const a_{_1};
					W  const w = method<N_lim, N_ord, N_kit>(XTAL_REF_(x), XTAL_REF_(o), a_);
					return get<0>(w);
				}
				XTAL_0IF (2 == N_ord) {
					auto const a1 = _2*a;
					A_ const a_{_1, a1};
					W  const w = method<N_lim, N_ord, N_kit>(XTAL_REF_(x), XTAL_REF_(o), a_);
					return get<0>(w);
				}
				XTAL_0IF (3 == N_ord) {
					auto const a1 = term_f(_1, _2, a);
					A_ const a_{_1, a1, a1};
					W  const w = method<N_lim, N_ord, N_kit>(XTAL_REF_(x), XTAL_REF_(o), a_);
					return get<0>(w);
				}
			}
			XTAL_0IF (1 == N_kit) {
				return XTAL_REF_(x);
			}
		}
		template <int N_lim=0, int N_ord=0, int N_kit=0> requires (1 <= N_ord and N_kit == 0)
		XTAL_DEF_(inline)
		XTAL_LET method(auto &&x
		,	real_number_q auto const &o
		,	algebra::sector_q auto const &a_
		)
		//\
		noexcept -> algebra::sector_t<XTAL_ALL_(x)[N_ord + 1]>
		noexcept -> auto
		{
			auto constexpr N  = N_ord + 1;
			auto constexpr N_ = N_ord + 0;
			auto constexpr M_ = N_ord - 1;

			using X  = XTAL_ALL_(x); using _op = bond::operate<X>;
			using W  = algebra::sector_t<X[N ]>;
			using W_ = algebra::sector_t<X[N_]>;

		//	static_assert(is_q<W , decltype(a )>);
			static_assert(is_q<W_, decltype(a_)>);

		//	algebra::series_t<X[N_]> z_(self().template head<zoom_type>());
		//	auto [s_, u_] = S_::template cache<W_, W_>();//NOTE: Can't access from lambda...
			auto cache = S_::template cache<W_, W_>();
			W w;
			auto &w_ = reinterpret_cast<W_ &>(w);
			auto &u_ = get<0>(cache);
			auto &s_ = get<1>(cache);

			auto &u0 = get<0>(u_), uN = _op::alpha_1;
			auto &w0 = get<0>(w), &wN = get<N_>(w);

		//	Initialize coefficients `u*` and voltages `w*`:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				auto &uI = get<I>(u_);
				uI += uI == X{};// Force non-zero... (annoying)
				uI *= get<I>(a_);
			});
		//	u_ *= a_;
			w0  = u0;
			bond::seek_forward_f<M_>([&] (auto I) XTAL_0FN {
				get<I + 1>(w_) = term_f(get<I + 1>(u_), get<I>(w_), o);
			});
			wN = root_f<-1, (4)>(term_f(uN, get<M_>(w_), o));

		//	Integrate states `s*` with coefficients/voltages:
			bond::seek_forward_f<1 + N_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {wN *= x - w_.product(s_);}
				XTAL_0IF (1 <= K) {wN  = x - w_.product(u_);}
				bond::seek_backward_f<N_>([&] (auto I) XTAL_0FN {
					auto const  &aI = get<I>(a_);
					auto const  &sI = get<I>(s_);
					auto const   wI = term_f(sI, get<I + 1>(w), o);
					auto const   uI = sigmoid_t<-1,-1>::template function<N_lim>(wI, aI);//, get<I>(z_));
					get<I>(w_) = wI;
					get<I>(u_) = uI;
				});
			});
		//	Finalize states and voltages/coefficients:
			bond::seek_forward_f<N_>([&] (auto I) XTAL_0FN {
				get<I>(s_) = term_f(-get<I>(s_), get<I>(w_), _op::alpha_2);
			});
		//	s_  = term_f(-s_, w_, _op::alpha_2);
			u_ /= a_;
		//	w_ *= u_;// TODO: Make this line optional?

			return w;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
