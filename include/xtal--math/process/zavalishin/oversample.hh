#pragma once
#include "./any.hh"

#include "./reuse.hh"
#include "./filter.hh"
#include "../taylor/cosy.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Oversamples the super-filter.
*/
template <auto ...As>	struct  oversample;
template <auto ...As>	using   oversample_t = confined_t<oversample<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <variable_q auto M_sup>
struct oversample<M_sup>
{
private:
	XTAL_VAL_(set) M_len = one << M_sup;

public:
	template <class S>
	using ectotype = typename S::data_type::value_type;

	template <class S>
	using endotype = atom::bucket_t<ectotype<S>[M_len]>;

	template <class S>
	using endokind = scheme::stashed<endotype<S>>;

	template <class S>
	class subtype : public bond::compose_s<S, endokind<S>>
	{
	//	static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, endokind<S>>;
		using T_ = typename S_::self_type;
		using D_ = typename S_::data_type;

	public:// CONSTRUCT
		using S_::S_;

	protected:// OPERATE
		XTAL_VAL_(return,inline,set)
		F_lop(auto &&u_)
		noexcept -> auto
		{
			using U_ = XTAL_ALL_(u_);
			int constexpr I_lim = U_::size();
			int constexpr I_end = I_lim -  1;
			return [u_=XTAL_REF_(u_)]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (U_{zero, u_.template element<I_end - I>()...})
					(bond::seek_to_t<U_::size() - 1>{});
		}
		XTAL_VAL_(return,inline,set)
		F_cos(auto &&...oo)
		noexcept -> auto
		{
			auto const w = taylor::cosy_t<1, 2>{}.template method<2>((one *...* XTAL_REF_(oo)));
			return term_f(-one, two, XTAL_MOV_(w));
		}
		XTAL_VAL_(return,inline,set)
		F_man(auto &&...oo)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(oo)...>;
			//\
			auto constexpr alpha = U_fit::ratio_f(4, 25);
			auto constexpr alpha = U_fit::ratio_f(2, 13);
			auto constexpr  beta = root_f<-1>(one - alpha);
			auto constexpr     u = beta*one;
			auto constexpr     w = beta*two;
			auto               x = F_cos(XTAL_REF_(oo)...);
			return term_f(term_f(u, w, x), -two, x + one)*(x - one);
		}
		template <int N_len, int N_off=1>
		XTAL_VAL_(return,inline,set)
		F_men(auto &&u)
		noexcept -> auto
		{
			using   U  = XTAL_ALL_(u);
			using   U_ = atom::quantity_t<U[N_len], xtd::plus_multiplies<>>;
			return [u=XTAL_REF_(u)]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (U_{(F_man(I + N_off, u)*u)...})
					(bond::seek_to_t<-N_len>{});
		}

	public:// OPERATE
		template <int N_sup=0, auto ...Ns> requires (0 == N_sup)
		XTAL_VAL_(return,inline,let)
		method(auto &&x
		,	auto &&...oo
		)	const
		noexcept -> auto
		{
			return S_::template method<Ns...>(XTAL_REF_(x), XTAL_REF_(oo)...);
		}
		template <int N_sup=0, auto ...Ns> requires (1 <= N_sup)
		XTAL_VAL_(return,inline,let)
		method(auto &&x, atom::math::dot_q auto &&a
		,	unstruct_t<decltype(x)> u
		,	auto &&...oo
		)	const
		noexcept -> auto
		{
			int constexpr K_div = two << N_sup;
			int constexpr K_len = one << N_sup;
			int constexpr K_end = K_len -  one;
			using U  = XTAL_ALL_(u);
			using X  = XTAL_ALL_(x);
			using Y  = XTAL_ALL_(XTAL_ANY_(S_).template method<Ns...>(x, a, u, oo...));
			using X_ = atom::quantity_t<X[K_len], xtd::plus_multiplies<>>;
			using Y_ = atom::quantity_t<Y[K_len], xtd::plus_multiplies<>>;

		//	TODO: Interpolate controls?
		//	TODO: Detect mutability/immutability for parallelization?
		//	TODO: Replace by signalling difference upstream,
		//		maybe with `occur::quartz_t<signed>(N_sup)`?
			u *= root_f<-1>((U) K_len);

			auto constexpr z_up = F_men<K_len, 1>(root_f<-1>((U) K_div));
			auto constexpr z_dn = F_lop(z_up);

		//	Upsampling
			auto const     x_up = x*U(K_div - one);
			auto const     x_dn = S_::stash1(x_up);
			auto           x_   = [&]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (X_{(get<I>(z_up)*(x_up) + get<I>(z_dn)*(x_dn))...})
					(bond::seek_to_t<K_len>{});

		//	Mapping/Downsampling:
			auto & y_dn  = S_::template stash1<Y_>();
			auto   y_up  = y_dn;
			get<0>(y_up) = S_::template method<Ns...>(get<0>(x_), a, u, oo...);
			bond::seek_to_e<K_len - 1, 1>([&, this] (auto I)
			XTAL_0FN {
				get<I>(y_up) +=
				get<I>(y_dn)  = S_::template method<Ns...>(get<I>(x_), a, u, oo...);
			});
			return dot_f(z_up, y_up);
		}

	};
};
template <constant_q auto M_sup>
struct oversample<M_sup>
{
	using superkind = oversample<(int) M_sup>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
	//	static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
		using D_ = typename S_::data_type;

	public:// CONSTRUCT
		using S_::S_;

		template <auto ...Ns>
		XTAL_VAL_(return,inline,let)
		method(auto &&...oo
		)	const
		noexcept -> auto
		{
			return S_::template method<M_sup, Ns...>(XTAL_REF_(oo)...);
		}

	};
};

////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
