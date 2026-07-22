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
\brief   Upsamples the super-filter.
*/
template <auto ...As>	struct  upsample;
template <auto ...As>	using   upsample_t = confined_t<upsample<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <variable_q auto M_sup>
struct upsample<M_sup>
{
private:
	XTAL_DEF_(set) M_len = one << M_sup;

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
		XTAL_DEF_(return,inline,set)
		F_lop(auto &&u_)
		noexcept -> auto
		{
			using U_ = XTAL_ALL_(u_);
			int constexpr I_lim = U_::size();
			int constexpr I_end = I_lim - 2;
			return [u_=XTAL_REF_(u_)]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (U_{u_.template element<I_end - I>()...})
					(bond::seek_to_t<U_::size()>{});
		}
		XTAL_DEF_(return,inline,set)
		F_cos(auto &&...oo)
		noexcept -> auto
		{
			auto const w = taylor::cosy_t<1, 2>{}.template method<2>((one *...* XTAL_REF_(oo)));
			return term_f(-one, two, XTAL_MOV_(w));
		}
		XTAL_DEF_(return,inline,set)
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
		XTAL_DEF_(return,inline,set)
		F_men(auto &&u)
		noexcept -> auto
		{
			using   U  = XTAL_ALL_(u);
			using   U_ = atom::quantity_plus_multiplies_t<U[N_len]>;
			return [u=XTAL_REF_(u)]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (U_{(F_man(I + N_off, u)*u)...})
					(bond::seek_to_t<N_len>{});
		}

	public:// OPERATE
		template <int N_sup=0, auto ...Ns> requires (0 == N_sup)
		XTAL_DEF_(return,inline,let)
		method(auto &&x
		,	auto &&...oo
		)	const
		noexcept -> auto
		{
			return S_::template method<Ns...>(XTAL_REF_(x), XTAL_REF_(oo)...);
		}
		template <int N_sup=0, auto ...Ns> requires (1 <= N_sup)
		XTAL_DEF_(return,inline,let)
		method(auto &&x, atom::math::dot_q auto &&a
		,	unstruct_t<decltype(x)> u
		,	auto &&...oo
		)	const
		noexcept -> auto
		{
			using U  = XTAL_ALL_(u);
			using X  = XTAL_ALL_(x);
			using Y  = XTAL_ALL_(XTAL_ANY_(S_).template method<Ns...>(x, a, u, oo...));
			int constexpr K_len =   one << N_sup; U constexpr k_len = K_len;
			int constexpr K_lim = K_len *    two; U constexpr k_lim = K_lim;
			int constexpr K_max = K_lim -    one; U constexpr k_max = K_max;

			u *= root_f<-1>(k_len);

			using X_ = atom::quantity_plus_multiplies_t<X[K_len]>;
			using Y_ = atom::quantity_plus_multiplies_t<Y[K_len]>;

			auto constexpr z_up = F_men<K_len, 1>(root_f<-1>(k_lim));
			auto         & y_dn = get<0>(S_::template stash<Y_>());
			auto         & x_dn = get<K_len - 1>(y_dn);

			auto constexpr z_dn = F_lop(z_up);
			auto           y_up = F_lop(y_dn);
			auto           x_up = x*k_max;
		//	auto           x_co = z_up*x_up + z_dn*x_dn;// TODO: Fix lifted scalar multiplication...
			auto           x_co = [&]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (X_{(get<I>(z_up)*(x_up) + get<I>(z_dn)*(x_dn))...})
					(bond::seek_to_t<K_len>{});

			bond::seek_to_e<-K_len>([&] (auto const I)
			XTAL_0FN_(do) (get<I>(y_up) +=
				get<I>(y_dn) = S_::template method<Ns...>(get<I>(x_co), a, u, oo...)
			));

			auto const y = dot_f(z_up, y_up); x_dn = x_up;
			return y;
		}

	};
};
template <constant_q auto M_sup>
struct upsample<M_sup>
{
	using superkind = upsample<(int) M_sup>;

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
		XTAL_DEF_(return,inline,let)
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
