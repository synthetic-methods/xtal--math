#pragma once
#include "./any.hh"

#include "../occur/dent.hh"
#include "../occur/dash.hh"
#include "./dot.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Rewires the interface to the superprocess
         via the routing/coefficient matrix `U_mat`.
*/
template <class ..._s>	struct  patch;
template <class ..._s>	using   patch_t = confined_t<patch<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <fixed_q U_mat>
struct patch<U_mat>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

		XTAL_TYP_(set) matrix_type = U_mat;
		XTAL_TYP_(set) vector_type = typename fixed<U_mat>::value_type;
		XTAL_DEF_(set) matrix_size = bond::pack_size<matrix_type>{};
		XTAL_DEF_(set) vector_size = bond::pack_size<vector_type>{};

		template <extent_type N_mask=1>
		struct attach
		{
			using superkind = typename occur::math::dent_s<matrix_type>::template attach<N_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				static_assert(process::any_q<R>);
				using R_ = bond::compose_s<R, superkind>;
				
			public:// CONSTRUCT
				using R_::R_;
			
			public:// ACCESS
				using R_::self;
				using R_::head;

				XTAL_FN1_(go) (XTAL_DEF_(return,inline,get) coefficients, [] (auto &&o, auto &&...oo)
				XTAL_0FN_(to) (XTAL_REF_(o).head(XTAL_REF_(oo)...)))

			};
			template <class R> requires complete_q<typename R::template head_t<U_mat>>
			class subtype<R> : public bond::compose_s<R>
			{
				static_assert(process::any_q<R>);
				using R_ = bond::compose_s<R>;
				
			public:// CONSTRUCT
				using R_::R_;

			};
		};
		template <extent_type N_mask=1>
		struct refix
		{
			using superkind = attach<N_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;
				
			public:// CONSTRUCT
				using R_::R_;
			
			public:// OPERATE
			
				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&...oo)
				noexcept -> auto
				{
					/**/
					return [this, o_=_std::tuple{XTAL_REF_(oo)...}]<auto ...I>(bond::seek_t<I...>)
					XTAL_0FN_(to) (R_::template method<Ns...>(dot_f(get<I>(R_::coefficients()), o_)...))
						(bond::seek_s<matrix_size> {});
					/*/
					auto const o_ = _std::tuple{XTAL_REF_(oo)...};
					return\
						([&, this]<auto ...I>(bond::seek_t<I...>)          XTAL_0FN_(to) (R_::template method<Ns...>((
						([&, this]<auto ...J>(bond::seek_t<J...>, auto _I) XTAL_0FN_(to)
							(zero +...+ (bond::pack_item_f<XTAL_ALL_(_I){}, J>(R_::coefficients())*get<J>(o_)))
							(bond::seek_s<vector_size> {}, constant_t<I>{})))...))
							(bond::seek_s<matrix_size> {}));
					/***/
				}

			};
		};
		template <extent_type N_mask=1>
		struct refit
		{
			using superkind = attach<N_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:// CONSTRUCT
				using R_::R_;
			
			public:// FLOW
			
				template <int N_ion>
				XTAL_DEF_(return,inline,let)
				flux(auto &&...oo)
				noexcept -> auto
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_DEF_(return,inline,let)
				flux(auto &&o_)
				noexcept -> auto
				requires same_q<vector_type, decltype(o_)> or occur::math::dash_p<U_mat, decltype(o_)>
				and      same_v<vector_size, bond::pack_size_v<decltype(o_)>>
				{
					return [&, this]<auto ...I>(bond::seek_t<I...>)
						XTAL_0FN_(to) (R_::template flux<N_ion>(
							occur::math::dash_f<U_mat>(dot_f(get<I>(R_::coefficients()), o_)...)))
							(bond::seek_s<matrix_size> {});
				}

			};
		};

		template <class ..._s>
		struct fix
		{
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;
			
			public:// FLOW
			
				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&...oo)
				noexcept -> auto
				{
					return R_::template method<Ns...>(bond::operate<_s>{}(XTAL_REF_(oo))...);
				}

			};
		};
		template <class ..._s>
		struct fit
		{
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;
			
			public:// FLOW
			
				template <int N_ion>
				XTAL_DEF_(return,inline,let)
				flux(auto &&...oo)
				noexcept -> auto
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_DEF_(return,inline,let)
				flux(occur::math::dash_q<U_mat> auto &&o_)
				noexcept -> auto
				requires same_v<vector_size, bond::pack_size_v<decltype(o_)>>
				{
					return XTAL_REF_(o_).apply([this] (auto &&...oo)
						XTAL_0FN_(to) (R_::template flux<N_ion>(bond::operate<_s>{}(XTAL_REF_(oo))...)));
				}

			};
		};

	};
};
template <fixed_q U_mat, class ..._s>
struct patch<U_mat, _s...>
{
	using superkind = patch<U_mat>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

		template <extent_type N_mask=1>
		struct refix
		:	bond::compose<void
			,	typename S_::template refix<N_mask>
			,	typename S_::template   fix<_s... >
			>
		{
		};
		template <extent_type N_mask=1>
		struct refit
		:	bond::compose<void
			,	typename S_::template refit<N_mask>
			,	typename S_::template   fit<_s... >
			>
		{
		};
	//	TODO: Remove `fi[tx]`?

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
