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

template <fixed_shaped_q U_mat>
struct patch<U_mat>
{
public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

		XTAL_TYP_(set) matrix_type = U_mat;
		XTAL_VAL_(set) outer_count = bond::pack_size<bond::pack_item_t<U_mat, 0>>{};
		XTAL_VAL_(set) inner_count = bond::pack_size<                  U_mat    >{};

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

				XTAL_FN1_(go) (XTAL_VAL_(return,inline,get) coefficients, [] (auto &&o, auto &&...oo)
				XTAL_0FN_(to) (XTAL_REF_(o).head(XTAL_REF_(oo)...)))

			protected:
				XTAL_VAL_(return,inline,let)
				dot(auto &&f_, bond::pack_q auto &&x_) const
				noexcept -> decltype(auto)
				{
					auto const &y_ = coefficients();
					return\
						([&]<auto ...O>(bond::seek_in_t<O...>)         XTAL_0FN_(to) (f_((
						([&]<auto ...I>(bond::seek_in_t<I...>, auto o) XTAL_0FN_(to)
							//\
							(zero +...+ (bond::pack_item_f<o, I>(y_)*get<I>(x_)))
							(zero +...+ (bond::pack_item_f<I, o>(y_)*get<I>(x_)))
							(bond::seek_to_t<inner_count> {}, constant_t<O>{})))...))
							(bond::seek_to_t<outer_count> {}));
				}

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
		struct rewire
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
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> auto
				{
					return R_::dot([&, this]
						XTAL_1FN_(call) (R_::template method<Ns...>),
						bond::pack_f(XTAL_REF_(oo)...));
				}
				
			};
		};
		template <extent_type N_mask=1>
		struct reflux
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
				XTAL_VAL_(return,inline,let)
				flux(auto &&...oo)
				noexcept -> signed
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion, class X_>
				XTAL_VAL_(return,inline,let)
				flux(X_ &&x_)
				noexcept -> signed
				requires same_v<inner_count, bond::pack_size_v<X_>>
				and in_v<1
				,	bond::tab_compatible_q<matrix_type, X_>
				,	   occur::math::dash_p<matrix_type, X_>
				>
				{
					return R_::dot([&, this]
						XTAL_1FN_(call) (flux_dash<N_ion>),
						XTAL_REF_(x_));
				}

			protected:
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux_dash(auto &&...oo)
				noexcept -> signed
				{
					using occur::math::dash_f;
					return R_::template flux<N_ion>(dash_f<U_mat>(XTAL_REF_(oo)...));
				}


			};
		};

		template <class ..._s>
		struct prewire
		{
		private:
			XTAL_VAL_(set) _N = sizeof...(_s);
			XTAL_VAL_(set) _M = inner_count;
			static_assert(1 <= _N);
			using U = bond::seek_front_t<_s...>;

		public:
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;
			
			public:// FLOW
			
				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> auto
				{
					XTAL_IF0
					XTAL_0IF (_N == _M) {return R_::template method<Ns...>(bond::operate<_s>{}(XTAL_REF_(oo))...);}
					XTAL_0IF (_N ==  1) {return R_::template method<Ns...>(bond::operate<U >{}(XTAL_REF_(oo))...);}
				}

			};
		};
		template <class ..._s>
		struct preflux
		{
		private:
			XTAL_VAL_(set) _N = sizeof...(_s);
			XTAL_VAL_(set) _M = inner_count;
			static_assert(1 <= _N);
			using U = bond::seek_front_t<_s...>;

		public:
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;
			
			public:// FLOW
			
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(auto &&...oo)
				noexcept -> auto
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(occur::math::dash_q<U_mat> auto &&o_)
				noexcept -> signed
				requires same_v<inner_count, bond::pack_size_v<decltype(o_)>>
				{
					return XTAL_REF_(o_).apply([this] (auto &&...oo)
					XTAL_0FN -> signed {
						XTAL_IF0
						XTAL_0IF (_N == _M) {return R_::template flux<N_ion>(bond::operate<_s>{}(XTAL_REF_(oo))...);}
						XTAL_0IF (_N ==  1) {return R_::template flux<N_ion>(bond::operate<U >{}(XTAL_REF_(oo))...);}
					});
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
		struct rewire
		:	bond::compose<void
			,	typename S_::template  rewire<N_mask>
			,	typename S_::template prewire<_s... >
			>
		{
		};
		template <extent_type N_mask=1>
		struct reflux
		:	bond::compose<void
			,	typename S_::template  reflux<N_mask>
			,	typename S_::template preflux<_s... >
			>
		{
		};
	//	TODO: Remove `fi[tx]`?

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
