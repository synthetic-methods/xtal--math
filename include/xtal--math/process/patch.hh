#pragma once
#include "./any.hh"

#include "../flow/dent.hh"
#include "./dot.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Rewires the interface to the superprocess
         via the routing/coefficient matrix `U_mtx`.
*/
template <class ..._s>	XTAL_TYP_(new) patch;
template <class ..._s>	XTAL_TYP_(let) patch_t = confined_t<patch<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <fixed_shaped_q U_mtx>
struct patch<U_mtx>
{
public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S>;

	protected:
		XTAL_VAL_(set)  outer_size = std::tuple_size<std::tuple_element_t<0, U_mtx>>{};
		XTAL_VAL_(set)  inner_size = std::tuple_size<                        U_mtx >{};

	public:// CONSTRUCT
		using S_::S_;

		template <extent_type N_mask=1>
		struct attach
		{
			template <class T>
			using endokind = bond::compose<void
			,	process::_detail::navigate<T, U_mtx>
			,	typename flow::math::dent_s<U_mtx>::template attach<N_mask>
			>;
			template <class R>
			class subtype : public bond::compose_s<R, endokind<subtype<R>>>
			{
				static_assert(process::any_q<R>);
				using R_ = bond::compose_s<R, endokind<subtype<R>>>;

			public:// CONSTRUCT
				using R_::R_;

			public:// ACCESS
				using R_::self;
				using R_::head;

				XTAL_FN4_(then) (subtype, XTAL_VAL_(return,inline,get) coefficients, [] (auto &&o, auto &&...oo)
				XTAL_0FN_(to) (XTAL_REF_(o).head(XTAL_REF_(oo)...)))

			protected:
				XTAL_VAL_(return,inline,let)
				dot     (bond::pack_q auto &&x_) const
				noexcept -> decltype(auto)
				{
					return dot_then(                XTAL_REF_(x_), bond::pack_f);
				}
				XTAL_VAL_(return,inline,let)
				dot_then(bond::pack_q auto &&x_, auto &&f_) const
				noexcept -> decltype(auto)
				{
					return dot_then(coefficients(), XTAL_REF_(x_), XTAL_REF_(f_));
				}

				XTAL_VAL_(return,inline,set)
				dot     (bond::pack_q auto const &y_, bond::pack_q auto &&x_)
				noexcept -> decltype(auto)
				{
					return dot_then( XTAL_REF_(y_), XTAL_REF_(x_), bond::pack_f);
				}
				XTAL_VAL_(return,inline,set)
				dot_then(bond::pack_q auto const &y_, bond::pack_q auto &&x_, auto &&f_)
				noexcept -> decltype(auto)
				{
					return\
						([&]<auto ...O>(bond::seek_in_t<O...>)         XTAL_0FN_(to) (f_((
						([&]<auto ...I>(bond::seek_in_t<I...>, auto o) XTAL_0FN_(to)
							//\
							(zero +...+ (bond::pack_item_f<o, I>(y_)*get<I>(x_)))
							(zero +...+ (bond::pack_item_f<I, o>(y_)*get<I>(x_)))
							(bond::seek_to_t<inner_size> {}, constant_t<O>{})))...))
							(bond::seek_to_t<outer_size> {}));
				}

			};
			template <class R> requires complete_q<typename R::template head_t<U_mtx>>
			class subtype<R> : public bond::compose_s<R>
			{
				static_assert(process::any_q<R>);
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;

			};
		};
		template <class ...Xs> XTAL_TYP_(let) deflux_t = atom::bundle_t<Xs..., union DEFLUX>;
		template <class ...Xs> XTAL_TYP_(let) reflux_t = bond::compose_s<occur::conferred_t<decltype(one)>, bond::tab<Xs..., union REFLUX>>;
		template <class ..._s> XTAL_TYP_(new) deflow;
		template <class ..._s> XTAL_TYP_(new) reflow;
		template <class ..._s> XTAL_TYP_(new) rewire;

		template <bond::pack_q X_> requires in_v<inner_size, bond::pack_size_v<X_>>
		struct deflow<X_>
		{
			using superkind = bond::tag<deflow>;
			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;
				using H_ = R_::template head_s<U_mtx>;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE

				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(auto &&...oo)
				noexcept -> signed
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(X_ x_)
				noexcept -> signed
				{
					return H_::dot_then(XTAL_MOV_(x_), [&, this]
						(auto &&...oo) XTAL_0FN_(to)
							(R_::template flux<N_ion>(reflux_t<>{}, XTAL_REF_(oo)...)));
				}

			};
		};
		template <class ...Xs> requires in_v<inner_size == sizeof...(Xs)>
		struct deflow<Xs...>
		{
			using superkind = deflow<deflux_t<Xs...>>;

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
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(Xs ...xs)
				noexcept -> signed
				requires (inner_size == sizeof...(xs))
				{
					return R_::template flux<N_ion>(deflux_t<Xs...>{XTAL_MOV_(xs)...});
				}

			};
		};

		template <class Y_> requires in_v<outer_size, std::tuple_size_v<Y_>>
		struct reflow<Y_>
		{
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
				noexcept -> signed
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(reflux_t<>, auto &&...oo)
				noexcept -> signed
				{
					return R_::template flux<N_ion>(bond::operate<Y_>{}(XTAL_REF_(oo)...));
				}

			};
		};
		template <class ...Ys>
		struct reflow
		{
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
				noexcept -> signed
				{
					return R_::template flux<N_ion>(XTAL_REF_(oo)...);
				}
				template <int N_ion>
				XTAL_VAL_(return,inline,let)
				flux(reflux_t<>, auto &&...oo)
				noexcept -> signed
				{
					XTAL_IF0
					XTAL_0IF (in_v<outer_size, sizeof...(Ys)>) {
						return R_::template flux<N_ion>(bond::operate<Ys>{}(XTAL_REF_(oo))...);
					}
					XTAL_0IF (in_v<      zero, sizeof...(Ys)>) {
						return R_::template flux<N_ion>(                    XTAL_REF_(oo) ...);
					}
				}

			};
		};

		template <class ...Xs>
		struct rewire
		{
			static_assert(0 == sizeof...(Xs));// For now...

			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;
				using H_ = R_::template head_s<U_mtx>;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> auto
				{
					return H_::dot_then(bond::pack_f(XTAL_REF_(oo)...), [&, this]
						XTAL_1FN_(call) (R_::template method<Ns...>));
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
