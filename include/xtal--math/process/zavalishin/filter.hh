#pragma once
#include "./any.hh"

#include "../../provision/shaper.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Integrating filter based on the Topology Preserving Transform (TPT), \
described in _The Art of VA Synthesis_ by Vadim Zavalishin. \

///\
The non-linearity is supplied as a `process`-`template` using `provision::shaper`. \
The process must conform to the signature `<M_ism, M_car>`, \
defining a stateless `static_method<N_var, ...>` within the `subtype`. \

///\
The parameters `M_ism` and `M_car` determine the type and return-value of curve, respectively. \
`M_ism` is expected to yield convex/concave curves for positive/negative values, \
and `M_car` is expected to return the slope when `== -1`. \

///\note\
The filter is preconfigured to support `std::complex<double>[4]` poles. \
With base-types of `double` and `std::complex<double>` respectively, \
the storage required is `16` and `32` bytes-per-pole. \

template <typename ...As> struct   filter;
template <typename ...As> using    filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class U_pole, int N_pole>
struct filter<U_pole[N_pole]>
{
	using   sample_type = occur::sample_t<>;
	using     zoom_type = occur::inferred_t<struct     ZOOM, typename bond::fixture<>::alpha_type>;
	using   select_type = occur::inferred_t<struct   SELECT, unsigned int, bond::word<2>>;
	using    order_type = occur::inferred_t<struct    ORDER, unsigned int, bond::seek_s<N_pole + 1>>;
	using topology_type = occur::inferred_t<struct TOPOLOGY, unsigned int, bond::word<2>>;
	using    limit_type = occur::inferred_t<struct    LIMIT, unsigned int, bond::word<2>>;

	using superkind = bond::compose<bond::tag<filter_t>
	,	provision::stowed<U_pole[N_pole << 1]>
	,	typename   sample_type::template   attach<>
	,	typename     zoom_type::template   attach<>
	,	typename   select_type::template dispatch<>
	,	typename    order_type::template dispatch<>
	,	typename topology_type::template dispatch<>
	,	typename    limit_type::template dispatch<>// TODO: Replace with `shape`'s parameters?
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			if constexpr (N_ion == +1) {
				if constexpr (in_q<decltype(o), order_type, topology_type>) {
					S_::stow(constant_n<0>);
				}
				if constexpr (occur::stage_q<decltype(o)>) {
					if (o == 0) {
						S_::stow(constant_n<0>);
					}
				}
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

		template <int N_sel=0, int N_ord=0, int N_top=0, auto ...Ns>
		XTAL_DEF_(inline,let)
		method(auto &&x_input
		,	absolve_u<decltype(x_input)> s_scale
		,	absolve_u<decltype(x_input)> s_damping
		,	absolve_u<decltype(x_input)> y_balance=0
		)
		noexcept -> auto
		{
			using X = XTAL_ALL_(x_input);
			using W = absolve_u<X>;

			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return XTAL_REF_(x_input);
			}
			XTAL_0IF (0 == N_top) {
				using U_outputs = atom::couple_t<X[N_ord + 1]>;
				using U_scalars = atom::couple_t<W[N_ord + 1]>;
				union U_io {U_outputs outputs; U_scalars scalars;};

				U_io io{[=] () XTAL_0FN -> U_scalars {
					auto const  &u = s_damping;
					auto const u02 =             two*u , u04 = u02*two;
					auto const u12 = term_f(one, two,u), w24 = u02*u12;
					XTAL_IF0
					XTAL_0IF (1 == N_ord) {return {one, one};}
					XTAL_0IF (2 == N_ord) {return {one, u02, one};}
					XTAL_0IF (3 == N_ord) {return {one, u12, u12, one};}
					XTAL_0IF (4 == N_ord) {return {one, u04, w24, u04, one};}
				}()};

				(void) static_edit<N_ord, N_top, Ns...>(XTAL_REF_(x_input), s_scale, io);// io.scalars -> io.outputs

				auto constexpr I_ =  static_cast<unsigned>(N_sel);
				auto constexpr I0 = _std::countr_one(I_ >>  0) +  0, J0 = I0 + 1;
				auto constexpr I1 = _std::countr_one(I_ >> J0) + J0, J1 = I1 + 1;
				W const &y1 = y_balance;
				//\
				W const  y0 = one;
				W const  y0 = term_f<-1, 2>(one, y1);
				if constexpr (I_ == 0) {
					static_assert(I0 == 0);
					static_assert(I1 == 1);
				}
				if constexpr (I1 < N_ord) {
					return term_f(y0*get<I0>(io.outputs), y1, get<I1>(io.outputs));
				}
				else {
					return y0;
				}
			}
			XTAL_0IF (1 == N_top) {
				return XTAL_REF_(x_input);
			}
		}
		template <int N_ord=0, int N_top=0, int N_lim=0, auto ...Ns> requires (1 <= N_ord and N_top == 0)
		XTAL_DEF_(inline,let)
		static_edit(auto const &x_input
		,	absolve_u<decltype(x_input)> s_scale
		,	auto &io
		)
		noexcept -> void
		{
			using X = XTAL_ALL_(x_input);
			using W = absolve_u<X>;

			using U_scalars_ = atom::couple_t<W[N_ord]>;
			using U_outputs_ = atom::couple_t<X[N_ord]>;
			using U_poles_   = atom::couple_t<X[N_ord]>;

			auto  scalars = io.scalars; auto scalars_ = scalars.self(constant_n<N_ord>);
			auto &outputs = io.outputs; auto outputs_ = outputs.self(constant_n<N_ord>);

			auto  stowed  = S_::template stow<U_poles_, U_poles_>();
			auto &states_ = get<0>(stowed);
			auto &slopes_ = get<1>(stowed);

		//	Initialize `scalars*`:
			slopes_.template blanket<1>();
			slopes_ *= scalars_;

		//	Initialize `outputs*`:
			get<0 >(outputs) = get<0>(slopes_);
			bond::seek_forward_f<N_ord - 1>([&] (auto I) XTAL_0FN {
				get<I + 1>(outputs_) = term_f(get<I + 1>(slopes_), get<I>(outputs_), s_scale);
			});
			get<N_ord>(outputs) = root_f<-1, (1)>(term_f(one, get<N_ord - 1>(outputs_), s_scale));

			auto constexpr K_lim = provision::shaper_q<S_> and above_n<0, N_lim>;

		//	Integrate `states*` with `scalars*` and `outputs*`:
			bond::seek_forward_f<1 + 1*K_lim>([&] (auto K) XTAL_0FN {
				XTAL_IF0
				XTAL_0IF (0 == K) {get<N_ord>(outputs) *= x_input - outputs_.product(states_);}
				XTAL_0IF (1 <= K) {get<N_ord>(outputs)  = x_input - outputs_.product(slopes_);}

				bond::seek_backward_f<N_ord>([&] (auto I) XTAL_0FN {
					XTAL_IF0
					XTAL_0IF (0 == K_lim) {
						get<I>(outputs_) = term_f(get<I>(states_), get<I + 1>(outputs), s_scale);
					}
					XTAL_0IF (1 == K_lim) {
						auto const exput = term_f(get<I>(states_), get<I + 1>(outputs), s_scale);
						auto const slope = S_::template shape_t<-1, -1>::template static_method<N_lim, Ns...>(exput) + (get<I>(scalars_) - one);
						get<I>(outputs_) = exput;
						get<I>( slopes_) = slope;
					}
				});
			});

		//	Finalize states and `outputs*`/`scalars*`:
			bond::seek_forward_f<N_ord>([&] (auto I) XTAL_0FN {
				get<I>(states_) = term_f(-get<I>(states_), two, get<I>(outputs_));
			});
			slopes_ /= scalars_;
		//	outputs_ *= slopes_;// TODO: Make this line optional?
		}

	};
};
template <>
struct filter<>
:	filter<typename bond::fixture<>::aphex_type[4]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
