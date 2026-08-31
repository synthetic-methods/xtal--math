#pragma once
#include "./any.hh"
#include "./serial.hh"

#include "../../process/imagine.hh"
#include "../../process/pade/unity.hh"


XTAL_ENV_(push)
namespace xtal::atom::math::fourier
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	XTAL_TYP_(new) series;
template <class   ..._s>	XTAL_TYP_(let) series_t = typename series<_s...>::type;
template <class   ...Ts>	XTAL_TYP_(ask) series_q = bond::tagged_as_p<series_t, Ts...>;

XTAL_VAL_(let) series_f = [] XTAL_1FN_(call) (_detail::factory<series_t>::make);


////////////////////////////////////////////////////////////////////////////////

template <class T> XTAL_TYP_(let) serious_t = typename std::remove_cvref_t<T>::transverse::type;
template <class T> XTAL_TYP_(ask) serious_q = requires {typename serious_t<T>;};

XTAL_VAL_(return,inline,let)
serious_f(serious_q auto &&t)
noexcept -> decltype(auto)
{
	using T = decltype(t);
	using F = xtd::qualify_cvref_t<T, serious_t<T>>;
	return reinterpret_cast<F>(XTAL_REF_(t));
};
XTAL_VAL_(return,inline,let)
serious_f(auto &&t)
noexcept -> decltype(auto)
{
	return XTAL_REF_(t);
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `serial` with multiplication via circular convolution.
*/
template <class ...Us> requires un_v<bond::devise_condensed_p<        Us...>>
struct series<Us...> :               bond::devise_condensed_s<series, Us...>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and un_q<xtd::plus<>, Us...>
struct series<Us...> : series<Us..., xtd::plus<>>
{};
template <class ...Us> requires in_v<bond::devise_condensed_p<Us...>> and in_q<xtd::plus<>, Us...>
struct series<Us...>
{
	template <class T>
	using endotype = typename serial<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<series_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_    = holotype<T>;
		using A_    = typename S_:: archetype;
		using I_    = typename S_::index_type;
		using U_    = typename S_::value_type;
		using V_    = typename S_::scale_type;
		using V_fit = bond::fit<V_>;
		static_assert(fixed_q<A_>);

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

	public:// CONSTRUCT
		using S_::S_;

		/*!
		\brief   Generates a section of the complex sinusoid determined by `std::pow(2, n)`.
		*/
		XTAL_VAL_(inline)
		XTAL_VAL_(new,explicit)
		homotype(constant_q auto const n, auto &&...oo)
		noexcept
		{
			generate<n>(XTAL_REF_(oo)...);
		}
		/*!
		\brief   Generates a section of the complex sinusoid determined by `std::pow(2, n)`.
		*/
		template <auto ...Ns>
		XTAL_VAL_(let)
		generate(constant_q auto const n, auto &&...oo)
		noexcept -> T &
		{
			generate<n, Ns...>(XTAL_REF_(oo)...);
		}

		/*!
		\brief   Populates `this` with powers of `u` from `N_dex` to `N_dex + N_count - 1`.
		\returns `this`.
		*/
		template <int N_dex=0, int N_ind=0, int N_step=1, int N_size=size>
		XTAL_VAL_(inline,let)
		generate(U_ const &u)
		noexcept -> T &
		{
			using process::math::cut_f;
			using process::math::monomial_f;
			using A_delta = typename S_::difference_type;

		//	Compute the start- and end-points for the required segment:
			A_delta constexpr N_lim = N_dex + N_size;
			A_delta constexpr _0 =  0*N_step;
			A_delta constexpr _1 =  1*N_step;
			A_delta constexpr _2 =  2*N_step;
			A_delta constexpr I0 = _1*N_dex + N_ind, J0 = _2*N_dex + N_ind;
			A_delta constexpr IZ = _1*N_lim + N_ind, JZ = _2*N_lim + N_ind;

			auto &s = self();
			XTAL_IF0
			XTAL_0IF (N_dex < 0) {
				generate<-N_dex, I0, N_step, N_lim>(u);
				if constexpr (N_lim&1) {
					get<N_lim>(s) = monomial_f<2>(get<N_lim/2>(s));
				}
			}
			XTAL_0IF (N_size < I0 + _1) {
			//	Populate the 0th-power only:
				get<I0>(s) = monomial_f<N_dex>(u);
			}
			XTAL_0IF_(else) {
			//	Populate the 0th and 1st powers:
				auto const o = monomial_f<N_dex>(u);
				get<I0 + _0>(s) = o;
				get<I0 + _1>(s) = o*u;

			//	Populate the remaining powers by squaring/multiplication:
				bond::seek_to_e<(N_size >> 1U)>([&] (auto M)
					XTAL_0FN -> void {
						auto constexpr UM = I0 + _1*M;
						auto constexpr WM = J0 + _2*M;

						auto const w = monomial_f<2>(cut_f<-3>(get<UM>(s)));
						get<WM + _0>(s) =   w;
						get<WM + _1>(s) = u*w;
					}
				);
			//	Compute the final value if `N_size` is odd:
				if constexpr (in_v<N_size&1> and in_v<N_size^1>) {
					get<IZ - _1>(s) = get<IZ - _2>(s)*(u);
				}
			}
			return s;
		}

		/*!
		\brief   Generates the complex sinusoid with length `2*PI*std::pow(2, N_shift)`.
		\note    To generate the FFT basis used by `transform` etc, use `N_shift == -1`.
		\returns `this`.
		\tparam  `N_shift >= -3`
		*/
		template <int N_shift=0> requires complex_field_q<U_>
		XTAL_VAL_(let)
		generate()
		noexcept -> T &
		{
			using std::span;
			using std::prev;
			using std::next;

		//	Initialize the forwards and backwards iterators:
			auto const i = S_::begin();
			auto const j = S_::rend() - 1;

		//	Compute the fractional sinusoid for this `size`:
			auto constexpr y = process::math::pade::unity_f<(+1)>(V_fit::ratio_f(-1, size << 1));

		//	Compute the initial `1/8`th then mirror the remaining segments:
			typename S_::difference_type constexpr M = size >> 2U;// `1/8`th
			static_assert(-4 <  N_shift);
			generate<0, 0, 1, M + (-3 <  N_shift)>(y);
			if constexpr (-2 <= N_shift) _detail::copy_to<[] (U_ const &v) XTAL_0FN_(to) (U_{-v.imag(), -v.real()})>(prev(j, 2*M), span(i, next(i, 1*M)));
			if constexpr (-1 <= N_shift) _detail::copy_to<[] (U_ const &v) XTAL_0FN_(to) (U_{ v.imag(), -v.real()})>(next(i, 2*M), span(i, next(i, 2*M)));
			if constexpr (-0 <= N_shift) _detail::copy_to<[] (U_ const &v) XTAL_0FN_(to) (U_{-v.real(), -v.imag()})>(next(i, 4*M), span(i, next(i, 4*M)));
			static_assert( 0 >= N_shift);// TODO: Extend to allow multiple copies using `bond::seek`.

			return self();
		}
		/*!
		\brief   Transforms `source` using the FFT, with `this` as the Fourier basis.
		\returns `source`.

		\note    The size of both `this` and `source` must be expressible as an integral power-of-two,
		         and `1 < source.size() <= this->size()`.

		\note    Around 10% to 20% faster than Eigen's default (based on KISSFFT).
		*/
		template <int N_dir=1> requires in_v<N_dir, 1, -1> and complex_field_q<U_>
		XTAL_VAL_(let)
		transform(isomorphic_q<T> auto &source) const
		noexcept -> decltype(auto)
		{
			using process::math::imagine_f;
			using    bond::math::bit_shift_f;
			using    bond::math::bit_count_f;
			using    bond::math::bit_reverse_f;
			unsigned constexpr N_conj = N_dir < 0;
			unsigned const     n_size = source.size();
			unsigned const     n_left = bit_shift_f(n_size);
			assert(2 <= n_size and 1 == bit_count_f(n_size));

		//	Move all inner entries to their bit-reversed locations:
			for (unsigned h{n_size >> 1}; 0 < --h;) {
				std::swap(source[h], source[bit_reverse_f(h, n_left)]);
			}
		//	Compute the transform of `source` using the precomputed half-period sinusoid in `this`:
			for (unsigned u_step{size}, n_step{1}; n_step < n_size;) {auto u = S_::data();
				for (unsigned m{   }; m < n_step;   u += u_step) {auto o = imagine_f<0, N_conj>(*u);
				for (unsigned n{m++}; n < n_size;   n += n_step) {
					auto &x = source[n], &y = source[n += n_step]; y *= o;
					auto  x_y = x - y; x += y; y = XTAL_MOV_(x_y);
				}}
				u_step >>= 1;
				n_step <<= 1;
			}
		//	Conjugate and scale the output if computing the inverse transform of the codomain:
			if constexpr (N_conj) {
				source /= V_fit::alpha_f(n_size);
			}
		//	Cast the output to the transformed domain:
			return serious_f(XTAL_REF_(source));
		}
		template <int N_dir=1> requires in_v<N_dir, 1, -1> and complex_field_q<U_>
		XTAL_VAL_(let)
		transform(isomorphic_q<T> auto &&source) const
		noexcept -> decltype(auto)
		{
			(void) transform<N_dir>(source);
			return serious_f(XTAL_MOV_(source));
		}
		/*!
		\returns a new `series` representing the FFT of `source`,
		using `this` as the Fourier basis.
		*/
		template <int N_dir=1> requires in_v<N_dir, 1, -1>
		XTAL_VAL_(return,inline,let)
		transformation(isomorphic_q<T> auto source) const
		noexcept -> auto
		{
			return transform<N_dir>(XTAL_MOV_(source));
		}

		/*!
		\returns `lhs` convolved with `rhs`, using `this` as the Fourier basis.
		*/
		XTAL_VAL_(let)
		convolve(isomorphic_q<T> auto &&y0, auto y1) const
		noexcept -> decltype(auto)
		{
			static_assert(same_q<decltype(y0), decltype(y1)>);
			return transform<-1>(transform<1>(XTAL_REF_(y0)) *= transform<1>(y1));
		}
		/*!
		\returns A new `series` representing the convolution of `y0` with `y1`,
		using `this` as the Fourier basis.
		*/
		XTAL_VAL_(return,inline,let)
		convolution(isomorphic_q<T> auto y0, auto const &y1) const
		noexcept -> auto
		{
			static_assert(same_q<decltype(y0), decltype(y1)>);
			return convolve(XTAL_MOV_(y0), y1);
		}

		/*!
		\brief   Multiplication by circular convolution.
		*/
		using S_::operator*=;

		XTAL_VAL_(return,inline,let)  operator * (auto const &              w) const noexcept -> auto   {return twin() *=   w ;}
		XTAL_VAL_(inline,let)         operator *=(std::initializer_list<U_> w)       noexcept -> auto & {return self() *= T(w);}

		XTAL_VAL_(let)
		operator *=(T const &t)
		noexcept -> T &
		{
			auto &s = self();
			if constexpr (complex_field_q<U_>) {
				T(constant_t<-1>{}).convolve(s, t);
			}
			else {
				using X = typename V_fit::aphex_type;
				using Y = typename series<X[size]>::type;
				Y s_(s);
				Y t_(t);
				Y(constant_t<-1>{}).convolve(s_, t_);
				_detail::move_to<[] XTAL_1FN_(call) (std::real)>(s.begin(), s_);
			}
			return s;
		}

		/*!
		\brief   The dual of `T`, defined by `hadamard::series`.
		*/
		struct transverse
		{
			template <class Y>
			using holotype = typename quantity<U_[size], xtd::multiplies<>>::template homotype<Y>;

			template <class Y>
			class homotype : public holotype<homotype<Y>>
			{
				using R_ = holotype<homotype<Y>>;

			public:
				using R_::R_;

				struct transverse {using type = T;};

			};
			using type = bond::derive_t<homotype>;

		};
	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
