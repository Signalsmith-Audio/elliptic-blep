#ifndef SIGNALSMITH_DSP_FFT2_H
#define SIGNALSMITH_DSP_FFT2_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <complex>
#include <vector>

namespace signalsmith { namespace fft2 {

/// A self-contained, reasonably fast power-of-2 FFT template
template<typename Sample>
struct SimpleFFT {
	using Complex = std::complex<Sample>;
	
	SimpleFFT(int maxSize=0) {
		resize(maxSize);
	}
	
	void resize(int maxSize) {
		twiddles.resize(maxSize/2);
		for (int i = 0; i < maxSize/2; ++i) {
			double twiddlePhase = -2*M_PI*i/maxSize;
			twiddles[i] = {
				Sample(std::cos(twiddlePhase)),
				Sample(std::sin(twiddlePhase))
			};
		}
		working.resize(maxSize);
	}
	
	void fft(int size, const Complex *time, Complex *freq) {
		if (size <= 1) {
			*freq = *time;
			return;
		}
		fftPass<false>(size, 1, time, freq, working.data());
	}

	void ifft(int size, const Complex *freq, Complex *time) {
		if (size <= 1) {
			*time = *freq;
			return;
		}
		fftPass<true>(size, 1, freq, time, working.data());
	}
private:
	std::vector<Complex> twiddles;
	std::vector<Complex> working;

	// Calculate a [size]-point FFT, where each element is a block of [stride] values
	template<bool inverse>
	void fftPass(int size, int stride, const Complex *input, Complex *output, Complex *working) const {
		if (size > 2) {
			// Calculate the two half-size FFTs (odd and even) by doubling the stride
			fftPass<inverse>(size/2, stride*2, input, working, output);
			combine2<inverse>(size, stride, working, output);
		} else {
			// The input can already be considered a 1-point FFT
			combine2<inverse>(size, stride, input, output);
		}
	}
	
	// Combine interleaved even/odd results into a single spectrum
	template<bool inverse>
	void combine2(int size, int stride, const Complex *input, Complex *output) const {
		auto twiddleStep = twiddles.size()*2/size;
		for (int i = 0; i < size/2; ++i) {
			Complex twiddle = twiddles[i*twiddleStep];
			
			const Complex *inputA = input + 2*i*stride;
			const Complex *inputB = input + (2*i + 1)*stride;
			Complex *outputA = output + i*stride;
			Complex *outputB = output + (i + size/2)*stride;
			for (int s = 0; s < stride; ++s) {
				Complex a = inputA[s];
				Complex b = inputB[s]*(inverse ? std::conj(twiddle) : twiddle);
				outputA[s] = a + b;
				outputB[s] = a - b;
			}
		}
	}
};

}} // namespace
#endif // include guard
