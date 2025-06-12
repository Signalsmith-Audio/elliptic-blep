"""
Designs a continous-time elliptic filter, shifts the results to make a Hilbert, and writing the results as a C++ struct
"""

import numpy
import scipy
import scipy.signal
import scipy.optimize
import sys

filter_order = 11
filter_ripple_db = 0.1
filter_stop_db = 97.5
lowpass_freq = 20000
highpass_freq = 20
ref_srate = 44100
blep_order = 3

print("\nDesigning BLEP filters", file=sys.stderr)

#------- Create the highpass filter

(hp_z, hp_p, hp_k) = scipy.signal.butter(blep_order, highpass_freq, btype="highpass", analog=True, output='zpk');
(hp_b, hp_a) = scipy.signal.zpk2tf(hp_z, hp_p, hp_k)
hp_coeffs, hp_poles, hp_direct = scipy.signal.residue(hp_b, hp_a)

#------- Create the impulse filter (elliptic lowpass)

(impulse_z, impulse_p, impulse_k) = scipy.signal.ellip(filter_order, filter_ripple_db, filter_stop_db, lowpass_freq, analog=True, output='zpk')

(impulse_b, impulse_a) = scipy.signal.zpk2tf(impulse_z, impulse_p, impulse_k)
impulse_coeffs, impulse_poles, impulse_direct = scipy.signal.residue(impulse_b, impulse_a)
if len(impulse_direct) > 0:
	exit(1)

#------- Create the direct-output filter

filter_z = numpy.append(impulse_z, hp_z)
filter_p = numpy.append(impulse_p, hp_p)
filter_k = impulse_k*hp_k

(filter_b, filter_a) = scipy.signal.zpk2tf(filter_z, filter_p, filter_k)
residue_coeffs, residue_poles, residue_direct = scipy.signal.residue(filter_b, filter_a)
if len(residue_direct) > 0:
	exit(1)

real_indices = numpy.where(numpy.abs(residue_poles.imag) < 1e-6)
real_poles = residue_poles[real_indices]
complex_indices = numpy.where(residue_poles.imag >= 1e-6)
complex_poles = residue_poles[complex_indices]

# Calculate the coefficients for direct synthesis of the output
real_coeffs_direct = residue_coeffs[real_indices]*1
complex_coeffs_direct = residue_coeffs[complex_indices]*2 # doubled since we're only using one and taking the real part

#------- Subtract direct-output filter from the highpass to get the BLEP residue filter

for i in range(len(hp_poles)):
	# Find the pole from the direct filter which matches this highpass one
	pole_distances = numpy.abs(residue_poles - hp_poles[i])
	pole_index = numpy.argmin(pole_distances)
	residue_coeffs[pole_index] -= hp_coeffs[i]

# Calculate the coefficients for the aliasing-cancellation signal
real_coeffs_blep = residue_coeffs[real_indices]*1
complex_coeffs_blep = residue_coeffs[complex_indices]*2 # doubled since we're only using one and taking the real part

#------- Add 0s into the impulse coeffs since it doesn't include the highpass

impulse_coeffs_padded = residue_coeffs*0;
for i in range(len(impulse_poles)):
	# Find the pole from the direct filter which matches this highpass one
	pole_distances = numpy.abs(residue_poles - impulse_poles[i])
	pole_index = numpy.argmin(pole_distances)
	impulse_coeffs_padded[pole_index] = impulse_coeffs[i]

# Calculate the coefficients for the aliasing-cancellation signal
real_coeffs_impulse = impulse_coeffs_padded[real_indices]*1
complex_coeffs_impulse = impulse_coeffs_padded[complex_indices]*2 # doubled since we're only using one and taking the real part

print("\tdone")

#------- Horrible brute-force search for a phase-compensation allpass

max_allpass_order = filter_order
max_delay_samples = max_allpass_order*2
allpass_design_density = 1024
allpass_design_peak_sensitivity = 1; # 1 = average error amplitude, 2 = average error energy, Inf = peak

def find_allpass_coeffs(opt_freq_density, max_allpass_order, max_delay_samples, error_metric=1):
	opt_freqs = numpy.linspace(0, numpy.pi, opt_freq_density)
	
	# Convert to digital filter using impulse-invariance (same as the C++ code will)
	hz_to_angular = 2*numpy.pi/ref_srate
	discrete_poles = numpy.exp(impulse_poles*hz_to_angular);
	discrete_coeffs = impulse_coeffs*hz_to_angular
	discrete_b, discrete_a = scipy.signal.invresz(discrete_coeffs, discrete_poles, [])
	discrete_z, discrete_p, discrete_k = scipy.signal.tf2zpk(discrete_b, discrete_a)

	def delay_allpass_coeffs(poles):
		(_, allpass_a) = scipy.signal.zpk2tf([], poles, 1)
		return allpass_a;
		
	def opt_allpass_response(poles):
		# Add the allpass poles/zeroes to the impulse filter
		combined_p = numpy.append(discrete_p, poles)
		combined_z = numpy.append(discrete_z, 1/poles);
		combined_k = discrete_k*abs(numpy.prod(poles))
		_, response = scipy.signal.freqz_zpk(combined_z, combined_p, combined_k, worN=opt_freqs);
		return response;

	# Decode real-valued optimisation coefficients into complex/real pole positions
	def opt_ri_to_poles(poles_ri):
		poles = []
		# conjugate pairs
		for i in range(0, len(poles_ri) - 1, 2):
			poles.append(poles_ri[i] + 1j*poles_ri[i + 1])
			poles.append(poles_ri[i] - 1j*poles_ri[i + 1])
		if len(poles_ri)%2 == 1:
			poles.append(poles_ri[-1])
		return numpy.array(poles)

	def opt_allpass_score(poles_ri, opt_ref_response):
		poles = opt_ri_to_poles(poles_ri)
		for pole in poles:
			if abs(pole) > 0.99: # inaccurate or unstable
				return 1e100

		response = opt_allpass_response(poles)
		diff = opt_ref_response - response
		
		score = sum(abs(diff)**error_metric) # 1 = error amplitude, 2 = energy, etc.
		return score

	best = None
	best_samples = 0
	_, impulse_response = scipy.signal.freqz_zpk(discrete_z, discrete_p, discrete_k, worN=opt_freqs)
	amp_response = numpy.abs(impulse_response);
	for allpass_order in range(1, max_allpass_order + 1):
		print("\tAllpass order: %i"%allpass_order, file=sys.stderr)
		for samples in range(1, max_delay_samples + 1):
			linear_response = numpy.exp(-1j*samples*opt_freqs)
			linear_response *= amp_response
			opt_start = numpy.zeros(allpass_order) + 0.1
			opt_result = scipy.optimize.minimize(opt_allpass_score, opt_start, args=(linear_response,), method='Nelder-Mead')
			if best == None or opt_result.fun < best.fun:
				best = opt_result
				best_samples = samples - 1
				poles = opt_ri_to_poles(best.x)
				average_error = (opt_result.fun/len(opt_freqs))**(1/error_metric);
				print("\t\t%i-sample delay -> %f"%(samples, average_error), file=sys.stderr)
	opt_allpass_poles = opt_ri_to_poles(best.x)
	return best_samples, delay_allpass_coeffs(opt_allpass_poles)

print("\nFinding phase-compensation allpass", file=sys.stderr)
linear_delay, linear_coeffs = find_allpass_coeffs(allpass_design_density, max_allpass_order, max_delay_samples, allpass_design_peak_sensitivity)

#------- Generate C++
print("\nWriting C++ code: out/coeffs.h")

cppCode = """template<typename Sample>
struct EllipticBlepCoeffs {
	static constexpr size_t maxIntegrals = %i;
	static constexpr size_t complexCount = %i;
	static constexpr int realCount = %i;
	std::array<std::complex<Sample>, complexCount> complexPoles{{
		%s
	}};
	std::array<Sample, realCount> realPoles{{
		%s
	}};
	// Coeffs for direct bandlimited synthesis of a polynomial-segment waveform
	std::array<std::complex<Sample>, complexCount> complexCoeffsDirect{{
		%s
	}};
	std::array<Sample, realCount> realCoeffsDirect{{
		%s
	}};
	// Coeffs for cancelling the aliasing from discontinuities in an existing waveform
	std::array<std::complex<Sample>, complexCount> complexCoeffsBlep{{
		%s
	}};
	std::array<Sample, realCount> realCoeffsBlep{{
		%s
	}};
	// Coeffs for adding an impulse in "residue" (not direct) mode
	std::array<std::complex<Sample>, complexCount> complexCoeffsImpulse{{
		%s
	}};
	std::array<Sample, realCount> realCoeffsImpulse{{
		%s
	}};
	// Allpass to make the phase approximately linear
	static constexpr int allpassLinearDelay = %i;
	static constexpr int allpassOrder = %i;
	std::array<Sample, allpassOrder> allpassCoeffs{{
		%s
	}};
};"""%(
	blep_order,
	len(complex_poles),
	len(real_poles),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_poles]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_poles]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs_direct]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs_direct]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs_blep]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs_blep]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs_impulse]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs_impulse]),
	linear_delay,
	len(linear_coeffs) - 1,
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in linear_coeffs[1:]]),
)
with open("out/coeffs.h", 'w') as file:
	file.write(cppCode)
