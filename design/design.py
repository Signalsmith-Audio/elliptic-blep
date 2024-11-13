"""
Designs a continous-time elliptic filter, shifts the results to make a Hilbert, and writing the results as a C++ struct
"""

import numpy
import scipy
import scipy.signal
import scipy.optimize

filter_order = 11
filter_ripple_db = 0.1
filter_stop_db = 97.5
lowpass_freq = 20000
highpass_freq = 20
blep_order = 3

#------- Create the highpass filter

(hp_z, hp_p, hp_k) = scipy.signal.butter(blep_order, highpass_freq, btype="highpass", analog=True, output='zpk');
(hp_b, hp_a) = scipy.signal.zpk2tf(hp_z, hp_p, hp_k)
hp_coeffs, hp_poles, hp_direct = scipy.signal.residue(hp_b, hp_a)

#------- Create the direct-output filter

(filter_z, filter_p, filter_k) = scipy.signal.ellip(filter_order, filter_ripple_db, filter_stop_db, lowpass_freq, analog=True, output='zpk')

filter_z = numpy.append(filter_z, hp_z)
filter_p = numpy.append(filter_p, hp_p)
filter_k *= hp_k

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

for i in range(len(hp_poles)):
	# Find the pole from the direct filter which matches this highpass one
	pole_distances = numpy.abs(residue_poles - hp_poles[i])
	pole_index = numpy.argmin(pole_distances)
	residue_coeffs[pole_index] -= hp_coeffs[i]

# Calculate the coefficients for the aliasing-cancellation signal
real_coeffs_blep = residue_coeffs[real_indices]*1
complex_coeffs_blep = residue_coeffs[complex_indices]*2 # doubled since we're only using one and taking the real part

# Write the results out as a C++ class
print("""template<typename Sample>
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
};"""%(
	blep_order,
	len(complex_poles),
	len(real_poles),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_poles]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_poles]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs_direct]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs_direct]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs_blep]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs_blep])
))
