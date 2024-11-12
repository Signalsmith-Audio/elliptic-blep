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

#------- Create the filter

(filter_z, filter_p, filter_k) = scipy.signal.ellip(filter_order, filter_ripple_db, filter_stop_db, lowpass_freq, analog=True, output='zpk')

(hp_z, hp_p, hp_k) = scipy.signal.butter(blep_order, highpass_freq, btype="highpass", analog=True, output='zpk');

filter_z = numpy.append(filter_z, hp_z)
filter_p = numpy.append(filter_p, hp_p)
filter_k *= hp_k

(filter_b, filter_a) = scipy.signal.zpk2tf(filter_z, filter_p, filter_k)
residue_coeffs, residue_poles, residue_direct = scipy.signal.residue(filter_b, filter_a)
if len(residue_direct) > 0:
	exit(1)

real_indices = numpy.where(numpy.abs(residue_poles.imag) < 1e-6)
real_poles = residue_poles[real_indices]
real_coeffs = residue_coeffs[real_indices]
complex_indices = numpy.where(residue_poles.imag >= 1e-6)
complex_poles = residue_poles[complex_indices]
complex_coeffs = residue_coeffs[complex_indices]*2 # doubled since we're only using one and taking the real part

# Write the results out as a C++ class
print("""template<typename Sample>
struct EllipticBlepCoeffs {
	static constexpr size_t maxIntegrals = %i;
	static constexpr size_t complexCount = %i;
	std::array<std::complex<Sample>, complexCount> complexPoles{{
		%s
	}};
	std::array<std::complex<Sample>, complexCount> complexCoeffs{{
		%s
	}};
	static constexpr int realCount = %i;
	std::array<Sample, realCount> realPoles{{
		%s
	}};
	std::array<Sample, realCount> realCoeffs{{
		%s
	}};
};"""%(
	blep_order,
	len(complex_poles),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_poles]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in complex_coeffs]),
	len(real_poles),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_poles]),
	",\n\t\t".join(["Sample(%s)"%(p.real.astype(str)) for p in real_coeffs])
))
