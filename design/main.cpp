#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << std::endl;

#include "../elliptic-blep.h"

#include "./simple-fft.h"
#include "plot/plot.h" // https://signalsmith-audio.co.uk/code/plot.git

#include <vector>

static constexpr double plotScaleX = 0.8, plotScaleY = 1;

template<typename Sample>
void plotBleps(std::string plotName, bool direct, const int oversample=16, const double sampleRate=48000) {
	using Complex = std::complex<Sample>;
	signalsmith::blep::EllipticBlep<Sample> blep(direct, sampleRate);
	
	constexpr double impulseSeconds = 10;
	size_t impulseLength = 1;
	while (impulseLength < sampleRate*impulseSeconds) impulseLength *= 2;
	std::vector<Complex> impulse(impulseLength);

	signalsmith::plot::Figure impulseFigure, spectrumFigure;
	auto &impulsePlot = impulseFigure(0, 0).plot(600*plotScaleX, 120*plotScaleY);
	impulsePlot.x.major(0).linear(0, 0.001).minor(0.001, "1ms").label("time");
	impulsePlot.y.major(0).label(direct ? "direct response" : "BLEP residue");
	if (direct) {
		impulsePlot.y.minor(1, " 1").linear(-0.5, 1.5);
	} else {
		impulsePlot.y.minor(-1, "-1").linear(-1.2, 0.8);
	}
	auto &impulsePlot2 = impulseFigure(0, 1).plot(600*plotScaleX, 120*plotScaleY);
	impulsePlot2.x.major(0).linear(0, 0.2).minor(0.1, "100ms").minor(0.2, "200ms");
	impulsePlot2.y.major(0).minor(100).minor(-100).linear(-150, 150);
	auto &spectrumPlot = spectrumFigure(0, 1).plot(600*plotScaleX, 250*plotScaleY);
	spectrumPlot.x.major(0).minor(20000, "20kHz").minor(24000, "24kHz").linear(0, 26000);
	spectrumPlot.y.major(0).minors(-30, -60, -90).label("dB").linear(-105, 5);
	auto &spectrumPlotZoom = spectrumFigure(0, 2).plot(600*plotScaleX, 80*plotScaleY);
	spectrumPlotZoom.x.copyFrom(spectrumPlot.x);
	spectrumPlotZoom.y.major(0).minors(-0.1).label("ripple").linear(-0.2, 0.1);
	auto &spectrumPlotZoom2 = spectrumFigure(0, 3).plot(600*plotScaleX, 80*plotScaleY);
	spectrumPlotZoom2.x.linear(0, 100).major(0).minors(20, 40, 60, 80).label("Hz");
	spectrumPlotZoom2.y.linear(-6, 0).major(0).minors(-3, -6).label("highpass");

	std::vector<Complex> spectrum(impulseLength);
	signalsmith::fft2::SimpleFFT<Sample> fft(impulseLength);
	
	auto &legend = impulsePlot.legend(1, 1);
	auto &legend2 = spectrumPlot.legend(0.05, 0);
	std::vector<const char *> labels = {"impulse", "1st-order BLEP", "2nd-order BLEP (BLAMP)", "3rd-order BLEP"};

	for (size_t o = 0; o <= blep.maxBlepOrder; ++o) {
		blep.reset();
		Sample impulseAmount = 1;
		blep.add(impulseAmount, o);
		for (size_t i = 0; i < impulseLength; ++i) {
			impulse[i] = blep.get();
			blep.step(1.0/oversample);
		}
		
		fft.fft(impulseLength, impulse.data(), spectrum.data());

		auto &impulseLine = impulsePlot.line();
		auto &impulseLine2 = impulsePlot2.line();
		auto &spectrumLine = spectrumPlot.line();
		auto &spectrumLineZoom = spectrumPlotZoom.line();
		auto &spectrumLineZoom2 = spectrumPlotZoom2.line();
		if (o < labels.size()) {
			legend.add(impulseLine, labels[o]);
			legend2.add(impulseLine, labels[o]);
		} else {
			legend.add(impulseLine, std::to_string(o) + "th-order BLEP");
			legend2.add(impulseLine, std::to_string(o) + "th-order BLEP");
		}
		for (size_t i = 0; i < impulseLength; ++i) {
			impulseLine.add(i/sampleRate/oversample, impulse[i].real());
			impulseLine2.add(i/sampleRate/oversample, impulse[i].real());
			Sample hz = i*sampleRate*oversample/impulseLength;
			Complex v = spectrum[i]*Sample(1.0/oversample);
			// Scale them so they cross over at 1000
			v *= Sample(std::pow(2*M_PI*1000/sampleRate, o));
			Sample db = 10*std::log10(std::norm(v) + 1e-30);
			spectrumLine.add(hz, db);
			if (o == 0) {
				spectrumLineZoom.add(hz, db);
				spectrumLineZoom2.add(hz, db);
			}
		}
	}
	
	impulseFigure.write(plotName + "-impulse.svg");
	impulseFigure(0, 0).write(plotName + "-impulse-top.svg");
	spectrumFigure.write(plotName + "-spectrum.svg");
}

template<typename Sample>
void plotPhase(std::string plotName) {
	using Complex = std::complex<Sample>;
	signalsmith::blep::EllipticBlep<Sample> blep;
	signalsmith::blep::EllipticBlepAllpass<Sample> allpass;

	size_t impulseLength = 16384;
	std::vector<Complex> impulse(impulseLength);
	std::vector<Complex> spectrum(impulseLength);
	signalsmith::fft2::SimpleFFT<Sample> fft(impulseLength);

	signalsmith::plot::Figure phaseFigure, timeFigure;
	auto &phasePlot = timeFigure(0, 1).plot(600*plotScaleX, 80*plotScaleY);
	phasePlot.y.linear(-M_PI, M_PI).major(0).minor(-M_PI, u8"-π").minor(M_PI, u8"π").label("phase");
	phasePlot.x.linear(0, 0.5).major(0).minor(20000/44100.0, "cutoff").minor(0.5, "Nyquist");
	auto &phaseAmpPlot = phaseFigure(0, -1).plot(600*plotScaleX, 120*plotScaleY);
	phaseAmpPlot.x.copyFrom(phasePlot.x).flip().blankLabels();
	phaseAmpPlot.y.linear(-80, 3).major(0).minors(-20, -40, -60, -80).label("dB");
	auto &groupDelayPlot = timeFigure(0, 2).plot(600*plotScaleX, 150*plotScaleY);
	groupDelayPlot.y.major(0).minor(int(allpass.linearDelay)).label("delay (samples)");
	groupDelayPlot.x.copyFrom(phasePlot.x);
	auto &impulsePlot = timeFigure(0, 0).plot(600*plotScaleX, 80*plotScaleY);
	impulsePlot.y.major(0).minor(1).label("impulse");
	impulsePlot.x.major(0).minor(int(allpass.linearDelay)).label("samples");
	impulsePlot.legend(1, 1).line(0, "with allpass").add(1, "uncorrected");

	phaseAmpPlot.legend(1, 0).add(0, "response").add(2, "difference from pure delay");
	double maxTime = 0, maxTimeAmp = 0;
	for (size_t skipAllpass = 0; skipAllpass < 2; ++skipAllpass) {
		auto &impulseLine = impulsePlot.line();

		allpass.reset();
		int oversampleFactor = 16;
		for (size_t i = 0; i < impulseLength; ++i) {
			Sample x = (i == 0);
			if (!skipAllpass) x = allpass(x);
			blep.add(x, 0);

			Sample v = blep.get();
			impulse[i] = v;
			if (i < allpass.linearDelay*3) {
				for (int o = 0; o < oversampleFactor; ++o) {
					Sample dt = double(o)/oversampleFactor;
					Sample v = blep.get(dt);
					if (skipAllpass && std::abs(v) > maxTimeAmp) {
						maxTimeAmp = std::abs(v);
						maxTime = i + dt;
					}
					impulseLine.add(i + dt, v);
				}
			}
			blep.step();
		}
		maxTime = std::round(maxTime*100)*0.01;
		impulsePlot.x.minor(maxTime);
		groupDelayPlot.y.minor(maxTime);
		
		fft.fft(impulseLength, impulse.data(), spectrum.data());

		auto *phaseLine = &phasePlot.line();
		auto &ampLine = phaseAmpPlot.line(skipAllpass ? 0 : 2);
		auto &delayLine = groupDelayPlot.line();
		Sample prevPhase = 0, prevFreq = -1;
		for (size_t i = 0; i <= impulseLength/2; ++i) {
			Sample freqNorm = Sample(i)/impulseLength;
			
			Sample phase = std::arg(spectrum[i]);
			if (std::abs(phase - prevPhase) > 6) {
				// Jump, by starting a new line
				phaseLine = &phasePlot.line(phaseLine->styleIndex);
			} else {
				Complex v = spectrum[i];
				if (!skipAllpass) {
					Complex expected = std::polar(Sample(1), Sample(-2*M_PI*freqNorm*allpass.linearDelay));
					v -= expected;
				}
				Sample db = 10*std::log10(std::norm(v) + 1e-30);
			
				phaseLine->add(freqNorm, phase);
				Sample groupDelay = (prevPhase - phase)/(freqNorm - prevFreq)/Sample(2*M_PI);
				if (i > 0) delayLine.add(freqNorm, groupDelay);
				ampLine.add(freqNorm, db);
			}
			
			prevPhase = phase;
			prevFreq = freqNorm;
		}
	}

	timeFigure.write(plotName + ".svg");
	phaseFigure.write(plotName + "-error.svg");
}

int main() {
	plotBleps<double>("double-direct", true);
	plotBleps<float>("float-direct", true);
	plotBleps<double>("double-residue", false);
	plotBleps<float>("float-residue", false);
	plotPhase<double>("double-phase");
}
