#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << std::endl;

#include "../elliptic-blep.h"

#include "./simple-fft.h"
#include "plot/plot.h" // https://signalsmith-audio.co.uk/code/plot.git

#include <vector>

template<typename Sample>
void plotBleps(std::string plotName, const int oversample=16, const double sampleRate=48000) {
	using Complex = std::complex<Sample>;
	signalsmith::blep::EllipticBlep<Sample> blep(sampleRate);
	
	constexpr double impulseSeconds = 10;
	size_t impulseLength = 1;
	while (impulseLength < sampleRate*impulseSeconds) impulseLength *= 2;
	std::vector<Complex> impulse(impulseLength);

	signalsmith::plot::Figure impulseFigure, spectrumFigure;
	auto &impulsePlot = impulseFigure(0, 0).plot(600, 150);
	impulsePlot.x.major(0).linear(0, 0.001).minor(0.001, "1ms").label("time");
	impulsePlot.y.major(0).minor(1, " 1").linear(-0.5, 1.5).label("response (48kHz)");
	auto &impulsePlot2 = impulseFigure(0, 1).plot(600, 250);
	impulsePlot2.x.major(0).linear(0, 0.2).minor(0.1, "100ms").minor(0.2, "200ms");
	impulsePlot2.y.major(0).minor(100).minor(-100).linear(-150, 150);
	auto &spectrumPlot = spectrumFigure(0, 1).plot(600, 250);
	spectrumPlot.x.major(0).minor(20000, "20kHz").minor(24000, "24kHz").linear(0, 26000);
	spectrumPlot.y.major(0).minors(-30, -60, -90).label("dB").linear(-105, 5);
	auto &spectrumPlotZoom = spectrumFigure(0, 2).plot(600, 80);
	spectrumPlotZoom.x.copyFrom(spectrumPlot.x);
	spectrumPlotZoom.y.major(0).minors(-0.1).label("ripple").linear(-0.2, 0.1);
	auto &spectrumPlotZoom2 = spectrumFigure(0, 3).plot(600, 80);
	spectrumPlotZoom2.x.linear(0, 100).major(0).minors(20, 40, 60, 80).label("Hz");
	spectrumPlotZoom2.y.linear(-6, 0).major(0).minors(-3, -6).label("highpass");

	std::vector<Complex> spectrum(impulseLength);
	signalsmith::fft2::SimpleFFT<Sample> fft(impulseLength);
	
	auto &legend = impulsePlot.legend(2, 1);
	auto &legend2 = spectrumPlot.legend(2, 1);
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

int main() {
	plotBleps<double>("double");
	plotBleps<float>("float");
}
