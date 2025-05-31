#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <vector>
#include <complex>

#include <string>
#include <algorithm>
#include <iostream>

#define _XOPEN_SOURCE 700
#include <matplot/matplot.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;
namespace mtp = matplot;

std::vector<double> generate_sin(double freq, double duration, double sampling_rate) {
    int num_samples = static_cast<int>(duration * sampling_rate);
    std::vector<double> signal(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        double t = i / sampling_rate;
        signal[i] = sin(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generate_cos(double freq, double duration, double sampling_rate) {
    int num_samples = static_cast<int>(duration * sampling_rate);
    std::vector<double> signal(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        double t = i / sampling_rate;
        signal[i] = cos(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generate_square(double freq, double duration, double sampling_rate) {
    int num_samples = static_cast<int>(duration * sampling_rate);
    std::vector<double> signal(num_samples);
    double period = 1.0 / freq;
    for (int i = 0; i < num_samples; ++i) {
        double t = i / sampling_rate;
        signal[i] = (std::fmod(t, period) < period / 2) ? 1.0 : -1.0;
    }
    return signal;
}

std::vector<double> generate_sawtooth(double freq, double duration, double sampling_rate) {
    int num_samples = static_cast<int>(duration * sampling_rate);
    std::vector<double> signal(num_samples);
    double period = 1.0 / freq;
    for (int i = 0; i < num_samples; ++i) {
        double t = i / sampling_rate;
        signal[i] = 2 * (t / period - std::floor(0.5 + t / period));
    }
    return signal;
}

std::vector<std::complex<double>> dft(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<std::complex<double>> result(N);
    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            sum += signal[n] * std::exp(std::complex<double>(0, angle));
        }
        result[k] = sum;
    }
    return result;
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& dft_coeff) {
    int N = dft_coeff.size();
    std::vector<std::complex<double>> result(N);
    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += dft_coeff[k] * std::exp(std::complex<double>(0, angle));
        }
        sum /= N;
        result[n] = sum;
    }
    return result;
}

void visualize_signal_and_dft(
    const std::vector<double>& signal,
    const std::vector<double>& time,
    const std::vector<std::complex<double>>& dft_result,
    const std::string& filename
) {
    auto f = mtp::figure(true);
    f->size(1200, 800);
    
    auto ax1 = f->add_axes({0.1, 0.55, 0.8, 0.4});
    mtp::plot(ax1, time, signal, "b-")->line_width(1.5);
    mtp::title(ax1, "Sygnal w dziedzinie czasu");
    mtp::xlabel(ax1, "Czas [s]");
    mtp::ylabel(ax1, "Amplituda");
    mtp::grid(ax1, true);
    
    std::vector<double> magnitude(dft_result.size());
    std::transform(dft_result.begin(), dft_result.end(), magnitude.begin(), [](const std::complex<double>& c) { return std::abs(c); });
    double fs = 1.0 / (time[1] - time[0]);
    std::vector<double> freqs(magnitude.size());
    for (size_t i = 0; i < freqs.size(); ++i) {
        freqs[i] = i * fs / magnitude.size();
    }
    
    size_t nyquist = magnitude.size() / 2;
    auto ax2 = f->add_axes({0.1, 0.1, 0.8, 0.4});
    mtp::stem(ax2, std::vector<double>(freqs.begin(), freqs.begin() + nyquist), std::vector<double>(magnitude.begin(), magnitude.begin() + nyquist));
    mtp::title(ax2, "Widmo czestotliwosciowe (DFT)");
    mtp::xlabel(ax2, "Czestotliwosc [Hz]");
    mtp::ylabel(ax2, "Amplituda");
    mtp::grid(ax2, true);
    mtp::xlim(ax2, {0, fs/2});
    
    mtp::save(filename);
}

void visualize_dft_idft(
    const std::vector<double>& original_signal,
    const std::vector<std::complex<double>>& dft_result,
    const std::vector<std::complex<double>>& reconstructed_signal,
    const std::string& filename
) {
    auto f = mtp::figure(true);
    f->size(1200, 800);
    f->title("DFT i Rekonstrukcja sygnalu");
    
    auto ax1 = f->add_axes({0.1, 0.7, 0.8, 0.2});
    mtp::plot(ax1, original_signal, "b-")->line_width(1.5).display_name("Oryginal");
    mtp::title(ax1, "Sygnal oryginalny");
    mtp::grid(ax1, true);
    
    std::vector<double> magnitude(dft_result.size());
    std::transform(dft_result.begin(), dft_result.end(), magnitude.begin(), [](const std::complex<double>& c) { return std::abs(c); });
    
    auto ax2 = f->add_axes({0.1, 0.4, 0.8, 0.2});
    mtp::stem(ax2, magnitude, "r-")->line_width(1).display_name("Amplituda DFT");
    mtp::title(ax2, "Widmo amplitudowe (DFT)");
    mtp::grid(ax2, true);
    
    std::vector<double> real_part(reconstructed_signal.size());
    std::transform(reconstructed_signal.begin(), reconstructed_signal.end(), real_part.begin(), [](const std::complex<double>& c) { return c.real(); });
    
    auto ax3 = f->add_axes({0.1, 0.1, 0.8, 0.2});
    mtp::plot(ax3, real_part, "g-")->line_width(1.5).display_name("Rekonstrukcja");
    mtp::title(ax3, "Sygnał zrekonstruowany (IDFT)");
    mtp::grid(ax3, true);
    
    mtp::save(filename);
}

std::vector<std::complex<double>> low_pass_filter_1d(
    const std::vector<std::complex<double>>& dft_data,
    double sampling_rate,
    double cutoff_freq)
{
    size_t N = dft_data.size();
    std::vector<std::complex<double>> filtered(dft_data);
    double freq_resolution = sampling_rate / N;

    for (size_t k = 0; k < N; ++k) {
        double freq = k * freq_resolution;
        if (freq > cutoff_freq) {
            filtered[k] = 0;
        }
    }
    return filtered;
}

void visualize_filtered_signal_1d(
    const std::vector<double>& original_signal,
    const std::vector<std::complex<double>>& filtered_dft,
    const std::vector<std::complex<double>>& reconstructed_signal,
    const std::string& filename)
{
    auto f = mtp::figure(true);
    f->size(1200, 800);
    f->title("Filtracja dolnoprzepustowa 1D");

    auto ax1 = f->add_axes({0.1, 0.7, 0.8, 0.2});
    mtp::plot(ax1, original_signal)->line_width(1.5);
    mtp::title(ax1, "Sygnal oryginalny");
    mtp::grid(ax1, true);

    std::vector<double> magnitude(filtered_dft.size());
    std::transform(filtered_dft.begin(), filtered_dft.end(), magnitude.begin(), [](const std::complex<double>& c) { return std::abs(c); });

    auto ax2 = f->add_axes({0.1, 0.4, 0.8, 0.2});
    mtp::stem(ax2, magnitude)->line_width(1.0);
    mtp::title(ax2, "Widmo po filtracji (DFT)");
    mtp::grid(ax2, true);

    std::vector<double> real_part(reconstructed_signal.size());
    std::transform(reconstructed_signal.begin(), reconstructed_signal.end(), real_part.begin(), [](const std::complex<double>& c) { return c.real(); });

    auto ax3 = f->add_axes({0.1, 0.1, 0.8, 0.2});
    mtp::plot(ax3, real_part)->line_width(1.5);
    mtp::title(ax3, "Sygnał po filtrze (IDFT)");
    mtp::grid(ax3, true);

    mtp::save(filename);
}

PYBIND11_MODULE(_core, m) {
    m.def("generate_sin", &generate_sin, py::arg("freq"), py::arg("duration"), py::arg("sampling_rate"));
    m.def("generate_cos", &generate_cos, py::arg("freq"), py::arg("duration"), py::arg("sampling_rate"));
    m.def("generate_square", &generate_square, py::arg("freq"), py::arg("duration"), py::arg("sampling_rate"));
    m.def("generate_sawtooth", &generate_sawtooth, py::arg("freq"), py::arg("duration"), py::arg("sampling_rate"));
    m.def("dft", &dft, py::arg("signal"));
    m.def("idft", &idft, py::arg("dft_coeff"));
    m.def("visualize_signal_and_dft", &visualize_signal_and_dft, py::arg("signal"), py::arg("time"), py::arg("dft_result"), py::arg("filename"));
    m.def("visualize_dft_idft", &visualize_dft_idft, py::arg("original_signal"), py::arg("dft_result"), py::arg("reconstructed_signal"), py::arg("filename"));
    m.def("visualize_filtered_signal_1d", &visualize_filtered_signal_1d, py::arg("original_signal"), py::arg("filtered_dft"), py::arg("reconstructed_signal"), py::arg("filename"));
    m.def("low_pass_filter_1d", &low_pass_filter_1d, py::arg("dft_data"), py::arg("sampling_rate"), py::arg("cutoff_freq"));
}