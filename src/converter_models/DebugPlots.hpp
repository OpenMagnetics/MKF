#pragma once

#ifdef DEBUG_PLOTS
#include <matplot/matplot.h>
#include <filesystem>
#include <string>
#include <vector>
#include <complex>

namespace debug_plots {
    using namespace matplot;

    static std::string base_folder = "debug_plots";

    inline void init_folder() {
        std::filesystem::create_directories(base_folder);
    }

    inline std::string path(const std::string& filename) {
        return base_folder + "/" + filename;
    }

    inline void plot_fft_db(const std::vector<double>& freq,
                        const std::vector<double>& spectrum,
                        const std::string& label,
                        const std::string& filename) {
    // Compute magnitude in dB
    std::vector<double> spectrum_dB(spectrum.size());
    std::transform(spectrum.begin(), spectrum.end(), spectrum_dB.begin(),
                   [](double v) {
                       double mag = std::abs(v);
                       return (mag > 1e-12) ? 20.0 * std::log10(mag) : -120.0; // floor
                   });

    figure(true);
    semilogx(freq, spectrum_dB)->display_name(label); // log-x axis
    xlabel("Frequency [Hz]");
    ylabel("Amplitude [dB]");
    title("FFT Spectrum (" + label + ")");
    grid(true);
    save(path(filename));
    }

    // === 1. Carrier Plot (only a few switching periods) ===
    inline void plot_carrier(const std::vector<double>& carrier,
                             double fsw) {
        auto t = linspace(0, 5.0 / fsw, carrier.size());
        figure(true);
        plot(t, carrier)->display_name("Carrier");
        xlabel("Time [s]");
        ylabel("Amplitude");
        title("Carrier over 5 switching periods");
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("carrier.png"));
    }

    // === 2. Carrier vs Reference (1 fundamental period) ===
    inline void plot_carrier_vs_ref(const std::vector<double>& carrier,
                                    const std::vector<double>& reference,
                                    double f1) {
        auto t = linspace(0, 1.0 / f1, carrier.size());
        figure(true);
        hold(true);
        plot(t, carrier)->display_name("Carrier");
        plot(t, reference)->display_name("Reference");
        legend();
        xlabel("Time [s]");
        ylabel("Amplitude");
        title("Carrier vs Reference");
        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("carrier_vs_ref.png"));
    }

        // === 2. Carrier vs References (1 fundamental period, all 3 phases) ===
    inline void plot_carrier_vs_refs(const std::vector<double>& carrier,
                                     const std::vector<double>& refA,
                                     const std::vector<double>& refB,
                                     const std::vector<double>& refC,
                                     double f1) {
        // Make sure all vectors are same length
        size_t N = std::min({carrier.size(), refA.size(), refB.size(), refC.size()});
        auto t = linspace(0, 1.0 / f1, N);

        figure(true);
        hold(true);

        plot(t, carrier)->display_name("Carrier");
        plot(t, refA)->display_name("Reference A");
        plot(t, refB)->display_name("Reference B");
        plot(t, refC)->display_name("Reference C");

        legend();
        xlabel("Time [s]");
        ylabel("Normalized amplitude");
        title("Carrier vs References (1 fundamental period)");
        grid(true);

        hold(false);
        gcf()->size(1920, 1080);   // 4K output
        save(path("carrier_vs_refs.png"));
    }

    // === 3. Comparator (PWM) outputs ===
    inline void plot_pwm_signals(const std::vector<int>& gateA,
                                const std::vector<int>& gateB,
                                const std::vector<int>& gateC,
                                double f1, double fs) {
        const double Tfund = 1.0 / f1;
        size_t Nfund = static_cast<size_t>(std::round(Tfund * fs));
        Nfund = std::min({Nfund, gateA.size(), gateB.size(), gateC.size()});

        std::vector<int> A(gateA.begin(), gateA.begin() + Nfund);
        std::vector<int> B(gateB.begin(), gateB.begin() + Nfund);
        std::vector<int> C(gateC.begin(), gateC.begin() + Nfund);

        auto t = linspace(0, Tfund, Nfund);
        figure(true);
        hold(true);
        stairs(t, A)->display_name("Gate A");
        stairs(t, B)->display_name("Gate B");
        stairs(t, C)->display_name("Gate C");
        legend();
        xlabel("Time [s]");
        ylabel("Gate signal (0/1)");
        title("PWM Comparator Outputs (1 fundamental period)");
        grid(true);
        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("pwm_signals.png"));
    }

    // === 3. Output comparison (va,vb,vc over 2 fsw) ===
    inline void plot_va_vb_vc_short(const std::vector<double>& va,
                                    const std::vector<double>& vb,
                                    const std::vector<double>& vc,
                                    double fsw,
                                    double fs) {
        const double Tshort = 2.0 / fsw;
        size_t Nshort = static_cast<size_t>(std::round(Tshort * fs));
        Nshort = std::min({Nshort, va.size(), vb.size(), vc.size()});

        std::vector<double> va_short(va.begin(), va.begin() + Nshort);
        std::vector<double> vb_short(vb.begin(), vb.begin() + Nshort);
        std::vector<double> vc_short(vc.begin(), vc.begin() + Nshort);

        auto t = linspace(0, Tshort, Nshort);
        figure(true);
        hold(true);
        plot(t, va_short)->display_name("va");
        plot(t, vb_short)->display_name("vb");
        plot(t, vc_short)->display_name("vc");
        legend();
        xlabel("Time [s]");
        ylabel("Voltage [V]");
        title("va, vb, vc (2 switching periods)");
        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("va_vb_vc_short.png"));
    }

    // === 4. Output comparison (va,vb,vc over 1 fundamental) ===
    inline void plot_va_vb_vc_fundamental(const std::vector<double>& va,
                                        const std::vector<double>& vb,
                                        const std::vector<double>& vc,
                                        double f1,
                                        double fs) {
        const double Tfund = 1.0 / f1;
        size_t Nfund = static_cast<size_t>(std::round(Tfund * fs));
        Nfund = std::min({Nfund, va.size(), vb.size(), vc.size()});

        std::vector<double> va_fund(va.begin(), va.begin() + Nfund);
        std::vector<double> vb_fund(vb.begin(), vb.begin() + Nfund);
        std::vector<double> vc_fund(vc.begin(), vc.begin() + Nfund);

        auto t = linspace(0, Tfund, Nfund);
        figure(true);
        hold(true);
        plot(t, va_fund)->display_name("va");
        plot(t, vb_fund)->display_name("vb");
        plot(t, vc_fund)->display_name("vc");
        legend();
        xlabel("Time [s]");
        ylabel("Voltage [V]");
        title("va, vb, vc (1 fundamental period)");
        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("va_vb_vc_fundamental.png"));
    }

    inline void plot_fft_vl1_il1(const std::vector<double>& freq,
                                const std::vector<double>& Vl1_fft,
                                const std::vector<double>& Il1_fft) {
        figure(true);
        hold(true);

        std::vector<double> Vl1_dB(Vl1_fft.size());
        std::transform(Vl1_fft.begin(), Vl1_fft.end(), Vl1_dB.begin(),
                    [](double v){ return 20.0*std::log10(std::max(v,1e-12)); });

        std::vector<double> Il1_dB(Il1_fft.size());
        std::transform(Il1_fft.begin(), Il1_fft.end(), Il1_dB.begin(),
                    [](double v){ return 20.0*std::log10(std::max(v,1e-12)); });

        semilogx(freq, Vl1_dB)->display_name("|V_L1(f)| [dB]");
        semilogx(freq, Il1_dB)->display_name("|I_L1(f)| [dB]");

        legend();
        xlabel("Frequency [Hz]");
        ylabel("Amplitude [dB]");
        title("FFT of vL1 and iL1");
        grid(true);

        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("fft_vl1_il1.png"));
    }

    // === 6. Instantaneous power p(t) over 1 fundamental + diagnostics ===
    inline void plot_power(const std::vector<double>& p,
                        double f1) {
        auto t = linspace(0, 1.0 / f1, p.size());
        figure(true);
        plot(t, p);
        xlabel("Time [s]");
        ylabel("Power [W]");
        title("Instantaneous Power p(t)");
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("power.png"));

        // --- Diagnostics ---
        double sum = std::accumulate(p.begin(), p.end(), 0.0);
        double p_avg = sum / static_cast<double>(p.size());

        double sq_sum = std::inner_product(p.begin(), p.end(), p.begin(), 0.0);
        double p_rms = std::sqrt(sq_sum / static_cast<double>(p.size()));

        auto [min_it, max_it] = std::minmax_element(p.begin(), p.end());
        double p_min = *min_it;
        double p_max = *max_it;

        std::cout << "[Power diagnostics] "
                << "avg = " << p_avg << " W, "
                << "RMS = " << p_rms << " W, "
                << "min = " << p_min << " W, "
                << "max = " << p_max << " W"
                << std::endl;
    }

    // === 7. vdc(t) ripple over fundamental ===
    inline void plot_vdc_ripple(const std::vector<double>& vdc,
                                double f1) {
        auto t = linspace(0, 1.0 / f1, vdc.size());
        figure(true);
        plot(t, vdc);
        xlabel("Time [s]");
        ylabel("Vdc [V]");
        title("DC Bus Voltage with Ripple");
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("vdc_ripple.png"));
    }

    // === 8. Final FFT result (vL1 & iL1 again after ripple) ===
    inline void plot_final_fft_vl1_il1(const std::vector<double>& freq,
                                    const std::vector<double>& Vl1_fft,
                                    const std::vector<double>& Il1_fft) {
        figure(true);
        hold(true);

        // Convert magnitude → dB (20*log10), clamp small values
        std::vector<double> Vl1_dB(Vl1_fft.size());
        std::transform(Vl1_fft.begin(), Vl1_fft.end(), Vl1_dB.begin(),
                    [](double v) {
                        double mag = std::abs(v);
                        return (mag > 1e-20) ? 20.0 * std::log10(mag) : -200.0; // dB floor
                    });

        std::vector<double> Il1_dB(Il1_fft.size());
        std::transform(Il1_fft.begin(), Il1_fft.end(), Il1_dB.begin(),
                    [](double v) {
                        double mag = std::abs(v);
                        return (mag > 1e-20) ? 20.0 * std::log10(mag) : -200.0;
                    });

        // Plot with log frequency axis
        semilogx(freq, Vl1_dB)->display_name("|V_L1(f)| final [dB]");
        semilogx(freq, Il1_dB)->display_name("|I_L1(f)| final [dB]");

        legend();
        xlabel("Frequency [Hz]");
        ylabel("Amplitude [dB]");
        title("Final FFT of vL1 and iL1");
        grid(true);

        hold(false);
        gcf()->size(1920, 1080);   // set figure size (in pixels)
        save(path("final_fft_vl1_il1.png"));
    }
} // namespace debug_plots
#endif