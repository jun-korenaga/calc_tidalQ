/*
 * calc_tidalQ.cc (radial tidal solver with the Sabadini et al. 2016 formulation)
 *
 * Jun Korenaga
 * January 2025
 */

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>
#include "tidal_sys_relaxS.h"

namespace {

void usage(std::ostream& out)
{
    out << "Usage: calc_tidalQ -W frequency -M order -S model "
           "[-I] [-R spacing] [-D dump_root] [-V]\n"
           "  -W  Angular forcing frequency (rad/s)\n"
           "  -M  Azimuthal order (0, 1, or 2) for degree 2\n"
           "  -S  Mantle model file (radius, density, bulk modulus, "
           "real and imaginary shear modulus in SI units)\n"
           "  -I  Disable mantle and core-boundary inertia terms\n"
           "  -R  Target radial grid spacing (m, default 1000)\n"
           "  -D  Write the radial solution to dump_root.01.dat\n"
           "  -V  Print grid information to standard error\n"
           "  -h, --help  Show this help\n";
}

bool readNumber(const char* text, double& value)
{
    char* end = nullptr;
    errno = 0;
    value = std::strtod(text, &end);
    return end != text && *end == '\0' && errno != ERANGE && std::isfinite(value);
}

} // namespace

int main(int argc, char** argv)
{
    bool got_frequency = false, inertia_off = false, verbose = false;
    const char* model = nullptr;
    const char* dump_root = nullptr;
    int order = -1;
    double frequency = 0.0;
    double spacing = 1000.0;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            usage(std::cout);
            return 0;
        }
        if (arg == "-I") {
            inertia_off = true;
            continue;
        }
        if (arg == "-V") {
            verbose = true;
            continue;
        }
        if (arg.size() < 2 || arg[0] != '-' ||
            std::string_view("WMSRD").find(arg[1]) == std::string_view::npos) {
            std::cerr << "Unknown option: " << arg << '\n';
            return 1;
        }
        const char* value = argv[i] + 2;
        if (*value == '\0') {
            if (i + 1 == argc || std::string_view(argv[i + 1]).empty() ||
                ((arg[1] == 'S' || arg[1] == 'D') && argv[i + 1][0] == '-')) {
                std::cerr << "Missing value for " << arg << '\n';
                return 1;
            }
            value = argv[++i];
        }

        switch (arg[1]) {
        case 'W':
            if (!readNumber(value, frequency) || !std::isfinite(frequency * frequency)) {
                std::cerr << "-W requires a finite angular frequency in rad/s.\n";
                return 1;
            }
            got_frequency = true;
            break;
        case 'M': {
            std::string_view order_text = value;
            if (order_text.size() != 1 || order_text[0] < '0' || order_text[0] > '2') {
                std::cerr << "-M requires order 0, 1, or 2.\n";
                return 1;
            }
            order = order_text[0] - '0';
            break;
        }
        case 'S':
            model = value;
            break;
        case 'R':
            if (!readNumber(value, spacing) || spacing <= 0.0) {
                std::cerr << "-R requires a positive grid spacing in meters.\n";
                return 1;
            }
            break;
        case 'D':
            dump_root = value;
            break;
        }
    }

    if (!got_frequency || order < 0 || model == nullptr) {
        usage(std::cerr);
        std::cerr << "-W, -M, and -S are required.\n";
        return 1;
    }
    if (!inertia_off && frequency * frequency == 0.0) {
        std::cerr << "Zero angular frequency requires -I.\n";
        return 1;
    }

    try {
        EarthProperties earth(spacing, spacing, model);
        if (verbose) {
            std::cerr << "nr = " << earth.nGrid() << ", dr = "
                      << (RE - Rc) / (earth.nGrid() - 1) << '\n';
        }

        TidalSystem sys(earth);
        if (dump_root != nullptr) sys.dumpSolutions(dump_root);
        if (inertia_off) sys.inertiaOff();

        std::complex<double> h2, k2;
        double delta_h, delta_k, tidalQ_h, tidalQ_k, residual;
        int solves = sys.calcLoveNumbers(frequency, order, h2, k2,
                                        delta_h, delta_k, tidalQ_h, tidalQ_k,
                                        residual);
        std::cout << frequency << " " << order << " "
                  << tidalQ_h << " " << tidalQ_k << " "
                  << h2.real() << " " << h2.imag() << " "
                  << k2.real() << " " << k2.imag() << " "
                  << delta_h << " " << delta_k << " "
                  << residual << " " << solves << '\n';
    } catch (const std::exception& e) {
        std::cerr << "calc_tidalQ: " << e.what() << '\n';
        return 1;
    }
}
