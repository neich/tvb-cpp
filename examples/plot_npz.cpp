


#include <vector>

#include <simulator/simulate.h>
#include <tools/npz_tools.h>
#include <matplotlibcpp.h>


const double gam = 0.15;

namespace plt = matplotlibcpp;

int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;

    const std::string filename = "Data_Raw/sim_test.npz";
    tvb::StateTrack result;
    result.m_states = tvb::npz2VecMatrixd(filename, "states");
    result.m_times = tvb::npz2VecDouble(filename, "times");

    int N = result.m_states[0].rows();

    size_t t_max = result.m_states.size();

    std::vector<std::vector<double>> y_plot(N, std::vector<double>(t_max));
    for (unsigned t = 0; t < t_max; ++t)
        for (unsigned r = 0; r < N; ++r)
            y_plot[r][t] = result.m_states[t](r, 3);

// Set the size of output image to 1200x780 pixels
    plt::figure_size(1200, 780);
    // Plot line from given x and y data. Color is selected automatically.
    for (unsigned i = 0; i < N; ++i)
        plt::plot(result.m_times, y_plot[i]);
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    plt::title("Sample figure");
    // Enable legend.
    plt::legend();
    // Save the image (file format is determined by the extension)
    plt::save("./sim_test.png");

}
