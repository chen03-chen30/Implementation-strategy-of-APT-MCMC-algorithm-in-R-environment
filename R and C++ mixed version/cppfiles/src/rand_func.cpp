#include <Rcpp.h> // 包含 Rcpp 以便使用 Rcpp::warning 或其他 Rcpp 功能（如果需要）
#include <Rmath.h>
#include <random> // C++11 随机数库
#include <cmath>  // std::sqrt, std::log, std::cos, std::exp
#include <omp.h>  // For omp_get_thread_num() to create unique seeds
#include <limits> // For std::numeric_limits
#include <chrono> // For an additional time-based seed component (optional)

// --- Thread-Safe Random Number Generation Setup ---

// 1. Define the Engine type we'll use
using Engine = std::mt19937;

// 2. Create a thread_local engine instance. Each thread gets its own copy.
//    Use a function static to ensure initialization happens correctly across threads.
Engine& get_thread_local_rng() {
    // thread_local ensures each thread has its own 'engine' variable.
    // It's initialized only once per thread.
    thread_local static Engine engine = []() {
        // Seeding logic: Executed once per thread when the engine is first created.
        // Combine multiple sources for a good seed:
        // - std::random_device: Non-deterministic source (if available)
        // - Thread ID: Ensures different threads get different starting points
        // - Time: Adds more variability, especially if threads start close together
        std::random_device rd;
        unsigned int seed = rd() +
                            static_cast<unsigned int>(omp_get_thread_num()) +
                            static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        return Engine(seed);
    }(); // Immediately invoked lambda to initialize
    return engine;
}

// --- Replacement Random Functions ---
// Note: Removed [[Rcpp::export]] from internal RNG functions,
// assuming they are only called from your C++ MCMC code.
// If you need to test them directly from R, you'd need a different approach
// or understand they won't be parallel in that R context.

/*
  Purpose: r8_normal_01() returns a unit pseudonormal R8 (thread-safe).
*/

/// 
// [[Rcpp::export]]
double r8_normal_01() {
    // Get the engine local to the current thread
    Engine& rng = get_thread_local_rng();
    // Create the appropriate distribution (mean 0, stddev 1)
    std::normal_distribution<double> dist(0.0, 1.0);
    // Generate and return the random number using the thread's engine
    return dist(rng);
}

/*
  Purpose: r8_normal_ab() returns a scaled pseudonormal R8 (thread-safe).
  Input: a = mean, b = standard deviation.
*/

/// 
// [[Rcpp::export]]
double r8_normal_ab(double a, double b) {
    // We can directly use the normal_distribution with specified mean and stddev
    Engine& rng = get_thread_local_rng();
    if (b < 0) {
         // Standard deviation cannot be negative
         Rcpp::warning("Negative standard deviation provided to r8_normal_ab: b=%f. Using abs(b).", b);
         b = std::abs(b);
    }
    std::normal_distribution<double> dist(a, b);
    return dist(rng);
    // Or, stick closer to the original implementation's logic:
    // return a + b * r8_normal_01(); // This also works and uses the thread-safe r8_normal_01
}

/*
  Purpose: r8_logunif_ab() returns a log uniform pseudo random number R8 (thread-safe).
  Input: a = lower limit, b = upper limit (must be > 0, a < b).
*/
/// 
// [[Rcpp::export]]
double r8_logunif_ab(double a, double b) {
    Engine& rng = get_thread_local_rng();
    if (a <= 0 || b <= 0 || a >= b) {
        Rcpp::warning("Invalid range for r8_logunif_ab: a=%f, b=%f. Requires 0 < a < b.", a, b);
        return std::numeric_limits<double>::quiet_NaN(); // Or throw error, or return default
    }
    // Generate a standard uniform number [0, 1) first
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double u01 = uniform_dist(rng);
    // Apply inverse transform sampling formula
    return a * std::exp(u01 * std::log(b / a));
}

/*
  Purpose: r8_unif_ab() returns a uniform pseudo random number R8 between a and b (thread-safe).
  Input: a = lower limit, b = upper limit. Range is [a, b).
*/
/// 
// [[Rcpp::export]]
double r8_unif_ab(double a, double b) {
    Engine& rng = get_thread_local_rng();
    if (a > b) { // Ensure correct order for distribution
        std::swap(a,b);
    } else if (a == b) {
        return a; // If a equals b, return that value
    }
    // std::uniform_real_distribution generates in [a, b)
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

/*
  Purpose: i4_unif_ab() returns a uniform pseudo random integer i4 between a and b (thread-safe).
  Input: a = lower limit, b = upper limit. Range is [a, b] (inclusive).
*/
/// 
// [[Rcpp::export]]
int i4_unif_ab(int a, int b) {
    Engine& rng = get_thread_local_rng();
    if (a > b) { // Ensure correct order
        std::swap(a, b);
    }
    // std::uniform_int_distribution generates in [a, b] (inclusive)
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

/*
  Purpose: i4_unif_0a() returns a rand int between i=0 to i=a-1 (thread-safe).
           Matches the floor() behavior of the original drand48 implementation.
  Input: a = upper limit (exclusive). Range is [0, a-1].
*/
/// 
// [[Rcpp::export]]
int i4_unif_0a(int a) {
    Engine& rng = get_thread_local_rng();
    if (a <= 0) {
        Rcpp::warning("Non-positive upper bound in i4_unif_0a: a=%d. Returning 0.", a);
        return 0; // Or throw error
    }
    // We want the range [0, a-1]
    std::uniform_int_distribution<int> dist(0, a - 1);
    return dist(rng);
}


// --- Potentially Exported Function (Already Thread-Safe) ---

/*
  Purpose: Calculates the Probability Density Function (PDF) of the normal distribution.
           This function is purely mathematical and was already thread-safe.
*/
/// 
// [[Rcpp::export]] // Keep export if you need to call this directly from R
double normal_pdf(double x, double mean, double sd) {
    static const double inv_sqrt_2pi = 0.3989422804014327; // 1 / sqrt(2*pi)
    if (sd <= 0) {
        // PDF is undefined or problematic for non-positive standard deviation
        Rcpp::warning("Non-positive standard deviation in normal_pdf: sd=%f", sd);
        if (x == mean) return std::numeric_limits<double>::infinity(); // Dirac delta at mean if sd=0
        else return 0.0; // Zero density elsewhere if sd=0
        // Or return NaN for sd < 0
        // return std::numeric_limits<double>::quiet_NaN();
    }
    double a = (x - mean) / sd;
    return inv_sqrt_2pi / sd * std::exp(-0.5 * a * a);
}
