#include "core/Distribution.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

namespace stochia {
namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kEpsilon = 1e-13;

double clampProbability(double value) {
    return std::clamp(value, 0.0, 1.0);
}

double regularizedGammaP(double a, double x) {
    if (x <= 0.0) return 0.0;
    if (x < a + 1.0) {
        double sum = 1.0 / a;
        double term = sum;
        double ap = a;
        for (int n = 1; n <= 200; ++n) {
            ap += 1.0;
            term *= x / ap;
            sum += term;
            if (std::abs(term) < std::abs(sum) * kEpsilon) break;
        }
        return clampProbability(sum * std::exp(-x + a * std::log(x) - std::lgamma(a)));
    }

    double b = x + 1.0 - a;
    double c = 1.0 / std::numeric_limits<double>::min();
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= 200; ++i) {
        const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = b + an / c;
        if (std::abs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::abs(delta - 1.0) < kEpsilon) break;
    }
    const double q = std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
    return clampProbability(1.0 - q);
}

double betaContinuedFraction(double a, double b, double x) {
    constexpr int maxIterations = 250;
    constexpr double fpMin = 1e-300;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < fpMin) d = fpMin;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= maxIterations; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpMin) d = fpMin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpMin) c = fpMin;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpMin) d = fpMin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpMin) c = fpMin;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::abs(delta - 1.0) < kEpsilon) break;
    }
    return h;
}

double regularizedBeta(double x, double a, double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double factor = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                                   + a * std::log(x) + b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return clampProbability(factor * betaContinuedFraction(a, b, x) / a);
    return clampProbability(1.0 - factor * betaContinuedFraction(b, a, 1.0 - x) / b);
}

double normalCdf(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

double p(const std::map<std::string, double>& values, const std::string& key) {
    const auto it = values.find(key);
    return it == values.end() ? 0.0 : it->second;
}

std::string number(double value) {
    std::ostringstream out;
    out.precision(6);
    out << value;
    return out.str();
}

class GenericDistribution final : public Distribution {
public:
    GenericDistribution(std::string id, std::map<std::string, double> values)
        : id_(std::move(id)), values_(std::move(values)) {}

    std::string id() const override { return id_; }
    std::map<std::string, double> parameters() const override { return values_; }

    std::string displayName() const override {
        for (const auto& spec : distributionCatalog())
            if (spec.id == id_) return spec.displayName;
        return id_;
    }

    DistributionType type() const override {
        for (const auto& spec : distributionCatalog())
            if (spec.id == id_) return spec.type;
        return DistributionType::Continuous;
    }

    std::pair<double, double> plotRange() const override {
        if (id_ == "bernoulli") return {-0.5, 1.5};
        if (id_ == "binomial") return {-0.5, p(values_, "n") + 0.5};
        if (id_ == "geometric") {
            const double q = std::max(0.001, p(values_, "p"));
            return {0.5, std::max(5.5, std::ceil(std::log(0.002) / std::log(1.0 - q)) + 0.5)};
        }
        if (id_ == "poisson") {
            const double lambda = p(values_, "lambda");
            return {-0.5, std::ceil(lambda + 5.0 * std::sqrt(lambda) + 2.0) + 0.5};
        }
        if (id_ == "negative_binomial") {
            const double m = mean(), s = std::sqrt(variance());
            return {-0.5, std::ceil(m + 5.0 * s + 2.0) + 0.5};
        }
        if (id_ == "uniform") return {p(values_, "a"), p(values_, "b")};
        if (id_ == "normal") {
            const double m = p(values_, "mu"), s = p(values_, "sigma");
            return {m - 4.0 * s, m + 4.0 * s};
        }
        if (id_ == "exponential") return {0.0, 7.0 / p(values_, "lambda")};
        if (id_ == "gamma") {
            const double m = mean(), s = std::sqrt(variance());
            return {0.0, m + 5.0 * s};
        }
        if (id_ == "beta") return {0.0, 1.0};
        if (id_ == "chi_squared") {
            const double m = mean(), s = std::sqrt(variance());
            return {0.0, m + 5.0 * s};
        }
        if (id_ == "student_t") return {-5.0, 5.0};
        if (id_ == "f") return {0.0, std::min(15.0, mean() + 5.0 * std::sqrt(std::max(variance(), 1.0)))};
        return {-5.0, 5.0};
    }

    double density(double x) const override {
        if (id_ == "bernoulli") {
            if (std::abs(x) < 1e-9) return 1.0 - p(values_, "p");
            if (std::abs(x - 1.0) < 1e-9) return p(values_, "p");
            return 0.0;
        }
        if (id_ == "binomial") {
            const int n = static_cast<int>(p(values_, "n"));
            const int k = static_cast<int>(std::llround(x));
            if (k < 0 || k > n || std::abs(x - k) > 1e-9) return 0.0;
            const double prob = p(values_, "p");
            return std::exp(std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0)
                            + k * std::log(prob) + (n - k) * std::log1p(-prob));
        }
        if (id_ == "geometric") {
            const int k = static_cast<int>(std::llround(x));
            if (k < 1 || std::abs(x - k) > 1e-9) return 0.0;
            const double prob = p(values_, "p");
            return prob * std::pow(1.0 - prob, k - 1);
        }
        if (id_ == "poisson") {
            const int k = static_cast<int>(std::llround(x));
            if (k < 0 || std::abs(x - k) > 1e-9) return 0.0;
            const double lambda = p(values_, "lambda");
            return std::exp(k * std::log(lambda) - lambda - std::lgamma(k + 1.0));
        }
        if (id_ == "negative_binomial") {
            const int k = static_cast<int>(std::llround(x));
            const int r = static_cast<int>(p(values_, "r"));
            if (k < 0 || std::abs(x - k) > 1e-9) return 0.0;
            const double prob = p(values_, "p");
            return std::exp(std::lgamma(k + r) - std::lgamma(k + 1.0) - std::lgamma(r)
                            + r * std::log(prob) + k * std::log1p(-prob));
        }
        if (id_ == "uniform") {
            const double a = p(values_, "a"), b = p(values_, "b");
            return x >= a && x <= b ? 1.0 / (b - a) : 0.0;
        }
        if (id_ == "normal") {
            const double mu = p(values_, "mu"), sigma = p(values_, "sigma");
            const double z = (x - mu) / sigma;
            return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * kPi));
        }
        if (id_ == "exponential") {
            const double lambda = p(values_, "lambda");
            return x < 0.0 ? 0.0 : lambda * std::exp(-lambda * x);
        }
        if (id_ == "gamma") {
            const double a = p(values_, "shape"), theta = p(values_, "scale");
            if (x <= 0.0) return 0.0;
            return std::exp((a - 1.0) * std::log(x) - x / theta - std::lgamma(a) - a * std::log(theta));
        }
        if (id_ == "beta") {
            const double a = p(values_, "alpha"), b = p(values_, "beta");
            if (x <= 0.0 || x >= 1.0) return 0.0;
            return std::exp((a - 1.0) * std::log(x) + (b - 1.0) * std::log1p(-x)
                            + std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b));
        }
        if (id_ == "chi_squared") {
            const double k = p(values_, "df");
            if (x <= 0.0) return 0.0;
            return std::exp((k / 2.0 - 1.0) * std::log(x) - x / 2.0
                            - (k / 2.0) * std::log(2.0) - std::lgamma(k / 2.0));
        }
        if (id_ == "student_t") {
            const double nu = p(values_, "df");
            return std::exp(std::lgamma((nu + 1.0) / 2.0) - std::lgamma(nu / 2.0))
                   / std::sqrt(nu * kPi) * std::pow(1.0 + x * x / nu, -(nu + 1.0) / 2.0);
        }
        if (id_ == "f") {
            const double d1 = p(values_, "d1"), d2 = p(values_, "d2");
            if (x <= 0.0) return 0.0;
            const double half1 = d1 / 2.0, half2 = d2 / 2.0;
            return std::exp(half1 * std::log(d1 / d2) + (half1 - 1.0) * std::log(x)
                            - (half1 + half2) * std::log1p(d1 * x / d2)
                            + std::lgamma(half1 + half2) - std::lgamma(half1) - std::lgamma(half2));
        }
        return 0.0;
    }

    double cdf(double x) const override {
        if (id_ == "bernoulli") {
            if (x < 0.0) return 0.0;
            if (x < 1.0) return 1.0 - p(values_, "p");
            return 1.0;
        }
        if (type() == DistributionType::Discrete) {
            if (x < 0.0) return 0.0;
            double sum = 0.0;
            const int upper = static_cast<int>(std::floor(x));
            const int start = id_ == "geometric" ? 1 : 0;
            for (int k = start; k <= upper; ++k) sum += density(k);
            return clampProbability(sum);
        }
        if (id_ == "uniform") {
            const double a = p(values_, "a"), b = p(values_, "b");
            return clampProbability((x - a) / (b - a));
        }
        if (id_ == "normal") return normalCdf((x - p(values_, "mu")) / p(values_, "sigma"));
        if (id_ == "exponential") return x <= 0.0 ? 0.0 : -std::expm1(-p(values_, "lambda") * x);
        if (id_ == "gamma")
            return x <= 0.0 ? 0.0 : regularizedGammaP(p(values_, "shape"), x / p(values_, "scale"));
        if (id_ == "beta")
            return regularizedBeta(x, p(values_, "alpha"), p(values_, "beta"));
        if (id_ == "chi_squared")
            return x <= 0.0 ? 0.0 : regularizedGammaP(p(values_, "df") / 2.0, x / 2.0);
        if (id_ == "student_t") {
            const double nu = p(values_, "df");
            const double ib = regularizedBeta(nu / (nu + x * x), nu / 2.0, 0.5);
            return x >= 0.0 ? 1.0 - ib / 2.0 : ib / 2.0;
        }
        if (id_ == "f") {
            if (x <= 0.0) return 0.0;
            const double d1 = p(values_, "d1"), d2 = p(values_, "d2");
            return regularizedBeta(d1 * x / (d1 * x + d2), d1 / 2.0, d2 / 2.0);
        }
        return 0.0;
    }

    double sample(std::mt19937_64& engine) const override {
        if (id_ == "bernoulli")
            return std::bernoulli_distribution(p(values_, "p"))(engine);
        if (id_ == "binomial")
            return std::binomial_distribution<int>(static_cast<int>(p(values_, "n")), p(values_, "p"))(engine);
        if (id_ == "geometric")
            return std::geometric_distribution<int>(p(values_, "p"))(engine) + 1.0;
        if (id_ == "poisson")
            return std::poisson_distribution<int>(p(values_, "lambda"))(engine);
        if (id_ == "negative_binomial")
            return std::negative_binomial_distribution<int>(
                static_cast<int>(p(values_, "r")), p(values_, "p"))(engine);
        if (id_ == "uniform")
            return std::uniform_real_distribution<double>(p(values_, "a"), p(values_, "b"))(engine);
        if (id_ == "normal")
            return std::normal_distribution<double>(p(values_, "mu"), p(values_, "sigma"))(engine);
        if (id_ == "exponential")
            return std::exponential_distribution<double>(p(values_, "lambda"))(engine);
        if (id_ == "gamma")
            return std::gamma_distribution<double>(p(values_, "shape"), p(values_, "scale"))(engine);
        if (id_ == "beta") {
            const double a = std::gamma_distribution<double>(p(values_, "alpha"), 1.0)(engine);
            const double b = std::gamma_distribution<double>(p(values_, "beta"), 1.0)(engine);
            return a / (a + b);
        }
        if (id_ == "chi_squared")
            return std::gamma_distribution<double>(p(values_, "df") / 2.0, 2.0)(engine);
        if (id_ == "student_t") {
            const double z = std::normal_distribution<double>()(engine);
            const double chi = std::gamma_distribution<double>(p(values_, "df") / 2.0, 2.0)(engine);
            return z / std::sqrt(chi / p(values_, "df"));
        }
        if (id_ == "f") {
            const double d1 = p(values_, "d1"), d2 = p(values_, "d2");
            const double x1 = std::gamma_distribution<double>(d1 / 2.0, 2.0)(engine) / d1;
            const double x2 = std::gamma_distribution<double>(d2 / 2.0, 2.0)(engine) / d2;
            return x1 / x2;
        }
        return 0.0;
    }

    double mean() const override {
        if (id_ == "bernoulli") return p(values_, "p");
        if (id_ == "binomial") return p(values_, "n") * p(values_, "p");
        if (id_ == "geometric") return 1.0 / p(values_, "p");
        if (id_ == "poisson") return p(values_, "lambda");
        if (id_ == "negative_binomial") return p(values_, "r") * (1.0 - p(values_, "p")) / p(values_, "p");
        if (id_ == "uniform") return (p(values_, "a") + p(values_, "b")) / 2.0;
        if (id_ == "normal") return p(values_, "mu");
        if (id_ == "exponential") return 1.0 / p(values_, "lambda");
        if (id_ == "gamma") return p(values_, "shape") * p(values_, "scale");
        if (id_ == "beta") return p(values_, "alpha") / (p(values_, "alpha") + p(values_, "beta"));
        if (id_ == "chi_squared") return p(values_, "df");
        if (id_ == "student_t") return p(values_, "df") > 1.0 ? 0.0 : std::numeric_limits<double>::quiet_NaN();
        if (id_ == "f") return p(values_, "d2") > 2.0 ? p(values_, "d2") / (p(values_, "d2") - 2.0)
                                                       : std::numeric_limits<double>::quiet_NaN();
        return 0.0;
    }

    double variance() const override {
        if (id_ == "bernoulli") return p(values_, "p") * (1.0 - p(values_, "p"));
        if (id_ == "binomial") return p(values_, "n") * p(values_, "p") * (1.0 - p(values_, "p"));
        if (id_ == "geometric") return (1.0 - p(values_, "p")) / std::pow(p(values_, "p"), 2.0);
        if (id_ == "poisson") return p(values_, "lambda");
        if (id_ == "negative_binomial")
            return p(values_, "r") * (1.0 - p(values_, "p")) / std::pow(p(values_, "p"), 2.0);
        if (id_ == "uniform") return std::pow(p(values_, "b") - p(values_, "a"), 2.0) / 12.0;
        if (id_ == "normal") return std::pow(p(values_, "sigma"), 2.0);
        if (id_ == "exponential") return 1.0 / std::pow(p(values_, "lambda"), 2.0);
        if (id_ == "gamma") return p(values_, "shape") * std::pow(p(values_, "scale"), 2.0);
        if (id_ == "beta") {
            const double a = p(values_, "alpha"), b = p(values_, "beta");
            return a * b / (std::pow(a + b, 2.0) * (a + b + 1.0));
        }
        if (id_ == "chi_squared") return 2.0 * p(values_, "df");
        if (id_ == "student_t") {
            const double nu = p(values_, "df");
            return nu > 2.0 ? nu / (nu - 2.0) : std::numeric_limits<double>::infinity();
        }
        if (id_ == "f") {
            const double d1 = p(values_, "d1"), d2 = p(values_, "d2");
            if (d2 <= 4.0) return std::numeric_limits<double>::infinity();
            return 2.0 * d2 * d2 * (d1 + d2 - 2.0)
                   / (d1 * std::pow(d2 - 2.0, 2.0) * (d2 - 4.0));
        }
        return 0.0;
    }

    std::string formula() const override {
        if (id_ == "bernoulli") return "P(X=x) = p^x(1-p)^(1-x),  x in {0,1}";
        if (id_ == "binomial") return "P(X=k) = C(n,k) p^k (1-p)^(n-k)";
        if (id_ == "geometric") return "P(X=k) = p(1-p)^(k-1),  k=1,2,...";
        if (id_ == "poisson") return "P(X=k) = exp(-lambda) lambda^k / k!";
        if (id_ == "negative_binomial") return "P(X=k) = C(k+r-1,k) p^r (1-p)^k";
        if (id_ == "uniform") return "f(x) = 1/(b-a),  a <= x <= b";
        if (id_ == "normal") return "f(x) = exp(-(x-mu)^2/(2 sigma^2)) / (sigma sqrt(2 pi))";
        if (id_ == "exponential") return "f(x) = lambda exp(-lambda x),  x >= 0";
        if (id_ == "gamma") return "f(x) = x^(alpha-1) exp(-x/theta) / (Gamma(alpha) theta^alpha)";
        if (id_ == "beta") return "f(x) = x^(alpha-1)(1-x)^(beta-1) / B(alpha,beta)";
        if (id_ == "chi_squared") return "f(x) = x^(k/2-1) exp(-x/2) / (2^(k/2) Gamma(k/2))";
        if (id_ == "student_t") return "f(x) = Gamma((nu+1)/2) / (sqrt(nu pi) Gamma(nu/2)) (1+x^2/nu)^(-(nu+1)/2)";
        if (id_ == "f") return "f(x) = (d1/d2)^(d1/2) x^(d1/2-1) / (B(d1/2,d2/2)(1+d1 x/d2)^((d1+d2)/2))";
        return {};
    }

private:
    std::string id_;
    std::map<std::string, double> values_;
};

bool validate(const DistributionSpec& spec,
              const std::map<std::string, double>& values,
              std::string& error) {
    for (const auto& parameter : spec.parameters) {
        const auto it = values.find(parameter.key);
        if (it == values.end()) {
            error = "Missing parameter: " + parameter.key;
            return false;
        }
        const double value = it->second;
        if (!std::isfinite(value) || value < parameter.minimum || value > parameter.maximum) {
            error = parameter.label + " is outside its valid range.";
            return false;
        }
        if (parameter.integral && std::abs(value - std::round(value)) > 1e-9) {
            error = parameter.label + " must be an integer.";
            return false;
        }
    }
    if (spec.id == "uniform" && p(values, "a") >= p(values, "b")) {
        error = "Uniform distribution requires a < b.";
        return false;
    }
    return true;
}

} // namespace

const std::vector<DistributionSpec>& distributionCatalog() {
    static const std::vector<DistributionSpec> catalog = {
        {"normal", "Normal · 正态", DistributionType::Continuous,
         {{"mu", "μ", 0.0, -1e6, 1e6}, {"sigma", "σ", 1.0, 1e-9, 1e6}}},
        {"uniform", "Uniform · 均匀", DistributionType::Continuous,
         {{"a", "a", 0.0, -1e6, 1e6}, {"b", "b", 1.0, -1e6, 1e6}}},
        {"exponential", "Exponential · 指数", DistributionType::Continuous,
         {{"lambda", "λ", 1.0, 1e-9, 1e6}}},
        {"gamma", "Gamma · 伽马", DistributionType::Continuous,
         {{"shape", "shape α", 2.0, 1e-9, 1e6}, {"scale", "scale θ", 1.0, 1e-9, 1e6}}},
        {"beta", "Beta · 贝塔", DistributionType::Continuous,
         {{"alpha", "α", 2.0, 1e-9, 1e6}, {"beta", "β", 5.0, 1e-9, 1e6}}},
        {"chi_squared", "Chi-square · 卡方", DistributionType::Continuous,
         {{"df", "degrees ν", 4.0, 1e-9, 1e6}}},
        {"student_t", "Student t · t 分布", DistributionType::Continuous,
         {{"df", "degrees ν", 5.0, 1e-9, 1e6}}},
        {"f", "F distribution · F 分布", DistributionType::Continuous,
         {{"d1", "d₁", 5.0, 1e-9, 1e6}, {"d2", "d₂", 10.0, 1e-9, 1e6}}},
        {"bernoulli", "Bernoulli · 伯努利", DistributionType::Discrete,
         {{"p", "p", 0.5, 1e-9, 1.0 - 1e-9}}},
        {"binomial", "Binomial · 二项", DistributionType::Discrete,
         {{"n", "n", 10.0, 1.0, 100000.0, true}, {"p", "p", 0.5, 1e-9, 1.0 - 1e-9}}},
        {"geometric", "Geometric · 几何", DistributionType::Discrete,
         {{"p", "p", 0.5, 1e-9, 1.0 - 1e-9}}},
        {"poisson", "Poisson · 泊松", DistributionType::Discrete,
         {{"lambda", "λ", 4.0, 1e-9, 1e6}}},
        {"negative_binomial", "Negative Binomial · 负二项", DistributionType::Discrete,
         {{"r", "r", 5.0, 1.0, 100000.0, true}, {"p", "p", 0.5, 1e-9, 1.0 - 1e-9}}},
    };
    return catalog;
}

std::shared_ptr<Distribution> createDistribution(
    const std::string& id,
    const std::map<std::string, double>& parameters,
    std::string* error) {
    const auto& catalog = distributionCatalog();
    const auto it = std::find_if(catalog.begin(), catalog.end(),
                                 [&id](const DistributionSpec& spec) { return spec.id == id; });
    if (it == catalog.end()) {
        if (error) *error = "Unknown distribution: " + id;
        return {};
    }
    std::string validationError;
    if (!validate(*it, parameters, validationError)) {
        if (error) *error = validationError;
        return {};
    }
    return std::make_shared<GenericDistribution>(id, parameters);
}

} // namespace stochia
