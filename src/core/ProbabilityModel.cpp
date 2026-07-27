#include "core/ProbabilityModel.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace stochia {
namespace {

bool validVariableName(const std::string& name) {
    if (name.empty() || !(std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_'))
        return false;
    return std::all_of(name.begin() + 1, name.end(), [](unsigned char c) {
        return std::isalnum(c) || c == '_';
    });
}

constexpr int kIntegrationSteps = 280;

double integrate(const std::function<double(double)>& function, double lower, double upper,
                 int steps = kIntegrationSteps) {
    if (!(upper > lower) || !std::isfinite(lower) || !std::isfinite(upper)) return 0.0;
    const double width = (upper - lower) / steps;
    double sum = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double x = lower + (i + 0.5) * width;
        const double value = function(x);
        if (std::isfinite(value)) sum += value;
    }
    return sum * width;
}

class AnalyticalDistribution final : public Distribution {
public:
    using Function = std::function<double(double)>;

    AnalyticalDistribution(std::string name, std::pair<double, double> range,
                           Function density, Function cdf, std::string formula,
                           double exactMean = std::numeric_limits<double>::quiet_NaN(),
                           double exactVariance = std::numeric_limits<double>::quiet_NaN())
        : name_(std::move(name)), range_(range), density_(std::move(density)),
          cdf_(std::move(cdf)), formula_(std::move(formula)),
          exactMean_(exactMean), exactVariance_(exactVariance) {}

    std::string id() const override { return "derived"; }
    std::string displayName() const override { return name_; }
    DistributionType type() const override { return DistributionType::Continuous; }
    std::map<std::string, double> parameters() const override { return {}; }
    std::pair<double, double> plotRange() const override { return range_; }
    double density(double x) const override {
        if (x < range_.first || x > range_.second) return 0.0;
        const double result = density_(x);
        return std::isfinite(result) ? std::max(0.0, result) : 0.0;
    }
    double cdf(double x) const override {
        if (x <= range_.first) return 0.0;
        if (x >= range_.second) return 1.0;
        return std::clamp(cdf_(x), 0.0, 1.0);
    }
    double sample(std::mt19937_64& engine) const override {
        const double target = std::uniform_real_distribution<double>(0.0, 1.0)(engine);
        double low = range_.first, high = range_.second;
        for (int i = 0; i < 70; ++i) {
            const double mid = (low + high) / 2.0;
            if (cdf(mid) < target) low = mid;
            else high = mid;
        }
        return (low + high) / 2.0;
    }
    double mean() const override {
        if (std::isfinite(exactMean_)) return exactMean_;
        return integrate([this](double x) { return x * density(x); }, range_.first, range_.second, 420);
    }
    double variance() const override {
        if (std::isfinite(exactVariance_)) return exactVariance_;
        const double m = mean();
        return integrate([this, m](double x) {
            const double delta = x - m;
            return delta * delta * density(x);
        }, range_.first, range_.second, 420);
    }
    std::string formula() const override { return formula_; }

private:
    std::string name_;
    std::pair<double, double> range_;
    Function density_;
    Function cdf_;
    std::string formula_;
    double exactMean_;
    double exactVariance_;
};

std::pair<double, double> affineRange(const Distribution& source, double scale, double shift) {
    const auto range = source.plotRange();
    const double a = scale * range.first + shift;
    const double b = scale * range.second + shift;
    return {std::min(a, b), std::max(a, b)};
}

std::shared_ptr<Distribution> affineDistribution(const std::shared_ptr<Distribution>& source,
                                                 double scale, double shift) {
    if (!source || std::abs(scale) < 1e-14) return {};
    if (source->id() == "normal") {
        const auto values = source->parameters();
        return createDistribution("normal",
            {{"mu", scale * values.at("mu") + shift},
             {"sigma", std::abs(scale) * values.at("sigma")}});
    }
    if (source->id() == "uniform") {
        const auto values = source->parameters();
        const double a = scale * values.at("a") + shift;
        const double b = scale * values.at("b") + shift;
        return createDistribution("uniform", {{"a", std::min(a, b)}, {"b", std::max(a, b)}});
    }
    const auto range = affineRange(*source, scale, shift);
    return std::make_shared<AnalyticalDistribution>(
        "Affine transform", range,
        [source, scale, shift](double z) {
            return source->density((z - shift) / scale) / std::abs(scale);
        },
        [source, scale, shift](double z) {
            const double value = source->cdf((z - shift) / scale);
            return scale > 0.0 ? value : 1.0 - value;
        },
        "f_Z(z) = f_X((z-b)/a) / |a|,  Z=aX+b");
}

std::pair<double, double> binaryRange(const Distribution& left, const Distribution& right,
                                      const std::string& op) {
    const auto a = left.plotRange();
    const auto b = right.plotRange();
    if (op == "+" || op == "-") {
        return op == "+" ? std::pair<double, double>{a.first + b.first, a.second + b.second}
                         : std::pair<double, double>{a.first - b.second, a.second - b.first};
    }
    std::vector<double> candidates;
    for (double x : {a.first, a.second}) {
        for (double y : {b.first, b.second}) {
            if (op == "*") candidates.push_back(x * y);
            else if (std::abs(y) > 1e-9) candidates.push_back(x / y);
        }
    }
    if (op == "/" && b.first < 0.0 && b.second > 0.0) {
        const double scale = std::max({std::abs(a.first), std::abs(a.second), 1.0})
                             / std::max(std::min(std::abs(b.first), std::abs(b.second)), 0.1);
        return {-10.0 * scale, 10.0 * scale};
    }
    const auto bounds = std::minmax_element(candidates.begin(), candidates.end());
    return {*bounds.first, *bounds.second};
}

std::shared_ptr<Distribution> closedSum(const std::shared_ptr<Distribution>& left,
                                        const std::shared_ptr<Distribution>& right,
                                        bool subtract = false) {
    if (!left || !right) return {};
    if (left->id() == "normal" && right->id() == "normal") {
        const auto a = left->parameters(), b = right->parameters();
        return createDistribution("normal",
            {{"mu", a.at("mu") + (subtract ? -b.at("mu") : b.at("mu"))},
             {"sigma", std::sqrt(a.at("sigma") * a.at("sigma")
                                  + b.at("sigma") * b.at("sigma"))}});
    }
    if (subtract) return {};

    const auto gammaShapeScale = [](const std::shared_ptr<Distribution>& distribution,
                                    double& shape, double& scale) {
        if (distribution->id() == "exponential") {
            shape = 1.0;
            scale = 1.0 / distribution->parameters().at("lambda");
            return true;
        }
        if (distribution->id() == "gamma") {
            shape = distribution->parameters().at("shape");
            scale = distribution->parameters().at("scale");
            return true;
        }
        return false;
    };
    double shapeA = 0.0, scaleA = 0.0, shapeB = 0.0, scaleB = 0.0;
    if (gammaShapeScale(left, shapeA, scaleA) && gammaShapeScale(right, shapeB, scaleB)
        && std::abs(scaleA - scaleB) < 1e-10) {
        return createDistribution("gamma", {{"shape", shapeA + shapeB}, {"scale", scaleA}});
    }
    if (left->id() == "bernoulli" && right->id() == "bernoulli") {
        const double p1 = left->parameters().at("p"), p2 = right->parameters().at("p");
        if (std::abs(p1 - p2) < 1e-10)
            return createDistribution("binomial", {{"n", 2.0}, {"p", p1}});
    }
    return {};
}

std::shared_ptr<Distribution> binaryDistribution(const std::shared_ptr<Distribution>& left,
                                                 const std::shared_ptr<Distribution>& right,
                                                 const std::string& op) {
    if (!left || !right || left->type() != DistributionType::Continuous
        || right->type() != DistributionType::Continuous) return {};
    if (op == "+" || op == "-") {
        if (auto closed = closedSum(left, right, op == "-")) return closed;
    }
    const auto leftRange = left->plotRange();
    const auto rightRange = right->plotRange();
    const auto range = binaryRange(*left, *right, op);

    if (op == "+") {
        return std::make_shared<AnalyticalDistribution>(
            "Numerical convolution", range,
            [left, right, leftRange](double z) {
                return integrate([&](double x) { return left->density(x) * right->density(z - x); },
                                 leftRange.first, leftRange.second);
            },
            [left, right, leftRange](double z) {
                return integrate([&](double x) { return left->density(x) * right->cdf(z - x); },
                                 leftRange.first, leftRange.second);
            },
            "f_Z(z) = integral f_X(x) f_Y(z-x) dx");
    }
    if (op == "-") {
        return std::make_shared<AnalyticalDistribution>(
            "Numerical convolution", range,
            [left, right, leftRange](double z) {
                return integrate([&](double x) { return left->density(x) * right->density(x - z); },
                                 leftRange.first, leftRange.second);
            },
            [left, right, leftRange](double z) {
                return integrate([&](double x) {
                    return left->density(x) * (1.0 - right->cdf(x - z));
                }, leftRange.first, leftRange.second);
            },
            "f_Z(z) = integral f_X(x) f_Y(x-z) dx");
    }
    if (op == "*") {
        return std::make_shared<AnalyticalDistribution>(
            "Product distribution", range,
            [left, right, leftRange](double z) {
                return integrate([&](double x) {
                    if (std::abs(x) < 1e-8) return 0.0;
                    return left->density(x) * right->density(z / x) / std::abs(x);
                }, leftRange.first, leftRange.second);
            },
            [left, right, leftRange](double z) {
                return integrate([&](double x) {
                    if (std::abs(x) < 1e-8) return 0.0;
                    const double conditional = x > 0.0 ? right->cdf(z / x)
                                                       : 1.0 - right->cdf(z / x);
                    return left->density(x) * conditional;
                }, leftRange.first, leftRange.second);
            },
            "f_Z(z) = integral f_X(x) f_Y(z/x) / |x| dx");
    }
    if (op == "/") {
        return std::make_shared<AnalyticalDistribution>(
            "Ratio distribution", range,
            [left, right, rightRange](double z) {
                return integrate([&](double y) {
                    return std::abs(y) * left->density(z * y) * right->density(y);
                }, rightRange.first, rightRange.second);
            },
            [left, right, rightRange](double z) {
                return integrate([&](double y) {
                    if (std::abs(y) < 1e-8) return 0.0;
                    const double conditional = y > 0.0 ? left->cdf(z * y)
                                                       : 1.0 - left->cdf(z * y);
                    return right->density(y) * conditional;
                }, rightRange.first, rightRange.second);
            },
            "f_Z(z) = integral |y| f_X(zy) f_Y(y) dy");
    }
    return {};
}

std::shared_ptr<Distribution> orderDistribution(const std::shared_ptr<Distribution>& left,
                                                const std::shared_ptr<Distribution>& right,
                                                bool maximum) {
    if (!left || !right || left->type() != DistributionType::Continuous
        || right->type() != DistributionType::Continuous) return {};
    const auto a = left->plotRange(), b = right->plotRange();
    const auto range = maximum ? std::pair<double, double>{std::max(a.first, b.first),
                                                           std::max(a.second, b.second)}
                               : std::pair<double, double>{std::min(a.first, b.first),
                                                           std::min(a.second, b.second)};
    if (maximum) {
        return std::make_shared<AnalyticalDistribution>(
            "Maximum distribution", range,
            [left, right](double z) {
                return left->density(z) * right->cdf(z) + right->density(z) * left->cdf(z);
            },
            [left, right](double z) { return left->cdf(z) * right->cdf(z); },
            "F_max(z) = F_X(z) F_Y(z)");
    }
    return std::make_shared<AnalyticalDistribution>(
        "Minimum distribution", range,
        [left, right](double z) {
            return left->density(z) * (1.0 - right->cdf(z))
                   + right->density(z) * (1.0 - left->cdf(z));
        },
        [left, right](double z) {
            return 1.0 - (1.0 - left->cdf(z)) * (1.0 - right->cdf(z));
        },
        "F_min(z) = 1 - (1-F_X(z))(1-F_Y(z))");
}

std::shared_ptr<Distribution> iidOrderDistribution(const std::shared_ptr<Distribution>& source,
                                                   long long count, bool maximum) {
    if (!source || source->type() != DistributionType::Continuous) return {};
    const auto range = source->plotRange();
    if (maximum) {
        return std::make_shared<AnalyticalDistribution>(
            "IID maximum", range,
            [source, count](double z) {
                return count * source->density(z) * std::pow(source->cdf(z), count - 1);
            },
            [source, count](double z) { return std::pow(source->cdf(z), count); },
            "F_max(z) = F_X(z)^n");
    }
    return std::make_shared<AnalyticalDistribution>(
        "IID minimum", range,
        [source, count](double z) {
            return count * source->density(z) * std::pow(1.0 - source->cdf(z), count - 1);
        },
        [source, count](double z) { return 1.0 - std::pow(1.0 - source->cdf(z), count); },
        "F_min(z) = 1 - (1-F_X(z))^n");
}

double irwinHallValue(double x, int count, bool cumulative) {
    if (cumulative) {
        if (x <= 0.0) return 0.0;
        if (x >= count) return 1.0;
    } else if (x <= 0.0 || x >= count) {
        return 0.0;
    }

    std::vector<long double> previous(static_cast<std::size_t>(count));
    for (int shift = 0; shift < count; ++shift) {
        const long double value = static_cast<long double>(x - shift);
        previous[static_cast<std::size_t>(shift)] = cumulative
            ? std::clamp(value, 0.0L, 1.0L)
            : (value >= 0.0L && value < 1.0L ? 1.0L : 0.0L);
    }
    for (int order = 2; order <= count; ++order) {
        const int outputSize = count - order + 1;
        const long double denominator = cumulative ? order : order - 1;
        for (int shift = 0; shift < outputSize; ++shift) {
            const long double value = static_cast<long double>(x - shift);
            previous[static_cast<std::size_t>(shift)] =
                value / denominator * previous[static_cast<std::size_t>(shift)]
                + (order - value) / denominator * previous[static_cast<std::size_t>(shift + 1)];
        }
    }
    return std::clamp(static_cast<double>(previous.front()), 0.0, cumulative ? 1.0
                                                                            : std::numeric_limits<double>::max());
}

std::shared_ptr<Distribution> uniformIidSumDistribution(double a, double b, long long count) {
    if (count < 1 || count > 128 || !(b > a)) return {};
    const int n = static_cast<int>(count);
    const double width = b - a;
    const double shift = count * a;
    const std::pair<double, double> range{count * a, count * b};
    return std::make_shared<AnalyticalDistribution>(
        "Scaled Irwin-Hall", range,
        [n, width, shift](double z) {
            return irwinHallValue((z - shift) / width, n, false) / width;
        },
        [n, width, shift](double z) {
            return irwinHallValue((z - shift) / width, n, true);
        },
        "F_S(s)=1/n! * sum_{k=0}^{floor(t)} (-1)^k C(n,k)(t-k)^n, "
        "t=(s-na)/(b-a)",
        count * (a + b) / 2.0,
        count * width * width / 12.0);
}

bool isNumber(const Expression::Structure& structure, double& value) {
    if (structure.kind != Expression::Structure::Kind::Number) return false;
    value = structure.number;
    return true;
}

} // namespace

bool ProbabilityModel::addDistribution(const std::string& name,
                                       std::shared_ptr<Distribution> distribution,
                                       std::string* error) {
    if (!validVariableName(name)) {
        if (error) *error = "Variable names must begin with a letter and contain only letters, numbers, or _.";
        return false;
    }
    if (!distribution) {
        if (error) *error = "Distribution is missing.";
        return false;
    }
    const bool isNew = variables_.find(name) == variables_.end();
    variables_[name] = RandomVariable{name, std::move(distribution), {}};
    if (isNew) order_.push_back(name);
    rebuildAnalytical();
    return true;
}

bool ProbabilityModel::addTransformation(const std::string& name,
                                         const std::string& source,
                                         std::string* error) {
    if (!validVariableName(name)) {
        if (error) *error = "Variable names must begin with a letter and contain only letters, numbers, or _.";
        return false;
    }
    std::string parseError;
    auto parsed = Expression::parse(source, &parseError);
    if (!parsed.isValid()) {
        if (error) *error = parseError;
        return false;
    }
    if (parsed.variables().empty()) {
        if (error) *error = "A transformation must reference at least one random variable.";
        return false;
    }
    for (const auto& dependency : parsed.variables()) {
        if (dependency == name) {
            if (error) *error = "A variable cannot depend on itself.";
            return false;
        }
        if (!find(dependency)) {
            if (error) *error = "Unknown random variable: " + dependency;
            return false;
        }
        if (hasPath(dependency, name)) {
            if (error) *error = "This definition would create a dependency cycle.";
            return false;
        }
    }
    const bool isNew = variables_.find(name) == variables_.end();
    auto expression = std::make_shared<Expression>(std::move(parsed));
    variables_[name] = RandomVariable{name, {}, std::move(expression)};
    if (isNew) order_.push_back(name);
    rebuildAnalytical();
    return true;
}

bool ProbabilityModel::remove(const std::string& name) {
    for (const auto& [otherName, variable] : variables_) {
        if (otherName != name && variable.expression && variable.expression->variables().count(name))
            return false;
    }
    if (!variables_.erase(name)) return false;
    order_.erase(std::remove(order_.begin(), order_.end(), name), order_.end());
    rebuildAnalytical();
    return true;
}

void ProbabilityModel::clear() {
    variables_.clear();
    order_.clear();
}

const RandomVariable* ProbabilityModel::find(const std::string& name) const {
    const auto it = variables_.find(name);
    return it == variables_.end() ? nullptr : &it->second;
}

std::vector<std::string> ProbabilityModel::orderedNames() const {
    return order_;
}

double ProbabilityModel::sampleOne(const std::string& name,
                                   std::mt19937_64& engine,
                                   std::unordered_map<std::string, double>& cache,
                                   std::set<std::string>& active,
                                   std::string* error) const {
    const auto cached = cache.find(name);
    if (cached != cache.end()) return cached->second;
    const auto* variable = find(name);
    if (!variable) {
        if (error) *error = "Unknown variable: " + name;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!active.insert(name).second) {
        if (error) *error = "Cyclic variable dependency.";
        return std::numeric_limits<double>::quiet_NaN();
    }

    double value = 0.0;
    if (variable->distribution) {
        value = variable->distribution->sample(engine);
    } else {
        std::unordered_map<std::string, double> inputs;
        for (const auto& dependency : variable->expression->ordinaryVariables()) {
            const double dependencyValue = sampleOne(dependency, engine, cache, active, error);
            if (!std::isfinite(dependencyValue)) {
                active.erase(name);
                return dependencyValue;
            }
            inputs.emplace(dependency, dependencyValue);
        }
        const auto independentSampler = [this, &engine, error](const std::string& dependency) {
            std::unordered_map<std::string, double> independentCache;
            std::set<std::string> independentActive;
            return sampleOne(dependency, engine, independentCache, independentActive, error);
        };
        const auto momentResolver = [this](const std::string& dependency, bool variance) {
            const auto* source = find(dependency);
            if (!source || !source->theory())
                throw std::runtime_error("E()/Var() requires a theoretical distribution for " + dependency);
            const double result = variance ? source->theory()->variance() : source->theory()->mean();
            if (!std::isfinite(result))
                throw std::runtime_error((variance ? "Variance" : "Expectation")
                                         + std::string(" is not finite for ") + dependency);
            return result;
        };
        value = variable->expression->evaluate(inputs, independentSampler, momentResolver, error);
    }
    active.erase(name);
    if (std::isfinite(value)) cache[name] = value;
    return value;
}

std::vector<double> ProbabilityModel::simulate(const std::string& name,
                                               std::size_t count,
                                               std::uint64_t seed,
                                               std::string* error) const {
    std::vector<double> values;
    values.reserve(count);
    std::mt19937_64 engine(seed);
    for (std::size_t i = 0; i < count; ++i) {
        std::unordered_map<std::string, double> cache;
        std::set<std::string> active;
        std::string sampleError;
        const double value = sampleOne(name, engine, cache, active, &sampleError);
        if (!std::isfinite(value)) {
            if (error) *error = sampleError;
            return {};
        }
        values.push_back(value);
    }
    return values;
}

SummaryStatistics ProbabilityModel::summarize(const std::vector<double>& samples) {
    SummaryStatistics result;
    result.count = samples.size();
    if (samples.empty()) return result;
    const double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
    result.mean = sum / static_cast<double>(samples.size());
    double squareSum = 0.0;
    for (double value : samples) squareSum += (value - result.mean) * (value - result.mean);
    result.variance = samples.size() > 1 ? squareSum / static_cast<double>(samples.size() - 1) : 0.0;
    const auto extremes = std::minmax_element(samples.begin(), samples.end());
    result.minimum = *extremes.first;
    result.maximum = *extremes.second;
    return result;
}

ProbabilityModel::AnalyticalResult ProbabilityModel::deriveAnalytical(
    const Expression::Structure& structure) const {
    struct Term {
        std::shared_ptr<Distribution> distribution;
        bool constant = false;
        double value = 0.0;
        std::set<std::string> roots;
        std::string derivation;
        bool closed = false;
    };

    std::function<void(const std::string&, std::set<std::string>&)> collectRoots =
        [this, &collectRoots](const std::string& name, std::set<std::string>& roots) {
        const auto* variable = find(name);
        if (!variable) return;
        if (!variable->expression) {
            roots.insert(name);
            return;
        }
        for (const auto& dependency : variable->expression->variables())
            collectRoots(dependency, roots);
    };

    std::function<Term(const Expression::Structure&)> derive =
        [this, &derive, &collectRoots](const Expression::Structure& node) -> Term {
        using Kind = Expression::Structure::Kind;
        if (node.kind == Kind::Number) return {{}, true, node.number, {}, {}, true};
        if (node.kind == Kind::Variable) {
            const auto* variable = find(node.text);
            if (!variable || !variable->theory()) return {};
            std::shared_ptr<Distribution> theory = variable->distribution
                                                   ? variable->distribution
                                                   : variable->analyticalDistribution;
            std::set<std::string> roots;
            collectRoots(node.text, roots);
            return {theory, false, 0.0, roots, variable->derivation,
                    !variable->isTransformation() || variable->closedForm};
        }
        if (node.kind == Kind::Unary && node.children.size() == 1) {
            auto source = derive(node.children[0]);
            if (source.constant) {
                source.value = node.text == "-" ? -source.value : source.value;
                return source;
            }
            if (source.distribution && node.text == "-") {
                source.distribution = affineDistribution(source.distribution, -1.0, 0.0);
                source.derivation = "Unary transform: Z=-X, f_Z(z)=f_X(-z)";
                return source;
            }
            return {};
        }
        if (node.kind == Kind::Binary && node.children.size() == 2) {
            auto left = derive(node.children[0]);
            auto right = derive(node.children[1]);
            if (left.constant && right.constant) {
                if (node.text == "+") left.value += right.value;
                else if (node.text == "-") left.value -= right.value;
                else if (node.text == "*") left.value *= right.value;
                else if (node.text == "/" && std::abs(right.value) > 1e-14) left.value /= right.value;
                else if (node.text == "^") left.value = std::pow(left.value, right.value);
                else return {};
                return left;
            }
            if (left.distribution && right.constant) {
                double scale = 1.0, shift = 0.0;
                if (node.text == "+") shift = right.value;
                else if (node.text == "-") shift = -right.value;
                else if (node.text == "*") scale = right.value;
                else if (node.text == "/" && std::abs(right.value) > 1e-14) scale = 1.0 / right.value;
                else return {};
                left.distribution = affineDistribution(left.distribution, scale, shift);
                left.derivation = "Affine transformation: Z=aX+b, f_Z(z)=f_X((z-b)/a)/|a|";
                left.closed = left.distribution && left.distribution->id() != "derived";
                return left;
            }
            if (left.constant && right.distribution) {
                double scale = 1.0, shift = 0.0;
                if (node.text == "+") shift = left.value;
                else if (node.text == "-") { scale = -1.0; shift = left.value; }
                else if (node.text == "*") scale = left.value;
                else return {};
                right.distribution = affineDistribution(right.distribution, scale, shift);
                right.derivation = "Affine transformation: Z=aX+b, f_Z(z)=f_X((z-b)/a)/|a|";
                right.closed = right.distribution && right.distribution->id() != "derived";
                return right;
            }
            if (!left.distribution || !right.distribution) return {};

            bool independent = true;
            for (const auto& root : left.roots)
                if (right.roots.count(root)) independent = false;
            if (!independent) return {};

            std::shared_ptr<Distribution> result;
            if (node.text == "+" || node.text == "-" || node.text == "*" || node.text == "/")
                result = binaryDistribution(left.distribution, right.distribution, node.text);
            if (!result) return {};
            left.roots.insert(right.roots.begin(), right.roots.end());
            left.distribution = result;
            left.closed = result->id() != "derived";
            if (node.text == "+")
                left.derivation = left.closed
                    ? "Closed-form convolution identified from the distribution family."
                    : "Convolution: f_Z(z)=integral f_X(x)f_Y(z-x)dx";
            else if (node.text == "-")
                left.derivation = left.closed
                    ? "Normal-family closure under independent subtraction."
                    : "Difference: f_Z(z)=integral f_X(x)f_Y(x-z)dx";
            else if (node.text == "*")
                left.derivation = "Product: f_Z(z)=integral f_X(x)f_Y(z/x)/|x| dx";
            else
                left.derivation = "Ratio: f_Z(z)=integral |y|f_X(zy)f_Y(y)dy";
            return left;
        }
        if (node.kind != Kind::Function || node.children.empty()) return {};

        if ((node.text == "e" || node.text == "mean" || node.text == "var")
            && node.children.size() == 1 && node.children[0].kind == Kind::Variable) {
            const auto* source = find(node.children[0].text);
            if (!source || !source->theory()) return {};
            const double value = node.text == "var" ? source->theory()->variance()
                                                    : source->theory()->mean();
            if (!std::isfinite(value)) return {};
            return {{}, true, value, {}, "Statistical moment resolved as an exact constant.", true};
        }

        if ((node.text == "sin" || node.text == "cos" || node.text == "log"
             || node.text == "exp" || node.text == "sqrt" || node.text == "abs")
            && node.children.size() == 1) {
            auto constant = derive(node.children[0]);
            if (constant.constant) {
                if (node.text == "sin") constant.value = std::sin(constant.value);
                else if (node.text == "cos") constant.value = std::cos(constant.value);
                else if (node.text == "log" && constant.value > 0.0)
                    constant.value = std::log(constant.value);
                else if (node.text == "exp") constant.value = std::exp(constant.value);
                else if (node.text == "sqrt" && constant.value >= 0.0)
                    constant.value = std::sqrt(constant.value);
                else if (node.text == "abs") constant.value = std::abs(constant.value);
                else return {};
                return constant;
            }
        }

        if (node.text.rfind("iid_", 0) == 0 && node.children.size() == 2
            && node.children[0].kind == Kind::Variable) {
            double rawCount = 0.0;
            if (!isNumber(node.children[1], rawCount)) return {};
            const auto count = static_cast<long long>(std::llround(rawCount));
            if (count < 1 || count > 1000000 || std::abs(rawCount - count) > 1e-9) return {};
            const auto* sourceVariable = find(node.children[0].text);
            if (!sourceVariable || !sourceVariable->theory()) return {};
            auto source = sourceVariable->distribution
                          ? sourceVariable->distribution : sourceVariable->analyticalDistribution;
            std::shared_ptr<Distribution> result;
            bool closed = false;
            std::string derivation;
            if (node.text == "iid_sum") {
                const auto parameters = source->parameters();
                if (source->id() == "normal") {
                    result = createDistribution("normal",
                        {{"mu", count * parameters.at("mu")},
                         {"sigma", std::sqrt(static_cast<double>(count)) * parameters.at("sigma")}});
                    closed = true;
                    derivation = "IID normal sum: sum X_i ~ Normal(n*mu, sqrt(n)*sigma)";
                } else if (source->id() == "exponential") {
                    result = createDistribution("gamma",
                        {{"shape", static_cast<double>(count)},
                         {"scale", 1.0 / parameters.at("lambda")}});
                    closed = true;
                    derivation = "IID exponential sum: sum X_i ~ Gamma(n, scale=1/lambda)";
                } else if (source->id() == "gamma") {
                    result = createDistribution("gamma",
                        {{"shape", count * parameters.at("shape")},
                         {"scale", parameters.at("scale")}});
                    closed = true;
                    derivation = "IID Gamma sum with common scale: Gamma(n*shape, scale)";
                } else if (source->id() == "bernoulli") {
                    result = createDistribution("binomial",
                        {{"n", static_cast<double>(count)}, {"p", parameters.at("p")}});
                    closed = true;
                    derivation = "IID Bernoulli sum: sum X_i ~ Binomial(n,p)";
                } else if (source->id() == "uniform") {
                    result = uniformIidSumDistribution(
                        parameters.at("a"), parameters.at("b"), count);
                    closed = static_cast<bool>(result);
                    derivation = result
                        ? "IID Uniform(a,b) sum: scaled Irwin-Hall distribution with exact finite-sum PDF/CDF."
                        : "Exact Irwin-Hall evaluation is limited to n<=128 for interactive stability.";
                } else if (source->id() == "poisson"
                           && count * parameters.at("lambda") <= 1e6) {
                    result = createDistribution("poisson",
                        {{"lambda", count * parameters.at("lambda")}});
                    closed = true;
                    derivation = "IID Poisson sum: sum X_i ~ Poisson(n*lambda)";
                } else if (source->id() == "binomial"
                           && count * parameters.at("n") <= 100000.0) {
                    result = createDistribution("binomial",
                        {{"n", count * parameters.at("n")}, {"p", parameters.at("p")}});
                    closed = true;
                    derivation = "IID Binomial sum with common p: Binomial(n*m,p)";
                } else if (source->id() == "negative_binomial"
                           && count * parameters.at("r") <= 100000.0) {
                    result = createDistribution("negative_binomial",
                        {{"r", count * parameters.at("r")}, {"p", parameters.at("p")}});
                    closed = true;
                    derivation = "IID negative-binomial sum with common p: NB(n*r,p)";
                } else if (source->id() == "chi_squared"
                           && count * parameters.at("df") <= 1e6) {
                    result = createDistribution("chi_squared",
                        {{"df", count * parameters.at("df")}});
                    closed = true;
                    derivation = "IID chi-square sum: Chi-square(n*df)";
                } else if (count <= 2 && source->type() == DistributionType::Continuous) {
                    result = source;
                    for (long long i = 1; i < count && result; ++i)
                        result = binaryDistribution(result, source, "+");
                    derivation = "IID sum computed by repeated numerical convolution.";
                }
            } else if (node.text == "iid_product") {
                if (count <= 2 && source->type() == DistributionType::Continuous) {
                    result = source;
                    for (long long i = 1; i < count && result; ++i)
                        result = binaryDistribution(result, source, "*");
                    derivation = "IID product computed from the product-density integral.";
                }
            } else if (node.text == "iid_min" || node.text == "iid_max") {
                result = iidOrderDistribution(source, count, node.text == "iid_max");
                derivation = node.text == "iid_max"
                    ? "IID maximum: F_M(z)=F_X(z)^n"
                    : "IID minimum: F_M(z)=1-(1-F_X(z))^n";
                closed = static_cast<bool>(result);
            }
            std::set<std::string> roots;
            collectRoots(node.children[0].text, roots);
            return {result, false, 0.0, roots, derivation, closed};
        }

        if ((node.text == "sum" || node.text == "product" || node.text == "prod"
             || node.text == "min" || node.text == "max") && node.children.size() >= 2) {
            auto result = derive(node.children[0]);
            if (!result.distribution) return {};
            for (std::size_t i = 1; i < node.children.size(); ++i) {
                auto next = derive(node.children[i]);
                if (!next.distribution) return {};
                if (i > 1 && result.distribution->id() == "derived"
                    && (node.text == "sum" || node.text == "product" || node.text == "prod"))
                    return {};
                for (const auto& root : result.roots)
                    if (next.roots.count(root)) return {};
                std::shared_ptr<Distribution> combined;
                if (node.text == "sum")
                    combined = binaryDistribution(result.distribution, next.distribution, "+");
                else if (node.text == "product" || node.text == "prod")
                    combined = binaryDistribution(result.distribution, next.distribution, "*");
                else
                    combined = orderDistribution(result.distribution, next.distribution, node.text == "max");
                if (!combined) return {};
                result.distribution = combined;
                result.closed = result.closed && next.closed && combined->id() != "derived";
                result.roots.insert(next.roots.begin(), next.roots.end());
            }
            if (node.text == "sum")
                result.derivation = result.distribution->id() != "derived"
                    ? "Known distribution family is closed under this independent sum."
                    : "Multi-variable sum evaluated by repeated convolution.";
            else if (node.text == "product" || node.text == "prod")
                result.derivation = "Multi-variable product evaluated by repeated product-density integrals.";
            else
                result.derivation = node.text == "max"
                    ? "Independent maximum: F_max(z)=product F_i(z)"
                    : "Independent minimum: F_min(z)=1-product(1-F_i(z))";
            return result;
        }

        auto source = derive(node.children[0]);
        if (!source.distribution || source.distribution->type() != DistributionType::Continuous) return {};
        if (node.text == "abs") {
            const auto old = source.distribution;
            const auto oldRange = old->plotRange();
            const std::pair<double, double> range{
                0.0, std::max(std::abs(oldRange.first), std::abs(oldRange.second))};
            source.distribution = std::make_shared<AnalyticalDistribution>(
                "Absolute-value distribution", range,
                [old](double z) { return z < 0.0 ? 0.0 : old->density(z) + old->density(-z); },
                [old](double z) { return z < 0.0 ? 0.0 : old->cdf(z) - old->cdf(-z); },
                "f_|X|(z)=f_X(z)+f_X(-z), z>=0");
            source.derivation = "Absolute-value transformation with two inverse branches.";
            source.closed = false;
            return source;
        }
        return {};
    };

    auto term = derive(structure);
    if (!term.distribution) return {};
    return {term.distribution, term.distribution->displayName(), term.derivation, term.closed};
}

void ProbabilityModel::rebuildAnalytical() {
    for (const auto& name : order_) {
        auto it = variables_.find(name);
        if (it == variables_.end() || !it->second.expression) continue;
        const auto result = deriveAnalytical(it->second.expression->structure());
        it->second.analyticalDistribution = result.distribution;
        it->second.analyticalName = result.name;
        it->second.derivation = result.derivation;
        it->second.closedForm = result.closedForm;
    }
}

bool ProbabilityModel::hasPath(const std::string& from, const std::string& to) const {
    if (from == to) return true;
    const auto* variable = find(from);
    if (!variable || !variable->expression) return false;
    for (const auto& dependency : variable->expression->variables())
        if (hasPath(dependency, to)) return true;
    return false;
}

} // namespace stochia
