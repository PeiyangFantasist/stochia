#pragma once

#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace stochia {

enum class DistributionType { Continuous, Discrete };

struct ParameterSpec {
    std::string key;
    std::string label;
    double defaultValue;
    double minimum;
    double maximum;
    bool integral = false;
};

struct DistributionSpec {
    std::string id;
    std::string displayName;
    DistributionType type;
    std::vector<ParameterSpec> parameters;
};

class Distribution {
public:
    virtual ~Distribution() = default;

    virtual std::string id() const = 0;
    virtual std::string displayName() const = 0;
    virtual DistributionType type() const = 0;
    virtual std::map<std::string, double> parameters() const = 0;
    virtual std::pair<double, double> plotRange() const = 0;
    virtual double density(double x) const = 0;
    virtual double cdf(double x) const = 0;
    virtual double sample(std::mt19937_64& engine) const = 0;
    virtual double mean() const = 0;
    virtual double variance() const = 0;
    virtual std::string formula() const = 0;
};

const std::vector<DistributionSpec>& distributionCatalog();
std::shared_ptr<Distribution> createDistribution(
    const std::string& id,
    const std::map<std::string, double>& parameters,
    std::string* error = nullptr);

} // namespace stochia
