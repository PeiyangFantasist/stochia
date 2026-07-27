#pragma once

#include "core/Distribution.h"
#include "core/Expression.h"

#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace stochia {

struct RandomVariable {
    std::string name;
    std::shared_ptr<Distribution> distribution;
    std::shared_ptr<Expression> expression;
    std::shared_ptr<Distribution> analyticalDistribution;
    std::string analyticalName;
    std::string derivation;
    bool closedForm = false;

    bool isTransformation() const { return static_cast<bool>(expression); }
    const Distribution* theory() const {
        return distribution ? distribution.get() : analyticalDistribution.get();
    }
};

struct SummaryStatistics {
    double mean = 0.0;
    double variance = 0.0;
    double minimum = 0.0;
    double maximum = 0.0;
    std::size_t count = 0;
};

class ProbabilityModel {
public:
    bool addDistribution(const std::string& name,
                         std::shared_ptr<Distribution> distribution,
                         std::string* error = nullptr);
    bool addTransformation(const std::string& name,
                           const std::string& expression,
                           std::string* error = nullptr);
    bool remove(const std::string& name);
    void clear();

    const RandomVariable* find(const std::string& name) const;
    const std::map<std::string, RandomVariable>& variables() const { return variables_; }
    std::vector<std::string> orderedNames() const;

    std::vector<double> simulate(const std::string& name,
                                 std::size_t count,
                                 std::uint64_t seed,
                                 std::string* error = nullptr) const;
    static SummaryStatistics summarize(const std::vector<double>& samples);

private:
    struct AnalyticalResult {
        std::shared_ptr<Distribution> distribution;
        std::string name;
        std::string derivation;
        bool closedForm = false;
    };

    AnalyticalResult deriveAnalytical(const Expression::Structure& structure) const;
    void rebuildAnalytical();
    double sampleOne(const std::string& name,
                     std::mt19937_64& engine,
                     std::unordered_map<std::string, double>& cache,
                     std::set<std::string>& active,
                     std::string* error) const;
    bool hasPath(const std::string& from, const std::string& to) const;

    std::map<std::string, RandomVariable> variables_;
    std::vector<std::string> order_;
};

} // namespace stochia
