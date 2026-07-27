#include "core/Distribution.h"
#include "core/Expression.h"
#include "core/ProbabilityModel.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unordered_map>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

bool closeTo(double a, double b, double tolerance) {
    return std::abs(a - b) <= tolerance;
}

} // namespace

int main() {
    std::string error;

    require(stochia::distributionCatalog().size() == 13, "all requested built-in distributions registered");
    for (const auto& spec : stochia::distributionCatalog()) {
        std::map<std::string, double> defaults;
        for (const auto& parameter : spec.parameters) defaults[parameter.key] = parameter.defaultValue;
        auto distribution = stochia::createDistribution(spec.id, defaults, &error);
        require(distribution != nullptr, "catalog distribution can be created from defaults");
        const auto range = distribution->plotRange();
        double previousCdf = -1.0;
        for (int i = 0; i <= 200; ++i) {
            const double x = range.first + (range.second - range.first) * i / 200.0;
            const double cdf = distribution->cdf(x);
            require(std::isfinite(cdf) && cdf >= -1e-12 && cdf <= 1.0 + 1e-12,
                    "CDF remains a probability");
            require(cdf + 1e-10 >= previousCdf, "CDF is monotone");
            previousCdf = cdf;
            const double density = distribution->density(x);
            require((std::isfinite(density) && density >= -1e-12) || std::isinf(density),
                    "density or PMF is non-negative");
        }
    }

    auto normal = stochia::createDistribution("normal", {{"mu", 0.0}, {"sigma", 1.0}}, &error);
    require(normal != nullptr, "normal distribution can be created");
    require(closeTo(normal->density(0.0), 0.39894228, 1e-7), "standard normal density");
    require(closeTo(normal->cdf(0.0), 0.5, 1e-12), "standard normal CDF");
    require(closeTo(normal->variance(), 1.0, 1e-12), "standard normal variance");

    auto binomial = stochia::createDistribution("binomial", {{"n", 10.0}, {"p", 0.5}}, &error);
    require(binomial != nullptr, "binomial distribution can be created");
    require(closeTo(binomial->density(5.0), 0.24609375, 1e-10), "binomial PMF");
    require(closeTo(binomial->cdf(10.0), 1.0, 1e-12), "binomial CDF");

    auto beta = stochia::createDistribution("beta", {{"alpha", 2.0}, {"beta", 2.0}}, &error);
    require(beta != nullptr, "beta distribution can be created");
    require(closeTo(beta->cdf(0.5), 0.5, 1e-10), "symmetric beta CDF");

    auto expression = stochia::Expression::parse("log(abs(X)) + max(Y, 2)^2", &error);
    require(expression.isValid(), "expression parses");
    const double evaluated = expression.evaluate({{"X", std::exp(1.0)}, {"Y", 1.0}}, &error);
    require(closeTo(evaluated, 5.0, 1e-12), "expression evaluates");
    require(expression.variables().size() == 2, "expression dependencies collected");

    auto variadic = stochia::Expression::parse(
        "sum(X,Y,3) + product(2,3,4) + min(X,Y,8) + max(X,Y,-1)", &error);
    require(variadic.isValid(), "variadic functions parse");
    require(closeTo(variadic.evaluate({{"X", 2.0}, {"Y", 5.0}}, &error), 41.0, 1e-12),
            "variadic sum, product, min and max evaluate");

    auto iidExpression = stochia::Expression::parse("iid_sum(X,4)", &error);
    int iidCalls = 0;
    const double iidValue = iidExpression.evaluate(
        {}, [&iidCalls](const std::string&) { return static_cast<double>(++iidCalls); }, &error);
    require(closeTo(iidValue, 10.0, 1e-12) && iidCalls == 4,
            "IID sum requests independent samples");

    auto momentsExpression = stochia::Expression::parse("E(X) + Var(Y)", &error);
    require(momentsExpression.isValid(), "E and Var operators parse");
    const double momentsValue = momentsExpression.evaluate(
        {}, {},
        [](const std::string& name, bool variance) {
            if (name == "X" && !variance) return 2.5;
            if (name == "Y" && variance) return 4.0;
            return 0.0;
        }, &error);
    require(closeTo(momentsValue, 6.5, 1e-12), "E and Var resolve theoretical moments");
    require(momentsExpression.ordinaryVariables().empty(),
            "moment operators do not consume random samples");

    stochia::ProbabilityModel model;
    require(model.addDistribution("X", normal, &error), "add X");
    require(model.addDistribution(
        "Y", stochia::createDistribution("exponential", {{"lambda", 2.0}}, &error), &error), "add Y");
    require(model.addTransformation("Z", "X + Y", &error), "add transformation");
    auto samples = model.simulate("Z", 50000, 42, &error);
    require(samples.size() == 50000, "simulate transformed variable");
    const auto summary = stochia::ProbabilityModel::summarize(samples);
    require(closeTo(summary.mean, 0.5, 0.03), "transformed sample mean");
    require(closeTo(summary.variance, 1.25, 0.06), "transformed sample variance");
    require(model.find("Z")->analyticalDistribution != nullptr,
            "continuous sum receives a theoretical distribution");

    stochia::ProbabilityModel iidModel;
    auto exponential = stochia::createDistribution("exponential", {{"lambda", 2.0}}, &error);
    require(iidModel.addDistribution("E", exponential, &error), "add exponential source");
    const bool addedIidSum = iidModel.addTransformation("S", "iid_sum(E,3)", &error);
    if (!addedIidSum) std::cerr << "IID parse error: " << error << '\n';
    require(addedIidSum, "add IID exponential sum");
    const auto* iidSum = iidModel.find("S");
    require(iidSum && iidSum->analyticalDistribution, "IID sum has analytical distribution");
    require(iidSum->closedForm && iidSum->analyticalDistribution->id() == "gamma",
            "IID exponential sum identified as Gamma");
    require(closeTo(iidSum->analyticalDistribution->parameters().at("shape"), 3.0, 1e-12),
            "Gamma shape equals IID count");
    require(closeTo(iidSum->analyticalDistribution->parameters().at("scale"), 0.5, 1e-12),
            "Gamma scale is inverse exponential rate");
    const auto iidSamples = iidModel.simulate("S", 30000, 99, &error);
    require(closeTo(stochia::ProbabilityModel::summarize(iidSamples).mean, 1.5, 0.04),
            "IID sampling uses exactly n independent copies");

    require(iidModel.addTransformation("Maximum", "iid_max(E,5)", &error), "add IID maximum");
    const auto* maximum = iidModel.find("Maximum");
    require(maximum && maximum->analyticalDistribution, "IID maximum theory available");
    const double expectedMaximumCdf = std::pow(exponential->cdf(1.0), 5.0);
    require(closeTo(maximum->analyticalDistribution->cdf(1.0), expectedMaximumCdf, 1e-10),
            "IID maximum CDF formula");

    stochia::ProbabilityModel closedModel;
    require(closedModel.addDistribution(
        "A", stochia::createDistribution("normal", {{"mu", 1.0}, {"sigma", 2.0}}, &error), &error),
        "add first normal");
    require(closedModel.addDistribution(
        "B", stochia::createDistribution("normal", {{"mu", 3.0}, {"sigma", 4.0}}, &error), &error),
        "add second normal");
    require(closedModel.addTransformation("N", "A-B", &error), "normal difference");
    const auto normalDifference = closedModel.find("N")->analyticalDistribution;
    require(normalDifference && normalDifference->id() == "normal", "normal difference stays normal");
    require(closeTo(normalDifference->mean(), -2.0, 1e-12), "normal difference mean");
    require(closeTo(normalDifference->variance(), 20.0, 1e-12), "normal difference variance");
    require(closedModel.addTransformation(
        "StandardA", "(A-E(A))/sqrt(Var(A))", &error), "standardization with E and Var");
    const auto standardized = closedModel.find("StandardA")->analyticalDistribution;
    require(standardized && standardized->id() == "normal", "standardization stays normal");
    require(closeTo(standardized->mean(), 0.0, 1e-12)
            && closeTo(standardized->variance(), 1.0, 1e-12),
            "E and Var participate in analytical simplification");
    require(closedModel.addTransformation("Order", "max(A,B)", &error), "analytical maximum");
    const auto orderTheory = closedModel.find("Order")->analyticalDistribution;
    require(orderTheory != nullptr, "maximum has theoretical PDF and CDF");
    require(closeTo(orderTheory->cdf(0.0),
                    closedModel.find("A")->distribution->cdf(0.0)
                    * closedModel.find("B")->distribution->cdf(0.0), 1e-10),
            "maximum CDF uses product formula");

    require(iidModel.addDistribution(
        "E", stochia::createDistribution("exponential", {{"lambda", 4.0}}, &error), &error),
        "edit existing source distribution");
    require(closeTo(iidModel.find("S")->analyticalDistribution->parameters().at("scale"), 0.25, 1e-12),
            "editing source rebuilds downstream analytical distributions");

    stochia::ProbabilityModel uniformModel;
    require(uniformModel.addDistribution(
        "U", stochia::createDistribution("uniform", {{"a", 0.0}, {"b", 1.0}}, &error), &error),
        "add standard uniform");
    require(uniformModel.addTransformation("T", "iid_sum(U,10)", &error),
            "add uniform IID sum");
    const auto* uniformSum = uniformModel.find("T");
    require(uniformSum && uniformSum->analyticalDistribution,
            "uniform IID sum receives exact theory");
    require(uniformSum->closedForm
            && uniformSum->analyticalDistribution->displayName() == "Scaled Irwin-Hall",
            "uniform IID sum identified as scaled Irwin-Hall");
    require(closeTo(uniformSum->analyticalDistribution->cdf(5.0), 0.5, 1e-10),
            "Irwin-Hall CDF is symmetric at its mean");
    require(closeTo(uniformSum->analyticalDistribution->mean(), 5.0, 1e-12),
            "Irwin-Hall exact mean");
    require(closeTo(uniformSum->analyticalDistribution->variance(), 10.0 / 12.0, 1e-12),
            "Irwin-Hall exact variance");
    require(uniformSum->analyticalDistribution->density(5.0) > 0.0,
            "Irwin-Hall PDF available");

    stochia::ProbabilityModel poissonModel;
    require(poissonModel.addDistribution(
        "P", stochia::createDistribution("poisson", {{"lambda", 3.0}}, &error), &error),
        "add Poisson source");
    require(poissonModel.addTransformation("PS", "iid_sum(P,4)", &error), "add Poisson IID sum");
    require(poissonModel.find("PS")->analyticalDistribution->id() == "poisson"
            && closeTo(poissonModel.find("PS")->analyticalDistribution->mean(), 12.0, 1e-12),
            "Poisson IID sum closure");

    require(!model.addTransformation("Bad", "Unknown + 1", &error), "unknown dependency rejected");
    require(!stochia::createDistribution("uniform", {{"a", 2.0}, {"b", 1.0}}, &error),
            "invalid uniform rejected");

    std::cout << "All Stochia core tests passed.\n";
    return 0;
}
