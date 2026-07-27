#pragma once

#include <memory>
#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace stochia {

class Expression {
public:
    struct Node;
    struct Structure {
        enum class Kind { Invalid, Number, Variable, Unary, Binary, Function };
        Kind kind = Kind::Invalid;
        double number = 0.0;
        std::string text;
        std::vector<Structure> children;
    };
    using IndependentSampler = std::function<double(const std::string&)>;
    using MomentResolver = std::function<double(const std::string&, bool variance)>;

    Expression();
    Expression(Expression&&) noexcept;
    Expression& operator=(Expression&&) noexcept;
    ~Expression();

    Expression(const Expression&) = delete;
    Expression& operator=(const Expression&) = delete;

    static Expression parse(const std::string& source, std::string* error = nullptr);

    bool isValid() const;
    double evaluate(const std::unordered_map<std::string, double>& variables,
                    std::string* error = nullptr) const;
    double evaluate(const std::unordered_map<std::string, double>& variables,
                    const IndependentSampler& independentSampler,
                    std::string* error = nullptr) const;
    double evaluate(const std::unordered_map<std::string, double>& variables,
                    const IndependentSampler& independentSampler,
                    const MomentResolver& momentResolver,
                    std::string* error = nullptr) const;
    const std::set<std::string>& variables() const;
    std::set<std::string> ordinaryVariables() const;
    const std::string& source() const;
    Structure structure() const;

private:
    std::unique_ptr<Node> root_;
    std::set<std::string> variables_;
    std::string source_;
};

} // namespace stochia
