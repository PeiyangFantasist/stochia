#include "core/Expression.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace stochia {

struct Expression::Node {
    enum class Kind { Number, Variable, Unary, Binary, Function };
    Kind kind = Kind::Number;
    double number = 0.0;
    std::string text;
    std::unique_ptr<Node> left;
    std::unique_ptr<Node> right;
    std::vector<std::unique_ptr<Node>> arguments;

    double evaluate(const std::unordered_map<std::string, double>& variables,
                    const Expression::IndependentSampler& independentSampler,
                    const Expression::MomentResolver& momentResolver) const {
        switch (kind) {
        case Kind::Number:
            return number;
        case Kind::Variable: {
            const auto it = variables.find(text);
            if (it == variables.end()) throw std::runtime_error("Unknown variable: " + text);
            return it->second;
        }
        case Kind::Unary: {
            const double value = left->evaluate(variables, independentSampler, momentResolver);
            return text == "-" ? -value : value;
        }
        case Kind::Binary: {
            const double a = left->evaluate(variables, independentSampler, momentResolver);
            const double b = right->evaluate(variables, independentSampler, momentResolver);
            if (text == "+") return a + b;
            if (text == "-") return a - b;
            if (text == "*") return a * b;
            if (text == "/") {
                if (std::abs(b) < 1e-15) throw std::runtime_error("Division by zero");
                return a / b;
            }
            if (text == "^") return std::pow(a, b);
            throw std::runtime_error("Unknown operator: " + text);
        }
        case Kind::Function: {
            if (text.rfind("iid_", 0) == 0) {
                if (!independentSampler)
                    throw std::runtime_error(text + "() is only available inside a probability model");
                if (arguments.size() != 2 || arguments[0]->kind != Kind::Variable)
                    throw std::runtime_error(text + "() requires a random-variable name as its first argument");
                const double rawCount = arguments[1]->evaluate(
                    variables, independentSampler, momentResolver);
                const auto count = static_cast<long long>(std::llround(rawCount));
                if (count < 1 || count > 1000000 || std::abs(rawCount - count) > 1e-9)
                    throw std::runtime_error(text + "() sample count must be an integer from 1 to 1000000");
                double result = (text == "iid_product") ? 1.0
                                : (text == "iid_min") ? std::numeric_limits<double>::infinity()
                                : (text == "iid_max") ? -std::numeric_limits<double>::infinity()
                                : 0.0;
                for (long long i = 0; i < count; ++i) {
                    const double value = independentSampler(arguments[0]->text);
                    if (text == "iid_sum") result += value;
                    else if (text == "iid_product") result *= value;
                    else if (text == "iid_min") result = std::min(result, value);
                    else if (text == "iid_max") result = std::max(result, value);
                }
                return result;
            }

            if (text == "e" || text == "mean" || text == "var") {
                if (!momentResolver)
                    throw std::runtime_error(text + "() requires a probability model");
                if (arguments.size() != 1 || arguments[0]->kind != Kind::Variable)
                    throw std::runtime_error(text + "() requires one random-variable name");
                return momentResolver(arguments[0]->text, text == "var");
            }

            if (arguments.empty()) throw std::runtime_error("Function has no arguments: " + text);
            const double a = arguments[0]->evaluate(variables, independentSampler, momentResolver);
            if (text == "sin") return std::sin(a);
            if (text == "cos") return std::cos(a);
            if (text == "log") {
                if (a <= 0.0) throw std::runtime_error("log() requires a positive value");
                return std::log(a);
            }
            if (text == "exp") return std::exp(a);
            if (text == "sqrt") {
                if (a < 0.0) throw std::runtime_error("sqrt() requires a non-negative value");
                return std::sqrt(a);
            }
            if (text == "abs") return std::abs(a);
            if (text == "sum") {
                double result = 0.0;
                for (const auto& argument : arguments)
                    result += argument->evaluate(variables, independentSampler, momentResolver);
                return result;
            }
            if (text == "product" || text == "prod") {
                double result = 1.0;
                for (const auto& argument : arguments)
                    result *= argument->evaluate(variables, independentSampler, momentResolver);
                return result;
            }
            if (text == "min") {
                double result = a;
                for (std::size_t i = 1; i < arguments.size(); ++i)
                    result = std::min(result, arguments[i]->evaluate(
                        variables, independentSampler, momentResolver));
                return result;
            }
            if (text == "max") {
                double result = a;
                for (std::size_t i = 1; i < arguments.size(); ++i)
                    result = std::max(result, arguments[i]->evaluate(
                        variables, independentSampler, momentResolver));
                return result;
            }
            throw std::runtime_error("Unknown function: " + text);
        }
        }
        return std::numeric_limits<double>::quiet_NaN();
    }
};

namespace {

enum class TokenType { End, Number, Identifier, Plus, Minus, Star, Slash, Caret, LeftParen, RightParen, Comma };

struct Token {
    TokenType type = TokenType::End;
    std::string text;
    double number = 0.0;
    std::size_t position = 0;
};

class Lexer {
public:
    explicit Lexer(const std::string& input) : input_(input) {}

    Token next() {
        while (position_ < input_.size() && std::isspace(static_cast<unsigned char>(input_[position_])))
            ++position_;
        if (position_ >= input_.size()) return {TokenType::End, {}, 0.0, position_};

        const std::size_t start = position_;
        const char c = input_[position_++];
        switch (c) {
        case '+': return {TokenType::Plus, "+", 0.0, start};
        case '-': return {TokenType::Minus, "-", 0.0, start};
        case '*': return {TokenType::Star, "*", 0.0, start};
        case '/': return {TokenType::Slash, "/", 0.0, start};
        case '^': return {TokenType::Caret, "^", 0.0, start};
        case '(': return {TokenType::LeftParen, "(", 0.0, start};
        case ')': return {TokenType::RightParen, ")", 0.0, start};
        case ',': return {TokenType::Comma, ",", 0.0, start};
        default: break;
        }

        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            while (position_ < input_.size()
                   && (std::isdigit(static_cast<unsigned char>(input_[position_]))
                       || input_[position_] == '.')) ++position_;
            if (position_ < input_.size() && (input_[position_] == 'e' || input_[position_] == 'E')) {
                ++position_;
                if (position_ < input_.size() && (input_[position_] == '+' || input_[position_] == '-'))
                    ++position_;
                while (position_ < input_.size()
                       && std::isdigit(static_cast<unsigned char>(input_[position_]))) ++position_;
            }
            const std::string text = input_.substr(start, position_ - start);
            std::size_t used = 0;
            const double number = std::stod(text, &used);
            if (used != text.size()) throw std::runtime_error("Invalid number at position " + std::to_string(start + 1));
            return {TokenType::Number, text, number, start};
        }

        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            while (position_ < input_.size()
                   && (std::isalnum(static_cast<unsigned char>(input_[position_]))
                       || input_[position_] == '_')) ++position_;
            return {TokenType::Identifier, input_.substr(start, position_ - start), 0.0, start};
        }
        throw std::runtime_error("Unexpected character at position " + std::to_string(start + 1));
    }

private:
    const std::string& input_;
    std::size_t position_ = 0;
};

class Parser {
public:
    Parser(const std::string& source, std::set<std::string>& variables)
        : lexer_(source), variables_(variables) {
        advance();
    }

    std::unique_ptr<Expression::Node> parse() {
        auto node = parseAdditive();
        if (current_.type != TokenType::End)
            fail("Unexpected token '" + current_.text + "'");
        return node;
    }

private:
    void advance() { current_ = lexer_.next(); }

    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error(message + " at position " + std::to_string(current_.position + 1));
    }

    std::unique_ptr<Expression::Node> makeBinary(std::string op,
                                                 std::unique_ptr<Expression::Node> left,
                                                 std::unique_ptr<Expression::Node> right) {
        auto node = std::make_unique<Expression::Node>();
        node->kind = Expression::Node::Kind::Binary;
        node->text = std::move(op);
        node->left = std::move(left);
        node->right = std::move(right);
        return node;
    }

    std::unique_ptr<Expression::Node> parseAdditive() {
        auto left = parseMultiplicative();
        while (current_.type == TokenType::Plus || current_.type == TokenType::Minus) {
            const std::string op = current_.text;
            advance();
            left = makeBinary(op, std::move(left), parseMultiplicative());
        }
        return left;
    }

    std::unique_ptr<Expression::Node> parseMultiplicative() {
        auto left = parseUnary();
        while (current_.type == TokenType::Star || current_.type == TokenType::Slash) {
            const std::string op = current_.text;
            advance();
            left = makeBinary(op, std::move(left), parseUnary());
        }
        return left;
    }

    std::unique_ptr<Expression::Node> parseUnary() {
        if (current_.type == TokenType::Plus || current_.type == TokenType::Minus) {
            const std::string op = current_.text;
            advance();
            auto node = std::make_unique<Expression::Node>();
            node->kind = Expression::Node::Kind::Unary;
            node->text = op;
            node->left = parseUnary();
            return node;
        }
        return parsePower();
    }

    std::unique_ptr<Expression::Node> parsePower() {
        auto left = parsePrimary();
        if (current_.type == TokenType::Caret) {
            advance();
            left = makeBinary("^", std::move(left), parseUnary());
        }
        return left;
    }

    std::unique_ptr<Expression::Node> parsePrimary() {
        if (current_.type == TokenType::Number) {
            auto node = std::make_unique<Expression::Node>();
            node->kind = Expression::Node::Kind::Number;
            node->number = current_.number;
            advance();
            return node;
        }

        if (current_.type == TokenType::Identifier) {
            std::string identifier = current_.text;
            std::transform(identifier.begin(), identifier.end(), identifier.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            const std::string original = current_.text;
            advance();
            if (current_.type != TokenType::LeftParen) {
                auto node = std::make_unique<Expression::Node>();
                node->kind = Expression::Node::Kind::Number;
                if (original == "pi") node->number = 3.14159265358979323846;
                else if (original == "e") node->number = 2.71828182845904523536;
                else {
                    node->kind = Expression::Node::Kind::Variable;
                    node->text = original;
                    variables_.insert(original);
                }
                return node;
            }

            static const std::set<std::string> oneArg = {"sin", "cos", "log", "exp", "sqrt", "abs"};
            static const std::set<std::string> variadic = {"sum", "product", "prod", "min", "max"};
            static const std::set<std::string> iid = {"iid_sum", "iid_product", "iid_min", "iid_max"};
            static const std::set<std::string> moments = {"e", "mean", "var"};
            if (!oneArg.count(identifier) && !variadic.count(identifier)
                && !iid.count(identifier) && !moments.count(identifier))
                fail("Unknown function '" + original + "'");
            advance();
            auto node = std::make_unique<Expression::Node>();
            node->kind = Expression::Node::Kind::Function;
            node->text = identifier;
            if (current_.type == TokenType::RightParen) fail("Expected a function argument");
            node->arguments.push_back(parseAdditive());
            while (current_.type == TokenType::Comma) {
                advance();
                node->arguments.push_back(parseAdditive());
            }
            if (current_.type != TokenType::RightParen) fail("Expected ')'");
            advance();
            if (oneArg.count(identifier) && node->arguments.size() != 1)
                fail(original + "() accepts exactly one argument");
            if (variadic.count(identifier) && node->arguments.size() < 2)
                fail(original + "() requires at least two arguments");
            if (iid.count(identifier)) {
                if (node->arguments.size() != 2)
                    fail(original + "() requires a variable and a count");
                if (node->arguments[0]->kind != Expression::Node::Kind::Variable)
                    fail(original + "() first argument must be a random-variable name");
            }
            if (moments.count(identifier)) {
                if (node->arguments.size() != 1
                    || node->arguments[0]->kind != Expression::Node::Kind::Variable)
                    fail(original + "() requires one random-variable name");
            }
            return node;
        }

        if (current_.type == TokenType::LeftParen) {
            advance();
            auto node = parseAdditive();
            if (current_.type != TokenType::RightParen) fail("Expected ')'");
            advance();
            return node;
        }
        fail("Expected a number, variable, or '('");
    }

    Lexer lexer_;
    Token current_;
    std::set<std::string>& variables_;
};

} // namespace

Expression::Expression() = default;
Expression::Expression(Expression&&) noexcept = default;
Expression& Expression::operator=(Expression&&) noexcept = default;
Expression::~Expression() = default;

Expression Expression::parse(const std::string& source, std::string* error) {
    Expression expression;
    expression.source_ = source;
    try {
        Parser parser(source, expression.variables_);
        expression.root_ = parser.parse();
    } catch (const std::exception& exception) {
        if (error) *error = exception.what();
        expression.root_.reset();
        expression.variables_.clear();
    }
    return expression;
}

bool Expression::isValid() const {
    return static_cast<bool>(root_);
}

double Expression::evaluate(const std::unordered_map<std::string, double>& variables,
                            std::string* error) const {
    return evaluate(variables, {}, {}, error);
}

double Expression::evaluate(const std::unordered_map<std::string, double>& variables,
                            const IndependentSampler& independentSampler,
                            std::string* error) const {
    return evaluate(variables, independentSampler, {}, error);
}

double Expression::evaluate(const std::unordered_map<std::string, double>& variables,
                            const IndependentSampler& independentSampler,
                            const MomentResolver& momentResolver,
                            std::string* error) const {
    if (!root_) {
        if (error) *error = "Expression is not valid.";
        return std::numeric_limits<double>::quiet_NaN();
    }
    try {
        const double value = root_->evaluate(variables, independentSampler, momentResolver);
        if (!std::isfinite(value)) throw std::runtime_error("Expression produced a non-finite value");
        return value;
    } catch (const std::exception& exception) {
        if (error) *error = exception.what();
        return std::numeric_limits<double>::quiet_NaN();
    }
}

const std::set<std::string>& Expression::variables() const {
    return variables_;
}

std::set<std::string> Expression::ordinaryVariables() const {
    std::set<std::string> result;
    std::function<void(const Node*, bool)> visit = [&visit, &result](const Node* node, bool iidSource) {
        if (!node) return;
        if (node->kind == Node::Kind::Variable && !iidSource) result.insert(node->text);
        if (node->kind == Node::Kind::Function
            && (node->text.rfind("iid_", 0) == 0
                || node->text == "e" || node->text == "mean" || node->text == "var")) {
            if (!node->arguments.empty()) visit(node->arguments[0].get(), true);
            for (std::size_t i = 1; i < node->arguments.size(); ++i)
                visit(node->arguments[i].get(), false);
            return;
        }
        visit(node->left.get(), iidSource);
        visit(node->right.get(), iidSource);
        for (const auto& argument : node->arguments) visit(argument.get(), iidSource);
    };
    visit(root_.get(), false);
    return result;
}

const std::string& Expression::source() const {
    return source_;
}

Expression::Structure Expression::structure() const {
    std::function<Structure(const Node*)> convert = [&convert](const Node* node) -> Structure {
        if (!node) return {};
        Structure result;
        result.kind = static_cast<Structure::Kind>(static_cast<int>(node->kind) + 1);
        result.number = node->number;
        result.text = node->text;
        if (node->left) result.children.push_back(convert(node->left.get()));
        if (node->right) result.children.push_back(convert(node->right.get()));
        for (const auto& argument : node->arguments)
            result.children.push_back(convert(argument.get()));
        return result;
    };
    return convert(root_.get());
}

} // namespace stochia
