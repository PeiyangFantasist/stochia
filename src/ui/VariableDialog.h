#pragma once

#include "core/Distribution.h"

#include <QDialog>

#include <map>

class QComboBox;
class QDoubleSpinBox;
class QFormLayout;
class QLineEdit;
class QLabel;

class VariableDialog final : public QDialog {
public:
    explicit VariableDialog(QWidget* parent = nullptr);
    void setDefinition(const QString& name, const std::string& distributionId,
                       const std::map<std::string, double>& parameters,
                       bool lockName = false);

    QString variableName() const;
    std::string distributionId() const;
    std::map<std::string, double> parameters() const;

private:
    void rebuildParameters();
    void validate();

    QLineEdit* nameEdit_ = nullptr;
    QComboBox* distributionCombo_ = nullptr;
    QFormLayout* parameterLayout_ = nullptr;
    QLabel* previewLabel_ = nullptr;
    std::map<std::string, QDoubleSpinBox*> parameterInputs_;
};
