#include "ui/VariableDialog.h"

#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QStyle>
#include <QVBoxLayout>

VariableDialog::VariableDialog(QWidget* parent) : QDialog(parent) {
    setWindowTitle(QStringLiteral("添加随机变量"));
    setMinimumWidth(420);

    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(24, 22, 24, 22);
    root->setSpacing(16);

    auto* heading = new QLabel(QStringLiteral("定义一个分布"));
    heading->setObjectName("dialogHeading");
    root->addWidget(heading);

    auto* form = new QFormLayout;
    form->setSpacing(12);
    nameEdit_ = new QLineEdit(QStringLiteral("X"));
    nameEdit_->setPlaceholderText(QStringLiteral("例如 X"));
    form->addRow(QStringLiteral("变量名"), nameEdit_);

    distributionCombo_ = new QComboBox;
    for (const auto& spec : stochia::distributionCatalog())
        distributionCombo_->addItem(QString::fromStdString(spec.displayName), QString::fromStdString(spec.id));
    form->addRow(QStringLiteral("分布"), distributionCombo_);
    root->addLayout(form);

    auto* parameterBox = new QWidget;
    parameterLayout_ = new QFormLayout(parameterBox);
    parameterLayout_->setContentsMargins(0, 0, 0, 0);
    parameterLayout_->setSpacing(12);
    root->addWidget(parameterBox);

    previewLabel_ = new QLabel;
    previewLabel_->setObjectName("formulaCard");
    previewLabel_->setWordWrap(true);
    root->addWidget(previewLabel_);

    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);
    buttons->button(QDialogButtonBox::Ok)->setText(QStringLiteral("添加变量"));
    buttons->button(QDialogButtonBox::Cancel)->setText(QStringLiteral("取消"));
    root->addWidget(buttons);

    connect(distributionCombo_, &QComboBox::currentIndexChanged, this,
            [this] { rebuildParameters(); });
    connect(nameEdit_, &QLineEdit::textChanged, this, [this] { validate(); });
    connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    rebuildParameters();
}

QString VariableDialog::variableName() const {
    return nameEdit_->text().trimmed();
}

std::string VariableDialog::distributionId() const {
    return distributionCombo_->currentData().toString().toStdString();
}

std::map<std::string, double> VariableDialog::parameters() const {
    std::map<std::string, double> result;
    for (const auto& [key, input] : parameterInputs_) result[key] = input->value();
    return result;
}

void VariableDialog::setDefinition(const QString& name, const std::string& distributionId,
                                   const std::map<std::string, double>& parameters,
                                   bool lockName) {
    nameEdit_->setText(name);
    nameEdit_->setReadOnly(lockName);
    const int index = distributionCombo_->findData(QString::fromStdString(distributionId));
    if (index >= 0) distributionCombo_->setCurrentIndex(index);
    rebuildParameters();
    for (const auto& [key, value] : parameters) {
        const auto it = parameterInputs_.find(key);
        if (it != parameterInputs_.end()) it->second->setValue(value);
    }
    validate();
}

void VariableDialog::rebuildParameters() {
    while (auto* item = parameterLayout_->takeAt(0)) {
        if (item->widget()) item->widget()->deleteLater();
        delete item;
    }
    parameterInputs_.clear();
    const auto id = distributionId();
    for (const auto& spec : stochia::distributionCatalog()) {
        if (spec.id != id) continue;
        for (const auto& parameter : spec.parameters) {
            auto* input = new QDoubleSpinBox;
            input->setDecimals(parameter.integral ? 0 : 5);
            input->setRange(parameter.minimum, parameter.maximum);
            input->setValue(parameter.defaultValue);
            input->setSingleStep(parameter.integral ? 1.0 : 0.1);
            input->setKeyboardTracking(false);
            parameterLayout_->addRow(QString::fromStdString(parameter.label), input);
            parameterInputs_[parameter.key] = input;
            connect(input, &QDoubleSpinBox::valueChanged, this, [this] { validate(); });
        }
        break;
    }
    validate();
}

void VariableDialog::validate() {
    std::string error;
    const auto distribution = stochia::createDistribution(distributionId(), parameters(), &error);
    const bool validName = !variableName().isEmpty();
    if (distribution) {
        previewLabel_->setText(QStringLiteral("f / p · %1").arg(QString::fromStdString(distribution->formula())));
        previewLabel_->setProperty("error", false);
    } else {
        previewLabel_->setText(QString::fromStdString(error));
        previewLabel_->setProperty("error", true);
    }
    previewLabel_->style()->unpolish(previewLabel_);
    previewLabel_->style()->polish(previewLabel_);
    if (auto* buttonBox = findChild<QDialogButtonBox*>())
        buttonBox->button(QDialogButtonBox::Ok)->setEnabled(distribution && validName);
}
