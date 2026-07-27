#include "ui/MainWindow.h"

#include "ui/PlotWidget.h"
#include "ui/VariableDialog.h"

#include <QAction>
#include <QApplication>
#include <QCloseEvent>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFrame>
#include <QHBoxLayout>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMenuBar>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTextBrowser>
#include <QToolBar>
#include <QVBoxLayout>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>

namespace {

QString formatValue(double value) {
    if (std::isnan(value)) return QStringLiteral("未定义");
    if (std::isinf(value)) return QStringLiteral("∞");
    return QString::number(value, 'g', 7);
}

QLabel* sectionLabel(const QString& text) {
    auto* label = new QLabel(text);
    label->setObjectName("sectionLabel");
    return label;
}

QWidget* metricCard(const QString& title, QLabel*& valueLabel) {
    auto* card = new QFrame;
    card->setObjectName("metricCard");
    auto* layout = new QVBoxLayout(card);
    layout->setContentsMargins(14, 12, 14, 12);
    layout->setSpacing(3);
    auto* titleLabel = new QLabel(title);
    titleLabel->setObjectName("metricTitle");
    valueLabel = new QLabel(QStringLiteral("—"));
    valueLabel->setObjectName("metricValue");
    layout->addWidget(titleLabel);
    layout->addWidget(valueLabel);
    return card;
}

std::vector<QPointF> empiricalCdf(std::vector<double> samples) {
    std::sort(samples.begin(), samples.end());
    std::vector<QPointF> result;
    if (samples.empty()) return result;
    const std::size_t stride = std::max<std::size_t>(1, samples.size() / 800);
    result.emplace_back(samples.front(), 0.0);
    for (std::size_t i = 0; i < samples.size(); i += stride)
        result.emplace_back(samples[i], static_cast<double>(i + 1) / samples.size());
    result.emplace_back(samples.back(), 1.0);
    return result;
}

std::vector<QPointF> histogram(const std::vector<double>& samples, bool discrete) {
    if (samples.empty()) return {};
    const auto extremes = std::minmax_element(samples.begin(), samples.end());
    if (discrete && *extremes.second - *extremes.first <= 160.0) {
        std::map<long long, std::size_t> counts;
        for (double value : samples) ++counts[std::llround(value)];
        std::vector<QPointF> result;
        result.reserve(counts.size());
        for (const auto& [value, count] : counts)
            result.emplace_back(static_cast<double>(value),
                                static_cast<double>(count) / samples.size());
        return result;
    }

    const double minimum = *extremes.first;
    const double maximum = *extremes.second;
    if (maximum <= minimum) return {{minimum, 1.0}};
    const int bins = std::clamp(static_cast<int>(std::sqrt(samples.size())), 20, 70);
    const double width = (maximum - minimum) / bins;
    std::vector<std::size_t> counts(static_cast<std::size_t>(bins));
    for (double value : samples) {
        const int index = std::clamp(static_cast<int>((value - minimum) / width), 0, bins - 1);
        ++counts[static_cast<std::size_t>(index)];
    }
    std::vector<QPointF> result;
    result.reserve(counts.size());
    for (int i = 0; i < bins; ++i)
        result.emplace_back(minimum + (i + 0.5) * width,
                            static_cast<double>(counts[static_cast<std::size_t>(i)])
                            / (samples.size() * width));
    return result;
}

} // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    buildUi();
    buildMenus();
    applyTheme();
    resize(1440, 860);
    setMinimumSize(1050, 680);
    setWindowTitle(QStringLiteral("Stochia · 随机空间"));

    std::string error;
    model_.addDistribution("X", stochia::createDistribution("normal", {{"mu", 0}, {"sigma", 1}}), &error);
    model_.addDistribution("Y", stochia::createDistribution("exponential", {{"lambda", 2}}), &error);
    model_.addTransformation("Z", "X + Y", &error);
    refreshVariableList();
    selectVariable(QStringLiteral("X"));
    markModified(false);
}

void MainWindow::buildUi() {
    auto* splitter = new QSplitter(Qt::Horizontal);
    splitter->setChildrenCollapsible(false);
    splitter->addWidget(buildVariablePanel());
    splitter->addWidget(buildVisualizationPanel());
    splitter->addWidget(buildInformationPanel());
    splitter->setStretchFactor(0, 0);
    splitter->setStretchFactor(1, 1);
    splitter->setStretchFactor(2, 0);
    splitter->setSizes({270, 850, 310});
    setCentralWidget(splitter);

    statusLabel_ = new QLabel(QStringLiteral("就绪 · Ready"));
    statusBar()->addWidget(statusLabel_, 1);
    statusBar()->setSizeGripEnabled(false);
}

QWidget* MainWindow::buildVariablePanel() {
    auto* panel = new QWidget;
    panel->setObjectName("sidePanel");
    panel->setMinimumWidth(230);
    panel->setMaximumWidth(360);
    auto* layout = new QVBoxLayout(panel);
    layout->setContentsMargins(18, 20, 18, 18);
    layout->setSpacing(12);

    auto* brandRow = new QHBoxLayout;
    auto* mark = new QLabel(QStringLiteral("S"));
    mark->setObjectName("brandMark");
    auto* brandText = new QVBoxLayout;
    auto* name = new QLabel(QStringLiteral("Stochia"));
    name->setObjectName("brandName");
    auto* subtitle = new QLabel(QStringLiteral("随机空间 · Probability Lab"));
    subtitle->setObjectName("brandSubtitle");
    brandText->addWidget(name);
    brandText->addWidget(subtitle);
    brandRow->addWidget(mark);
    brandRow->addLayout(brandText);
    brandRow->addStretch();
    layout->addLayout(brandRow);

    projectLabel_ = new QLabel(QStringLiteral("未命名实验"));
    projectLabel_->setObjectName("projectLabel");
    layout->addWidget(projectLabel_);
    layout->addWidget(sectionLabel(QStringLiteral("随机变量 · VARIABLES")));

    variableList_ = new QListWidget;
    variableList_->setObjectName("variableList");
    variableList_->setSpacing(4);
    variableList_->setContextMenuPolicy(Qt::CustomContextMenu);
    layout->addWidget(variableList_, 1);
    connect(variableList_, &QListWidget::currentTextChanged, this, [this](const QString& text) {
        const int separator = text.indexOf('\n');
        if (!text.isEmpty()) selectVariable(separator < 0 ? text : text.left(separator));
    });

    auto* addButton = new QPushButton(QStringLiteral("＋  添加分布"));
    addButton->setObjectName("primaryButton");
    auto* editButton = new QPushButton(QStringLiteral("编辑所选"));
    editButton->setObjectName("subtleButton");
    auto* removeButton = new QPushButton(QStringLiteral("移除所选"));
    removeButton->setObjectName("subtleButton");
    layout->addWidget(addButton);
    layout->addWidget(editButton);
    layout->addWidget(removeButton);
    connect(addButton, &QPushButton::clicked, this, [this] { addDistribution(); });
    connect(editButton, &QPushButton::clicked, this, [this] { editSelected(); });
    connect(removeButton, &QPushButton::clicked, this, [this] { removeSelected(); });
    return panel;
}

QWidget* MainWindow::buildVisualizationPanel() {
    auto* panel = new QWidget;
    panel->setObjectName("mainPanel");
    auto* layout = new QVBoxLayout(panel);
    layout->setContentsMargins(26, 20, 26, 18);
    layout->setSpacing(14);

    auto* titleRow = new QHBoxLayout;
    auto* titleColumn = new QVBoxLayout;
    plotTitle_ = new QLabel(QStringLiteral("选择一个随机变量"));
    plotTitle_->setObjectName("plotTitle");
    typeBadge_ = new QLabel(QStringLiteral("—"));
    typeBadge_->setObjectName("typeBadge");
    titleColumn->addWidget(plotTitle_);
    titleColumn->addWidget(typeBadge_, 0, Qt::AlignLeft);
    titleRow->addLayout(titleColumn);
    titleRow->addStretch();

    viewCombo_ = new QComboBox;
    viewCombo_->addItem(QStringLiteral("密度 / PMF"), QStringLiteral("density"));
    viewCombo_->addItem(QStringLiteral("累积分布 CDF"), QStringLiteral("cdf"));
    viewCombo_->addItem(QStringLiteral("模拟对照"), QStringLiteral("simulation"));
    viewCombo_->setMinimumWidth(160);
    titleRow->addWidget(viewCombo_);
    layout->addLayout(titleRow);
    connect(viewCombo_, &QComboBox::currentIndexChanged, this, [this] { refreshPlot(); });

    plotWidget_ = new PlotWidget;
    plotWidget_->setObjectName("plotCard");
    layout->addWidget(plotWidget_, 1);

    auto* transformCard = new QFrame;
    transformCard->setObjectName("transformCard");
    auto* transformLayout = new QHBoxLayout(transformCard);
    transformLayout->setContentsMargins(16, 13, 16, 13);
    transformLayout->setSpacing(10);
    auto* equals = new QLabel(QStringLiteral("构造"));
    equals->setObjectName("transformLabel");
    transformName_ = new QLineEdit(QStringLiteral("W"));
    transformName_->setMaximumWidth(70);
    transformName_->setAlignment(Qt::AlignCenter);
    auto* equalSign = new QLabel(QStringLiteral("="));
    transformExpression_ = new QLineEdit;
    transformExpression_->setPlaceholderText(
        QStringLiteral("例如 (X-E(X))/sqrt(Var(X)) 或 iid_sum(X,10)"));
    auto* defineButton = new QPushButton(QStringLiteral("定义变量"));
    defineButton->setObjectName("primaryButton");
    transformLayout->addWidget(equals);
    transformLayout->addWidget(transformName_);
    transformLayout->addWidget(equalSign);
    transformLayout->addWidget(transformExpression_, 1);
    transformLayout->addWidget(defineButton);
    layout->addWidget(transformCard);
    connect(defineButton, &QPushButton::clicked, this, [this] { addTransformation(); });
    connect(transformExpression_, &QLineEdit::returnPressed, this, [this] { addTransformation(); });
    return panel;
}

QWidget* MainWindow::buildInformationPanel() {
    auto* scroll = new QScrollArea;
    scroll->setObjectName("infoScroll");
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setMinimumWidth(270);
    scroll->setMaximumWidth(390);
    auto* panel = new QWidget;
    panel->setObjectName("infoPanel");
    auto* layout = new QVBoxLayout(panel);
    layout->setContentsMargins(18, 20, 18, 20);
    layout->setSpacing(12);

    layout->addWidget(sectionLabel(QStringLiteral("数学信息 · INSPECTOR")));
    layout->addWidget(metricCard(QStringLiteral("期望  E(X)"), expectationValue_));
    layout->addWidget(metricCard(QStringLiteral("方差  Var(X)"), varianceValue_));
    layout->addWidget(metricCard(QStringLiteral("模拟均值  x̄"), sampleValue_));

    layout->addWidget(sectionLabel(QStringLiteral("分布公式")));
    formulaView_ = new QTextBrowser;
    formulaView_->setObjectName("mathCard");
    formulaView_->setMinimumHeight(118);
    formulaView_->setMaximumHeight(170);
    layout->addWidget(formulaView_);

    layout->addWidget(sectionLabel(QStringLiteral("推导 / 计算路径")));
    derivationView_ = new QTextBrowser;
    derivationView_->setObjectName("mathCard");
    derivationView_->setMinimumHeight(150);
    layout->addWidget(derivationView_);

    layout->addWidget(sectionLabel(QStringLiteral("蒙特卡洛")));
    auto* simulationRow = new QHBoxLayout;
    sampleCount_ = new QSpinBox;
    sampleCount_->setRange(100, 1000000);
    sampleCount_->setSingleStep(1000);
    sampleCount_->setValue(20000);
    sampleCount_->setSuffix(QStringLiteral(" 次"));
    auto* runButton = new QPushButton(QStringLiteral("运行模拟"));
    runButton->setObjectName("primaryButton");
    simulationRow->addWidget(sampleCount_, 1);
    simulationRow->addWidget(runButton);
    layout->addLayout(simulationRow);
    layout->addStretch();
    connect(runButton, &QPushButton::clicked, this, [this] { simulate(); });

    scroll->setWidget(panel);
    return scroll;
}

void MainWindow::buildMenus() {
    auto* toolbar = addToolBar(QStringLiteral("项目工具"));
    toolbar->setMovable(false);
    toolbar->setToolButtonStyle(Qt::ToolButtonTextOnly);

    auto* newAction = new QAction(QStringLiteral("新建"), this);
    newAction->setShortcut(QKeySequence::New);
    auto* openAction = new QAction(QStringLiteral("打开"), this);
    openAction->setShortcut(QKeySequence::Open);
    auto* saveAction = new QAction(QStringLiteral("保存"), this);
    saveAction->setShortcut(QKeySequence::Save);
    auto* saveAsAction = new QAction(QStringLiteral("另存为"), this);
    saveAsAction->setShortcut(QKeySequence::SaveAs);
    auto* addAction = new QAction(QStringLiteral("添加分布"), this);
    addAction->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_D));
    auto* simulateAction = new QAction(QStringLiteral("运行模拟"), this);
    simulateAction->setShortcut(Qt::Key_F5);

    toolbar->addAction(newAction);
    toolbar->addAction(openAction);
    toolbar->addAction(saveAction);
    toolbar->addSeparator();
    toolbar->addAction(addAction);
    toolbar->addAction(simulateAction);

    auto* fileMenu = menuBar()->addMenu(QStringLiteral("文件"));
    fileMenu->addActions({newAction, openAction, saveAction, saveAsAction});
    auto* modelMenu = menuBar()->addMenu(QStringLiteral("模型"));
    modelMenu->addActions({addAction, simulateAction});
    auto* editAction = modelMenu->addAction(QStringLiteral("编辑所选变量"));
    editAction->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_E));

    auto* helpMenu = menuBar()->addMenu(QStringLiteral("帮助"));
    auto* guideAction = helpMenu->addAction(QStringLiteral("使用指南"));
    guideAction->setShortcut(QKeySequence::HelpContents);
    auto* syntaxAction = helpMenu->addAction(QStringLiteral("表达式与 IID 语法"));
    helpMenu->addSeparator();
    auto* aboutAction = helpMenu->addAction(QStringLiteral("关于 Stochia"));

    connect(newAction, &QAction::triggered, this, [this] { newProject(); });
    connect(openAction, &QAction::triggered, this, [this] { openProject(); });
    connect(saveAction, &QAction::triggered, this, [this] { saveProject(); });
    connect(saveAsAction, &QAction::triggered, this, [this] { saveProjectAs(); });
    connect(addAction, &QAction::triggered, this, [this] { addDistribution(); });
    connect(simulateAction, &QAction::triggered, this, [this] { simulate(); });
    connect(editAction, &QAction::triggered, this, [this] { editSelected(); });
    connect(guideAction, &QAction::triggered, this, [this] { showHelp(QStringLiteral("guide")); });
    connect(syntaxAction, &QAction::triggered, this, [this] { showHelp(QStringLiteral("syntax")); });
    connect(aboutAction, &QAction::triggered, this, [this] { showHelp(QStringLiteral("about")); });
}

void MainWindow::applyTheme() {
    qApp->setStyleSheet(QStringLiteral(R"(
        * { font-family: "Segoe UI", "Microsoft YaHei UI", sans-serif; font-size: 13px; color: #283044; }
        QMainWindow, #mainPanel { background: #f4f6fb; }
        QMenuBar { background: #ffffff; border-bottom: 1px solid #e4e7ef; padding: 2px; }
        QMenuBar::item { padding: 6px 10px; border-radius: 5px; }
        QMenuBar::item:selected { background: #eef1ff; color: #5366d8; }
        QMenu { background: white; border: 1px solid #dde1eb; padding: 5px; }
        QMenu::item { padding: 7px 28px; border-radius: 4px; }
        QMenu::item:selected { background: #eef1ff; }
        QToolBar { background: white; border: 0; border-bottom: 1px solid #e3e6ee; spacing: 4px; padding: 5px 10px; }
        QToolButton { padding: 6px 11px; border-radius: 6px; }
        QToolButton:hover { background: #eef1ff; color: #5366d8; }
        #sidePanel, #infoPanel, #infoScroll { background: #ffffff; }
        #sidePanel { border-right: 1px solid #e3e6ee; }
        #infoPanel { border-left: 1px solid #e3e6ee; }
        #brandMark { background: #596be3; color: white; border-radius: 10px; font-size: 20px; font-weight: 700;
                     min-width: 42px; max-width: 42px; min-height: 42px; max-height: 42px; qproperty-alignment: AlignCenter; }
        #brandName { font-size: 20px; font-weight: 700; color: #20263a; }
        #brandSubtitle { color: #8b92a4; font-size: 11px; }
        #projectLabel { background: #f5f6fb; border-radius: 7px; padding: 8px 10px; color: #646d82; }
        #sectionLabel { color: #8b92a4; font-size: 10px; font-weight: 700; padding-top: 6px; }
        #variableList { border: 0; background: transparent; outline: none; }
        #variableList::item { padding: 10px 11px; border-radius: 8px; }
        #variableList::item:selected { background: #edf0ff; color: #485bd2; }
        #variableList::item:hover:!selected { background: #f5f6fa; }
        QPushButton { border: 0; border-radius: 7px; padding: 8px 13px; }
        #primaryButton { background: #596be3; color: white; font-weight: 600; }
        #primaryButton:hover { background: #495bd1; }
        #primaryButton:pressed { background: #3f4fbc; }
        #subtleButton { background: #f1f3f8; color: #687087; }
        #subtleButton:hover { background: #e8eaf2; }
        #plotTitle { font-size: 24px; font-weight: 700; color: #20263a; }
        #typeBadge { color: #5366d8; background: #e9edff; border-radius: 5px; padding: 3px 8px; font-size: 11px; }
        #plotCard, #transformCard, #metricCard, #mathCard { background: white; border: 1px solid #e1e5ee; border-radius: 10px; }
        #transformLabel { font-weight: 700; color: #596273; }
        #metricTitle { color: #8b92a4; font-size: 11px; }
        #metricValue { color: #242b42; font-size: 20px; font-weight: 650; }
        #formulaCard { background: #f5f7fc; border-radius: 8px; padding: 12px; color: #596273; }
        #formulaCard[error="true"] { background: #fff0f1; color: #b54752; }
        #dialogHeading { font-size: 20px; font-weight: 700; }
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
            background: white; border: 1px solid #d9dde8; border-radius: 7px; padding: 7px 9px; min-height: 20px;
        }
        QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border: 1px solid #6677e5; }
        QComboBox::drop-down { border: 0; width: 24px; }
        QTextBrowser { padding: 9px; }
        QSplitter::handle { background: #e3e6ee; width: 1px; }
        QStatusBar { background: white; border-top: 1px solid #e3e6ee; }
        QStatusBar QLabel { padding-left: 8px; color: #747d91; }
        QScrollBar:vertical { width: 8px; background: transparent; }
        QScrollBar::handle:vertical { background: #cdd2df; border-radius: 4px; min-height: 28px; }
    )"));
}

QString MainWindow::nextVariableName() const {
    for (QChar letter = 'X'; letter <= 'Z'; letter = QChar(letter.unicode() + 1))
        if (!model_.find(QString(letter).toStdString())) return QString(letter);
    for (int index = 1; index < 1000; ++index) {
        const QString candidate = QStringLiteral("X%1").arg(index);
        if (!model_.find(candidate.toStdString())) return candidate;
    }
    return QStringLiteral("V");
}

void MainWindow::addDistribution() {
    VariableDialog dialog(this);
    if (dialog.exec() != QDialog::Accepted) return;
    std::string error;
    auto distribution = stochia::createDistribution(dialog.distributionId(), dialog.parameters(), &error);
    if (!model_.addDistribution(dialog.variableName().toStdString(), distribution, &error)) {
        QMessageBox::warning(this, QStringLiteral("无法添加变量"), QString::fromStdString(error));
        return;
    }
    refreshVariableList();
    selectVariable(dialog.variableName());
    markModified();
    showStatus(QStringLiteral("已添加随机变量 %1").arg(dialog.variableName()));
}

void MainWindow::addTransformation() {
    const QString name = transformName_->text().trimmed();
    const QString source = transformExpression_->text().trimmed();
    std::string error;
    if (!model_.addTransformation(name.toStdString(), source.toStdString(), &error)) {
        showStatus(QString::fromStdString(error), true);
        QMessageBox::warning(this, QStringLiteral("表达式无效"), QString::fromStdString(error));
        return;
    }
    refreshVariableList();
    selectVariable(name);
    transformExpression_->clear();
    transformName_->setText(nextVariableName());
    markModified();
    showStatus(QStringLiteral("已构造 %1 = %2").arg(name, source));
}

void MainWindow::editSelected() {
    const auto* variable = model_.find(selectedName_.toStdString());
    if (!variable) return;
    std::string error;
    if (variable->distribution) {
        VariableDialog dialog(this);
        dialog.setWindowTitle(QStringLiteral("编辑随机变量"));
        dialog.setDefinition(selectedName_, variable->distribution->id(),
                             variable->distribution->parameters(), true);
        if (dialog.exec() != QDialog::Accepted) return;
        auto distribution = stochia::createDistribution(
            dialog.distributionId(), dialog.parameters(), &error);
        if (!model_.addDistribution(selectedName_.toStdString(), distribution, &error)) {
            QMessageBox::warning(this, QStringLiteral("无法更新变量"), QString::fromStdString(error));
            return;
        }
    } else {
        bool accepted = false;
        const QString expression = QInputDialog::getText(
            this, QStringLiteral("编辑变换变量"),
            QStringLiteral("%1 =").arg(selectedName_),
            QLineEdit::Normal, QString::fromStdString(variable->expression->source()), &accepted);
        if (!accepted) return;
        if (!model_.addTransformation(selectedName_.toStdString(),
                                      expression.trimmed().toStdString(), &error)) {
            QMessageBox::warning(this, QStringLiteral("无法更新变量"), QString::fromStdString(error));
            return;
        }
    }
    refreshVariableList();
    selectVariable(selectedName_);
    markModified();
    showStatus(QStringLiteral("已更新随机变量 %1，并重新计算所有下游解析分布").arg(selectedName_));
}

void MainWindow::removeSelected() {
    if (selectedName_.isEmpty()) return;
    if (!model_.remove(selectedName_.toStdString())) {
        QMessageBox::information(this, QStringLiteral("无法移除"),
                                 QStringLiteral("其他变量仍依赖 %1，请先移除相关变换。").arg(selectedName_));
        return;
    }
    selectedName_.clear();
    refreshVariableList();
    const auto names = model_.orderedNames();
    if (!names.empty()) selectVariable(QString::fromStdString(names.front()));
    else {
        plotWidget_->setSeries({});
        refreshInformation();
    }
    markModified();
}

QString MainWindow::variableDescription(const stochia::RandomVariable& variable) const {
    if (variable.expression)
        return QStringLiteral("%1\n  = %2").arg(QString::fromStdString(variable.name),
                                                QString::fromStdString(variable.expression->source()));
    QStringList arguments;
    for (const auto& [key, value] : variable.distribution->parameters())
        arguments << QStringLiteral("%1=%2").arg(QString::fromStdString(key), formatValue(value));
    return QStringLiteral("%1\n  ~ %2(%3)")
        .arg(QString::fromStdString(variable.name),
             QString::fromStdString(variable.distribution->id()),
             arguments.join(QStringLiteral(", ")));
}

void MainWindow::refreshVariableList() {
    const QString keep = selectedName_;
    variableList_->blockSignals(true);
    variableList_->clear();
    for (const auto& name : model_.orderedNames()) {
        const auto* variable = model_.find(name);
        auto* item = new QListWidgetItem(variableDescription(*variable));
        item->setData(Qt::UserRole, QString::fromStdString(name));
        item->setSizeHint(QSize(0, 54));
        variableList_->addItem(item);
        if (QString::fromStdString(name) == keep) variableList_->setCurrentItem(item);
    }
    variableList_->blockSignals(false);
    transformName_->setText(nextVariableName());
}

void MainWindow::selectVariable(const QString& name) {
    const auto* variable = model_.find(name.toStdString());
    if (!variable) return;
    selectedName_ = name;
    for (int i = 0; i < variableList_->count(); ++i) {
        auto* item = variableList_->item(i);
        if (item->data(Qt::UserRole).toString() == name) {
            variableList_->blockSignals(true);
            variableList_->setCurrentItem(item);
            variableList_->blockSignals(false);
            break;
        }
    }
    plotTitle_->setText(variableDescription(*variable).replace('\n', QStringLiteral("  ")));
    if (variable->distribution) {
        typeBadge_->setText(variable->distribution->type() == stochia::DistributionType::Continuous
                                ? QStringLiteral("连续型 · CONTINUOUS")
                                : QStringLiteral("离散型 · DISCRETE"));
    } else {
        if (variable->analyticalDistribution) {
            typeBadge_->setText(variable->closedForm
                ? QStringLiteral("变换变量 · 解析闭式")
                : QStringLiteral("变换变量 · 理论数值解"));
        } else {
            typeBadge_->setText(QStringLiteral("变换变量 · MONTE CARLO"));
        }
    }
    refreshPlot();
    refreshInformation();
}

void MainWindow::refreshPlot() {
    const auto* variable = model_.find(selectedName_.toStdString());
    if (!variable) {
        plotWidget_->setSeries({});
        return;
    }
    const QString mode = viewCombo_->currentData().toString();
    std::vector<PlotWidget::Series> series;

    const auto* theoryDistribution = variable->theory();
    if (theoryDistribution && mode != QStringLiteral("simulation")) {
        const auto range = theoryDistribution->plotRange();
        const bool discrete = theoryDistribution->type() == stochia::DistributionType::Discrete;
        PlotWidget::Series theory;
        theory.name = mode == QStringLiteral("cdf") ? QStringLiteral("理论 CDF") :
                                                     (discrete ? QStringLiteral("理论 PMF") : QStringLiteral("理论 PDF"));
        theory.color = QColor("#596be3");
        theory.style = discrete ? PlotWidget::Style::Bars
                                : (mode == QStringLiteral("cdf") ? PlotWidget::Style::Curve : PlotWidget::Style::Curve);
        theory.filled = !discrete && mode == QStringLiteral("density");
        if (discrete) {
            const int start = static_cast<int>(std::ceil(range.first));
            const int end = static_cast<int>(std::floor(range.second));
            for (int k = start; k <= end; ++k)
                theory.points.emplace_back(k, mode == QStringLiteral("cdf")
                                                 ? theoryDistribution->cdf(k)
                                                 : theoryDistribution->density(k));
            if (mode == QStringLiteral("cdf")) theory.style = PlotWidget::Style::Steps;
        } else {
            constexpr int points = 500;
            for (int i = 0; i <= points; ++i) {
                const double x = range.first + (range.second - range.first) * i / points;
                theory.points.emplace_back(x, mode == QStringLiteral("cdf")
                                                  ? theoryDistribution->cdf(x)
                                                  : theoryDistribution->density(x));
            }
        }
        series.push_back(std::move(theory));
        plotWidget_->setAxisLabels(QStringLiteral("x"),
                                   mode == QStringLiteral("cdf") ? QStringLiteral("F(x)")
                                                                 : (discrete ? QStringLiteral("P(X=x)")
                                                                             : QStringLiteral("f(x)")));
    } else {
        std::string error;
        const std::size_t quickCount = std::min(20000, sampleCount_->value());
        auto samples = model_.simulate(selectedName_.toStdString(), quickCount, 0x57A0C41AULL, &error);
        if (samples.empty()) {
            showStatus(QString::fromStdString(error), true);
            plotWidget_->setSeries({});
            return;
        }
        const bool discrete = theoryDistribution
                              && theoryDistribution->type() == stochia::DistributionType::Discrete;
        PlotWidget::Series empirical;
        empirical.name = QStringLiteral("模拟估计");
        empirical.color = QColor("#ef8b4c");
        empirical.style = PlotWidget::Style::Bars;
        empirical.points = histogram(samples, discrete);
        series.push_back(std::move(empirical));

        if (theoryDistribution) {
            PlotWidget::Series theory;
            theory.name = discrete ? QStringLiteral("理论 PMF") : QStringLiteral("理论 PDF");
            theory.color = QColor("#596be3");
            theory.style = discrete ? PlotWidget::Style::Bars : PlotWidget::Style::Curve;
            const auto range = theoryDistribution->plotRange();
            if (discrete) {
                for (int k = static_cast<int>(std::ceil(range.first));
                     k <= static_cast<int>(std::floor(range.second)); ++k)
                    theory.points.emplace_back(k, theoryDistribution->density(k));
            } else {
                for (int i = 0; i <= 500; ++i) {
                    const double x = range.first + (range.second - range.first) * i / 500.0;
                    theory.points.emplace_back(x, theoryDistribution->density(x));
                }
            }
            series.insert(series.begin(), std::move(theory));
        }
        plotWidget_->setAxisLabels(QStringLiteral("x"), discrete ? QStringLiteral("概率") : QStringLiteral("密度"));
        refreshInformation(&samples);
    }
    plotWidget_->setSeries(std::move(series));
}

void MainWindow::refreshInformation(const std::vector<double>* samples) {
    const auto* variable = model_.find(selectedName_.toStdString());
    if (!variable) {
        expectationValue_->setText(QStringLiteral("—"));
        varianceValue_->setText(QStringLiteral("—"));
        sampleValue_->setText(QStringLiteral("—"));
        formulaView_->clear();
        derivationView_->clear();
        return;
    }
    if (variable->distribution) {
        expectationValue_->setText(formatValue(variable->distribution->mean()));
        varianceValue_->setText(formatValue(variable->distribution->variance()));
        formulaView_->setText(QStringLiteral("<p style='font-family:Cambria Math;font-size:15px'>%1</p>")
                                  .arg(QString::fromStdString(variable->distribution->formula()).toHtmlEscaped()));
        derivationView_->setHtml(QStringLiteral(
            "<p><b>模型路径</b></p>"
            "<p>参数 → 理论 PDF / PMF → CDF</p>"
            "<p>采样器：C++ &lt;random&gt;</p>"
            "<p style='color:#7d8598'>模拟结果可与理论曲线直接比较。</p>"));
    } else {
        if (variable->analyticalDistribution && variable->closedForm) {
            expectationValue_->setText(formatValue(variable->analyticalDistribution->mean()));
            varianceValue_->setText(formatValue(variable->analyticalDistribution->variance()));
        } else {
            expectationValue_->setText(variable->analyticalDistribution
                                            ? QStringLiteral("由理论密度积分")
                                            : QStringLiteral("数值估计"));
            varianceValue_->setText(variable->analyticalDistribution
                                         ? QStringLiteral("由理论密度积分")
                                         : QStringLiteral("数值估计"));
        }
        QString theoryHtml;
        if (variable->analyticalDistribution) {
            theoryHtml = QStringLiteral(
                "<p style='color:#5366d8'><b>%1</b></p>"
                "<p style='font-family:Cambria Math'>%2</p>")
                .arg(QString::fromStdString(variable->analyticalName).toHtmlEscaped(),
                     QString::fromStdString(variable->analyticalDistribution->formula()).toHtmlEscaped());
        }
        formulaView_->setHtml(QStringLiteral(
            "<p style='font-family:Cambria Math;font-size:16px'><b>%1</b> = %2</p>%3")
            .arg(selectedName_.toHtmlEscaped(),
                 QString::fromStdString(variable->expression->source()).toHtmlEscaped(),
                 theoryHtml));
        const auto& deps = variable->expression->variables();
        QString dependencyText;
        for (const auto& dep : deps) {
            if (!dependencyText.isEmpty()) dependencyText += QStringLiteral(", ");
            dependencyText += QString::fromStdString(dep);
        }
        QString derivation;
        if (variable->analyticalDistribution) {
            derivation = QStringLiteral(
                "<p><b>%1</b></p><p>%2</p>"
                "<p style='color:#7d8598'>“密度 / PMF”和“CDF”显示理论结果；“模拟对照”显示蒙特卡洛验证。</p>")
                .arg(variable->closedForm ? QStringLiteral("已识别解析闭式")
                                          : QStringLiteral("已建立理论积分"),
                     QString::fromStdString(variable->derivation).toHtmlEscaped());
        } else {
            derivation = QStringLiteral(
                "<p><b>变量变换</b></p><p>输入：%1</p>"
                "<p>逐样本计算 g(%1)，再估计经验分布。</p>"
                "<p style='color:#7d8598'>当不存在稳定解析形式时自动采用蒙特卡洛。</p>").arg(dependencyText);
        }
        derivationView_->setHtml(derivation);
    }
    if (samples && !samples->empty()) {
        const auto stats = stochia::ProbabilityModel::summarize(*samples);
        sampleValue_->setText(formatValue(stats.mean));
        if (variable->isTransformation()) {
            expectationValue_->setText(formatValue(stats.mean));
            varianceValue_->setText(formatValue(stats.variance));
        }
    } else {
        sampleValue_->setText(QStringLiteral("尚未运行"));
    }
}

void MainWindow::simulate() {
    if (selectedName_.isEmpty()) return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const auto seed = static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::string error;
    auto samples = model_.simulate(selectedName_.toStdString(),
                                   static_cast<std::size_t>(sampleCount_->value()), seed, &error);
    QApplication::restoreOverrideCursor();
    if (samples.empty()) {
        QMessageBox::warning(this, QStringLiteral("模拟失败"), QString::fromStdString(error));
        return;
    }
    viewCombo_->setCurrentIndex(2);
    const auto* variable = model_.find(selectedName_.toStdString());
    const auto* theoryDistribution = variable->theory();
    const bool discrete = theoryDistribution
                          && theoryDistribution->type() == stochia::DistributionType::Discrete;
    std::vector<PlotWidget::Series> series;
    if (theoryDistribution) {
        PlotWidget::Series theory;
        theory.name = discrete ? QStringLiteral("理论 PMF") : QStringLiteral("理论 PDF");
        theory.color = QColor("#596be3");
        theory.style = discrete ? PlotWidget::Style::Bars : PlotWidget::Style::Curve;
        const auto range = theoryDistribution->plotRange();
        if (discrete) {
            for (int k = static_cast<int>(std::ceil(range.first));
                 k <= static_cast<int>(std::floor(range.second)); ++k)
                theory.points.emplace_back(k, theoryDistribution->density(k));
        } else {
            for (int i = 0; i <= 500; ++i) {
                const double x = range.first + (range.second - range.first) * i / 500.0;
                theory.points.emplace_back(x, theoryDistribution->density(x));
            }
        }
        series.push_back(std::move(theory));
    }
    series.push_back({QStringLiteral("模拟估计"), QColor("#ef8b4c"), PlotWidget::Style::Bars,
                      histogram(samples, discrete), false});
    plotWidget_->setSeries(std::move(series));
    refreshInformation(&samples);
    const auto stats = stochia::ProbabilityModel::summarize(samples);
    showStatus(QStringLiteral("完成 %1 次模拟 · 均值 %2 · 方差 %3")
                   .arg(samples.size()).arg(formatValue(stats.mean), formatValue(stats.variance)));
}

void MainWindow::newProject() {
    if (!maybeSave()) return;
    model_.clear();
    currentPath_.clear();
    selectedName_.clear();
    projectLabel_->setText(QStringLiteral("未命名实验"));
    refreshVariableList();
    plotWidget_->setSeries({});
    refreshInformation();
    markModified(false);
}

void MainWindow::openProject() {
    if (!maybeSave()) return;
    const QString path = QFileDialog::getOpenFileName(this, QStringLiteral("打开 Stochia 项目"),
                                                      QString(), QStringLiteral("Stochia 项目 (*.pvis)"));
    if (!path.isEmpty()) readProject(path);
}

bool MainWindow::saveProject() {
    return currentPath_.isEmpty() ? saveProjectAs() : writeProject(currentPath_);
}

bool MainWindow::saveProjectAs() {
    QString path = QFileDialog::getSaveFileName(this, QStringLiteral("保存 Stochia 项目"),
                                                currentPath_.isEmpty() ? QStringLiteral("MyProbabilityModel.pvis")
                                                                       : currentPath_,
                                                QStringLiteral("Stochia 项目 (*.pvis)"));
    if (path.isEmpty()) return false;
    if (!path.endsWith(QStringLiteral(".pvis"), Qt::CaseInsensitive)) path += QStringLiteral(".pvis");
    return writeProject(path);
}

bool MainWindow::writeProject(const QString& path) {
    QJsonArray variables;
    for (const auto& name : model_.orderedNames()) {
        const auto* variable = model_.find(name);
        QJsonObject object;
        object["name"] = QString::fromStdString(name);
        if (variable->distribution) {
            object["kind"] = QStringLiteral("distribution");
            object["distribution"] = QString::fromStdString(variable->distribution->id());
            QJsonObject parameters;
            for (const auto& [key, value] : variable->distribution->parameters())
                parameters[QString::fromStdString(key)] = value;
            object["parameters"] = parameters;
        } else {
            object["kind"] = QStringLiteral("transformation");
            object["expression"] = QString::fromStdString(variable->expression->source());
        }
        variables.append(object);
    }
    QJsonObject root;
    root["format"] = QStringLiteral("stochia-project");
    root["version"] = 1;
    root["variables"] = variables;
    root["selected"] = selectedName_;
    root["view"] = viewCombo_->currentData().toString();
    root["sampleCount"] = sampleCount_->value();

    QFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        QMessageBox::critical(this, QStringLiteral("保存失败"), file.errorString());
        return false;
    }
    if (file.write(QJsonDocument(root).toJson(QJsonDocument::Indented)) < 0) {
        QMessageBox::critical(this, QStringLiteral("保存失败"), file.errorString());
        return false;
    }
    currentPath_ = path;
    projectLabel_->setText(QFileInfo(path).baseName());
    markModified(false);
    showStatus(QStringLiteral("项目已保存到 %1").arg(path));
    return true;
}

bool MainWindow::readProject(const QString& path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(this, QStringLiteral("打开失败"), file.errorString());
        return false;
    }
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
        QMessageBox::critical(this, QStringLiteral("项目格式错误"), parseError.errorString());
        return false;
    }
    const auto root = document.object();
    if (root.value("format").toString() != QStringLiteral("stochia-project")) {
        QMessageBox::critical(this, QStringLiteral("项目格式错误"), QStringLiteral("这不是有效的 Stochia 项目。"));
        return false;
    }

    stochia::ProbabilityModel loaded;
    for (const auto& value : root.value("variables").toArray()) {
        const auto object = value.toObject();
        const std::string name = object.value("name").toString().toStdString();
        std::string error;
        bool ok = false;
        if (object.value("kind").toString() == QStringLiteral("distribution")) {
            std::map<std::string, double> parameters;
            const auto jsonParameters = object.value("parameters").toObject();
            for (auto it = jsonParameters.begin(); it != jsonParameters.end(); ++it)
                parameters[it.key().toStdString()] = it.value().toDouble();
            auto distribution = stochia::createDistribution(
                object.value("distribution").toString().toStdString(), parameters, &error);
            ok = loaded.addDistribution(name, distribution, &error);
        } else {
            ok = loaded.addTransformation(name, object.value("expression").toString().toStdString(), &error);
        }
        if (!ok) {
            QMessageBox::critical(this, QStringLiteral("项目内容错误"), QString::fromStdString(error));
            return false;
        }
    }

    model_ = std::move(loaded);
    currentPath_ = path;
    projectLabel_->setText(QFileInfo(path).baseName());
    sampleCount_->setValue(root.value("sampleCount").toInt(20000));
    const QString view = root.value("view").toString();
    const int viewIndex = viewCombo_->findData(view);
    if (viewIndex >= 0) viewCombo_->setCurrentIndex(viewIndex);
    refreshVariableList();
    QString selected = root.value("selected").toString();
    if (!model_.find(selected.toStdString())) {
        const auto names = model_.orderedNames();
        selected = names.empty() ? QString() : QString::fromStdString(names.front());
    }
    if (!selected.isEmpty()) selectVariable(selected);
    else {
        selectedName_.clear();
        plotWidget_->setSeries({});
        refreshInformation();
    }
    markModified(false);
    showStatus(QStringLiteral("已打开 %1").arg(path));
    return true;
}

bool MainWindow::maybeSave() {
    if (!modified_) return true;
    const auto result = QMessageBox::question(
        this, QStringLiteral("保存更改"),
        QStringLiteral("当前实验有尚未保存的更改。是否保存？"),
        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel,
        QMessageBox::Save);
    if (result == QMessageBox::Cancel) return false;
    return result == QMessageBox::Discard || saveProject();
}

void MainWindow::markModified(bool modified) {
    modified_ = modified;
    QString title = QStringLiteral("Stochia · 随机空间");
    if (!currentPath_.isEmpty()) title += QStringLiteral(" — %1").arg(QFileInfo(currentPath_).baseName());
    if (modified_) title += QStringLiteral(" *");
    setWindowTitle(title);
}

void MainWindow::showStatus(const QString& text, bool error) {
    statusLabel_->setText(text);
    statusLabel_->setStyleSheet(error ? QStringLiteral("color:#b54752;") : QString());
}

void MainWindow::showHelp(const QString& topic) {
    QDialog dialog(this);
    dialog.setWindowTitle(topic == QStringLiteral("about")
                              ? QStringLiteral("关于 Stochia")
                              : QStringLiteral("Stochia 帮助中心"));
    dialog.resize(760, 620);
    auto* layout = new QVBoxLayout(&dialog);
    layout->setContentsMargins(22, 20, 22, 18);
    auto* browser = new QTextBrowser;
    browser->setOpenExternalLinks(true);

    if (topic == QStringLiteral("about")) {
        browser->setHtml(QStringLiteral(
            "<h1>Stochia · 随机空间</h1>"
            "<p><b>让随机变量成为可以被构造、操作和观察的数学对象。</b></p>"
            "<p>版本 0.2.1 · Qt 6 概率论实验室</p>"
            "<p>支持理论分布、表达式变换、解析推断、数值积分与蒙特卡洛验证。</p>"));
    } else if (topic == QStringLiteral("syntax")) {
        browser->setHtml(QStringLiteral(
            "<h1>表达式与 IID 语法</h1>"
            "<h2>普通运算</h2>"
            "<p><code>+ - * / ^</code>，例如 "
            "<code>X+Y</code>、<code>X/Y</code>、<code>(X+Y)^2</code>。</p>"
            "<h2>初等函数</h2>"
            "<p><code>sin(X)</code>、<code>cos(X)</code>、<code>log(X)</code>、"
            "<code>exp(X)</code>、<code>sqrt(X)</code>、<code>abs(X)</code>。</p>"
            "<h2>可变参数批量运算</h2>"
            "<ul>"
            "<li><code>sum(X,Y,Z)</code>：连加</li>"
            "<li><code>product(X,Y,Z)</code> 或 <code>prod(X,Y,Z)</code>：连乘</li>"
            "<li><code>min(X,Y,Z)</code>：批量最小值</li>"
            "<li><code>max(X,Y,Z)</code>：批量最大值</li>"
            "</ul>"
            "<h2>数字特征</h2>"
            "<ul>"
            "<li><code>E(X)</code> 或 <code>mean(X)</code>：理论期望</li>"
            "<li><code>Var(X)</code>：理论方差</li>"
            "</ul>"
            "<p>例如 <code>(X-E(X))/sqrt(Var(X))</code> 会构造 X 的标准化变量。"
            "统计算子返回常数，不会额外消耗随机样本。</p>"
            "<h2>独立同分布构造</h2>"
            "<p>第一个参数是已存在的随机变量，第二个参数是正整数：</p>"
            "<ul>"
            "<li><code>iid_sum(X,10)</code>：10 个 X 的独立副本之和</li>"
            "<li><code>iid_product(X,5)</code>：5 个独立副本之积</li>"
            "<li><code>iid_min(X,20)</code>：20 个独立副本的最小值</li>"
            "<li><code>iid_max(X,20)</code>：20 个独立副本的最大值</li>"
            "</ul>"
            "<p><b>区别：</b><code>X+X</code> 使用同一个样本，等于 <code>2*X</code>；"
            "<code>iid_sum(X,2)</code> 使用两个相互独立的 X 样本。</p>"
            "<h2>解析解示例</h2>"
            "<ul>"
            "<li>若 <code>X~Exponential(lambda)</code>，则 "
            "<code>iid_sum(X,n)~Gamma(n,1/lambda)</code>。</li>"
            "<li>若 <code>U~Uniform(a,b)</code>，则 "
            "<code>iid_sum(U,n)</code> 为缩放 Irwin–Hall 分布，并使用精确有限和 PDF/CDF。</li>"
            "<li>正态变量的独立加减仍为正态分布。</li>"
            "<li><code>iid_max(X,n)</code> 的 CDF 为 <code>F_X(x)^n</code>。</li>"
            "<li><code>iid_min(X,n)</code> 的 CDF 为 <code>1-(1-F_X(x))^n</code>。</li>"
            "</ul>"
            "<p>变量名区分大小写；函数名和逗号应使用英文字符。</p>"));
    } else {
        browser->setHtml(QStringLiteral(
            "<h1>Stochia 使用指南</h1>"
            "<h2>1. 创建基础随机变量</h2>"
            "<p>点击左侧“添加分布”，输入变量名，选择分布并填写参数。</p>"
            "<h2>2. 构造变换变量</h2>"
            "<p>在中央底部填写新变量名和表达式。例如 "
            "<code>Z = sum(X,Y)</code> 或 <code>M = iid_max(X,10)</code>。</p>"
            "<h2>3. 选择计算结果</h2>"
            "<p>“密度 / PMF”和“累积分布 CDF”显示理论结果。若识别到已知闭式，"
            "右侧会标注具体分布；一般连续加减乘除和最值则使用理论积分公式。"
            "“模拟对照”用于用蒙特卡洛样本验证理论曲线。</p>"
            "<h2>4. 修改变量</h2>"
            "<p>选中左侧变量后点击“编辑所选”，可直接修改分布参数或变换表达式。"
            "下游变量的解析结果会自动重新计算。</p>"
            "<h2>5. 保存实验</h2>"
            "<p>使用“文件 → 保存”生成 <code>.pvis</code> 项目；稍后可以继续打开。</p>"
            "<h2>常见错误</h2>"
            "<ul>"
            "<li><code>log()</code> 的输入必须为正数。</li>"
            "<li>除数不能为零。</li>"
            "<li>表达式引用的变量必须已经存在。</li>"
            "<li>不允许循环依赖，例如先定义 Z=X+1，再把 X 改为 Z+1。</li>"
            "</ul>"
            "<p>完整语法请打开“帮助 → 表达式与 IID 语法”。</p>"));
    }
    layout->addWidget(browser);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    buttons->button(QDialogButtonBox::Close)->setText(QStringLiteral("关闭"));
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    layout->addWidget(buttons);
    dialog.exec();
}

void MainWindow::closeEvent(QCloseEvent* event) {
    maybeSave() ? event->accept() : event->ignore();
}
